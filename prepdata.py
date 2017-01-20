#!/usr/bin/env python
#RMS Jan 2017

import glob
import sys
import obspy
import os
import numpy as np

#SAC prep scripts
SACnormstack='/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/SAC_normstack.sh'
SACnormstack2='/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/SAC_normstack2.sh'
normalizemacro='/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/normalize.m'
weightingsscript='/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/GMT_weightings.sh'


def InitialStackCorrelate(badfileslist,event,infile):

	'''Extracts the relative arrival time residuals from AIMBAT (or dbpick) output, appends relative arrival times to the SAc headers
	and runs the Boyce SAC code, creating an initial stack and correlation files'''

	indata = open(infile,'r')
	lines = indata.readlines()
	indata.close()

	assessed_stations = []

	os.system('rm *trim*') 
	os.system('rm *stk*')
	os.system('rm *p0*')
	os.system('rm *p1*')

	for line in lines:
		#get the datafile associated with each record - just BHZ for now, but can be changed
		vals = line.split()
		stationcode = vals[3].strip()
		relativeresid = float(vals[-3])

		#Read the files and input the MCCC arrival time into t4 location in the SAC header

		sacfiles=glob.glob('*%s*.BHZ' %stationcode)

		if stationcode not in assessed_stations:

			#read the file, crop and add to the stream
			trace = obspy.read(sacfiles[0],format='SAC')
			Ppred = float(trace[0].stats.sac.user1) #predicted P travel time
			Prel = Ppred + relativeresid #the alignment time from MCCCS
			trace[0].stats.sac.t4 = Prel #assign the algnment time to a SAC header for further use

			station = trace[0].stats.station
			network = trace[0].stats.network
			channel = trace[0].stats.channel

			newname = sacfiles[0][4:]+'.p0' #removes the 'vel' flag from the filename

			trace[0].write(newname,format='SAC')

			assessed_stations.append(stationcode)

	#sys.exit(1)

 	#----------------------------------------
	#Method from Boyce et al. (2017) starts here 
	#----------------------------------------

	#Generate the stack, called stack.sac and compute the X correlation between 
	#each sacfile and the stack 

	os.system('%s %s' %(SACnormstack,normalizemacro))


def generate_weightings(badfileslist,event,scheme='RMS',function=7,SNR_cutoff=3.5,XC_cutoff=0.90,XC_time_cuttoff=0.25):

	''' Scheme can be 'RMS' or 'XC'. Function is one of the 9 weight functions provided in the code '''


	if scheme == 'RMS':

		if os.path.isfile('rms_pre_arrival_signal.out'):

			print '--------------------------------------'
			print 'Generating weightings for RMS noise'
			print '--------------------------------------'

			outfile1 = open('SNR_norm_values.txt','w')
			outfile2 = open('SNR_norm_files.txt','w')
			infile = open('rms_pre_arrival_signal.out','r')

			lines = infile.readlines()
			infile.close()

			fnames = []
			snrs = []

			for line in lines:
				vals = line.split()
				filename = vals[0]
				ampsnoise = float(vals[1])
				ampssignal = float(vals[2])
				SNR = ampssignal/ampsnoise

				if SNR > SNR_cutoff:
					fnames.append(filename)
					snrs.append(SNR)
				else:
					print 'file %s has SNR of %g. Less than the cutoff of %g' %(filename,SNR,SNR_cutoff)
					print 'file will not be included in subsequent stacking'
					badfileslist.write('%s %s\n' %(event,filename))

			#normalize the SNRs and write to file

			snrs = np.array(snrs)
			snrs_normed = snrs/np.max(snrs)

			for element in zip(fnames,snrs_normed):
				outfile1.write('%g\n' %(element[1]))
				outfile2.write('%s\n' %(element[0]))

			outfile1.close()
			outfile2.close()

			#Generate weighting values associated with each file and produce output file contining filename and weights, for stacking
			#Also create a new SAC macro that will stack the selected files according to their new weights

			os.system('%s SNR_norm_values.txt %s %i' %(weightingsscript,'SNR_weighted_values.txt',function))
			os.system('paste SNR_norm_files.txt SNR_norm_values.txt SNR_weighted_values.txt > SNR_weights.txt')
			os.system('rm SNR_norm_files.txt SNR_norm_values.txt SNR_weighted_values.txt')

			outfile = open('addstack.m','w')
			infile = open('SNR_weights.txt','r')
			lines = infile.readlines()
			infile.close()

			for line in lines:
				vals = line.split()
				filename = vals[0]
				weight = float(vals[-1])

				outfile.write('addstack %s weight %g\n' %(filename[:-4]+'.stk2',weight))

			outfile.close()

			os.system('echo "read *.BHZ.p0.stk1; write append .stk2" >> mk_stk2.m') #start writing a macro to stack the traces again
			os.system('rm *.rms')

		else:

			print 'No outfile from RMS determination found!'
			sys.exit(1)

	elif scheme == 'XC':

		if os.path.isfile('correlation.out'):

			print '--------------------------------------'
			print 'Generating weightings for XC case'
			print '--------------------------------------'

			outfile1 = open('XC_norm_values.txt','w')
			outfile2 = open('XC_norm_files.txt','w')
			infile = open('correlation.out','r')

			lines = infile.readlines()
			infile.close()

			fnames = []
			xcs = []
			atimes = []

			for line in lines:
				vals = line.split()
				filename = vals[0]
				XCval = float(vals[1])
				adjusttime = float(vals[2])

				if (XCval > XC_cutoff) and (adjusttime < XC_time_cuttoff):
					fnames.append(filename)
					xcs.append(XCval)
					atimes.append(adjusttime)
				else:
					print 'file %s has XC of %g. Less than the cutoff of %g' %(filename,XCval,XC_cutoff)
					print 'file will not be included in subsequent stacking'
					badfileslist.write('%s %s\n' %(event,filename))

			#normalize the SNRs and write to file

			xcs = np.array(xcs)
			xcs_normed = xcs/np.max(xcs)

			for element in zip(fnames,xcs_normed,atimes):

				#Adjust the SAC header to that the relative arrival time value is shifted according 
				#to the cross correlation value

				#we want to read the untrimmed version of the file and write a new, adjusted version
				fname_read = element[0][:-12]+'p0'
				trace = obspy.read(fname_read,format='SAC')
				newname = fname_read[:-3]+'.p1'
				tnew = trace[0].stats.sac.t4 + element[2]
				trace[0].stats.sac.t4 = tnew
				trace.write(newname,format='SAC')

				os.system('echo "m %s %s" >> norm_stk2.m' %(normalizemacro,newname))

				outfile1.write('%g\n' %(element[1]))
				outfile2.write('%s\n' %(element[0]))

			outfile1.close()
			outfile2.close()

			#Generate weighting values associated with each file and produce output file contining filename and weights, for stacking
			#Also create a new SAC macro that will stack the selected files according to their new weights

			os.system('%s XC_norm_values.txt %s %i' %(weightingsscript,'XC_weighted_values.txt',function))
			os.system('paste XC_norm_files.txt XC_norm_values.txt XC_weighted_values.txt > XC_weights.txt')
			os.system('rm XC_norm_files.txt XC_norm_values.txt XC_weighted_values.txt')

			outfile = open('addstack.m','w')
			infile = open('XC_weights.txt','r')
			lines = infile.readlines()
			infile.close()

			for line in lines:
				vals = line.split()
				filename = vals[0]
				weight = float(vals[-1])

				outfile.write('addstack %s weight %g\n' %(filename[:-12]+'p1.stk2',weight))

			outfile.close()

			#write required commands to the second stacking macro

			os.system('echo "read *.BHZ.p1" >> mk_stk2.m')
			os.system('echo "cuterr fillz" >> mk_stk2.m')
			os.system('echo "cut t4 -10 10" >> mk_stk2.m')
			os.system('echo "read" >> mk_stk2.m')
			os.system('echo "synchronize" >> mk_stk2.m')
			os.system('echo "interpolate delta 0.025" >> mk_stk2.m')
			os.system('echo "int" >> mk_stk2.m')
			os.system('echo "taper width 0.3" >> mk_stk2.m')
			os.system('echo "chnhdr b -10" >> mk_stk2.m')
			os.system('echo "write append .stk2" >> mk_stk2.m')
			os.system('echo "cut off" >> mk_stk2.m')
			os.system('echo "m norm_stk2.m" >> mk_stk2.m')
			os.system('rm *.corr')



	else:

		print 'Weighting scheme not recongized!'
		sys.exit(1)


def SecondStackCorrelate(event):

	print '--------------------------------------'
	print 'Begin second round of stacking/XC'
	print '--------------------------------------'

	#----------------------------------------
	#Method from Boyce et al. (2017)
	#----------------------------------------

	#Generate the second stack, called stack2.sac and compute the X correlation between 
	#each sacfile and the stack 

	os.system('%s %s' %(SACnormstack2,normalizemacro) )






def main():

	datadir = '/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/test_data'

	os.chdir(datadir)
	events = glob.glob('20*')

	badfileslist = open('bad_files.txt','w')

	for event in events:

		os.chdir(event)
		os.chdir('BH_VEL')

		targetfile = 'P_0.02.0.1_AIMBAT.out'

		if os.path.isfile(targetfile):
			traces = InitialStackCorrelate(badfileslist,event,targetfile)
			generate_weightings(badfileslist,event,scheme='XC')
			SecondStackCorrelate(event)



		os.chdir(datadir)

	badfileslist.close()




if __name__ == '__main__':

	main()