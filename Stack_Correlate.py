#!/usr/bin/env python 
#RMS Jan 2017

#Main stage of the Boyce et al (2017) absolute arrival times method
#load data, use the two stacking methods to produce a series of stacked SAC files
#The next stage involves picking on those files

import boyce2017functions as AT
import os
import glob



def main():


	datadir = '/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/test_alaska/2016-08-01_2016-09-25'
	targetfile = 'P_0.02.0.1_AIMBAT.out'


	SCHEME = 'RMS'
	SNR_cutoff = 3.5
	XC_cutoff = 0.9
	XC_time_cuttoff = 0.25
	Weighting_function = 7

	if not os.path.exists('event_stacks'):
	 	os.system('mkdir event_stacks')

	os.chdir('event_stacks')
	eventstackdir = os.getcwd()
	os.chdir(datadir)

	events = glob.glob('20*')

	badfileslist = open('bad_files.txt','wa')

	for event in events:

		print '\n=======================================\n'
		print 'IN EVENT %s' %event
		print '\n=======================================\n'

		os.chdir(event)
		os.chdir('BH_VEL')

		if os.path.isfile(targetfile):

			#prelim check to see if there is data in the file
			infile = open(targetfile,'r')
			lines = infile.readlines()
			infile.close()
			if len(lines) < 5:
				print 'No data in targetfile for event %s' %event

			else:

				#Append relative arrival times to the file headers, window the unfiltered data, stack and cross correlate all the
				#windowed traces with the stack

				AT.InitialStackCorrelate(badfileslist,event,targetfile)

				#QC stage - generate weightings for the second stack, either by using the SNR or the XC values

				AT.GenerateWeightings(badfileslist,event,scheme=SCHEME, function=Weighting_function, SNR_cutoff=SNR_cutoff, XC_cutoff=XC_time_cuttoff, XC_time_cuttoff=XC_time_cuttoff)

				#Do the second stack

				AT.SecondStackCorrelate(event,scheme=SCHEME)

				#Write a list of all the files whose XC/SNR values are sufficiently high that we can be confident in finding their absolute delay times

				AT.CheckErrors(badfileslist,event, XC_cutoff=XC_cutoff, XC_time_cuttoff=XC_time_cuttoff)

				#Clean up (testing)
				os.system('rm *stk*')

				#Put the stacked files in another directory, ready for picking
				os.system('mv *STACK* %s' %eventstackdir)



		os.chdir(datadir)

	badfileslist.close()


if __name__ == '__main__':

	main()