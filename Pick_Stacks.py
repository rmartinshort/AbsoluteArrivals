#!/usr/bin/env python

#Runs plot/pick on each of the absolute traveltime stacks, allowing the user to choose an absolute arrival time

import glob
import argparse
import os
import sys

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('-inpath',action='store',dest='inpath',help='The full path to the folder where the stacked sac files are stored')
	results = parser.parse_args()

	inpath = results.inpath

	try:

		if os.path.exists(inpath):

			sacfiles = glob.glob(inpath+'/*.sac')

			evts = []
			for sacfile in sacfiles:
				evt = sacfile.split('.')[0]
				if evt not in evts:
					evts.append(evt)

			for event in evts:

				print '------------------------------------\n'
				print 'Event %s' %event 
				print '\n------------------------------------'

				sacfiles = glob.glob(event+'*.sac')
				if len(sacfiles) == 2:
					os.system('./Plot_Pick.py -infile1 %s -infile2 %s' %(sacfiles[0],sacfiles[1]))
				else:
					print 'Number of sacfiles in %s != 2!' %(sacfiles)
					os.system('./Plot_Pick.py -infile1 %s' %(sacfiles[0]))

		else:

			print 'provided path does not exist!'
			sys.exit(1)

	except:

		print 'No vaid path entered'
		print 'Usage: Pick_Stats.py -inpath [full path]'
		sys.exit(1)

if __name__ == '__main__':

	main()



