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

			for sacfile in sacfiles:
				os.system('./Plot_Pick.py -infile %s' %sacfile)

		else:

			print 'provided path does not exist!'
			sys.exit(1)

	except:

		print 'No vaid path entered'
		print 'Usage: Pick_Stats.py -inpath [full path]'
		sys.exit(1)

if __name__ == '__main__':

	main()



