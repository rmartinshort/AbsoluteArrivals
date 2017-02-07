#!/usr/bin/env python 
#RMS Jan 2017

#Last stage of the Boyce et al (2017) absolute arrival times method

#Take the arrival corrections from stacks in the event_stacks directory and paste them into the headers of the original sacfiles
#Also write a file of 'absolute arrival times' for each event. This can then be used in future tomogrraphy etc

import boyce2017functions as AT


def main():

	datadir = '/Users/rmartinshort/Documents/Berkeley/Alaska/Tomography/Joint_surface/python_abstimes/test_from_lower48/2006-07-01_2006-09-01'

	AT.AppendAbsoluteArrivals(datadir)


if __name__ == '__main__':

	main()