#!/usr/bin/env python

#Use PyQTgraph to plot a seismogram stack, pick the absolute arrival and write to the SAC header

from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
import sys
import os
import obspy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-infile1',action='store',dest='infile1',help='The full file path to the data you want to plot/pick. Should be a stacked SACfile')
parser.add_argument('-infile2',action='store',dest='infile2',help='The full file path to the data you want to plot/pick. Should be a stacked SACfile')
results = parser.parse_args()


if os.path.isfile(results.infile1):

	insac1 = results.infile1

	trace1 = obspy.read(insac1,format='SAC')

	if results.infile2:
		insac2 = results.infile2
		tracesecond = obspy.read(insac2,format='SAC')
	else:
		print 'No second trace found/provided'
		tracesecond = None


	#Set up the QTGui ap
	app = QtGui.QApplication([])
	win = pg.GraphicsWindow(title="SeismoViewer plot")
	pg.setConfigOptions(antialias=True)

	#Set up cursor location label
	label = pg.LabelItem(justify='right')
	win.addItem(label)

	#Set up x axis
	npts = trace1[0].stats.npts
	x = np.linspace(-10,10,npts)
	DATA = trace1[0].data

	#Plot the trace 
	p1 = win.addPlot(title='%s' %insac1)
	p1.plot(x,DATA,pen=(255,0,0))
	p1.setLabel('bottom', "Time", units='s')
	p1.showGrid(x=True,y=True,alpha=0.5)

	#plot the second trace if its present
	if tracesecond:
		DATA2 = tracesecond[0].data
		p1.plot(x,DATA2,pen=(0,255,0))

	#plot the position of the relative arrival (always ay zero)
	vZeroline = pg.InfiniteLine(angle=90,movable=False,pen=(255,255,255))
	p1.addItem(vZeroline, ignoreBounds=True)
	vZeroline.setPos(0.0)

	#pick scroll pointer position
	vLinescroll = pg.InfiniteLine(angle=90, movable=False)
	p1.addItem(vLinescroll, ignoreBounds=True)

	#pick indicator position
	vLinepick = pg.InfiniteLine(angle=90, movable=False)
	p1.addItem(vLinepick, ignoreBounds=True)

	vb = p1.vb

	def mouseMoved(evt):

		''' Move the posotion of the pick scroller around on the plot'''

		pos = evt[0]  ## using signal proxy turns original arguments into a tuple
		if p1.sceneBoundingRect().contains(pos):
			mousePoint = vb.mapSceneToView(pos)
			index = int(mousePoint.x())

			if index >= min(x) and index <= max(x):
				label.setText("<span style='font-size: 16pt'>x = %0.4f s" % (mousePoint.x()))

			vLinescroll.setPos(mousePoint.x())

	def mouseClicked(evt):

		'''On mouse click, select a pick and change the sac header accordingly'''

		pos = evt.scenePos() #location of the moise click event within the scene 
		if p1.sceneBoundingRect().contains(pos):
			mousePoint = vb.mapSceneToView(pos)

			#The location of the pick
			vLinepick.setPos(mousePoint.x())

            #Set the sac header and write it to file
			trace1[0].stats.sac.t0 = mousePoint.x()
			trace1.write(insac1,format='SAC')

			#Write the same time to trace2 if present

			if tracesecond:
				tracesecond[0].stats.sac.t0 = mousePoint.x()
				tracesecond.write(insac2,format='SAC')


	#Add the mouse moved and mouse clicked functions to the plot object 
	proxy = pg.SignalProxy(p1.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)
	p1.scene().sigMouseClicked.connect(mouseClicked)

else:

	print 'Provided file %s not found' %results.infile
	sys.exit(1)






## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

