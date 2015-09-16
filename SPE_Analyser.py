#!/usr/bin/python2
#-------------------------------------------------------------------------------
# Name:		GUI Test
# Purpose:
#
# Author:	  scholi
#
# Created:	 06.05.2012
# Copyright:   (c) scholi 2012
# Licence:	 <your licence>
#-------------------------------------------------------------------------------

import sys
import os
from PyQt4 import QtCore, QtGui
from Spectro import Ui_MainWindow

from SPE import SPEopen

import numpy as np

import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.mlab as mlab

class SpectroGUI(QtGui.QMainWindow):
	def __init__(self, filename=None, arg=None, parent=None):
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		self.dpi = 100
		self.fig = Figure((5.0, 4.0), dpi=self.dpi)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.ui.widget)

		# Since we have only one plot, we can use add_axes
		# instead of add_subplot, but then the subplot
		# configuration tool in the navigation toolbar wouldn't
		# work.
		#
		self.axes = self.fig.add_subplot(411)
		self.raw = self.fig.add_subplot(412)
		self.plt = self.fig.add_subplot(413)
		self.ref = self.fig.add_subplot(414)
#		self.vari = self.fig.add_subplot(515)


		# Bind the 'pick' event for clicking on one of the bars
		#
		self.canvas.mpl_connect('pick_event', self.on_pick)

		# Create the navigation toolbar, tied to the canvas
		#
		self.mpl_toolbar = NavigationToolbar(self.canvas, self.ui.widget)

		vbox = QtGui.QVBoxLayout()
		vbox.addWidget(self.canvas)
		vbox.addWidget(self.mpl_toolbar)
		self.ui.widget.setLayout(vbox)

		# Signals/Slots
		self.connect(self.ui.actionOpen, QtCore.SIGNAL('triggered()'), self.loadSPE)
		self.connect(self.ui.sigTop, QtCore.SIGNAL('valueChanged(int)'), self.on_draw)
		self.connect(self.ui.sigBottom, QtCore.SIGNAL('valueChanged(int)'), self.on_draw)
		self.connect(self.ui.ref1Top, QtCore.SIGNAL('valueChanged(int)'), self.on_draw)
		self.connect(self.ui.ref1Bottom, QtCore.SIGNAL('valueChanged(int)'), self.on_draw)
		self.connect(self.ui.rawsigStdDev, QtCore.SIGNAL('stateChanged(int)'), self.on_draw)
		self.connect(self.ui.sigStdDev, QtCore.SIGNAL('stateChanged(int)'), self.on_draw)
		self.connect(self.ui.refStdDev, QtCore.SIGNAL('stateChanged(int)'), self.on_draw)
		self.connect(self.ui.NumSignals, QtCore.SIGNAL('valueChanged(int)'), self.on_draw)
		self.connect(self.ui.SigSpacing, QtCore.SIGNAL('valueChanged(double)'), self.on_draw)
		self.connect(self.ui.actionClose, QtCore.SIGNAL('triggered()'), self.closeEvent)
		self.connect(self, QtCore.SIGNAL('triggered()'), self.closeEvent)
		self.connect(self.ui.multisig, QtCore.SIGNAL("toggled(bool)"), self.on_draw)

		if filename!=None:
			self.loadSPE(filename)

	def loadSPE(self, filename=None):
		if filename==None:
			fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file')
		else:
			fname=filename
		if not os.path.exists(fname):
			 QtGui.QMessageBox.warning(self, "Warning", "The file \"%s\" cannot be accessed!"%(fname))
		self.data=SPEopen(fname)
		self.axes.imshow(self.data["rawdata"])

		self.lup = [mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="red")]
		self.ldown = [mpl.lines.Line2D([0,self.data["var"]['nx']],[self.data["var"]['ny'],self.data["var"]['ny']],color="red")]
		self.axes.add_line(self.lup[0])
		self.axes.add_line(self.ldown[0])

		self.lref1up = mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="green")
		self.lref1down = mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="green")
		self.lref2up = mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="green")
		self.lref2down = mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="green")
		self.axes.add_line(self.lref1up)
		self.axes.add_line(self.lref1down)
		self.axes.add_line(self.lref2up)
		self.axes.add_line(self.lref2down)

		self.ui.sigBottom.setRange(1,self.data["var"]["ny"])
		self.ui.sigTop.setRange(1,self.data["var"]["ny"])
		self.ui.ref1Bottom.setRange(1,self.data["var"]["ny"])
		self.ui.ref1Top.setRange(1,self.data["var"]["ny"])
		self.ui.ref2Top.setRange(1,self.data["var"]["ny"])
		self.ui.ref2Bottom.setRange(1,self.data["var"]["ny"])
		trd={
			"SignalUp":self.ui.sigTop,
			"SignalDown":self.ui.sigBottom,
			"Ref1Up":self.ui.ref1Top,
			"Ref1Down":self.ui.ref1Bottom,
			"Ref2Up":self.ui.ref2Top,
			"Ref2Down":self.ui.ref2Bottom,
			"NumSignals":self.ui.NumSignals}
			
		if os.path.exists(fname[:-4]+".txt"):
			ff=open(fname[:-4]+".txt","r")
			inf={}
			for x in ff.readlines():
				c=x.split()
				inf[c[0]]=c[1]
			for x in trd:
				if x in inf:
					trd[x].setValue(int(inf[x]))
			if "SignalSpacing" in inf: self.ui.SigSpacing.setValue(float(inf["SignalSpacing"]))
		elif fname[-6:-4]=="90" and os.path.exists(fname[:-6]+".txt"):
			print("Bingo!")
			ff=open(fname[:-6]+".txt","r")
			inf={}
			for x in ff.readlines():
				c=x.split()
				inf[c[0]]=c[1]
			for x in trd:
				if x in inf:
					trd[x].setValue(int(inf[x]))
			if "SignalSpacing" in inf: self.ui.SigSpacing.setValue(float(inf["SignalSpacing"]))

	def on_pick(self):
		pass
	def closeEvent(self, event=None):
		print("Closing...")
		yup = self.ui.sigTop.value()
		ydown = self.ui.sigBottom.value()
		yrefup= self.ui.ref1Top.value()
		yrefdown = self.ui.ref1Bottom.value()
		yref2up = self.ui.ref2Top.value()
		yref2down = self.ui.ref2Bottom.value()
		n=self.ui.NumSignals.value()
		if n==0: n=1
		dh=self.ui.SigSpacing.value()
		fname=self.data["var"]["filename"]
		# Write Infos
		ff=open(fname[:-4]+".txt","w")	
		ff.write("SignalUp\t%i\nSignalDown\t%i\nRef1Up\t%i\nRef1Down\t%i\nRef2Up\t%i\nRef2Down\t%i\nNumSignals\t%i\nSignalSpacing\t%f"%(yup,ydown,yrefup,yrefdown,yref2up,yref2down,n,dh))
		ff.close()
		# Write Data
		# prepare Data
		#ff=open(fname[:-4]+".dat","w")
		ffr=open(fname[:-4]+"_raw.dat","w")	
		print("Saving...")
		if yrefup!=yrefdown:
			ref=np.mean(self.data["rawdata"][yrefup:yrefdown],0)
		w=self.data["wavelength"]
		#data=[np.mean(self.data["rawdata"][yup-i*dh:ydown-i*dh],0) for i in range(n)]
		rawdata=[np.mean(self.data["rawdata"][yup-i*dh:ydown-i*dh],0) for i in range(n)]
		#ff.write("# wavelength[nm]")
		ffr.write("# wavelength[nm]")
		for i in range(n):
			#ff.write("\tIntensity_antenna_%i[-]"%(i))
			ffr.write("\tIntensity_antenna_%i[-]"%(i))
		if yrefup!=yrefdown: ffr.write("\tReference[-]")
		ffr.write("\n")
		#if yrefup!=yrefdown:
			#for i in range(n):
				#data[i]=data[i]-ref
		for x in range(len(rawdata[0])):
			#ff.write("%f"%(w[x]))
			ffr.write("%f"%(w[x]))
			for i in range(n):
				#ff.write("\t%f"%(data[i][x]))
				ffr.write("\t%f"%(rawdata[i][x]))
			if yrefup!=yrefdown: ffr.write("\t%f"%(ref[x]))
			#ff.write("\n")
			ffr.write("\n")
		#ff.close()
		ffr.close()
		self.destroy()

	def on_draw(self):
		yup = self.ui.sigTop.value()
		ydown = self.ui.sigBottom.value()
		yrefup= self.ui.ref1Top.value()
		yrefdown = self.ui.ref1Bottom.value()
		yref2up = self.ui.ref2Top.value()
		yref2down = self.ui.ref2Bottom.value()
		n=self.ui.NumSignals.value()
		if n==0: n=1
		dh=self.ui.SigSpacing.value()

		self.lup[0].set_data([0,self.data["var"]['nx']],[yup,yup])
		self.ldown[0].set_data([0,self.data["var"]['nx']],[ydown,ydown])
		self.lref1up.set_data([0,self.data["var"]['nx']],[yrefup,yrefup])
		self.lref1down.set_data([0,self.data["var"]['nx']],[yrefdown,yrefdown])
		self.lref2up.set_data([0,self.data["var"]['nx']],[yref2up,yref2up])
		self.lref2down.set_data([0,self.data["var"]['nx']],[yref2down,yref2down])
			

		self.plt.clear()
		self.ref.clear()
		self.raw.clear()
		self.ref.clear()

	# Signals
		data=[np.mean(self.data["rawdata"][yup:ydown],0)]
		self.raw.plot(self.data["wavelength"],data[0],color="black")
		if self.ui.multisig.checkState()==QtCore.Qt.Checked:
			if len(self.lup)>n:
				for i in range(n+1,len(self.lup)):
					self.axes.remove(self.lup[i])
					self.axes.remove(self.ldown[i])
				self.lup=self.lup[:n]
				self.ldown=self.ldown[:n]
			if n>1:
				if len(self.lup)<n:
					for i in range(len(self.lup),n):
						self.lup.append(mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="red"))
						self.ldown.append(mpl.lines.Line2D([0,self.data["var"]['nx']],[0,0],color="red"))
						self.axes.add_line(self.lup[i])
						self.axes.add_line(self.ldown[i])
				for i in range(1,n):
					data.append(np.mean(self.data["rawdata"][int(yup-dh*i):int(ydown-dh*i)],0))
					self.raw.plot(self.data["wavelength"],data[i])
					self.lup[i].set_data([0,self.data["var"]['nx']],[int(yup-dh*i),int(yup-dh*i)])
					self.ldown[i].set_data([0,self.data["var"]['nx']],[int(ydown-dh*i),int(ydown-dh*i)])
	
		# Sig StdDev
		drs=np.sqrt(np.var(self.data["rawdata"][yup:ydown],0))
		if self.ui.rawsigStdDev.checkState()==QtCore.Qt.Checked:
				vx,vy = mlab.poly_between(self.data["wavelength"],data[0]-drs,data[0]+drs)
				self.raw.fill(vx,vy,color="red")

		# Reference Signal
		if yrefup!=yrefdown:
			ref=np.mean(self.data["rawdata"][yrefup:yrefdown],0)
			self.ref.plot(self.data["wavelength"],ref,color="black")
			# Std Dev
			delta=np.sqrt(np.var(self.data["rawdata"][yrefup:yrefdown],0))
			if self.ui.refStdDev.checkState()==QtCore.Qt.Checked:
				vx,vy = mlab.poly_between(self.data["wavelength"],ref-delta,ref+delta)
				self.ref.fill(vx,vy,color="blue")
				DS=np.sqrt((drs/ref)**2+((data[0]*delta)/(ref**2))**2)
				data[0]=data[0]/ref
				self.plt.plot(self.data["wavelength"],data[0],color="black")
			if self.ui.multisig.checkState()==QtCore.Qt.Checked:
				if n>1:
					for i in range(1,n):
						data[i]=data[i]-ref
						self.plt.plot(self.data["wavelength"],data[i])
			if self.ui.sigStdDev.checkState()==QtCore.Qt.Checked:
					vx,vy = mlab.poly_between(self.data["wavelength"],data[0]-DS,data[0]+DS)
					self.plt.fill(vx,vy)
		self.canvas.draw()

if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	arg=None
	if len(sys.argv)>2:
		arg=sys.argv[2]
	if len(sys.argv)>1:
		myapp = SpectroGUI(sys.argv[1],arg)
	else:	
		myapp = SpectroGUI(None,arg)
	myapp.show()
	sys.exit(app.exec_())
