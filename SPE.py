#!/usr/bin/python2
# -*- coding: utf-8 -*-

import struct
#import sys
import numpy as np

def SPEopen(path):
	# Read the SPE file and load it's content into RAM
	f=open(path,"rb")
	var={"filename":path}
	t={"exp":('f',10,4),"nx":('H',42,2),"ny":('H',656,2),"nframes":('I',1446,4),"datatype":('H',108,2),"date":('10s',20,10),"time":('6s',172,6),"usercomment1":('80s',200,80),"usercomment2":('80s',280,80),"usercomment3":('80s',360,80),"usercomment4":('80s',440,80),"usercomment5":('80s',520,80),'cw':('f',72,4)}
	for x in t:
		f.seek(t[x][1])
		var[x]=struct.unpack(t[x][0],f.read(t[x][2]))[0]
#	print(x,var[x])
	length=var['nx']*var['ny']*var['nframes']
	f.seek(4100)
	if var['datatype']==0:
		rawdata=np.array(struct.unpack("%if"%(length),f.read(length*4)))
	elif var['datatype']==1:
		rawdata=np.array(struct.unpack("%ii"%(length),f.read(length*4)))
	elif var['datatype']==2:
		rawdata=np.array(struct.unpack("%iH"%(length),f.read(length*4)))
	elif var['datatype']==3:
		rawdata=np.array(struct.unpack("%ih"%(length),f.read(length*4)))
	rawdata=np.reshape(rawdata,(var['ny'],var['nx']))
	x0 = var['nx']/2+0.5
	x=np.arange(var['nx'])
	if var['cw']==600.0:
		a2=-1.54334e-5
		a1=0.382039
		wavelength=var['cw']+a1*(x-x0)+a2*(x-x0)**2
	elif var['cw']==700.0:
		a2=-1.58478e-5
		a1=0.378275
		wavelength=var['cw']+a1*(x-x0)+a2*(x-x0)**2
	elif var['cw']==800.0:
		a2=-1.62586e-5
		a1=0.374399
		wavelength=var['cw']+a1*(x-x0)+a2*(x-x0)**2
	elif var['cw']==900.0:
		a2=-1.66658e-5
		a1=0.37041
		wavelength=var['cw']+a1*(x-x0)+a2*(x-x0)**2
	else:
		print("pixel -> wavelength function unknown! Saving pixel value for x")
		wavelength=x
	
	return {"var":var,"wavelength":wavelength,"rawdata":rawdata}
