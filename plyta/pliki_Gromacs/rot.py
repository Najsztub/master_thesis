# -*- coding: utf-8 -*-
import sys
import os
from math import *
from numpy import array


home=os.getcwd()+"/"
#home='/home/mateusz/klustek/tit_s/new/'
#plik=open(home+sys.argv[1])
plik=open(home+'npt_rot.gro')
wect1=array([0.,0.,0.])
wect2=array([0.,0.,0.])
FILE=plik.readlines()
input=array([1,0,0])

#def sphe(co):
	#r = sqrt(co[0]**2 + co[1]**2 + co[2]**2)
	#phi = atan2(co[1], co[0])
	#if r != 0.0:
		#theta = acos( co[2] / r)
	#else:
		#theta = pi / 2 * sgn(co[2])
	#return array([r, theta*180/pi, phi*180/pi])

def sphe(co):
	r = sqrt(co[0]**2 + co[1]**2 + co[2]**2)
	if co[0] >= 0: 
		theta= asin(co[1]/sqrt(co[0]**2+co[1]**2))
	
	else: 
		theta= pi - asin(co[1]/sqrt(co[0]**2+co[1]**2))
	
	if r != 0.0:
		phi = acos( co[2] / r)
	else:
		phi = pi / 2 * sgn(co[2])
	return array([r, phi*180/pi, theta*180/pi])


def cart(co):
	x=co[0]*cos(co[2]*pi/180)*sin(co[1]*pi/180)
	y=co[0]*sin(co[2]*pi/180)*sin(co[1]*pi/180)
	z=co[0]*cos(co[1]*pi/180)
	return array([x,y,z])

def angZ(co):
	len=sqrt(pow(co[0],2)+pow(co[1],2))
	if len<=0: 
		ang="Err"
		return ang
	x_r=co[0]/len
	y_r=co[1]/len
	ang=asin(y_r)*180/pi
	if co[0]<0: ang= 180 - ang
	return ang

def angX(co):
	len=sqrt(pow(co[1],2)+pow(co[2],2))
	if len<=0: 
		ang="Err"
		return ang
	x_r=co[1]/len
	y_r=co[2]/len
	ang=asin(y_r)*180/pi
	if co[1]<0: ang= 180 - ang
	return ang

def angY(co):
	len=sqrt(pow(co[2],2)+pow(co[0],2))
	if len<=0: 
		ang="Err"
		return ang
	x_r=co[2]/len
	y_r=co[0]/len
	ang=asin(y_r)*180/pi
	if co[2]<0: ang= 180 - ang
	return ang

def angle(co):
	ang=array([angX(co), angY(co), angZ(co)])
	return ang

def convert(co, ile):
	temp=array([float(co[21:28]), float(co[29:36]), float(co[37:44])])
	#print temp
	wsp_sphe=sphe(temp)
	#print wsp_sphe
	theta=wsp_sphe[1]-ile[1]
	phi=wsp_sphe[2]-ile[2]
	if theta > 180:
		theta=360 - theta
		phi=phi+180
	wsp_sphe[1]=theta
	wsp_sphe[2]=phi
	#print wsp_sphe
	return cart(wsp_sphe)

### Wczyt. danych ###
index=0
while 1:
	try:
		temp=FILE[index].rsplit()
		if temp[0]=='338ILE' and temp[1]=='CA': wect2=array([float(temp[3]),float(temp[4]),float(temp[5])])
		if temp[0]=='18LYS' and temp[1]=='CA': wect1=array([float(temp[3]),float(temp[4]),float(temp[5])])
	except IndexError:
		break
	index=index+1

#####################
wynik=wect1-wect2
print wynik
print angle(wynik)

#obrot=sphe(input) - sphe(wynik)
#obrot =array([1.,9.,0.])

##'{0:8.3f}{1:8.3f}{2:8.3f}'.format(co[0], co[1], co[2])
#FILE_NEW=[]
#index=0
#lf=len(FILE)
#for index in range(lf):
	#if index==0 or index==1 or index==lf-1: 
		#FILE_NEW.append(FILE[index])
		#continue
	#try:
		#wsp=convert(FILE[index], obrot)
		#FILE_NEW.append(FILE[index][0:20]+'{0:8.3f}{1:8.3f}{2:8.3f}'.format(wsp[0], wsp[1], wsp[2])+'\n')
	#except IndexError:
		#continue

#out=open(home+'tit_tmp2.gro','w')
#for index in range(lf):
	#out.write(FILE_NEW[index])
