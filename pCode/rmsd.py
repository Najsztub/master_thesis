# -*- coding: utf-8 -*-
#import matplotlib
#matplotlib.use('GTK')

from pylab import *
import numpy


def conv(inpu, n1):
	resu=[0, 0]
	temp = inpu.rsplit()
	try:
		resu=[float(temp[0]), float(temp[n1])]
	except ValueError:
		resu=[0, 0]
	except IndexError:
		resu=[0, 0]
	return resu


def read_data(nazwa,n1):
	zwrot=[]
	plik=open(nazwa)
	while 1:
		temp=plik.readline()
		if not temp: break
		temp1=conv(temp,n1)
		if temp1 != [0, 0]: zwrot.append(temp1)		
	plik.close()
	return zwrot

def smooth(list,degree=5):  
	window=degree*2-1  
	weight=numpy.array([1.0]*window)  
	weightGauss=[]  
	for i in range(window):
		i=i-degree+1  
		frac=i/float(window)  
		gauss=1/(numpy.exp((4*(frac))**2))  
		weightGauss.append(gauss)  
	weight=numpy.array(weightGauss)*weight  
	smoothed=[0.0]*(len(list)-window)  
	for i in range(len(smoothed)):  
		smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
	return smoothed  

def max_s(input):
	x, y = 0, 0
	for k in range(len(input)):
		if input[k] > x:
			x=input[k]
			y=k
	ret=[x,y]
	return ret

#dist = read_data("/home/mateusz/klustek/1TKI/p0.4/pullx.xvg", 2, 3)
rmsd = read_data("/home/mateusz/klustek/1TKI/npt_rmsd.xvg", 1)
pres = read_data("/home/mateusz/klustek/1TKI/npt_p.xvg", 1)
vol = read_data("/home/mateusz/klustek/1TKI/npt_v.xvg", 1)
tempe = read_data("/home/mateusz/klustek/1TKI/npt_t.xvg", 1)

dr, dp, dv, dt =[], [], [], []
r1, p1, v1, t1 = [], [], [], []
v=0.8
i=0
while 1:
		try:
			r1.append(rmsd[i][1])
			dr.append(rmsd[i][0]/1000)
			i=i+1	
		except IndexError:
			break
			
i=0
while 1:
		try:
			p1.append(pres[i][1])
			dp.append(pres[i][0]/1000)
			i=i+1	
		except IndexError:
			break
			
i=0
while 1:
		try:
			v1.append(vol[i][1])
			dv.append(vol[i][0]/1000)
			i=i+1	
		except IndexError:
			break
			
i=0
while 1:
		try:
			t1.append(tempe[i][1])
			dt.append(tempe[i][0]/1000)
			i=i+1	
		except IndexError:
			break
			

r1_s=smooth(r1, 20)
p1_s=smooth(p1, 20)
v1_s=smooth(v1, 20)
t1_s=smooth(t1, 20)

subplot(221)
plot(dr,r1, color='#269823', aa=1)
#plot(dr[(len(r1)-len(r1_s))/2:len(r1_s)+(len(r1)-len(r1_s))/2], r1_s, linewidth=3, aa=1, color='#269823')
yticks( size=30)
xticks( size=30)
xlabel('Czas [ns]', size=40)
ylabel('RMSD [nm]',size=40)

subplot(222)
plot(dp,p1, color='#C9CBF0', aa=1)
plot(dp[(len(p1)-len(p1_s))/2:len(p1_s)+(len(p1)-len(p1_s))/2], p1_s, linewidth=3, aa=1, color='#2A32D4')
yticks( size=30)
xticks( size=30)
xlabel('Czas [ns]', size=40)
ylabel('Ci'+u'ś'+'nienie [bar]',size=40)

subplot(223)
plot(dv,v1, color='#EECB74', aa=1)
plot(dv[(len(v1)-len(v1_s))/2:len(v1_s)+(len(v1)-len(v1_s))/2], v1_s, linewidth=3, aa=1, color='#8B6914')
yticks( size=30)
xticks( size=30)
xlabel('Czas [ns]', size=40)
ylabel('Obj'+u'ę'+'to'+u'ść'+' [mn^3]',size=40)

subplot(224)
plot(dt,t1, color='#FFAEC1', aa=1)
plot(dt[(len(t1)-len(t1_s))/2:len(t1_s)+(len(t1)-len(t1_s))/2], t1_s, linewidth=3, aa=1, color='#F21B4E')
yticks( size=30)
xticks( size=30)
xlabel('Czas [ns]', size=40)
ylabel('Temperatura [K]',size=40)


show()			
#plot(d1[(len(f1)-len(f1_s)-1)/2:len(f1_s)+(len(f1)-len(f1_s)-1)/2], f1_s, d1, f1)
#plot(d1[(len(f1)-len(f1_s))/2:len(f1_s)+(len(f1)-len(f1_s))/2], f1_s, d1[(len(f2)-len(f2_s))/2:len(f2_s)+(len(f2)-len(f2_s))/2], f2_s)
#xscale('log')
			
