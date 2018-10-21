# -*- coding: utf-8 -*-
#import matplotlib
#matplotlib.use('GTK')

from pylab import *
import numpy


def conv(inpu, n1, n2):
	resu=[0, 0, 0]
	temp = inpu.rsplit()
	try:
		resu=[float(temp[0]), float(temp[n1]), float(temp[n2])]
	except ValueError:
		resu=[0, 0, 0]
	except IndexError:
		resu=[0, 0, 0]
	return resu


def read_data(nazwa,n1,n2):
	zwrot=[]
	plik=open(nazwa)
	while 1:
		temp=plik.readline()
		if not temp: break
		temp1=conv(temp,n1,n2)
		if temp1 != [0, 0, 0]: zwrot.append(temp1)		
	plik.close()
	return zwrot

def smooth(list,degree):  
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
force = read_data("/home/mateusz/klustek/1TKI/2p0.8/pullf.xvg", 1, 2)

d1 =[]
d2 =[]
f1, f2 = [], []
v=0.8
i=0
#d1_0, d2_0 = dist[0][1],dist[0][2]
while 1:
		try:
			#d1.append(d1_0-dist[i][1])
			#d2.append(dist[i][2]-d2_0)
			#f1.append(force[i*10][1])
			#f2.append(force[i*10][2])
			f1.append(force[i*100][1]/0.602)
			f2.append(force[i*100][2]/0.602)
			d1.append(force[i*100][0]*v*1e-3)
			i=i+1	
		except IndexError:
			break
f1_s=smooth(f1, int(0.16/v/10e-3))
f2_s=smooth(f2, int(0.16/v/10e-3))
max1 = max_s(f1_s)
max2 = max_s(f2_s)
#max2=[742, 36924]
print max1
print max2
plot(d1,f1, color='#CAF0C9', aa=1)
plot(d1,f2, color='#C9E6F0')
plot(d1[(len(f1)-len(f1_s))/2:len(f1_s)+(len(f1)-len(f1_s))/2], f1_s, linewidth=3, aa=1, color='#269823')
plot(d1[(len(f2)-len(f2_s))/2:len(f2_s)+(len(f2)-len(f2_s))/2], f2_s, linewidth=3, aa=1, color='#227797')
yticks( size=30)
xticks( size=30)
xlabel('Rozci'+u'ą'+'gni'+u'ę'+'cie [nm]', size=40)

ylabel('Si'+u'ł'+'a [pN]',size=40)
annotate(str(int(max1[0]))+' pN', xy=( d1[max1[1]+int(0.16/v/10e-3)], max1[0]),  xycoords='data', xytext=(-50, 30), textcoords='offset points', arrowprops=dict(arrowstyle="->"),  size=25)
annotate(str(int(max2[0]))+' pN', xy=( d1[max2[1]+int(0.16/v/10e-3)], max2[0]),  xycoords='data', xytext=(-50, 30), textcoords='offset points', arrowprops=dict(arrowstyle="->"),  size=25)
show()			
#plot(d1[(len(f1)-len(f1_s)-1)/2:len(f1_s)+(len(f1)-len(f1_s)-1)/2], f1_s, d1, f1)
#plot(d1[(len(f1)-len(f1_s))/2:len(f1_s)+(len(f1)-len(f1_s))/2], f1_s, d1[(len(f2)-len(f2_s))/2:len(f2_s)+(len(f2)-len(f2_s))/2], f2_s)
#xscale('log')
			
