# -*- coding: utf-8 -*-
from pylab import *
import numpy
import matplotlib


def conv(inpu):
	resu=[0, 0]
	temp = inpu.rsplit()
	try:
		resu=[float(temp[0]), float(temp[1])]
	except ValueError:
		resu=[0, 0]
	except IndexError:
		resu=[0, 0]
	return resu


def read_data(nazwa):
	zwrot=[]
	plik=open(nazwa)
	while 1:
		temp=plik.readline()
		if not temp: break
		temp1=conv(temp)
		if temp1 != [0, 0]: zwrot.append(temp1)		
	plik.close()
	return zwrot

force=[]
force = read_data("/home/mateusz/klustek/1TKI/npt_sm.xvg")
E_pot, step = [], []
i=0
while 1:
		try:
			E_pot.append(force[i][1])
			step.append(force[i][0]/1000)
			i=i+1	
		except IndexError:
			break

fig_width_pt = 300.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

golden_mean = (sqrt(5)-1.0)/2.0
#rcParams.update({'backend': 'ps',
   #'ps.usedistiller': 'xpdf',
   #'axes.labelsize': 10,
   #'text.fontsize': 10,
   #'xtick.labelsize': 10,
   #'ytick.labelsize': 10,
   #'figure.figsize': fig_size,
   #'text.usetex': True })


clf()
plot(step, E_pot, linewidth=1)
yticks( size=30)
xticks( size=30)
xlabel('Czas [ns]', size=40)
ylim(0.05, 0.3)
ylabel('RMSD [nm]',size=40)
#savefig('fog1.pdf')
show()


