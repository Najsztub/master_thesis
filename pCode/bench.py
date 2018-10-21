# -*- coding: utf-8 -*-
#import matplotlib
#matplotlib.use('GTK')

from pylab import *
import numpy

sc=2

v_PME=[0.535 , 1.214 , 1.785, 2.332, 2.789, 3.310, 3.780, 4.258, 3.874, 4.356, 4.594 , 4.509 , 4.515, 5.495, 5.218, 5.530]
x=linspace(1, 16, 100)
subplot(221)
plot(range(1,17), v_PME, 'x',  markersize=14, mew=3)
plot( x, x*0.535, linewidth=2)
yscale('log')
xscale('log')
xticks([1, 2, 4, 8, 16], ('1', '2', '4', '8', '16'))
xlim(1, 17)
ylim(0.4, 16)
yticks( size=int(14*sc))
xticks( size=int(14*sc))
#xlabel('Ilo'+u'ś'+u'ć'+' rdzeni (z PME)', size=int(16*sc))
ylabel('Wydajno'+u'ś'+u'ć'+'\n[ns/dzie'+u'ń'+']', size=int(16*sc))
v_noPME=[0.535, 1.216 , 1.784 , 2.341, 2.793, 3.321, 3.786, 4.285, 3.900, 4.366, 4.020,  4.380 , 4.082, 4.611, 4.800, 5.087]
subplot(222)
plot(range(1,17), v_noPME, 'x',  markersize=14, mew=3)
plot(x, x*0.535, linewidth=2)
yscale('log')
xscale('log')
xticks([1, 2, 4, 8, 16], ('1', '2', '4', '8', '16'))
xlim(1, 17)
ylim(0.4, 16)
yticks( size=int(14*sc))
xticks( size=int(14*sc))
#xlabel('Ilo'+u'ś'+u'ć'+' rdzeni (bez PME)', size=int(16*sc))
ylabel('Wydajno'+u'ś'+u'ć'+'\n[ns/dzie'+u'ń'+']', size=int(16*sc))

subplot(223)
diff_PME=[]
for i in range(1,17):
	diff_PME.append((v_PME[i-1]-i*0.535)/(i*0.535))

plot(range(1,17), diff_PME, 'x',  markersize=14, mew=3, mec='#FF0000')
plot(x,x*0, linewidth=2)
xscale('log')
xticks([1, 2, 4, 8, 16], ('1', '2', '4', '8', '16'))
xlim(1, 17)
#ylim(0.4, 16)
yticks( size=int(14*sc))
xticks( size=int(14*sc))
xlabel('Ilo'+u'ś'+u'ć'+' rdzeni (z PME)', size=int(16*sc))
ylabel('Odchylenie \nod wydajno'+u'ś'+'ci \nteoretycznej', size=int(16*sc))

subplot(224)
diff_PME=[]
for i in range(1,17):
	diff_PME.append((v_noPME[i-1]-i*0.535)/(i*0.535))

plot(range(1,17), diff_PME, 'x',  markersize=14, mew=3, mec='#FF0000')
plot(x,x*0, linewidth=2)
xscale('log')
xticks([1, 2, 4, 8, 16], ('1', '2', '4', '8', '16'))
xlim(1, 17)
#ylim(0.4, 16)
yticks( size=int(14*sc))
xticks( size=int(14*sc))
xlabel('Ilo'+u'ś'+u'ć'+' rdzeni (bez PME)', size=int(16*sc))
ylabel('Odchylenie \nod wydajno'+u'ś'+'ci \nteoretycznej', size=int(16*sc))

show()
