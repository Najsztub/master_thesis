# -*- coding: utf-8 -*-
from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.integrate import quad
from numpy import exp
from scipy.special import lambertw

#from mpmath import lambertw
#from mpmath import mp
#mp.dps = 25; mp.pretty = True

kb=1.3807e-23	# J/K
ks=5.031e20	# N/m /(kb*T)c -> 1/m^2
T=300.		#K


###		Prawdopodobieństwo, że cząsteczka nie uległa rozpleceniu w czasie t; Wzór nr 6 w Mat Sup- Schulten.pdf

s = lambda t, vel_l, p_l: exp(-p_l[0] * exp(-0.5*ks*p_l[1]**2)/(vel_l*ks*p_l[1]*(p_l[2]/(p_l[2]+ks))**(1.5))*(exp(ks*vel_l*p_l[1]*t - (0.5*(ks*vel_l*t)**2)/(ks+p_l[2]))-1))

###		Obliczanie wartości teoretycznych sił rozplatania w zależności od parametrów i szybkości

def F(p_l, x_l):
	k0_l, xb_l, km_l = p_l[0], p_l[1], p_l[2]
	D_l=k0_l * sqrt(2*pi) / (km_l**(1.5)*xb_l*exp(-0.5*km_l*xb_l**2)) #Rozwiązanie równiania nr 3 w Mat Sup- Schulten.pdf
	b=D_l*(km_l+ks)
	out=[]
	for i in range(len(x_l)):
		vel_l = x_l [i]
		a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
		tau=(lambertw(-exp(-a))+a)/b #Rozwiązanie równania nr 4 w Supp Mat - Szulten
		out.append(-ks*T*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p_l), epsabs=1e-16)[0])) #Wzór nr 5 w Mat Sup- Schulten.pdf
	return array(out)

###		Funkcja błędu; suma kwadratu różnic między siłą teoretyczną dla prędkości V i siłą f z symulacji
def err(p_l):
	out=0
	fs=F(p_l,V)
	out=((fs - f)**2).sum()
	#print 'k0, xb, km ', p_l, out
	return out

###		Obliczanie R^2

def R_sq(par, x_l, f_sim):
	f=F(par,x_l)
	f_sim_m=0
	l=len(x_l)
	for i in range(l):
		f_sim_m+=f_sim[i]
	f_sim_m=f_sim_m/l
	s_l, s_m =0. , 0.
	for j in range(l):
		s_l+=(f_sim[j]-f[j])**2
		s_m+=(f_sim[j]-f_sim_m)**2
		#print s_l, f_sim[j]-f[j], s_m
	return 1-s_l/s_m


f=array([1581.2e-12, 1209.1e-12, 978.1e-12, 863.2e-12, 795.4e-12, 592.4e-12, 491e-12]) # siły rozplatania N
V=array([28., 8., 2.8, 0.8, 0.28, 0.028, 0.0028]) # v rozplatania m/s

p=[1.6e-1, 2e-10, 9e20] # Parametry startowe [k0, xb, km]
p2=[1.6e-11, 3.85e-10, 6.905e20] # Parametry Schultena
p3=[3.3e-4, 2.50e-10, 1.096e21] # Parametry eksperymantalne
p_full=[  4.24307789e-04 ,  2.36799911e-10  , 1.21931408e+21]

p1=fmin(err, p, ftol=1e-1, maxiter=1000000, maxfun=1000000) #Minimalizacja 
#p1=[  7.66024933e-01,   1.71966849e-10 ,  1.81616113e+21]
#p1=[  2.36061090e-01  , 1.82607521e-10 ,  1.67820881e+21]

f_exp=array([125e-12,  166e-12, 230e-12, 283e-12 ])
V_exp=array([0.01e-6, 0.1e-6, 1e-6, 10e-6])

#print 'Po optymalizacji: ', b1
#print 'Schultena: ', b2

subplot(121)
x_log=logspace(-8,2, 200)
#plot(x_log, F(p_full, x_log)/1e-12)
plot(x_log, F(p1, x_log)/1e-12)
plot(x_log, F(p2, x_log)/1e-12)
#plot(x_log, F(p3, x_log)/1e-12)

plot(V_exp, f_exp/1e-12, 'xr')
plot(V, f/1e-12, 'og')

xscale('log')
xlabel('v [m/s]')
ylabel('F [pN]')
legend( ('Optymalizacja', 'Par. Schultena', 'Z AFM', 'Z symulacji'), loc=2 )


####################################################### kinaza 1TKI


ks=2.00e20		# N/m /(kb*T)c -> 1/m^2
subplot(122)
f1=array([1323e-12, 796e-12, 645e-12, 366e-12, 495e-12, 452e-12]) #f dla rozplatania N-końca Grubmiller.pdf
#f1=array([1376e-12, 710e-12, 526e-12, 559e-12, 667e-12, 538e-12]) #f dla rozplatania C-końca Grubmiller.pdf

#f=array([1382e-12, 861e-12, 746e-12, 609e-12, 366e-12]) #f dla rozplatania N-końca mojej sym.
f=array([1618e-12, 822e-12, 745e-12, 650e-12, 459e-12, 429e-12]) # f dla rozplatania C-końca mojej sym.
V=array([50, 10, 5, 2, 0.8, 0.4]) # V dla rozplatania - takie same

p=[1.6e2, 4e-10, 6e20] # Parametry startowe [k0, xb, km]

p3=fmin(err, p, xtol=1e-5, maxiter=1000000, maxfun=1000000) #Minimalizacja 
#p3=[2.78690714e+03  , 1.96046057e-10  , 8.85721077e+20]
#p2=[  3.96721306e-01  , 3.47040930e-10 ,  4.35590375e+20]


x_log=logspace(-3,2, 300)
plot(x_log, F(p3, x_log)/1e-12)
plot(V, f/1e-12, 'o')

f=f1
p4=fmin(err, p, ftol=1e-1, maxiter=1000000, maxfun=1000000)

print 'Parametry po optymalizacji', p1
print 'Parametry z pracy Schultena', p2
print 'Optymalizacjia z symulacji', p3
print 'Optymalizacjia z pracy', p4

plot(x_log, F(p4, x_log)/1e-12)
plot(V, f/1e-12, 'x')
xscale('log')
xlabel('v [m/s]')
ylabel('F [pN]')
legend( ('Optymalizacja, symulacja', 'Wyniki symulacji', 'Optymalizacja publikacja', 'Wyniki publikacji'), loc=2)
show()



