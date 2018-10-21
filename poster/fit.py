# -*- coding: utf-8 -*-
from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.integrate import quad
from numpy import exp

kb=1.3807e-23	# J/K
k0=1.6e-11		# 1/s
ks=2.00e20		# N/m /(kb*T)c -> 1/m^2
xb=3.85e-10		# m

###		Funkcja Lambera omega

def lambertW(x_l, prec = 1E-6, maxiters = 1000000):
    w = 0
    for i in range(maxiters):
        we = w * pow(math.e,w)
        w1e = (w + 1) * pow(math.e,w)
        if prec > abs((x_l - we) / w1e):
            return w
        w -= (we - x_l) / (w1e - (w+2) * (we-x_l) / (2*w+2))
    raise ValueError("Brak zbieżności dla abs(z) = %f" % abs(x_l))


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
		tau=(lambertW(-exp(-a))+a)/b #Rozwiązanie równania nr 4 w Supp Mat - Szulten
		out.append(-ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p_l), epsabs=1e-16)[0])) #Wzór nr 5 w Mat Sup- Schulten.pdf
	return out

###		Funkcja błędu; suma kwadratu różnic między siłą teoretyczną dla prędkości V i siłą f z symulacji
def err(p_l):
	out=0
	fs=F(p_l,V)
	out=((fs - f)**2).sum()
	print 'k0, xb, km ', p_l, out
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

#f=array([1323e-12, 796e-12, 645e-12, 366e-12, 495e-12, 452e-12]) #f dla rozplatania C-końca Grubmiller.pdf
f1=array([1376e-12, 710e-12, 526e-12, 559e-12, 667e-12, 538e-12]) #f dla rozplatania N-końca Grubmiller.pdf

#f=array([1382e-12, 861e-12, 746e-12, 609e-12, 366e-12]) #f dla rozplatania C-końca mojej sym.
f=array([1618e-12, 822e-12, 745e-12, 650e-12, 459e-12, 429e-12]) # f dla rozplatania N-końca mojej sym.
V=array([50, 10, 5, 2, 0.8, 0.4]) # V dla rozplatania - takie same

p=[1.6e3, 2e-10, 4e20] # Parametry startowe [k0, xb, km]

p3=fmin(err, p, xtol=1e-5, maxiter=1000000, maxfun=1000000) #Minimalizacja 


#p1=[2.78690714e+03  , 1.96046057e-10  , 8.85721077e+20]
#p2=[7.94429056e+07  , 4.28672569e-11  , 9.25515245e+21]
#R2=0.99680229230609763
#R1=0.99197968889279697

x_log=logspace(-3,2, 300)
plot(x_log, F(p3, x_log))
plot(V, f, 'o')

f=f1
p2=fmin(err, p, xtol=1e-5, maxiter=1000000, maxfun=1000000)
plot(x_log, F(p2, x_log))
plot(V, f, 'x')
xscale('log')
print p3, p2
show()



