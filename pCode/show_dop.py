# -*- coding: utf-8 -*-
from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.integrate import quad
from numpy import exp
from scipy.special import lambertw

kb=1.3807e-23
k0=1.6e-11
ks=2.0e20
xb=3.85e-10
km=6.9e20
v=28.
T=300.

params = {'axes.labelsize': 26,
          'text.fontsize': 26,
          'legend.fontsize': 20,
          'xtick.labelsize': 22,
          'ytick.labelsize': 22}
rcParams.update(params)


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


#f=[1382e-12, 861e-12, 746e-12, 609e-12, 366e-12]
V=[50, 10, 5, 2, 0.8, 0.4]
V2=[50, 10, 5, 2, 0.8, 0.4]
x_log=logspace(-7,2, 200)
f12=array([1618e-12, 822e-12, 745e-12, 650e-12, 459e-12, 429e-12])
f11=array([1372e-12, 815e-12, 720e-12, 576e-12, 443e-12, 334e-12])
p12=[2.78690714e+03  , 1.96046057e-10  , 8.85721077e+20]
p11=[7.94429056e+07  , 4.28672569e-11  , 9.25515245e+21]
p11_50=[  6.76379499e+00 ,  2.82525955e-10  , 5.77520154e+20]
p12_50=[  5.67230187e+01  , 2.56388211e-10 ,  6.47123456e+20]


f31=array([1323e-12, 796e-12, 645e-12, 366e-12, 495e-12, 452e-12]) # <-N
f32=array([1376e-12, 710e-12, 526e-12, 559e-12, 667e-12, 538e-12]) # <-C
p31=[3.97035645e-01 ,  3.47028735e-10 ,  4.35607856e+20]
p32= [ 1.21944193e-26 ,  1.10098013e-09  , 1.40503037e+20]
p31_50=[  3.42099632e-01 ,  3.49305043e-10  , 4.32387983e+20]
p32_50=[  4.14731885e-05 ,  4.42650488e-10  , 3.62924718e+20]

subplots_adjust(left=0.06, bottom=0.11, right=0.97, hspace=0.30, top=0.94, wspace=0.17)
subplot(211)
plot(x_log, F(p11, x_log)/1e-12, color='#E4010C')
plot(x_log, F(p12, x_log)/1e-12, color='#4042CB')
plot(V, f11/1e-12, 'x', markersize=21, mew=3, mec='#E4010C', alpha=0.8)
plot(V, f12/1e-12, '+', markersize=21, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12/1e-12, '+', markersize=21, mew=3, mec='#000000', alpha=0.8)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]')
ylabel('Si'+u'ł'+'a [pN]')
xscale('log')
ylim(ymin=-0.5e-9/1e-12, ymax=2.0e-9/1e-12)
legend( ('Dopasowanie N-koniec', 'Dopasowanie C-koniec', 'Wyniki N-koniec', 'Wyniki C-koniec', 'Wynik AFM'), loc=2, markerscale=0.8)

subplot(212)
plot(x_log, F(p31, x_log)/1e-12, color='#E4010C')
plot(x_log, F(p32, x_log)/1e-12, color='#4042CB')
plot(V2, f31/1e-12, 'x', markersize=21, mew=3, mec='#E4010C', alpha=0.8)
plot(V2, f32/1e-12, '+', markersize=21, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12/1e-12, '+', markersize=21, mew=3, mec='#000000', alpha=0.8)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]',)
ylabel('Si'+u'ł'+'a [pN]')
xscale('log')
ylim(ymin=-0.5e-9/1e-12, ymax=2.0e-9/1e-12)
legend( ('Dopasowanie N-koniec', 'Dopasowanie C-koniec', 'Wyniki N-koniec', 'Wyniki C-koniec', 'Wynik AFM'), loc=2, markerscale=0.8 )

show()
#text=(0.2,1.9e-9, 'cus')


