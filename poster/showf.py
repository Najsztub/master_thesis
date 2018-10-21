# -*- coding: utf-8 -*-
from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.integrate import quad

kb=1.3807e-23
k0=1.6e-11
ks=2.0e20
xb=3.85e-10
km=6.9e20
v=28.

D=k0 * sqrt(2*pi) / (km**(1.5)*xb*exp(-0.5*km*xb**2))

def s(t,vel_l, p_l):
	return exp(-p_l[0] * exp(-0.5*ks*p_l[1]**2)/(vel_l*ks*p_l[1]*(p_l[2]/(p_l[2]+ks))**(1.5))*(exp(ks*vel_l*p_l[1]*t - (0.5*(ks*vel_l*t)**2)/(ks+p_l[2]))-1))

def lambertW(x_l, prec = 1E-6, maxiters = 1000000):
    w = 0
    for i in range(maxiters):
        we = w * pow(math.e,w)
        w1e = (w + 1) * pow(math.e,w)
        if prec > abs((x_l - we) / w1e):
            return w
        w -= (we - x_l) / (w1e - (w+2) * (we-x_l) / (2*w+2))
    raise ValueError("W doesn't converge fast enough for abs(z) = %f" % abs(x_l))

def F(p_l, x_l):
	k0_l, xb_l, km_l = p_l[0], p_l[1], p_l[2]
	D_l=k0_l * sqrt(2*pi) / (km_l**(1.5)*xb_l*exp(-0.5*km_l*xb_l**2))
	b=D_l*(km_l+ks)
	out=[]
	for i in range(len(x_l)):
		vel_l = x_l [i]
		a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
		tau=(lambertW(-exp(-a))+a)/b
		out.append(-ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p_l))[0]))
	return out
	#vel_l = x
	#a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
	#tau=(lambertW(-exp(-a))+a)/b
	#return -ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p))[0])


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
f11=[1618e-12, 822e-12, 745e-12, 650e-12, 459e-12, 429e-12]
f12=[1372e-12, 815e-12, 720e-12, 576e-12, 443e-12, 334e-12]
p11=[2.78690714e+03  , 1.96046057e-10  , 8.85721077e+20]
p12=[7.94429056e+07  , 4.28672569e-11  , 9.25515245e+21]
p11_50=[  6.76379499e+00 ,  2.82525955e-10  , 5.77520154e+20]
p12_50=[  5.67230187e+01  , 2.56388211e-10 ,  6.47123456e+20]

f21=[1393e-12, 816e-12, 797e-12, 630e-12, 467e-12]
f22=[1559e-12, 872e-12, 769e-12, 650e-12, 459e-12]
p21=[5.93962667e+04 ,  1.35022476e-10  , 1.60558846e+21]
p22=[  2.32222859e+05 ,  1.24004147e-10 ,  1.66997383e+21]
p21_50=[  6.40501180e+01 ,  2.21665449e-10 ,  8.71650659e+20]
p22_50=[  3.86013440e+01 , 2.39675429e-10  , 7.48977113e+20]

f32=[1323e-12, 796e-12, 645e-12, 366e-12, 495e-12, 452e-12]
f31=[1376e-12, 710e-12, 526e-12, 559e-12, 667e-12, 538e-12]
p32=[3.97035645e-01 ,  3.47028735e-10 ,  4.35607856e+20]
p31= [ 1.21944193e-26 ,  1.10098013e-09  , 1.40503037e+20]
p32_50=[  3.42099632e-01 ,  3.49305043e-10  , 4.32387983e+20]
p31_50=[  4.14731885e-05 ,  4.42650488e-10  , 3.62924718e+20]

subplots_adjust(left=0.06, bottom=0.11, right=0.97, hspace=0.30, top=0.94, wspace=0.17)
subplot(221)
plot(x_log, F(p11, x_log), color='#E4010C')
plot(x_log, F(p12, x_log), color='#4042CB')
plot(V, f11, 'x', markersize=14, mew=3, mec='#E4010C', alpha=0.8)
plot(V, f12, '+', markersize=14, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12, '+', markersize=14, mew=3, mec='#000000', alpha=0.8)
yticks( size=14)
xticks( size=14)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]', size=16)
ylabel('Si'+u'ł'+'a [N]', size=16)
xscale('log')
ylim(ymin=-0.5e-9, ymax=2.0e-9)

x_log2=logspace(-7,2, 200)
subplot(222)
plot(x_log2, F(p11_50, x_log2), color='#E4010C')
plot(x_log2, F(p12_50, x_log2), color='#4042CB')
plot(V, f11, 'x', markersize=14, mew=3, mec='#E4010C', alpha=0.8)
plot(V, f12, '+', markersize=14, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12, '+', markersize=14, mew=3, mec='#000000', alpha=0.8)
yticks( size=14)
xticks( size=14)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]', size=16)
ylabel('Si'+u'ł'+'a [N]', size=16)
xscale('log')
ylim(ymin=-0.5e-9, ymax=2.0e-9)

subplot(223)
plot(x_log, F(p31, x_log), color='#E4010C')
plot(x_log, F(p32, x_log), color='#4042CB')
plot(V2, f31, 'x', markersize=14, mew=3, mec='#E4010C', alpha=0.8)
plot(V2, f32, '+', markersize=14, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12, '+', markersize=14, mew=3, mec='#000000', alpha=0.8)
yticks( size=14)
xticks( size=14)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]', size=16)
ylabel('Si'+u'ł'+'a [N]', size=16)
xscale('log')
ylim(ymin=-0.5e-9, ymax=2.0e-9)

subplot(224)
plot(x_log, F(p31_50, x_log), color='#E4010C')
plot(x_log, F(p32_50, x_log), color='#4042CB')
plot(V2, f31, 'x', markersize=14, mew=3, mec='#E4010C', alpha=0.8)
plot(V2, f32, '+', markersize=14, mew=3, mec='#4042CB', alpha=0.8)
plot(1e-6, 50e-12, '+', markersize=14, mew=3, mec='#000000', alpha=0.8)

yticks( size=14)
xticks( size=14)
xlabel('Szybko'+u'ś'+u'ć'+' rozci'+u'ą'+'gania [nm/ns]', size=16)
ylabel('Si'+u'ł'+'a [N]', size=16)
xscale('log')
ylim(ymin=-0.5e-9, ymax=2.0e-9)

show()
#text=(0.2,1.9e-9, 'cus')


