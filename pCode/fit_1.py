from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.integrate import quad

kb=1.3807e-23
k0=1.6e-11
ks=5.031e20
xb=3.85e-10
km=6.9e20
v=28.

D=k0 * sqrt(2*pi) / (km**(1.5)*xb*exp(-0.5*km*xb**2))

def s(t,vel_l, p_l):
	return exp(-p_l[0] * exp(-0.5*ks*p_l[1]**2)/(vel_l*ks*p_l[1]*(p_l[2]/(p_l[2]+ks))**(1.5))*(exp(ks*vel_l*p_l[1]*t - (0.5*(ks*vel_l*t)**2)/(ks+p_l[2]))-1))

def lambertW(x_l, prec = 1E-8, maxiters = 1000000):
    w = 0
    for i in range(maxiters):
        we = w * pow(math.e,w)
        w1e = (w + 1) * pow(math.e,w)
        if prec > abs((x_l - we) / w1e):
            return w
        w -= (we - x_l) / (w1e - (w+2) * (we-x_l) / (2*w+2))
    raise ValueError("W doesn't converge fast enough for abs(z) = %f" % abs(x_l))

def F(p_l, x_l):
	k0_l = 3.3e-4
	xb_l = 2.5e-10
	km_l = p_l[2]
	p_l[0]=k0_l
	p_l[1]=xb_l
	D_l=k0_l * sqrt(2*pi) / (km_l**(1.5)*xb_l*exp(-0.5*km_l*xb_l**2))
	#k_l=1/( sqrt(2*pi))*D_l*km_l**(1.5)*xb_l*exp(-0.5*km_l*xb_l**2)
	b=D_l*(km_l+ks)
	out=[]
	for i in range(len(x_l)):
		vel_l = x_l [i]
		a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
		tau=(lambertW(-exp(-a))+a)/b
		out.append(-ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p_l), epsabs=1e-16)[0]))
	return out
	#vel_l = x
	#a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
	#tau=(lambertW(-exp(-a))+a)/b
	#return -ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p))[0])

def err(p_l):
	y=array(f)
	x=array(V)
	out=0
	fs=F(p_l,x)
	for i in range(len(x)):
		out+=(fs[i] - y[i])**2
	print 'out ', p_l, out
	return out

def R_sq(par, x_l, f_sim):
	f_l=F(par,x_l)
	f_sim_m=0.
	l=len(x_l)
	for i in range(l):
		f_sim_m+=f_sim[i]
	f_sim_m=f_sim_m/l
	s_l, s_m =0. , 0.
	for j in range(l):
		s_l+=(f_sim[j]-f_l[j])**2
		s_m+=(f_sim[j]-f_sim_m)**2
		#print s_l, f_sim[j]-f[j], s_m
	return 1-s_l/s_m
def R_sum(par, x_l, f_sim):
	f_l = F(par,x_l)
	s_l = 0.
	for j in range(len(x_l)):
		s_l+=(f_sim[j]-f_l[j])**2
	return s_l

f=[1581e-12, 1208e-12, 997e-12, 863e-12, 794e-12, 592e-12, 491e-12]
V=[28., 8., 2.8, 0.8, 0.28, 0.028, 0.0028]
fmin(err, [1.6e-11, 3.85e-10, 6.9e20], xtol=1e-6, maxiter=1000000, maxfun=1000000)
p=[1.6e-11, 3.85e-10, 6.9e20]
plot(V, F([p[0], p[1], p[2]], V))

x_log=logspace(-3,2, 300)
plot(x_log, F(p1, x_log))
plot(V, f, 'o')
xscale('log')
show()



