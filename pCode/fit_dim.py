from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import fmin
from scipy.optimize import leastsq
from scipy.integrate import quad

T=300
kb=1.3807e-23
k0=1.6e-11
ks=2.0e20
xb=3.85e-10
km=6.9e20
#v=28.

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
	print k0_l, xb_l, km_l
	for i in range(len(x_l)):
		vel_l = x_l [i]
		a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
		tau=(lambertW(-exp(-a))+a)/b
		out.append(-ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p_l), epsabs=1e-16)[0]))
	return array(out)
	#vel_l = x
	#a=(xb_l*D_l*(km_l+ks)**2/(ks*vel_l)+1)
	#tau=(lambertW(-exp(-a))+a)/b
	#return -ks*300*kb*(xb_l-vel_l*quad(s,0,tau, args=(vel_l, p))[0])


def residuals(p_l, y, x):
	err = y - F(p_l,x)
	print 'out ', p_l, err
	return err


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
f=[1618, 822, 745, 650, 459, 429]
V=[50, 10, 5, 2, 0.8, 0.4]
p=[1.6, 3.85, 280]
x_log=logspace(-3,2, 300)
p1=fmin(err, p, xtol=1e-5, maxiter=1000000, maxfun=1000000)
p1=[2.78690714e+03  , 1.96046057e-10  , 8.85721077e+20]
p2=[7.94429056e+07  , 4.28672569e-11  , 9.25515245e+21]
p1_50=[  6.76379499e+00 ,  2.82525955e-10  , 5.77520154e+20]
p2_50=[  5.67230187e+01  , 2.56388211e-10 ,  6.47123456e+20]
R2=0.99680229230609763
R1=0.99197968889279697
R1_50=0.99420449831835578
R2_50=0.98624415872029281
plot(x_log, F(p1, x_log)
plot(V, f, 'o')
show()



