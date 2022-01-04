from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def a_paraca (t, x, v):
	return g-gamma*abs(v)*v

# g, gamma son variables globales
def Euler (acel, x, v, tv, n=10):
	nv, xv, vv = tv.size, zeros_like(tv), zeros_like(tv)
	xv[0], vv[0] = x, v
	for k in range(1, nv):
		t, Dt = tv[k-1], (tv[k]-tv[k-1])/float(n)
		for i in range(n):
			a = acel(t, x, v)
			t, v, x = t+Dt, v+a*Dt, x+v*Dt
		xv[k], vv[k] = x, v
	return (xv, vv)

def EulerCromer (acel, x, v, tv, n=10):
	nv, xv, vv = tv.size, zeros_like(tv), zeros_like(tv)
	xv[0], vv[0] = x, v
	for k in range(1, nv):
		t, Dt = tv[k-1], (tv[k]-tv[k-1])/float(n)
		for i in range(n):
			a = acel(t, x, v)
			t, v = t+Dt, v+a*Dt
			x=x+v*Dt
		xv[k], vv[k] = x, v
	return (xv, vv)

def RK2(acel,x,v,tv,n=10):
	nv, xv, vv = tv.size, zeros_like(tv), zeros_like(tv)
	xv[0], vv[0] = x, v
	for k in range(1, nv):
		t, Dt = tv[k-1], (tv[k]-tv[k-1])/float(n)
		for i in range(n):
			a = acel(t, x, v)
			vast=v+a/2*Dt
			xast=x+v/2*Dt
			tast=t+Dt/2
			a=acel(tast,xast,vast)
			v=v+a*Dt
			x=x+vast*Dt
			t=t+Dt			
		xv[k], vv[k] = x, v
	return (xv, vv)


g=9.80 # m/s^2
gamma =0.10 # m^-1

tv = linspace(0, 3, 30)
x, v = 0, 0 # condiciones iniciales
xv, vv = Euler(a_paraca, x, v, tv)
xv3,vv3=RK2(a_paraca,x,v,tv)
xv4,vv4=EulerCromer(a_paraca,x,v,tv)

def funder (t, y):
	x, v = y[0], y[1]
	return array([v, g-gamma*abs(v)*v])

yv2 = solve_ivp(funder, [tv[0], tv[-1]], [x, v], t_eval=tv).y.T
xv2, vv2 = yv2[:,0], yv2[:,1]

difEulerx,difEulerv=abs(xv-xv2),abs(vv-vv2)
difRK2x,difRK2v=abs(xv3-xv2),abs(vv3-vv2)
difEuCrx,difEuCrv=abs(xv4-xv2),abs(vv4-vv2)

plt.figure(figsize=(15,9))
plt.suptitle('Comparación métodos con solve_ivp')

plt.subplots_adjust(wspace=0.25,hspace=0.40)
plt.subplot(2,2,1,title='Euler')
plt.plot(tv,xv,'b-',label='x(t) Euler')
plt.plot(tv,vv,'r-',label='v(t) Euler')
plt.xlabel("t (s)")

plt.subplot(2,2,2,title='solve_ivp')
plt.plot(tv,xv2,'b-',label='x(t) solve_ivp')
plt.plot(tv,vv2,'r-',label='v(t) solve_ivp')
plt.xlabel("t (s)")

plt.subplot(2,2,3,title='RK2')
plt.plot(tv,xv3,'b-',label='x(t) RK2')
plt.plot(tv,vv3,'r-',label='v(t) RK2')
plt.xlabel("t (s)")

plt.subplot(2,2,4,title='EulerCromer')
plt.plot(tv,xv4,'b-',label='x(t)')
plt.plot(tv,vv4,'r-',label='v(t)')
plt.xlabel("t (s)")

plt.legend()
plt.show()

plt.figure(figsize=(15,9))
plt.suptitle('Comparación métodos con solve_ivp')

plt1=plt.subplot(1,3,1,title='Euler')
plt.plot(tv,difEulerx,'b-',label='x(t)')
plt.plot(tv,difEulerv,'r-',label='v(t)')
plt.xlabel("t (s)")
plt1.set_yscale('log')

plt2=plt.subplot(1,3,2,title='RK2')
plt.plot(tv,difRK2x,'b-',label='x(t)')
plt.plot(tv,difRK2v,'r-',label='v(t)')
plt.xlabel("t (s)")
plt2.set_yscale('log')

plt3=plt.subplot(1,3,3,title='EulerCromer')
plt.plot(tv,difEuCrx,'b-',label='x(t)')
plt.plot(tv,difEuCrv,'r-',label='v(t)')
plt.xlabel("t (s)")
plt3.set_yscale('log')

plt.legend()
plt.show()
