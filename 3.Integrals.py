from numpy import *
from scipy.integrate import quad
import matplotlib.pyplot as plt
from numpy.random import uniform

Amp=arange(0,180)
a,b=0.0,pi/2

def f(ang):
	return 1/sqrt(1-k**2*sin(ang)**2)

def integ(f,a,b,n=100):
	h=(b-a)/n
	s=0
	for i in range(n):
		xi=a+(i+0.5)*h
		s+=f(xi)

	return h*s

def AGM(x,y,eps=1e-15):
	for i in range(100):
		if abs(x-y)< eps:break
		x=(x+y)/2
		y=sqrt(x*y)
	return x

def Montecarlo(f,a,b,n=int(1e5)): 
	x=uniform(a,b,n) #creamos un array con n numeros aleatorios entre a y b
	sumatorio=sum(f(x)) #hacemos un sumatorio con todos los puntos de la función elegidos
	integral=(b-a)/n*sumatorio #hacemos la media del area de todos los rectangulos con altura f(xi) y base (b-a)
	return integral

I,II,III,IIII=[],[],[],[]
print("%8s %16s %16s %16s %16s" % ("A_max", "rectángulos", "quad", "AGM","Montecarlo"))
for A in Amp:
	k=sin(deg2rad(A)/2)
	I1=2/pi*integ(f,a,b)
	I2=2/pi*quad(f,a,b,epsabs=1e-15)[0]
	I3=1/AGM(1,sqrt(1-k**2))
	I4=2/pi*Montecarlo(f,a,b)
	print("%8.1f %16.12f %16.12f %16.12f %16.12f" %(A, I1, I2, I3, I4))
	I.append(I1)
	II.append(I2)
	III.append(I3)
	IIII.append(I4)

plt.figure('Período de un péndulo')
plt.plot(Amp,I,'r-',label='Cuadratura')
plt.plot(Amp,II,'g--',label='quad')
plt.plot(Amp,III,'b-',label='AGM')
plt.plot(Amp,IIII,'y-.',label='Montecarlo') 
plt.legend()
plt.title('Comparación métodos integración')
plt.xlabel('Amplitud (º)')
plt.ylabel('Período (s)')

plt.show() #se puede ver que los métodos quad, cuadratura y Montecarlo coinciden casi totalmente