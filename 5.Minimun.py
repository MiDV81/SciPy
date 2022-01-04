from numpy import *
from scipy.optimize import fmin
import matplotlib.pyplot as plt

g, h, sigma, rho=9.8 , 1.0, 0.00728, 1000

def deriv(f, FAC=2.0e-6):
    def g(x):
        dx = where(x != 0.0, FAC*x, FAC) # para usar vectorizada
        return (f(x+dx)-f(x-dx))/(2.0*dx)
    return g


def fmax (f):
    return lambda x: -f(x)


def biseccion (f, a, b, xtol=1.0e-8):   
    fa, fb = f(a), f(b)    
    if fa == 0.0: return a
    if fb == 0.0: return b
    assert fa*fb < 0, "Error al aplicar la bisección. No hay cambio de signo"   
    n = 1+int(log(abs(b-a)/xtol)/log(2.0))

    for i in range(n):        
        x = 0.5*(a+b)
        fx = f(x)             
        if fx == 0.0:   return x      
        if fb*fx < 0.0:
            a, fa = x, fx       
        else:
            b, fb = x, fx                
    return x


def vel(lam):
    global nfun
    nfun +=1
    return sqrt((g*lam/(2*pi)+2*pi*sigma/(rho*lam))*tanh(2*pi*h/lam))


def barrido(f,a,b,n=1000):
    global nfun
    x=linspace(a,b,n)
    fx=f(x)
    ia=argmin(fx)
    nfun=n
    return x[ia],fx[ia]

def SeccionAurea(f,a,b,err=1e-8):
    aur=(sqrt(5)-1)/2
    for i in range(100): 
        x1=a+(b-a)*aur #dos puntos entre a y b, de manera que se va acotando el intervalo donde se busca
        x2=b-(b-a)*aur       
        if f(x1)<f(x2): 
            a=x1            
        else:
            b=x2
        if abs(b-a)<err:
            break            
    return a,f(a)


a=1.0e-3
b=0.2

#biseccion
print('Bisección de la derivada')
nfun=0
der=deriv(vel)
xmin1=biseccion(der,a,b)
vmin1=vel(xmin1)
print('x=%20.16f; vmin=%20.16f; nfun=%d \n '%(xmin1,vmin1,nfun))


#barrido
print('Barrido')

nfun=0
xmin2,vmin2=barrido(vel,a,b)
print('x=%20.16f; vmin=%20.16f; nfun=%d \n '%(xmin2,vmin2,nfun))


#scipy
print('fmin de scipy')
nfun=0
xmin3=fmin(vel,a,xtol=1e-10,disp=False)
vmin3=vel(xmin3)
print('x=%20.16f; vmin=%20.16f; nfun=%d \n ' %(xmin3,vmin3,nfun))

#sección aurea
print('Sección Aurea')

nfun=0
xmin4,vmin4=SeccionAurea(vel,a,b)
print('x=%20.16f; vmin=%20.16f; nfun=%d \n '%(xmin4,vmin4,nfun))

x=linspace(a,b,1000)
fx=vel(x)
plt.figure()
plt.plot(x,fx)
plt.show()