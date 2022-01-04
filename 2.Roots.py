from numpy import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

R = 8.314472
a, b = 0.3640, 4.267e-5
T = arange(100.0, 1000.1, 100.0)
P = logspace(0.0, 2.0, 13)


def f(v):
    
    global nfun
    nfun += 1    
    
    return (pbar*v**2+a)*(v-b)-R*t*v**2

der = lambda v: 3*pbar*v**2-2*b*pbar*v+a-2*R*t*v #definimos la derivada de f para el método Newton-Raphson


def biseccion (f, a, b, xtol=1.0e-8):
   
    fa, fb = f(a), f(b)
    
    if fa == 0.0: return a
    if fb == 0.0: return b

    assert fa*fb < 0, "Error al aplicar la bisección. No hay cambio de signo"
   
    n = 1+int(log(abs(b-a)/xtol)/log(2.0))
    
    for i in range(n):
        
        x = 0.5*(a+b)
        fx = f(x)        
        
        if fx == 0.0:
            return x      
        if fb*fx < 0.0:
            a, fa = x, fx       
        else:
            b, fb = x, fx
                
    return x



def falsi(f,a,b,xtol=1e-10):

    fa,fb=f(a),f(b)
    
    if fa==0: return a
    if fb==0: return b
    assert fa*fb<0, 'Error en regula'

    while True:

        c=b-fb*(a-b)/(fa-fb)#se igualan las pendientes entre la recta que pasa por (a,f(a)) y (c,f(c)) con la recta que pasa por (b,f(b)) y (c,f(c)) y se saca c (punto donde f(c)=0). 
        fc=f(c)

        if abs(fc)<=xtol: break

        if fc*fb < 0.0:
            a,fa = c,fc
        else: 
            b,fb = c,fc

    return c

def NewRap(f,a,xtol=1e-8):
    fa=f(a)
    da=der(a)

    if fa==0: return a

    while True:
        a = a-(fa/da) #se encuentra la recta tangente a fa y se saca el punto donde se cruza con el eje de abscisas, el cual será nuestro nuevo 'a'
        fa,da = f(a),der(a)
        if abs(fa) <= xtol: break

    return a


print("===============================================================")
print("Ecuación de van der Waals comparada con la de los gases ideales")
print("Gas: CO2, constantes de van der Waals (SI):")
print("a = %e, b = %e" % (a, b))
print("===============================================================")

plt.figure('Interacciones de distintos métodos',figsize=(16,9))
plt.suptitle('Nº interacciones vs. Presiones (bar)')

for t in T: 
    Nbis=[]
    Nfsol=[]
    Nfal=[]
    NNew=[]

    print('')
    print('%4s %9s %14s %12s %7s %12s %7s %12s %9s %14s %7s'%('T','P','v ideal','v bisec','nfun','v fsolve','nfun','v falsi', 'nfun','v New-Rap','nfun'))
    print('%5s %10s %12s %12s %7s %12s %7s %12s %9s %14s %7s'%('(K)','(atm)','(L/mol)','(L/mol)','','(L/mol)','','(L/mol)','','(L/mol)',''))
    plt.subplot(2,5,t/100)
    plt.title(str(int(t))+' K')
    
    for e in P: 

        nfun=0
        pbar = 101325.0*e
        vid = R*t/pbar 
        sol1 = biseccion(f, 1.0e-6*vid, 10*vid)      
        nbis = nfun 
        Nbis.append(nbis)
        
        nfun = 0        
        sol2 = fsolve(f, vid)       
        nfsol = nfun
        Nfsol.append(nfsol)

        nfun=0
        sol3=falsi(f,1.0e-3*vid,1.2*vid)
        nfal=nfun
        Nfal.append(nfal)

        nfun=0
        sol4=NewRap(f,vid)
        nNew=nfun
        NNew.append(nNew)
        
        print(" %5d %8.2f %15.8f %12.8f %4d %15.8f %4d %15.8f %8d %15.8f %4d" %(t, e, 1.0e3*vid, 1.0e3*sol1, nbis, 1.0e3*sol2,nfsol,1.0e3*sol3,nfal,1.0e3*sol4,nNew))

    plt.plot(P,Nbis, 'y-',label='Bisección')
    plt.plot(P,Nfsol, 'g-',label='f solve')
    plt.plot(P,NNew, 'r-',label='Newton-Raphson')
    plt.plot(P,Nfal,'b-', label='Regula falsi')
    if t==600:
        plt.legend()

plt.show()

#el método regula falsi realiza demasiadas interacciones, cuando en teoría debería realizar aproximadamente las de la bisección. El método de Newton-Raphson es el que menos interacciones realiza de todos los probados. 