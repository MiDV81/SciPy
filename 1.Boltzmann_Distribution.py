# -*- coding: utf-8 -*-
# Distribución de energía de Boltzmann
from numpy import *
from numpy.random import randint, uniform
import matplotlib.pyplot as plt

# Parametros del problema
np = 100000 # número total de partículas
ei = 5   # energía inicial de cada partícula
ni = 800000 # número de interacciones
nis = 10000 # número de interacciones entre salida de resultados

# Energia inicial de cada particula... 
ep = ei*ones(np) # Todas las particulas con Em

# Funcion que pinta un histograma con los datos
# * ep: Distribucion de energias
# * emax: Valor maximo del eje de energias
# * e,nt: Distribucion teorica
def pinta(ep, emax, e, nt):
    # Clear the current figure
    plt.clf() 

    # Dibujamos el histograma
    plt.hist(ep, bins=emax, range=[0,emax], rwidth=0.6, color='b')

    # Dibujamos la distribucion teorica
    plt.plot(e, nt, '-r')

    # Lo adornamos todo un poco
    plt.xlabel(u'Energía')
    plt.ylabel(u'Fracción de partículas con energía dada')
    plt.title(u'interacción número '+str(i))

    # Pintar
    plt.draw()

    # Mostrar y seguir
    plt.pause(0.01)
    



# Para el dibujo
# + Energía máxima para dibujar
emax = 5*ei 
# + Distribución teórica
e = arange(emax+1)
nt = np/ei*exp(-e/ei) # la energía media es 'ei'


# Iteramos 
for i in range(ni+1):
    # Pintamos los resultados cuando ’i’ sea multiplo de ’nis’
    if mod(i, nis) == 0:
        pinta(ep, emax, e, nt)
        
    # Se sortean dos partículas diferentes
    ip1 = randint(np)
    ip2 = randint(np)
    while ip1==ip2:
        ip2 = randint(np)

    # Calculamos la energía total de ambas particulas
    et = ep[ip1] + ep[ip2] 
    if et == 0: continue


    # la energia ’et’ se reparte entre ambas particulas aleatoriamente
    e1 = uniform(0,et)
    e2 = et-e1
    ep[ip1] = e1
    ep[ip2] = e2

# El final
plt.show()
