from numpy import *
import argparse
import matplotlib.pyplot as plt

def funplot(numfunc,f,numsub,subconf,ran,lab,ttl,legend,log,save,line):
	sumafunc=0 #nos permite contar en que posición de la lista de funciones se encuentra
	plt.figure(figsize=(12,9)) #cambia el tamaño de la figura creada
	for e in range(0,numsub):
		x=linspace(ran[0+e*3],ran[1+e*3],ran[2+e*3])
		plt.subplot(subconf[0],subconf[1],e+1)
		plt.subplots_adjust(wspace=0.25,hspace=0.40) #modificar la distancia entre las subplots
		if lab!=0: #si no se especifica, no se pondran titulos en los ejes
			plt.xlabel(lab[0+2*e],weight='bold') #poner los titulos en negrita
			plt.ylabel(lab[1+2*e],weight='bold')
		if ttl!=0: #si no se especifica, no se pondra titulos a las subplots
			plt.title(ttl[e],weight='bold',fontsize='large')
		if legend!=0:
			if log[e]==0:
				for i in range(0,numfunc[e]):				
					plt.plot(x,eval(f[sumafunc+i])(x), line[sumafunc+i], label=legend[sumafunc+i])
				
			else:
				for i in range(0,numfunc[e]):
					plt.loglog(x,eval(f[sumafunc+i])(x),line[sumafunc+i],label=legend[sumafunc+i])
			plt.legend()
		else: 
			if log[e]==0:
				for i in range(0,numfunc[e]):
					plt.plot(x,eval(f[sumafunc+i])(x),line[sumafunc+i])
			else:
				for i in range(0,numfunc[e]):					
					plt.loglog(x,eval(f[sumafunc+i])(x),line[sumafunc+i]) 
		sumafunc+=numfunc[e]
	if save!=0:
		plt.savefig(save)
	plt.show()
	
	
#Algunas funciones para representar:
def x2(x):
	return x**3	
def x1(x):
	return x**2+x+1	
def x3(x):
	return x**3+x**2+x+1	
def x4(x):
	return x**3+x**2
def x5(x):
	return x+1
def x6(x):
	return x**3-x*+2-x		
def x7(x):
	return x**4
	

#Argumentos de parser:

parser = argparse.ArgumentParser(description='Crear una función')
parser.add_argument('--numsub',type=int,default=1,help='número de subplots que se quieren')
parser.add_argument('--subconf',type=int,default=[1,1],nargs=2,help='configuración de las subplots (filas columnas)')
parser.add_argument('--func',required=True,type=str,nargs='*',help='funciones que se quieren representar')
parser.add_argument('--range',required=True,type=int,nargs='*',help= 'determina los rangos que representa la función, funciona en grupos de 3 (xmin,xmax y número de puntos)')
parser.add_argument('--labels',type=str,default=0,nargs='*',help= 'Permite añadir títulos a los ejes')
parser.add_argument('--title',type=str,default=0,nargs='*',help='Permite añadir un título a cada subplot')
parser.add_argument('--legend',type=str,default=0,nargs='*',help='Permite añadir una leyenda a las gráficas, se puede poner en todas o no poner en ninguna')
parser.add_argument('--log',type=int,default=[0,0,0,0],nargs='*',help='Permite cambiar a escala logarítmica (0-no, 1-si), y se requieren tantos como subplots se quieran') #si no queremos especificar permite como máximo 4 subplots
parser.add_argument('--save',type=str,default='0',help='especificar el archivo donde se quiere guardar el resultado')
parser.add_argument('--line', type=str, default='b-''g-''y-''k-''m-', nargs='*',help='formato de las líneas que representan cada función, tantas como funciones se quieran') #permite como máximo 10 funciones si no se especifica
parser.add_argument('--numfunc',type=int,default=[1,1,1,1],nargs='*',help='número de funciones que se quieren por subplot')


#Llamamos a la función:

arg=parser.parse_args()
funplot(arg.numfunc,arg.func,arg.numsub,arg.subconf,arg.range,arg.labels,arg.title,arg.legend,arg.log,arg.save,arg.line)

