# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 03:52:01 2019

@author: wabng
"""

#  Este programa realiza la solucion de manera numerica
#  por el metodo de diferencias finitas, del sistema acoplado de ecuaciones
#  parciales diferenciales, una ecuacion nos determina la evolucion de la parte imaginaria de la
#  funcion de onda, mientras que otra nos determina la evolución de la parte real de la funcion de onda

#    INICIALIZANDO RECURSOS NECESARIOS
import matplotlib
matplotlib.use('TkAgg')
from pylab import *

#  PARAMETROS 
n = 500 # tamñaño de iteraciones de la malla
Dh = 1. / n # tamaño de paso espacial por cada iteracion punto
Dt = 0.002 # incremento temporal 

k=0.02


def initialize():
    global u, v, nextu, nextv
    u = zeros([n, n])
    v = zeros([n, n])
    for x in range(n):
        for y in range(n):
            u[x, y] =k*exp(-(((10*x*Dh)**2 + ((10*y*Dh)**2))/2 ))*(4*(((10*Dh*x)**2))-2)*(8*(((10*y*Dh)**3))-12*(10*y*Dh)) 
            v[x, y] = 0 
            #v[x, y] =exp(-((x*Dh)**2 + (y*Dh)**2)/2 )*(4*((Dh*x)**2)-2)*(8*(((y*Dh)**3)-12*(y*Dh))) 
    nextu = zeros([n, n])
    nextv = zeros([n, n])


def observe():
    global u, v, nextu, nextv
    subplot(1, 2, 1)
    cla()
    imshow(u, vmin = -0.5, vmax =0.5, cmap = cm.cool)
    title('Parte real')
    subplot(1, 2, 2)
    cla()
    imshow(v, vmin = -0.5, vmax = 0.5, cmap = cm.cool)
    title('Parte imaginaria')
    
    
def update():
    
    global u, v, nextu, nextv
    for x in range(n):
        for y in range(n):
            # state-transition function
            uC, uR, uL, uU, uD = u[x,y], u[(x+1)%n,y], u[(x-1)%n,y], \
            u[x,(y+1)%n], u[x,(y-1)%n]
            vC, vR, vL, vU, vD = v[x,y], v[(x+1)%n,y], v[(x-1)%n,y], \
            v[x,(y+1)%n], v[x,(y-1)%n]
            d_re = (uR + uL + uU + uD - 4 * uC) / (Dh**2)
            d_im = (vR + vL + vU + vD - 4 * vC) / (Dh**2)
            nextu[x,y] = uC - (d_im-(((x*Dh)**2)+((y*Dh)**2))*vC) * Dt/2
            nextv[x,y] = vC+((d_re-((((x*Dh)**2)+((y*Dh)**2))*uC)) * Dt/2)
      
    
    u, nextu = nextu, u
    v, nextv = nextv, v
          
initialize()
    
#for i in range(100):
 #   update()
#observe()
#show()
    
    
    
import pycxsimulator
pycxsimulator.GUI(stepSize = 1).start(func=[initialize, observe, update])
