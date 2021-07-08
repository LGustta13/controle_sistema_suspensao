# -*- coding: utf-8 -*-

## Saturação

%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from control.matlab import *

###############################################
# DADOS DO PROBLEMA
###############################################

## Variáveis para o modelo GOL
kt = 178140
Cs = 1425
ks = 15817
m = 34.2
M = 225.4

##################################################
# Modelo do sistema
##################################################

# Função de transferência de malha aberta
s = tf('s')
g = (1/s)*(kt*(Cs*s + ks))/(m*M*(s**4) + (m+M)*Cs*(s**3) + (kt*M+(m+M)*ks)*(s**2) + kt*Cs*s + kt*ks)
t = np.linspace(0,7,10000)

d = np.zeros(len(t))
d[4000:5000] = 1

##################################################
# LGR
##################################################

gl = (s**2 + 5.52*s + 67.22)/(s**3 + 22.83*s**2 + 130.3*s)

# Saída de saturação
GL = feedback(minreal(gl*g),1)

yl,t = step(GL,t,0,d)

plt.figure()
plt.plot(t,yl)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')

##################################################
# FREQUÊNCIA
##################################################

gf = 2*(s + 4)/(s + 3)

# Saída de saturação
GF = feedback(minreal(gf*g),1)

yf,t = step(GF,t)

plt.figure()
plt.plot(t,yf)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')

##################################################
# ESTADOS
##################################################

Am = np.block([[-8.53853634e+01,  2.57621579e+01,  3.40143976e+01,
         -1.61738210e+01, -1.59776924e+00],
        [-1.00000000e+02,  9.34698797e-15, -2.23173714e-14,
          1.39192692e-14,  0.00000000e+00],
        [ 0.00000000e+00,  1.00000000e+01,  7.53564761e-15,
         -1.65687226e-14,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00, -1.00000000e+01,
         -2.83680058e-15,  0.00000000e+00],
        [-0.00000000e+00, -5.12236868e-15, -3.29303461e+01,
          3.65515287e+01,  0.00000000e+00]])

Bm = np.block([[0.],
       [0.],
       [0.],
       [0.],
       [1.]])

Cm = np.block([[ 0.00000000e+00,  5.12236868e-15,  3.29303461e+01,
         -3.65515287e+01,  0.00000000e+00],
               [37.39660268,  31.65214547,  -1.08405153, 
        -20.37770771, 1.5977692437822202]]) 

Dm = np.block([[0],[0]])

GE = ss(Am,Bm,Cm,Dm)

t = np.linspace(0,5,10000)
ye,t = step(GE,t)

plt.figure()
plt.plot(t,ye)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')
plt.legend(['Saída do sistema','Sinal de Controle'])

plt.figure()
plt.plot(t,ye, 'r-', t, yl, 'b-', t, yf, 'g-')
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')
plt.legend(['Saída do sistema','Sinal de Controle ESTADOS', 'Sinal de Controle LGR', 'Sinal de Controle FREQUENCIA'])

