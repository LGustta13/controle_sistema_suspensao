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
g = (kt*(Cs*s + ks))/(m*M*(s**4) + (m+M)*Cs*(s**3) + (kt*M+(m+M)*ks)*(s**2) + kt*Cs*s + kt*ks)
t = np.linspace(0,5,10000)

##################################################
# LGR
##################################################

gl = 8.79*((s - complex(-2.76,7.72))*(s - complex(-2.76, -7.72)))/(s*(s + 11.417)**2)

# Saída de saturação
GL = feedback(g*gl,1)

yl,t = step(GL,t)

plt.figure()
plt.plot(t,yl)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')

##################################################
# FREQUÊNCIA
##################################################

gf = 2*(s + 4)/(s*(s + 3))

# Saída de saturação
GF = feedback(g*gf,1)

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
          
               ]) 
#[37.39660268,  31.65214547,  -1.08405153, -20.37770771, 1.5977692437822202]
     

Dm = np.block([[0]])

GE = ss(Am,Bm,Cm,Dm)

t = np.linspace(0,5,10000)
ye,t = step(GE,t)

plt.figure()
plt.plot(t,ye)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')
plt.legend(['Saída do sistema','Sinal de Controle'])

#################

num = -3.285*(10**(-14))*s**3 - 2.846*(10**(-12))*s**2 + 5.262*(10**(4))*s + 5.84*(10**5) 
dennum = s**5 + 85.39*s**4 + 2576*s**3 + 3.401*(10**4)*s**2 + 2.144*(10**5)*s + 5.84*(10**5)
den = dennum-num

numg = 2.538*(10**8)*s + 2.818*10**9
deng = 7709*s**4 + 3.699*(10**5)*s**3 + 4.426*(10**7)*s**2 + 2.538*(10**8)*s + 2.818*10**9

numc = num/numg
denc = den/deng

CC = numc/denc

t = np.linspace(0,5,10000)
F = feedback(CC,g)
y,t = step(F,t)

plt.figure()
plt.plot(t,y)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Sinal de Controle')



plt.figure()
plt.plot(t,yl, 'r-', t, yf, 'b-', t, ye, 'g-')
plt.xlabel('Tempo [segundos]')
plt.ylabel('Resposta ao degrau')
plt.legend(['LGR', 'FREQUENCIA','ESTADOS'])

