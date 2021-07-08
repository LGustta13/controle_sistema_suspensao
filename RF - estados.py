# -*- coding: utf-8 -*-

## Relatório 3 - Simulações

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

# Requisitos de desempenho
UP = 0.02 # ultrapassagem percentual
e = 0     # constante de erro de velocidade
Ts = 0.7  # tempo de assentamento

gss = tf2ss(g)
A,B,C,D = ssdata(gss)

##################################################
# Polos desejados
##################################################
zeta = -np.log(UP)/np.sqrt(np.pi**2+np.log(UP)**2) # fator de amortecimento requerido
print('zeta: ',np.round(zeta,2))

# Polos dominantes
wn = 4/(Ts*zeta)

sd1 = complex(-zeta*wn,wn*np.sqrt(1-zeta**2))
sd2 = complex(-zeta*wn,-wn*np.sqrt(1-zeta**2))

# Polos substituindo zero
sd3 = -11.09964912;

# Polos adicionais 
sd4 = -5*wn*zeta # distante do polo dominante
sd5 = -6*wn*zeta # distante do polo dominante

# Polos desejados
sd = np.block([sd1,sd2,sd3,sd4,sd5])


##################################################
# Projeto do controle integral
##################################################

# Planta modificada
Ab = np.block([[A, np.zeros((4,1))], [-C,0]])
Bb = np.block([[B],[0]])

Kb = place(Ab,Bb,sd)

# Ganhos da realimentação de estado e do termo integral
K = Kb[0,0:4]
Ki = -Kb[0,4]

###################################################
# Malha fechada
###################################################
Am = np.block([[A-B*K,B*Ki],[-C,0]])
Bm = np.block([[np.zeros((4,1))],[1]])
Cm = np.block([C,0])
Dm = np.block([0])

T = ss(Am,Bm,Cm,Dm)

###################################################
# Simulação da resposta ao degrau
###################################################
t = np.linspace(0,5,10000)
y,t = step(T,t)

plt.figure()
plt.plot(t,y)
plt.xlabel('Tempo [segundos]')
plt.ylabel('Posição do sistema de amortecimento')

mfinfo = stepinfo(T,t)

print('Ultrapassagem percentual da malha fechada:',np.round(mfinfo['Overshoot'],2))
print('Tempo de acomodação da malha fechada:',np.round(mfinfo['SettlingTime'],2))


###################################################
# Descobrindo a função de transferência do controlador
###################################################

Tt = ss2tf(T)
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
plt.ylabel('Posição do sistema de amortecimento')

mfinfo = stepinfo(T,t)

print('Ultrapassagem percentual da malha fechada:',np.round(mfinfo['Overshoot'],2))
print('Tempo de acomodação da malha fechada:',np.round(mfinfo['SettlingTime'],2))
