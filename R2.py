%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from control.matlab import *

## Variáveis para o modelo GOL
kt = 178140
Cs = 1425
ks = 15817
m = 34.2
M = 225.4

## Função de transferência para a PLANTA (Gp)
numGp = [kt*Cs, kt*ks]
denGp = [m*M, (m+M)*Cs, (kt*M+(m+M)*ks), kt*Cs, kt*ks]
Gp = tf(numGp, denGp)

## Função de transferência para Gd e Gi (derivativo e integral)
numGd = [1,0.31]
denGd = [1]
Gd = tf(numGd,denGd)

numGi = [1,0.019]
denGi = [1,0]
Gi = tf(numGi,denGi)

## Lugar das raízes
rlocus(Gp, plotstr=True, Plot=True, PrintGain=True, grid=True);

## Ganho K
k = 5

## Malha fechada
T1 = feedback(k*Gp,1)
T2 = feedback(Gp,1)

## Simulação da resposta ao degrau
y1,t1 = step(T1)
y2,t2 = step(T2)
fig = plt.figure()
plt.plot(t1,y1,'r-s', t2,y2,'b-o')
plt.yticks(np.linspace(0,1.1,10))
plt.xticks(np.linspace(0,max(t1),5))
plt.grid()
plt.title('Resposta ao degrau')
plt.xlabel('Tempo[s]')
plt.ylabel('Amplitude')
plt.legend('Com compensador','Sem compensador')
plt.show()

## Simulação da resposta ao degrau em malha aberta
y,t = step(Gp)
fig = plt.figure()
plt.plot(t,y,'r')
plt.yticks(np.linspace(0,1.1,10))
plt.xticks(np.linspace(0,max(t1),5))
plt.grid()
plt.title('Resposta ao degrau')
plt.xlabel('Tempo[s]')
plt.ylabel('Amplitude')
plt.show()


