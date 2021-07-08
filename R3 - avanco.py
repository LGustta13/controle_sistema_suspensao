# -*- coding: utf-8 -*-

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

# Função de transferência de malha aberta
s = tf('s')
g = (kt*(Cs*s + ks))/(m*M*(s**4) + (m+M)*Cs*(s**3) + (kt*M+(m+M)*ks)*(s**2) + kt*Cs*s + kt*ks)

# Requisitos de desempenho
UP = 0.02 # ultrapassagem percentual
e = 0.01     # constante de erro de velocidade
Ts = 0.7  # tempo de assentamento

###############################################
# AJUSTE DE GANHO E CÁLCULO DE MARGEM DE FASE
###############################################

# Ganho K para requisito de erro estacionario
kpi = dcgain(g) # Kp inicial
kp_novo = (1/e)-1 
print('Kp novo:',np.round(kp_novo,2)) 
k = kp_novo/kpi 
print('Ganho K necessário para condição de erro:',np.round(k,2))

# Fator de amortecimento requerido
zeta = -np.log(UP)/np.sqrt(np.pi**2+np.log(UP)**2)

# Margem de fase requerida
pm_req = np.rad2deg(np.arctan((2*zeta)/(np.sqrt(-2*(zeta**2)+np.sqrt(1+4*(zeta**4))))))
print('Margem de fase requerida: ',round(pm_req,1))

# Diagrama de Bode
fig = plt.figure()
bode(g,k*g,dB=True,Hz=False)

# Margens de estabilidade inserindo ganho K
gm,pm,wgm,wpm = margin(k*g)
print('Margem de fase sem compensação: ',round(pm,1))

###############################################
# PROJETO DO COMPENSADOR EM AVANÇO
###############################################

# Contribuição de fase do compensador => PHImax
pm_seg = 5 # fase de segurança
phi = (pm_req + pm_seg) - pm # compensacao de fase do compensador - PHImax
print('Fase máxima do compensador: ',round(phi,1))

# Parâmetro beta do compensador
beta = (1-np.sin(np.deg2rad(phi)))/(1+np.sin(np.deg2rad(phi)))
print('beta: ',beta)

# Magnitude do compensador na fase maxima
A = k/np.sqrt(beta)
AdB = 20*np.log10(A)
print('Magnitude do compensador na fase máxima: ',AdB)

# Nova frequencia de margem de fase => wmax
fig = plt.figure()
bode(g,dB=True,Hz=False) # procurar a amplitude -AdB para obter wmax
wmax = 190 # valor encontrado no diagrama de bode de magnitude

# T e polos/zeros do compensador
T = 1/(wmax*np.sqrt(beta))
zc = -1/T
pc = -1/(beta*T)

# Função de trasferência do compensador
kc = k/beta # compensador deve ter ganho estatico unitario
gc = kc*(s-zc)/(s-pc)

###############################################
# Verificação de resultados
###############################################

gmc,pmc,wgmc,wpmc = margin(gc*g) # verificação das margens 
print('Margem de fase do sistema compensado: ',round(pmc,1))

fig = plt.figure()
bode(gc*g,k*g,dB=True,Hz=False) 
plt.legend(['KG_c(s)G(s)','KG(s)'])

# Malha fechada
mf1 = feedback(k*g,1);
mf2 = feedback(minreal(k*gc*g),1); 

# Resposta ao degrau
t = np.linspace(0,3,10000)
fig = plt.figure()
y1, t = step(mf1,t)
y2, t = step(mf2,t)
plt.plot(t,y1,'r--',t,y2,'k--')
plt.legend(['ajuste de ganho','sistema compensado'])

st1 = stepinfo(mf1,t)
st2 = stepinfo(mf2,t)
print('Ultrapassagem percentual com ajuste de ganho: ',st1['Overshoot'])
print('Ultrapassagem percentual com compensador avanço: ',st2['Overshoot'])
print('Tempo de pico com ajuste de ganho: ',st1['PeakTime'])
print('Tempo de pico percentual com compensador avanço: ',st2['PeakTime'])


# Resposta à rampa
N = 10000
t = np.linspace(0,100,N)
mf1 = feedback(k*g,1);
mf2 = feedback(k*gc*g,1); 
fig = plt.figure()
y1, t = step(mf1/s,t) # resposta ao degrau de 1/s*T é a resposta à rampa de T
y2, t = step(mf2/s,t)
plt.plot(t,t,'b',t,y1,'r--',t,y2,'k--')
plt.legend(['rampa','ajuste de ganho','sistema compensado'])

plt.figure()
bode(mf2)
