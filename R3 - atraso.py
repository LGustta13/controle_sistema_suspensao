# -*- coding: utf-8 -*-

%matplotlib
import matplotlib.pyplot as plt
import numpy as np
from control.matlab import *

## Relatório 3 - Simulações

## última atualização
# Achei aquele gráfico zuado UP0 e Ts 2.3s, onde já terminei uma simulação e agora tentar manipular a fase do compensador para diminuir a margem

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
g = (1/s)*(kt*(Cs*s + ks))/(m*M*(s**4) + (m+M)*Cs*(s**3) + (kt*M+(m+M)*ks)*(s**2) + kt*Cs*s + kt*ks)

# Requisitos de desempenho
UP = 0.02 # ultrapassagem percentual
e = 0.01     # constante de erro de velocidade
Ts = 0.7  # tempo de assentamento

###############################################
# AJUSTE DE GANHO E CÁLCULO DE MARGEM DE FASE
###############################################


fig2 = plt.figure()
bode(g,cs*g,dB=True,Hz=False) 
###########################################

gmc,pmc,wgmc,wpmc = margin(g) # verificação das margens 

# Margem de fase requerida
zeta = -np.log(UP)/np.sqrt(np.pi**2+np.log(UP)**2) # fator de amortecimento requerido
print('zeta: ',np.round(zeta,2))
PM = np.rad2deg(np.arctan((2*zeta)/(np.sqrt(-2*(zeta**2)+np.sqrt(1+4*(zeta**4))))))  # margem de fase requerida
print('Margem de fase requerida: ',np.round(PM,2))
PHI_seg = 12 # valor de segurança => entre 5 e 12 graus
PHI = -180 + (PM+PHI_seg) # fase na nova freq. de margem de fase 
print('Fase na nova frequência de margem de fase ',np.round(PHI))

# Ajuste de ganho para garantir a margem de fase requerida
fig1 = plt.figure()
bode(g,dB=True,Hz=False)
wpm = 5.9 # freq. de margem de fase retirada do diagrama de Bode em PHI_PM
AdB = 11; # ganho a ser adicionado na frequencia wpm
k1 = np.power(10,AdB/20) # valor absoluto do ganho
print('Valor de K para requisito de ultrapassagem percentual: ',np.round(k1,1))

# verificação das margens 
bode(g, g*k1,dB=True,Hz=False)

###############################################
# PROJETO DO COMPENSADOR EM ATRASO
###############################################

# Nova frequencia de margem de fase requerida => wpm_nova
PHI_seg = 12 # valor de segurança => entre 5 e 12 graus
PHI = -180 + (PM+PHI_seg) # fase na nova freq. de margem de fase 
print('Fase na nova frequência de margem de fase ',np.round(PHI))

fig = plt.figure()
bode(g,'r-',g,'b--',dB=True,Hz=False) # procurar wpm_nova => fase PHI em k*g
plt.legend(['Ajuste de ganho para overshoot','sem compensação'])
wpm_nova = 2.7 # frequencia na fase phi
print('Nova frequência de margem de fase ',wpm_nova)


# Ganho do compensador para garantir magnitude 0dB na frequencia wpm_nova
kpm_dB = 8 # magnitude (com valor negativo) em k*G(s) na frequencia wpm_nova

# Zero do compensador (zc)
zc = -wpm_nova/10 # aproximadamente uma década abaixo da freq. wpm_nova

# Polo do compensador (pc)
# f(w) = -20log(w) + B => eq. da reta 
# f(w) = kpm_dB => B = f(w)+20log(w)
# f(w) = 0 => w = 10^(B/20)
B = kpm_dB +20*np.log10(wpm_nova/10)
wpc = np.power(10,B/20)
pc = -wpc

# Funcao de transferencia compensador Gc = Kc(s-zc)/(s-pc)
kc = pc/zc # ganho para impor compensador com ganho estatico unitario
gc = kc*((s-zc))/(s-pc) # funcao de transferencia do compensador


###############################################
# Verificação de resultados
###############################################

gmc,pmc,wgmc,wpmc = margin(gc*g) # verificação das margens 
print('Margem de fase do sistema compensado: ',round(pmc,1))

fig = plt.figure()
bode(gc*g,g,dB=True,Hz=False) 
plt.legend(['G_c(s)G(s)','G(s)'])

# Malha fechada
mf1 = feedback(g,1);
mf2 = feedback(minreal(gc*g),1); 

# Resposta ao degrau
t = np.linspace(0,3,10000)
st1 = stepinfo(mf1,t)
st2 = stepinfo(mf2,t)
print('Ultrapassagem percentual com ajuste de ganho inicial: ',st1['Overshoot'])
print('Ultrapassagem percentual com compensador atraso: ',st2['Overshoot'])

fig = plt.figure()
y1, t = step(mf1,t)
y2, t = step(mf2,t)
plt.plot(t,y1,'r--',t,y2,'k--')
plt.legend(['sem compensador','sistema compensado'])

# Resposta à rampa
N = 10000
t = np.linspace(0,100,N)
mf1 = feedback(k1*g,1);
mf2 = feedback(k*gc*g,1); 
fig = plt.figure()
y1, t = step(mf1/s,t) # resposta ao degrau de 1/s*T é a resposta à rampa de T
y2, t = step(mf2/s,t)
plt.plot(t,t,'b',t,y1,'r--',t,y2,'k--')
plt.legend(['rampa','ajuste de ganho','sistema compensado'])
e1 = y1[N-1]-t[N-1] # erro_estacionario = y(inf)-r(inf), sendo r = t
e2 = y2[N-1]-t[N-1] # quanto maior N, mais o erro se aproxima do erro estacionário
print('Erro à rampa do sistema não compensado: ',round(e1,4))
print('Erro à rampa do sistema compensado: ',round(e2,4))