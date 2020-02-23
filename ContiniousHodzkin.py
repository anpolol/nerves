import numpy as np
import time
import matplotlib.pyplot as plt
import math
#from scipy import signal
start_time = time.time()

mu=10
sigma=1

Time = 10
l = 50
dx = 0.1
dt = dx*dx/6
K = int(Time / dt)+1
I =int(round(l/dx))+1
L = 1
gNa = 120
phiNa = -115
gK = 36
phiK = 12
gl = 0.3
phil = -10.613
a = 238e-6
R = 35.4

g_m = gNa/gK
V_Na = phiNa/phiK
g = gl/gK
V_l=phil/phiK


V1 = np.zeros((I,K))
V2 = np.zeros((I,K))
V3 = np.zeros((I,K))
n_array_1 = np.zeros((I,K))
m_array_1 = np.zeros((I,K))
h_array_1 = np.zeros((I,K))
n_array_2 = np.zeros((I,K))
m_array_2 = np.zeros((I,K))
h_array_2 = np.zeros((I,K))
n_array_3 = np.zeros((I,K))
m_array_3 = np.zeros((I,K))
h_array_3 = np.zeros((I,K))

for i in range(I):
    n_array_1[i][0] = 0.5
    h_array_1 [i][0] = 0.5
    n_array_2[i][0] = 0.5
    h_array_2[i][0] = 0.5
    n_array_3[i][0] = 0.5
    h_array_3[i][0] = 0.5


def alpha_h(V):
    return 0.07*math.exp(V/20)
def alpha_m(V):
    return 0.1*(V+25)/(math.exp((V+25)/10)-1)
def alpha_n(V):
    return 0.01*(V+10)/(math.exp((V+10)/10)-1)
def betta_h(V):
    return 1/(math.exp((V+30)/10)+1)
def betta_m(V):
    return 4*math.exp(V/18)
def betta_n(V):
    return 0.125*math.exp(V/80)


def n(tt,v):
    return (L/gK)*alpha_n(v*phiK)*(1-tt) - (L/gK)*betta_n(v*phiK)*tt
def m(tt,v):
    return (L/gK)*alpha_m(v*phiK)*(1-tt) - (L/gK)*betta_m(v*phiK)*tt
def h(tt,v):
    return (L/gK)*alpha_h(v*phiK)*(1-tt) - (L/gK)*betta_h(v*phiK)*tt
def f(y1,n_n,h_h,m_m):
    return (y1 - 1)*(pow(n_n,4)) + (g_m)*(y1 - V_Na)*pow(m_m,3)*(h_h) + (g)*(y1 - V_l)

##ачальные условия
for j in range(I):
    V1[j][0] = math.exp(-((j) * dx - mu) * ((j) * dx - mu) / (2 * sigma)) * (
                1 / (sigma * sigma * math.sqrt(2 * math.pi)))  # распределение по всему волокну в момент времени 0
    V2[j][0] = math.exp(-((j) * dx - mu) * ((j) * dx - mu) / (2 * sigma)) * (
                1 / (sigma * sigma * math.sqrt(2 * math.pi)))  # распределение по всему волокну в момент времени 0
    V3[j][0] = math.exp(-((j) * dx - mu) * ((j) * dx - mu) / (2 * sigma)) * (
                1 / (sigma * sigma * math.sqrt(2 * math.pi)))  # распределение по всему волокну в момент времени 0
#signal.unit_impulse(I, 'mid')[i]
for k in range(K):
    V1[0][k] =  math.exp(-(0 - mu)*(0- mu)/(2*sigma)) *(1/(sigma*sigma*math.sqrt(2*math.pi)))
    V1[I-1][k] =math.exp(-((l) - mu)*((l) - mu)/(2*sigma)) *(1/(sigma*sigma*math.sqrt(2*math.pi)))
    V2[0][k] = math.exp(-(0 - mu) * (0 - mu) / (2 * sigma)) * (1 / (sigma * sigma * math.sqrt(2 * math.pi)))
    V2[I - 1][k] = math.exp(-((l) - mu) * ((l) - mu) / (2 * sigma)) * (1 / (sigma * sigma * math.sqrt(2 * math.pi)))
    V3[0][k] =  math.exp(-(0 - mu)*(0- mu)/(2*sigma)) *(1/(sigma*sigma*math.sqrt(2*math.pi)))
    V3[I-1][k] =math.exp(-((l) - mu)*((l) - mu)/(2*sigma)) *(1/(sigma*sigma*math.sqrt(2*math.pi)))

for j in range(1,K-1):
    for i in range(2,I-1):
        n1 = n(n_array_1[i][j-1], V1[i][j-1])  # u1
        m1 = m(m_array_1[i][j-1], V1[i][j-1])  # v1
        h1 = h(h_array_1[i][j-1], V1[i][j-1])  # u2

        n2 = n(n_array_1[i][j-1] + (dt / 2) * n1, V1[i][j-1])
        m2 = m(m_array_1[i][j-1] + (dt / 2) * m1, V1[i][j-1])
        h2 = h(h_array_1[i][j-1] + (dt / 2) * h1, V1[i][j-1])

        n3 = n(n_array_1[i][j-1] + (dt / 2) * n2, V1[i][j-1])
        m3 = m(m_array_1[i][j-1] + (dt / 2) * m2, V1[i][j-1])
        h3 = h(h_array_1[i][j-1] + (dt / 2) * h2, V1[i][j-1])

        n4 = n(n_array_1[i][j-1] + (dt) * n3, V1[i][j-1])
        m4 = m(m_array_1[i][j-1] + (dt) * m3, V1[i][j-1])
        h4 = h(h_array_1[i][j-1] + (dt) * h3, V1[i][j-1])

        n_array_1[i][j] = (n_array_1[i][j-1] + (1 / 6) * (n1 + 2 * n2 + 2 * n3 + n4))
        m_array_1[i][j] = (m_array_1[i][j-1] + (1 / 6) * (m1 + 2 * m2 + 2 * m3 + m4))
        h_array_1[i][j] = (h_array_1[i][j-1] + (1 / 6) * (h1 + 2 * h2 + 2 * h3 + h4))
        V1[i][j+1] = (dt/(dx*dx))* V1[i-1][j] + (1-2*dt/(dx*dx))*V1[i][j] -dt*(f(V1[i][j],n_array_1[i][j],h_array_1[i][j],m_array_1[i][j]))+ (dt/(dx*dx))*V1[i+1][j]

j = []
y1 = []
y2 = []
y3 =[]
y4=[]
y5=[]
y6=[]
for i in range(I):
    j.append(i*dx)
    y1.append((V1[i][0]))
    #y2.append((V1[i][int(K/5)]))
    #y3.append((V1[i][int(2*K/5)]))
    #y4.append((V1[i][int(3*K/5)]))
    #y5.append((V1[i][int(4*K/5)]))
    y6.append((V1[i][K-1]))

fig, ax = plt.subplots()                        # будет 1 график, на нем:
ax.plot(j, y1, color="blue", label="u1(x)")      # функция y1(x), синий, надпись y(x)3ax.plot(j, y2, color="orange", label="u1(x)")
ax.plot(j, y3, color="black", label="u1(x)")
#ax.plot(j, y4, color="yellow", label="u1(x)")
#ax.plot(j, y5, color="red", label="u1(x)")
#ax.plot(j, y6, color="brown", label="u1(x)")   # функция y2(x), красный, надпись y'(x)
plt.show()                                      # показать рисунок
#fig.savefig('1.png')                            # сохранить в файл 1.png
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))
