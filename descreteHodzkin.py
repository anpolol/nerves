import numpy as np
import time
import matplotlib.pyplot as plt
import math
import random
start_time = time.time()
import scipy

Time =1000
dt = 1
tau =0
c = -0.01
amount = int(Time / dt)
gNa = 120
phiNa = -115
gK = 36
phiK = 12
gl = 0.3
phil = -10.613
a = 1E-3
R = 10E8
L=1 #электроемкость мембраны


u1 = []
u2 = []
n1 = []
n2 = []
m1 = []
m2 = []
h1 = []
h2 = []
l1 = 0.59
l2 = 0.32
tau=194
#l1 = random.random()
#l2 = random.random()
print(l1,l2)
for i in range(int(tau/dt)+1):
    u1.append(0)
    u2.append(0)

    h1.append(l1)
    h2.append(l1)
    n1.append(l2)
    n2.append(l2)
    m1.append(0.05)
    m2.append(0.05)
u1[int(tau/dt)]=0
u2[int(tau/dt)]=0.5

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
    return (L/gK)*alpha_n(v*phiK)*(1-tt)- (L/gK)*betta_n(v*phiK)*tt
def m(tt,v):
    return (L/gK)*alpha_m(v*phiK)*(1-tt) - (L/gK)*betta_m(v*phiK)*tt
def h(tt,v):
    return (L/gK)*alpha_h(v*phiK)*(1-tt) - (L/gK)*betta_h(v*phiK)*tt


def f(y1,y2,n_n,h_h,m_m):
    return (-(y1 - 1)*(pow(n_n,4)) - (gNa/gK)*(y1 - phiNa/phiK)*pow(m_m,3)*(h_h) - (gl/gK)*(y1 - phil/phiK)) #+ (c)*math.tanh(y2))

for i in range(int(tau/dt),amount+int(tau/dt)-1):
   # print(i,n1[i],u1[i],n2[i],u2[i])
    t = (i-int(tau/dt))*dt
    k1 = f(u1[i], u2[i - int(tau / dt)], n1[i], h1[i], m1[i])  # u1
    q1 = n(n1[i], u1[i])
    o1 = m(m1[i], u1[i])
    r1 = h(h1[i], u1[i])
    w1 = f(u2[i], u1[i - int(tau / dt)], n2[i], h2[i], m2[i])
    e1 = n(n2[i], u2[i])
    t1 = m(m2[i], u2[i])
    i1 = h(h2[i], u2[i])

    k2 = f(u1[i] + dt * k1 / 2, u2[i - int(tau / dt)] + dt * w1 / 2, n1[i] + q1 * dt / 2, h1[i] + dt * r1 / 2,m1[i] + dt * o1 / 2)
    q2 = n(n1[i] + q1 * dt / 2, u1[i] + dt * k1 / 2)
    o2 = m(m1[i] + dt * o1 / 2, u1[i] + dt * k1 / 2)
    r2 = h(h1[i] + dt * r1 / 2, u1[i] + dt * k1 / 2)
    w2 = f(u2[i] + w1 * dt / 2, u1[i - int(tau / dt)] + dt * k1 / 2, n2[i] + e1 * dt / 2, h2[i] + i1 * dt / 2,m2[i] + t1 * dt / 2)
    e2 = n(n2[i] + e1 * dt / 2, u2[i] + dt * w1 / 2)
    t2 = m(m2[i] + t1 * dt / 2, u2[i] + dt * w1 / 2)
    i2 = h(h2[i] + i1 * dt / 2, u2[i] + dt * w1 / 2)

    k3 = f(u1[i] + dt * k2 / 2, u2[i - int(tau / dt)] + dt * w2 / 2, n1[i] + q2 * dt / 2, h1[i] + dt * r2 / 2, m1[i] + dt * o2 / 2)
    q3 = n(n1[i] + q2 * dt / 2, u1[i] + dt * k2 / 2)
    o3 = m(m1[i] + dt * o2 / 2, u1[i] + dt * k2 / 2)
    r3 = h(h1[i] + dt * r2 / 2, u1[i] + dt * k2 / 2)
    w3 = f(u2[i] + w2 * dt / 2, u1[i - int(tau / dt)] + dt * k2 / 2, n2[i] + e2 * dt / 2, h2[i] + i2 * dt / 2,m2[i] + t2 * dt / 2)
    e3 = n(n2[i] + e2 * dt / 2, u2[i] + dt * w2 / 2)
    t3 = m(m2[i] + t2 * dt / 2, u2[i] + dt * w2 / 2)
    i3 = h(h2[i] + i2 * dt / 2, u2[i] + dt * w2 / 2)

    k4 = f(u1[i] + dt * k3, u2[i - int(tau / dt)] + dt * w3, n1[i] + q3 * dt, h1[i] + dt * r3, m1[i] + dt * o3)
    q4 = n(n1[i] + q3 * dt, u1[i] + dt * k3)
    o4 = m(m1[i] + dt * o3, u1[i] + dt * k3)
    r4 = h(h1[i] + dt * r3, u1[i] + dt * k3)
    w4 = f(u2[i] + w3 * dt, u1[i - int(tau / dt)] + dt * k3, n2[i] + e3 * dt, h2[i] + i3 * dt, m2[i] + t3 * dt)
    e4 = n(n2[i] + e3 * dt, u2[i] + dt * w3)
    t4 = m(m2[i] + t3 * dt, u2[i] + dt * w3)
    i4 = h(h2[i] + i3 * dt, u2[i] + dt * w3)

    u1.append(u1[i] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4))
    n1.append(n1[i] + (dt / 6) * (q1 + 2 * q2 + 2 * q3 + q4))
    m1.append(m1[i] + (dt / 6) * (o1 + 2 * o2 + 2 * o3 + o4))
    h1.append(h1[i] + (dt / 6) * (r1 + 2 * r2 + 2 * r3 + r4))
    u2.append(u2[i] + (dt / 6) * (w1 + 2 * w2 + 2 * w3 + w4))
    n2.append(n2[i] + (dt / 6) * (e1 + 2 * e2 + 2 * e3 + e4))
    m2.append(m2[i] + (dt / 6) * (t1 + 2 * t2 + 2 * t3 + t4))
    h2.append(h2[i] + (dt / 6) * (i1 + 2 * i2 + 2 * i3 + i4))

    #print(m1[i],h1[i],n1[i])



j=[]

for i in range(amount-1):
    j.append((i)*dt)

fig, ax = plt.subplots()

b1= []
b2= []
b3=[]
for i in range(amount-1):
    b1.append(u1[i+int(tau/dt)])
    b2.append(u2[i + int(tau / dt)])
    b3.append(h1[i + int(tau / dt)])
ax.plot(j, b3, color="blue", label="m(t)")
#ax.plot(j, b1, color="orange", label="u1(t)") # функция y1(x), синий, надпись y(x)
#ax.plot(j, b2, color="black", label="u2(t)")      # функция y2(x), красный, надпись y'(x)
plt.legend()

plt.show()                                      # показать рисунок
#fig.savefig('1.png')                            # сохранить в файл 1.png
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))
