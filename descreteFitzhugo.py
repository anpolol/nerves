import numpy as np
import time
import matplotlib.pyplot as plt
import math
start_time = time.time()

Time = 500
dt = 1
amount = int(Time / dt)
tau = 30#14.9#14.94974
a = 0.25
b = 0.02
c = 0.2
gamma = 0.02
u1 = []
u2 = []
v1 = []
v2 = []
for i in range(int(tau/dt)):
    u1.append(0)
    u2.append(0)
    v1.append(0)
    v2.append(0)
u1.append(0.5)
u2.append(0)
v1.append(0)
v2.append(0)
print(u1[0],u1[int(tau/dt)-1],u1[int(tau/dt)])

def f1(y1,w1,y2,t):
    return -a*y1 + (a+1)*pow(y1,2)-pow(y1,3) - w1 + c*math.tanh(y2)
def f2(y1,w1):
    return b*y1 - gamma*w1

for i in range(int(tau/dt),amount+int(tau/dt)-1):
    t = (i-int(tau/dt))*dt
    k1 = f1(u1[i], v1[i], u2[i-int(tau/dt)],t) # u1
    m1 = f2(u1[i],v1[i]) #v1
    o1 = f1(u2[i], v2[i], u1[i-int(tau/dt)],t) #u2
    r1 = f2(u2[i],v2[i]) #v2

    k2 = f1(u1[i]+dt*k1/2, v1[i]+dt*m1/2, u2[i-int(tau/dt)]+dt*o1/2, t+dt/2)
    m2 = f2(u1[i]+dt*k1/2, v1[i]+dt*m1/2)
    o2 = f1(u2[i]+dt*o1/2,v2[i]+dt*r1/2,u1[i-int(tau/dt)]+dt*k1/2,   t+dt/2)
    r2 = f2(u2[i]+dt*o1/2, v2[i]+dt*r1/2)

    k3 = f1(u1[i]+dt*k2/2, v1[i]+dt*m2/2, u2[i-int(tau/dt)]+dt*o2/2, t+dt/2)
    m3 = f2(u1[i]+dt*k2/2, v1[i]+dt*m2/2)
    o3 = f1( u2[i]+dt*o2/2, v2[i]+dt*r2/2,u1[i-int(tau/dt)]+dt*k2/2, t+dt/2)
    r3 = f2(u2[i]+dt*o2/2, v2[i]+dt*r2/2)

    k4 = f1(u1[i]+dt*k3, v1[i]+dt*m3, u2[i-int(tau/dt)]+dt*o3, t+dt)
    m4 = f2(u1[i]+dt*k3, v1[i]+dt*m3)
    o4 = f1( u2[i]+dt*o3, v2[i]+dt*r3 ,u1[i-int(tau/dt)]+dt*k3,t+dt)
    r4 = f2(u2[i]+dt*o3, v2[i]+dt*r3)

    u1.append(u1[i] + (1/6) *(k1 + 2*k2 + 2*k3 + k4))
    v1.append(v1[i] + (1/6) *( m1 + 2*m2 + 2*m3 + m4))
    u2.append(u2[i] + (1/6) *(o1 + 2*o2 + 2*o3 + o4))
    v2.append(v2[i] + (1/6) *(r1 + 2*r2 + 2*r3 + r4))


j=[]

for i in range(amount-1):
    j.append((i)*dt)

fig, ax = plt.subplots()

w1= []
w2= []
for i in range(amount-1):
    w1.append(u1[i+int(tau/dt)])
    w2.append(u2[i + int(tau / dt)])
#ax.plot(j, w1, color="blue", label="u1(x)")      # функция y1(x), синий, надпись y(x)
ax.plot(j, w2, color="orange", label="u2(x)")      # функция y2(x), красный, надпись y'(x)

plt.show()                                      # показать рисунок
#fig.savefig('1.png')                            # сохранить в файл 1.png
end_time = time.time()
print("--- %s seconds ---" % (end_time - start_time))

