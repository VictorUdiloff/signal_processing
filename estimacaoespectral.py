import numpy as np
import matplotlib.pyplot as plt

th = np.random.rand(1)
t = np.linspace(0,999,1000)
x = np.cos(0.3*t+th)

L = 100
K = 20

v  = np.zeros(1000)
n = np.random.rand(1000) - 0.5
for i in range(0,999):
    v[i+1] = 0.95*v[i] + n[i]

y = x+v
In = np.zeros(1000)
U = 1


for k in range(0,K):
    w = np.zeros(1000)
    w[50*(k-1):50*k] = 1
    I = np.abs(np.fft.fft(y*w))
    In += (I**2)/(L*U)
In /= 20


plt.plot(np.abs(np.fft.fft(x+v)))
plt.plot(10*np.log(In))
plt.show()