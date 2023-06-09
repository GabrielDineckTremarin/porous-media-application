import matplotlib.pyplot as plt
#import plotly.plotly as py
import numpy as np
from scipy.signal import argrelextrema
import scipy.fftpack
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

Lc_factor = 1.0

Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)

Nx = int(Nx)
Ny = int(Ny)
y0 = int(y0)

xc,yc = np.loadtxt('transient/data/grid.dat', skiprows=0, unpack=True)

xc = np.array(xc[0: Nx])* Lc_factor
yc = np.array(yc[0::Nx])* Lc_factor


N = input('enter the final iteration:  ')
O = input('enter the step:  ')

M=102
N=N+102
NM=int(N-M)/O
u =     ['']*NM
v =     ['']*NM
speed = ['']*NM
Z =     ['']*NM
P =     ['']*NM
T =     ['']*NM
Xf =    ['']*NM
Yf =    ['']*NM
print 'data size:', NM
t = np.zeros(NM)# ['']*N
y = np.zeros(NM)# ['']*N

height = input('enter the y location of the probe:  ')

for k in range (0,NM*O,O):
        i = k/O
        print i
	u[i],v[i],Z[i],P[i],T[i] = np.loadtxt('transient/data/plot_field'+str(k+M)+'.dat', skiprows=0, unpack=True)
	t[i] = np.loadtxt('transient/data/time'+str(k+M)+'.dat', skiprows=0, unpack=True)
        v[i] = zip(*[iter(v[i])]*Nx)
        T[i] = zip(*[iter(T[i])]*Nx)
        P[i] = zip(*[iter(P[i])]*Nx)
        v[i] = np.array(v[i])
        P[i] = np.array(P[i])
        T[i] = np.array(P[i])
        t[i] = np.array(t[i])
        for k in range (1,len(yc),1):
            if yc[k-1] <= height:
                if yc[k] > height:
                    y[i] = v[i][k,0]

col_format = "{:<16}" * 2 + "\n"   # 2 left-justfied columns with 5 character width

with open("data/fft_probe.dat", 'w') as of:
    for x in zip(t, y):
        of.write(col_format.format(*x))
