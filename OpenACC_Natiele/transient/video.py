import matplotlib.pyplot as plt
import numpy as np
#from scipy import stats
#import pylab
import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#from strp2 import *
#from spf import *

Lc_factor = 1.0
S = 19.0


Xstream = np.random.rand(200)*20
Ystream = np.random.rand(200)*20-10

Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)

Nx = int(Nx)
Ny = int(Ny)
y0 = int(y0)

Mi = 1 #input('enter initial iteration:')
Mi = int(Mi)
M=Mi+200
N=M+300
O=1
#x = ['']*N
#y = ['']*N
u = ['']*N
v = ['']*N
speed = ['']*N
Z = ['']*N
P = ['']*N
T = ['']*N
Xf = ['']*N
Yf = ['']*N
t = ['']*N

x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)

k = input('enter number:1 for temperature, 2 for pressure, 3 for velocity =  ')
k = int(k)
cont = 100


for i in range (M,N,O):
        print (i-100)
        cont = cont + 1
        u[i] = np.loadtxt('transient/data/u'+str(i)+'.dat', skiprows=0, unpack=True)
        v[i] = np.loadtxt('transient/data/v'+str(i)+'.dat', skiprows=0, unpack=True)
        P[i] = np.loadtxt('transient/data/P'+str(i)+'.dat', skiprows=0, unpack=True)
        T[i] = np.loadtxt('transient/data/T'+str(i)+'.dat', skiprows=0, unpack=True)
        Z[i] = np.loadtxt('transient/data/Z'+str(i)+'.dat', skiprows=0, unpack=True)
        t[i] = np.loadtxt('transient/data/time'+str(i)+'.dat', skiprows=0, unpack=True)
        u[i] = np.array(u[i])
        v[i] = np.array(v[i])
        T[i] = np.array(T[i])
        speed[i] = np.sqrt(u[i]*u[i] + v[i]*v[i])
        #################
        #k=1, temperature
        #k=2, pressure
        #k=3, velocity
        #################

        if k == 1:
        ############ TEMPERATURE
                plt.figure(figsize=(8, 8))
                CS = plt.contourf(x,y,T[i],50,cmap='RdBu_r') 
                plt.contour(x,y,Z[i],colors=('k',),linestyles=('--',),linewidths=(2,), levels=[1.0/(S+1)])
                cbar = plt.colorbar(CS) 
                plt.axis('scaled')
                plt.savefig('transient/t'+str(cont)+'.png', bbox_inches='tight')
                plt.close()
        
        if k == 2:
        ############ PRESSURE
                plt.figure(figsize=(8, 8))
                CS = plt.contourf(x,y,P[i],50,cmap='RdBu_r') 
                cbar = plt.colorbar(CS) 
                plt.axis('scaled')
                plt.savefig('transient/t'+str(cont)+'.png', bbox_inches='tight')
                plt.close()
        
        if k == 3:
        ############# VELOCITY 
                plt.figure(figsize=(8, 8))
                CS = plt.contourf(x,y,speed[i],50,cmap='RdBu_r') 
                cbar = plt.colorbar(CS) 
                plt.axis('scaled')
                plt.savefig('transient/t'+str(cont)+'.png', bbox_inches='tight')
                plt.close()




#        A = plt.contour(x,y,Z[i], levels=[1]) #"Plota" o contorno para Z=1
#        V = A.collections[0].get_paths()[0].vertices #Faz um lance que eu nao entendo, copiei do Stack Overflow
#        Xf[i] = V[:,0] #Separa as componentes
#        Yf[i] = V[:,1]
#        print np.max(Yf[i]), np.max(Xf[i]), i

#        plt.figure(figsize=(11, 7))
#        plt.plot( u[i][:,len(x)/2],y  ,'k.-', linewidth=1.0, label='t')
#        plt.xlabel('u')
#        plt.ylabel('y')
#        plt.xlim(-1,1)
#        plt.legend( loc= 'best')
##################
###        circle = plt.Circle((0, 0), 1, color='k', linewidth=2.0, fill=False)
#        fig0, ax0 = plt.subplots(figsize=(12, 12))
##        ax0.add_artist(circle)
##        strmplt(x,y,u[j],v[j],Xstream,Ystream,[-10,10,-5,20])
##        strmplt(x,y,u[i],v[i],Xstream,Ystream,[0,20,-10,10])
##        plt.contourf(x,y,np.transpose(curlnum),20,cmap='jet'), plt.colorbar()
#        plt.contour(x,y,Z[i],colors=('k',),linestyles=('--',),linewidths=(2,), levels=[1])
#        T[i] = 2 * 300.0 * T[i]
#        #CS2 = plt.contourf(-x,y,T[i],50,cmap='jet') 
#        CS = plt.contourf(x,y,speed[i],50,cmap='RdBu_r') 
##        plt.streamplot(x, y, u[i], v[i],density=2, linewidth=0.8, color='k',arrowstyle='-')
###        for c in CS.collections:
###            c.set_edgecolor("face")
###        #for c in CS2.collections:
###        #    c.set_edgecolor("face")
#        cbar = plt.colorbar(CS) 
###        cbar = plt.colorbar(CS, orientation="vertical",ticks=[300,650,1000,1350,1700, 2050])
###        cbar.ax.set_ylabel('$\hat{T}$ [K]')
###        plt.clim(np.min(T[i]), np.max(T[i]))
####        cbar.set_ticklabels([300,650,1000,1350,1700, 2100])
###        cbar.set_label(r'$\hat{T}$ [K]', x=-0.12, labelpad = 4)
#
#
#        plt.axis('scaled')
#        #SPF()
#        #plt.ylim(-1,6)
#        #plt.xlim(-2,2) 

#v = np.array(v)


#f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharex='col',figsize=(5.7, 4))
#i = 220
#CS = ax1.contourf(x,y,T[i],50,cmap='jet') 
#CS2 = ax1.contourf(-x,y,T[i],50,cmap='jet') 
#for c in CS.collections:
#    c.set_edgecolor("face")
#for c in CS2.collections:
#    c.set_edgecolor("face")
#ax1.axis('scaled')
#ax1.set_title('1.57 s')
#ax1.set_ylim(-1,8)
#ax1.set_xlim(-1.5,1.5) 
#i = 222
#CS = ax2.contourf(x,y,T[i],50,cmap='jet') 
#CS2 = ax2.contourf(-x,y,T[i],50,cmap='jet') 
#for c in CS.collections:
#    c.set_edgecolor("face")
#for c in CS2.collections:
#    c.set_edgecolor("face")
#ax2.axis('scaled')
#ax2.set_title('1.59 s')
#ax2.set_ylim(-1,8)
#ax2.set_xlim(-1.5,1.5) 
#i = 224
#CS = ax3.contourf(x,y,T[i],50,cmap='jet') 
#CS2 = ax3.contourf(-x,y,T[i],50,cmap='jet') 
#for c in CS.collections:
#    c.set_edgecolor("face")
#for c in CS2.collections:
#    c.set_edgecolor("face")
#ax3.axis('scaled')
#ax3.set_title('1.62 s')
#ax3.set_ylim(-1,8)
#ax3.set_xlim(-1.5,1.5) 
#cbar = plt.colorbar(CS, orientation="vertical",ticks=[300,650,1000,1350,1700, 2050])
#cbar.ax.set_ylabel('$\hat{T}$ [K]')
##ax3.clim(np.min(T[i]), np.max(T[i]))
##        cbar.set_ticklabels([300,650,1000,1350,1700, 2100])
#cbar.set_label(r'$\hat{T}$ [K]', x=-0.12, labelpad = 4)
##SPF()
#plt.savefig('transient/multiplot.eps', bbox_inches='tight')
#plt.show()
