import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
from strp2 import *

S = 4.65

rad = 1.0

Nx ,Ny , y0 = np.loadtxt('data/mesh.dat', skiprows=0, unpack=True)

x = np.loadtxt('data/x.dat', skiprows=0, unpack=True)
y = np.loadtxt('data/y.dat', skiprows=0, unpack=True)

u = np.loadtxt('data/u.dat', skiprows=0, unpack=True)
v = np.loadtxt('data/v.dat', skiprows=0, unpack=True)

P   = np.loadtxt('data/P.dat', skiprows=0, unpack=True)
T   = np.loadtxt('data/T.dat', skiprows=0, unpack=True)
Z   = np.loadtxt('data/Z.dat', skiprows=0, unpack=True)

#itc, error_u, error_v, error_p, mass, residual_p = np.loadtxt('data/error.dat', skiprows=0, unpack=True)

x_grid,y_grid= np.loadtxt('data/grid.dat', skiprows=0, unpack=True)
#x_grid_u,y_grid_u= np.loadtxt('data/grid_u.dat', skiprows=0, unpack=True)
#x_grid_v,y_grid_v= np.loadtxt('data/grid_v.dat', skiprows=0, unpack=True)
#x_grid2,y_grid2= np.loadtxt('data/grid_boundary.dat', skiprows=0, unpack=True)
#x_grid3,y_grid3= np.loadtxt('data/grid_boundary_side.dat', skiprows=0, unpack=True)
#x_grid4,y_grid4= np.loadtxt('data/grid_droplet.dat', skiprows=0, unpack=True)

Nx = int(Nx)
Ny = int(Ny)
y0 = int(y0)

x = np.array(x)
y = np.array(y)
u = np.array(u)
v = np.array(v)
T = np.array(T)

#transform to uniform mesh to be able to plot using pythons streamplot function
x_int = x
y_int = y
Nx = int(len(x))
Ny = int(len(y))

x_int = np.linspace(min(x), max(x), num=Nx,endpoint='True')
y_int = np.linspace(min(y), max(y), num=Ny,endpoint='True')

X,Y = np.meshgrid(x,y)
X2,Y2 = np.meshgrid(x_int,y_int)
interp_spline = RectBivariateSpline(y,x,v)
v_int = interp_spline(y_int,x_int)

X,Y = np.meshgrid(x,y)
X2,Y2 = np.meshgrid(x_int,y_int)
interp_spline = RectBivariateSpline(y,x,u)
u_int = interp_spline(y_int,x_int)

speed = np.sqrt(u*u + v*v)


############# GRID
plt.figure(figsize=(8, 8))
plt.plot( x_grid  ,y_grid  ,'k-', linewidth=1.0, label='Grid')
#plt.plot( x_grid_u  ,y_grid_u  ,'ks', linewidth=1.0, label='u Grid')
#plt.plot( x_grid_v  ,y_grid_v  ,'kP', linewidth=1.0, label='v Grid')
#plt.plot( x_grid2  ,y_grid2  ,'bo', linewidth=0.50, label='Boundary')
#plt.plot( x_grid3  ,y_grid3  ,'ro', linewidth=0.50, label='Boundary_side')
#plt.plot( x_grid4  ,y_grid4  ,'go', linewidth=0.50, label='Boundary_side')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
#plt.xlim(0,0.4)
#plt.ylim(-0.4,0.4)
plt.legend( loc= 'upper right')
plt.savefig('output/mesh.png', bbox_inches='tight')
plt.close()


########### PRESSURE
plt.figure(figsize=(8, 8))
plt.contourf(x,y,P,50, cmap=plt.cm.get_cmap('RdBu_r'))
#plt.plot( x_grid  ,y_grid  ,'k-', linewidth=0.2, label='Grid')
plt.colorbar()
plt.axis('scaled')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Pressure - P')
circle = plt.Circle((0, 0), rad, color='k', linewidth=2.0, fill=False)
plt.gca().add_patch(circle)
plt.savefig('output/pressure.png', bbox_inches='tight')
plt.show()


########### TEMPERATURE
plt.figure(figsize=(8, 8))
plt.contour(x,y,Z,colors=('k',),linestyles=('--',),linewidths=(2,), levels=[1/(S+1)])
plt.contourf(x,y,T,50, cmap=plt.cm.get_cmap('RdBu_r'))
#plt.plot( x_grid  ,y_grid  ,'k-', linewidth=0.2, label='Grid')
plt.colorbar()
plt.axis('scaled')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Temperature - T')
circle = plt.Circle((0, 0), rad, color='k', linewidth=2.0, fill=False)
plt.gca().add_patch(circle)
plt.savefig('output/temperature.png', bbox_inches='tight')
plt.show()

#plt.figure(figsize=(8, 8))
#plt.contourf(x,y,Z,50, cmap=plt.cm.get_cmap('RdBu_r'))
#plt.plot( x_grid  ,y_grid  ,'k-', linewidth=0.2, label='Grid')
#plt.colorbar()
#plt.axis('scaled')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.show()

############# STREAMLINES

X_stream = np.linspace(0,40,100)
Y_stream = X_stream*+9 -10 #y[len(y)/2]


plt.figure(figsize=(8, 8))
plt.streamplot(x_int, y_int, u_int, v_int,density=3, linewidth=0.8, color='k',arrowstyle='-')
#strmplt(x,y,u,v,X_stream,Y_stream)
#plt.plot( x_grid  ,y_grid  ,'k-', linewidth=0.2, label='Grid')
#plt.plot(X_stream,Y_stream, 'k-.',linewidth=0.4)
plt.contourf( x,y,speed,50, cmap=plt.cm.get_cmap('RdBu_r'))
plt.colorbar()
plt.title('Velocity Magnitude')
plt.axis('scaled')
plt.ylim(min(y),max(y))
plt.xlim(min(x),max(x))
circle = plt.Circle((0, 0), rad, color='k', linewidth=2.0, fill=False)
plt.gca().add_patch(circle)
plt.savefig('output/streamlines.png', bbox_inches='tight')
plt.show()



############ ERROR
#plt.figure(figsize=(8, 8))
#plt.plot( itc,error_u,'r-', linewidth=1.0, label='Error u')
#plt.plot( itc,error_v,'k-', linewidth=1.0, label='Error v')
#plt.plot( itc,error_p,'b-', linewidth=1.0, label='Error p')
#plt.plot( itc,residual_p,'g-', linewidth=1.0, label='Residual p')
#plt.yscale('log')
##plt.xscale('log')
#plt.xlabel('iterations')
#plt.ylabel('$Error$')
#plt.legend( loc= 'best')
#plt.savefig('output/error.png', bbox_inches='tight')
#plt.show()
#



#zero_line = np.linspace(0,10,10)
#plt.figure(figsize=(11, 7))
#plt.plot( u[:,10],y  ,'k.-', linewidth=1.0, label='Symm Axis')
#plt.plot(zero_line, 0*zero_line,'r-', linewidth=0.5)
#plt.xlabel('u')
#plt.ylabel('y')
#plt.xlim(min(y),max(y))
#plt.legend( loc= 'best')
##circle = plt.Circle((0, 0), 1, color='r', linewidth=2.0, fill=True)
##plt.gca().add_patch(circle)
#plt.savefig('output/u.png', bbox_inches='tight')
#plt.show()
