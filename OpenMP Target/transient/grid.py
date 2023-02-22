import matplotlib.pyplot as plt
import numpy as np
x_grid,y_grid= np.loadtxt('data/grid.dat', skiprows=0, unpack=True)
x_grid_u,y_grid_u= np.loadtxt('data/grid_u.dat', skiprows=0, unpack=True)
x_grid_v,y_grid_v= np.loadtxt('data/grid_v.dat', skiprows=0, unpack=True)
x_grid2,y_grid2= np.loadtxt('data/grid_boundary.dat', skiprows=0, unpack=True)
x_grid3,y_grid3= np.loadtxt('data/grid_boundary_side.dat', skiprows=0, unpack=True)
x_grid4,y_grid4= np.loadtxt('data/grid_droplet.dat', skiprows=0, unpack=True)


#

############# GRID
plt.figure(figsize=(8, 8))
plt.plot( x_grid  ,y_grid  ,'k-', linewidth=1.0, label='Grid')
#plt.plot( x_grid_u  ,y_grid_u  ,'ks', linewidth=1.0, label='u Grid')
#plt.plot( x_grid_v  ,y_grid_v  ,'kP', linewidth=1.0, label='v Grid')
plt.plot( x_grid2  ,y_grid2  ,'bo', linewidth=0.50, label='Boundary')
plt.plot( x_grid3  ,y_grid3  ,'ro', linewidth=0.50, label='Boundary_side')
plt.plot( x_grid4  ,y_grid4  ,'go', linewidth=0.50, label='Boundary_side')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
#plt.xlim(0,0.4)
#plt.ylim(-0.4,0.4)
plt.legend( loc= 'upper right')
plt.savefig('output/mesh.png', bbox_inches='tight')
plt.show()

