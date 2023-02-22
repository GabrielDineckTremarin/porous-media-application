import matplotlib.pyplot as plt
import numpy as np


itc, error_u, error_v, error_p = np.loadtxt('data/error.dat', skiprows=0, unpack=True)

########### ERROR
plt.figure(figsize=(8, 8))
plt.plot( itc,error_u,'r-', linewidth=1.0, label='Error u')
plt.plot( itc,error_v,'k-', linewidth=1.0, label='Error v')
plt.plot( itc,error_p,'b-', linewidth=1.0, label='Error p')
plt.yscale('log')
#plt.xscale('log')
plt.xlabel('iterations')
plt.ylabel('$Error$')
plt.legend( loc= 'best')
plt.savefig('output/error.png', bbox_inches='tight')
plt.show()

