import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as si

raio = 1/10.5

def strmplt(x,y,u,v,X0,Y0):
	Fu = si.RectBivariateSpline(x,y,u.T) # Funcoes interpoladoras
	Fv = si.RectBivariateSpline(x,y,v.T)

	for k in range(2):
		if k==0: # Linhas de corrente a jusante
			t = np.linspace(0,400,5000) # Vetor tempo
		if k==1: # Linhas de corrente a montante
			t = -np.linspace(0,400,5000) # Vetor tempo

		X = np.zeros([len(X0),len(t)]) # Matrizes de istrimilainis
		Y = np.zeros([len(X0),len(t)])

		X[:,0] = X0 # Alimenta as seeds
		Y[:,0] = Y0

		for i in range(len(t)-1): # Integracao no tempo
			# Metodo de Heun, eu acho
			X[:,i+1] = X[:,i] + Fu.ev(X[:,i],Y[:,i])*(t[i+1]-t[i]) # Preditor
			Y[:,i+1] = Y[:,i] + Fv.ev(X[:,i],Y[:,i])*(t[i+1]-t[i])
			X[:,i+1] = X[:,i] + .5*(Fu.ev(X[:,i+1],Y[:,i+1])+Fu.ev(X[:,i],Y[:,i]))*(t[i+1]-t[i]) # Corretor
			Y[:,i+1] = Y[:,i] + .5*(Fv.ev(X[:,i+1],Y[:,i+1])+Fv.ev(X[:,i],Y[:,i]))*(t[i+1]-t[i])

		for i in range(len(X0)): # Plota as linhas de corrente
		#	X[i,(X[i,:]**2+Y[i,:]**2)<raio] = np.nan # Interior do cilindro
			plt.plot(X[i,:],Y[i,:], 'k-', linewidth=0.4) # Plota i-esima linha de corrente


def strmplt_inv(x,y,u,v,X0,Y0):
	Fu = si.RectBivariateSpline(x,y,u.T) # Funcoes interpoladoras
	Fv = si.RectBivariateSpline(x,y,v.T)

	for k in range(2):
		if k==0: # Linhas de corrente a jusante
			t = np.linspace(0,400,5000) # Vetor tempo
		if k==1: # Linhas de corrente a montante
			t = -np.linspace(0,400,5000) # Vetor tempo

		X = np.zeros([len(X0),len(t)]) # Matrizes de istrimilainis
		Y = np.zeros([len(X0),len(t)])

		X[:,0] = X0 # Alimenta as seeds
		Y[:,0] = Y0

		for i in range(len(t)-1): # Integracao no tempo
			# Metodo de Heun, eu acho
			X[:,i+1] = X[:,i] + Fu.ev(X[:,i],Y[:,i])*(t[i+1]-t[i]) # Preditor
			Y[:,i+1] = Y[:,i] + Fv.ev(X[:,i],Y[:,i])*(t[i+1]-t[i])
			X[:,i+1] = X[:,i] + .5*(Fu.ev(X[:,i+1],Y[:,i+1])+Fu.ev(X[:,i],Y[:,i]))*(t[i+1]-t[i]) # Corretor
			Y[:,i+1] = Y[:,i] + .5*(Fv.ev(X[:,i+1],Y[:,i+1])+Fv.ev(X[:,i],Y[:,i]))*(t[i+1]-t[i])

		for i in range(len(X0)): # Plota as linhas de corrente
			X[i,(X[i,:]**2+Y[i,:]**2)<raio] = np.nan # Interior do cilindro
			plt.plot(-X[i,:],Y[i,:], 'k-', linewidth=0.4) # Plota i-esima linha de corrente


#	plt.plot(np.sin(np.linspace(0,2*np.pi)),np.cos(np.linspace(0,2*np.pi)),'k') # Plota o cilindro
