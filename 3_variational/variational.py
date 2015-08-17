import numpy as np 
import matplotlib.pyplot as plt
from scipy import misc
from scipy.integrate import odeint
import pickle
import time

""" Some constants """
N=4 #image dimension


""" Load "EL function" """
fobj=open("sample.pickle","rb")
F=pickle.load(fobj)
fobj.close()



# solve the system dc/dt = f(c, t)
def f(c,t):
	retval=[F[i](*c) for i in range(N**2)]
	return retval


# initial condition vector
c0=np.random.rand(N**2)

# time grid
t  = np.linspace(0, 5., 1000)

# solve that shit
soln = odeint(f, c0, t)