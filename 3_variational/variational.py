import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt
from scipy import misc
from scipy.integrate import odeint
import pickle
import cloudpickle
import time

fobj=open("sample.pickle","rb")
F=pickle.load(fobj)
fobj.close()
