import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt
from scipy import misc
from scipy.integrate import odeint
import cPickle


def lamb(func):
    return sp.lambdify(c, func, modules='numpy')

""" Defining some functions """
def psi1(x):
	return np.square(x)

def psi2(x):
	return np.square(x)


""" Taking discrete values of f function """
fv=misc.imread('/home/martin/Documents/utfsm/CC3/3_variational/Circle-Missing-data-reduced-gray.png').astype(float)


""" Some constants and parameter values """
N=10 # Dimensions of image, assuming square form
eps=0.5 #shape Parameters
alpha=1.
beta=1.


""" Defining Lagrangian funtion L(x,y,u,ux,uy,uxx,uxy,uyy) and calculating it's derivatives """

#Variables of lagrangian
x1=sp.Symbol('x1')
x2=sp.Symbol('x2')
x3=sp.Symbol('x3')
x4=sp.Symbol('x4')
x5=sp.Symbol('x5')

L=(x1)**2 + alpha*psi1(x2**2+x3**2) + beta*psi2((x4+x5)**2)
Lu=sp.diff(L,x1)
Lux=sp.diff(L,x2)
Luy=sp.diff(L,x3)
Luxx=sp.diff(L,x4)
Luyy=sp.diff(L,x5)


""" Defining RBF and calculating it's derivatives """

#Mean variables (x,y)
x,y,xi,yi=sp.symbols('x y,xi,yi')


#RBF and it's derivatives
phi=sp.exp(-eps**2.*((x-xi)**2+(y-yi)**2))
phix=sp.diff(phi,x).simplify()
phiy=sp.diff(phi,y).simplify()
phixx=sp.diff(phix,x).simplify()
phiyy=sp.diff(phiy,y).simplify()


""" Computing collotacion points """
xv=np.linspace(0.,1.,N+1)[1::]-(1./(2.*N))
yv=np.linspace(0.,1.,N+1)[1::]-(1./(2.*N))


""" Creating symbols for c_i of RBFs """
carray=["c"+str(i) for i in range(N**2)]
c=sp.symbols(carray)


""" Computing u, ux, uy, uxx, uyy with N points"""
u=0.
ux=0.
uy=0.
uxx=0.
uyy=0.
points=[]


for i in range(N):
	tmpx=xv[i]
	for j in range(N):
		tmpy=yv[j]
		u+=c[10*i+j]*phi.subs(xi,tmpx).subs(yi,tmpy)
		ux+=c[10*i+j]*phix.subs(xi,tmpx).subs(yi,tmpy)
		uy+=c[10*i+j]*phiy.subs(xi,tmpx).subs(yi,tmpy)
		uxx+=c[10*i+j]*phixx.subs(xi,tmpx).subs(yi,tmpy)
		uyy+=c[10*i+j]*phiyy.subs(xi,tmpx).subs(yi,tmpy)
		points.append((tmpx,tmpy))


""" Building the EL equation """
Lu=Lu.subs(x1,u)
Lux=Lux.subs(x2,ux).subs(x3,uy)
Luy=Luy.subs(x2,ux).subs(x3,uy)
Luxx=Luxx.subs(x4,uxx).subs(x5,uyy)
Luyy=Luyy.subs(x4,uxx).subs(x5,uyy)
EL=Lu-sp.diff(Lux,x)-sp.diff(Luy,y)+sp.diff(Luxx,x,2)+sp.diff(Luyy,y,2)


""" Evaluating EL equation at domain points [0,1]x[0,1] (falta restar 2f) """
evEL=[EL.subs(x,points[i][0]).subs(y,points[i][1]) for i in range(N**2)]


""" Coneverting into NumPy evaluated functions """
F=map(lamb,evEL)

# #Now we have EL equation evaluated on the domain points (only dependent on ci coef)


""" Serializing it and storing with Tickle... """
f=file("sample.pickle","w")
pickle.dump(F,f) # Dump F object into f file
f.close()









