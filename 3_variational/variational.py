import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def lamb(func):
    return sp.lambdify((x, y), func, modules='numpy')

""" Defining some functions """
def psi1(x):
	return np.square(x)

def psi2(x):
	return np.square(x)

def f(x,y):
	return x+y



""" Some important constants """
eps=0.5 #shape parameter
Nx=10 #dimensions x on image
Ny=10 #dimensions y on image
N=Nx*Ny


""" Defining Lagrangian funtion F(x,y,u,ux,uy,uxx,uxy,uyy) and calculating it's derivatives """

#Parameters
a=sp.Symbol('a')
b=sp.Symbol('b')

#Variables of lagrangian
x1=sp.Symbol('x1')
x2=sp.Symbol('x2')
x3=sp.Symbol('x3')
x4=sp.Symbol('x4')
x5=sp.Symbol('x5')
f=sp.Symbol('f')

L=(x1-f)**2 + a*psi1(x2**2+x3**2) + b*psi2((x4+x5)**2)
Lu=sp.diff(L,x1)
Lux=sp.diff(L,x2)
Luy=sp.diff(L,x3)
Luxx=sp.diff(L,x4)
Luyy=sp.diff(L,x5)


""" Defining RBF and calculating it's derivatives """

#Mean variables (x,y)
x,y,xi,yi=sp.symbols('x y,xi,yi')

#Shape parameter
e=sp.Symbol('e')

#RBF and it's derivatives
phi=sp.exp(-e**2*((x-xi)**2+(y-yi)**2))
phix=sp.diff(phi,x).simplify()
phiy=sp.diff(phi,y).simplify()
phixx=sp.diff(phix,x).simplify()
phiyy=sp.diff(phiy,y).simplify()


""" Computing collotacion points """
xv=np.linspace(0.,1.,Nx+1)[1::]-(1./(2.*Nx))
yv=np.linspace(0.,1.,Ny+1)[1::]-(1./(2.*Ny))


""" Creating symbols for c_i of RBFs """
carray=["c"+str(i) for i in range(N)]
c=sp.symbols(carray)


""" Computing u, ux, uy, uxx, uyy with N points"""
u=0.
ux=0.
uy=0.
uxx=0.
uyy=0.


for i in range(Nx):
	tmpx=xv[i]
	for j in range(Ny):
		tmpy=yv[j]
		u+=c[10*i+j]*phi.subs(xi,tmpx).subs(yi,tmpy)
		ux+=c[10*i+j]*phix.subs(xi,tmpx).subs(yi,tmpy)
		uy+=c[10*i+j]*phiy.subs(xi,tmpx).subs(yi,tmpy)
		uxx+=c[10*i+j]*phixx.subs(xi,tmpx).subs(yi,tmpy)
		uyy+=c[10*i+j]*phiyy.subs(xi,tmpx).subs(yi,tmpy)






""" Building the EL equation """
Lu=Lu.subs(x1,u)
Lux=Lux.subs(x2,ux).subs(x3,uy)
Luy=Luy.subs(x2,ux).subs(x3,uy)
Luxx=Luxx.subs(x4,uxx).subs(x4,uyy)
Luyy=Luyy.subs(x5,uxx).subs(x5,uyy)
EL=Lu - sp.diff(Lux,x) - sp.diff(Luy,y) + sp.diff(Luxx,x,2) + sp.diff(Luyy,y,2)
