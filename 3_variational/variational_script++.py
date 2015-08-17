import numpy as np 
import sympy as sp
from scipy import misc
import pickle
import cloudpickle
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

t0=time.clock()

def lamb(func):
    return sp.lambdify(c, func, modules='numpy')

""" Defining some functions """
def psi1(x):
	return x**2

def psi2(x):
	return x**2

""" Taking discrete values of f function """
fv=misc.imread('Circle-Missing-data-reduced-gray.png').astype(float)


""" Some constants and parameter values """
N=4 # Dimensions of image, assuming square form
eps=sp.Float(0.5) #shape Parameters
alpha=sp.Float(1.)
beta=sp.Float(1.)	


""" Defining Lagrangian funtion L(x,y,u,ux,uy,uxx,uxy,uyy) and calculating it's derivatives """

#Variables of lagrangian
x1,x2,x3,x4,x5=sp.symbols('x1 x2 x3 x4 x5')

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
	tmpx=sp.Float(xv[i])
	for j in range(N):
		tmpy=sp.Float(yv[j])
		u+=c[N*i+j]*phi.subs(xi,tmpx).subs(yi,tmpy)
		ux+=c[N*i+j]*phix.subs(xi,tmpx).subs(yi,tmpy)
		uy+=c[N*i+j]*phiy.subs(xi,tmpx).subs(yi,tmpy)
		uxx+=c[N*i+j]*phixx.subs(xi,tmpx).subs(yi,tmpy)
		uyy+=c[N*i+j]*phiyy.subs(xi,tmpx).subs(yi,tmpy)
		points.append((tmpx,tmpy))


""" Building the EL equation """

Lu=Lu.subs(x1,u)
Lux=Lux.subs(x2,ux).subs(x3,uy)
Luy=Luy.subs(x2,ux).subs(x3,uy)
Luxx=Luxx.subs(x4,uxx).subs(x5,uyy)
Luyy=Luyy.subs(x4,uxx).subs(x5,uyy)

def dx(expr):
	return sp.diff(expr,x)

def dy(expr):
	return sp.diff(expr,y)

def d2x(expr):
	return sp.diff(expr,x,2)

def d2y(expr):
	return sp.diff(expr,y,2)

result=[]

with ProcessPoolExecutor(max_workers=4) as executor:
	result.append(executor.submit(dx,Lux))
	result.append(executor.submit(dy,Luy))
	result.append(executor.submit(d2x,Luxx))
	result.append(executor.submit(d2y,Luyy))

EL=Lu-sp.diff(Lux,x)-sp.diff(Luy,y)+sp.diff(Luxx,x,2)+sp.diff(Luyy,y,2)


""" Evaluating EL equation at domain points [0,1]x[0,1] (falta restar 2f) """
evEL=[EL.subs(x,points[i][0]).subs(y,points[i][1]) for i in range(N**2)]

#Now we have EL equation evaluated on the domain points (only dependent on ci coef)


""" Coneverting into NumPy evaluated functions """
F=map(lamb,evEL)


""" test """
params=np.linspace(0,1,N**2)
for i in range(N**2):
	F[i](*params)



""" Serializing it and storing with Tickle... """
fobj=open("sample.pickle","wb")
cloudpickle.dump(F,fobj) # Dump F object into fobj file
fobj.close()

tf=time.clock()

print("execution time:",tf-t0)