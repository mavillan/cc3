import numpy as np 
import sympy as sp

#some constants
N = 50
eps = 0.5

#symbols definition
x,y = sp.symbols('x y')

#RBF and it's derivatives
phi = sp.exp(-eps**2.*(x**2+y**2))

#first order
phix = sp.diff(phi,x).simplify()
phiy = sp.diff(phi,y).simplify()

#second order
phixx = sp.diff(phix,x).simplify()
phixy = sp.diff(phix,y).simplify()
phiyx = sp.diff(phiy,x).simplify()
phiyy = sp.diff(phiy,y).simplify()


#third order
phixxx = sp.diff(phixx,x).simplify()
phixxy = sp.diff(phixx,y).simplify()
phiyyx = sp.diff(phiyy,x).simplify()
phiyyy = sp.diff(phiyy,y).simplify()


#fourth order
phixxxx = sp.diff(phixxx,x).simplify()
phixxyy = sp.diff(phixxy,y).simplify()
phiyyxx = sp.diff(phiyyx,x).simplify()
phiyyyy = sp.diff(phiyyy,y).simplify()


#lambdify all previous symbolic expressions
phi = sp.lambdify((x,y), phi, modules='numpy')
phix = sp.lambdify((x,y), phix, modules='numpy')
phiy = sp.lambdify((x,y), phiy, modules='numpy')
phixx = sp.lambdify((x,y), phixx, modules='numpy')
phixy = sp.lambdify((x,y), phixy, modules='numpy')
phiyx = sp.lambdify((x,y), phiyx, modules='numpy')
phiyy = sp.lambdify((x,y), phiyy, modules='numpy')
phixxx = sp.lambdify((x,y), phixxx, modules='numpy')
phixxy = sp.lambdify((x,y), phixxy, modules='numpy')
phiyyx = sp.lambdify((x,y), phiyyx, modules='numpy')
phiyyy = sp.lambdify((x,y), phiyyy, modules='numpy')
phixxxx = sp.lambdify((x,y), phixxxx, modules='numpy')
phixxyy = sp.lambdify((x,y), phixxyy, modules='numpy')
phiyyxx = sp.lambdify((x,y), phiyyxx, modules='numpy')
phiyyyy = sp.lambdify((x,y), phiyyyy, modules='numpy')



""" Computing collocation points """
xv=np.linspace(0.,1.,N+1)[1::]-(1./(2.*N))
yv=np.linspace(0.,1.,N+1)[1::]-(1./(2.*N))

""" Computing distances matrix """
Dx = np.empty((N,N))
Dy = np.empty((N,N))

for k in range(N):
	Dx[k,:] = (xv[k] - xv)
	Dy[k,:] = (yv[k] - yv)


""" Computing matrices """
phi_m = np.empty((N**2,N**2))
phix_m = np.empty((N**2,N**2))
phiy_m = np.empty((N**2,N**2))
phixx_m = np.empty((N**2,N**2))
phixy_m = np.empty((N**2,N**2))
phiyx_m = np.empty((N**2,N**2))
phiyy_m = np.empty((N**2,N**2))
phixxx_m = np.empty((N**2,N**2))
phixxy_m = np.empty((N**2,N**2))
phiyyx_m = np.empty((N**2,N**2))
phiyyy_m = np.empty((N**2,N**2))
phixxxx_m = np.empty((N**2,N**2))
phixxyy_m = np.empty((N**2,N**2))
phiyyxx_m = np.empty((N**2,N**2))
phiyyyy_m = np.empty((N**2,N**2))


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phi_m[r_index,:] = phi(dx,dy).ravel()
		r_index += 1
outfile = open('phi_m.npy','wb')
np.save(outfile, phi_m)
outfile.close()
del phi_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phix_m[r_index,:] = phix(dx,dy).ravel()
		r_index += 1
outfile = open('phix_m.npy','wb')
np.save(outfile, phix_m)
outfile.close()
del phix_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiy_m[r_index,:] = phiy(dx,dy).ravel()
		r_index += 1
outfile = open('phiy_m.npy','wb')
np.save(outfile, phiy_m)
outfile.close()
del phiy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixx_m[r_index,:] = phixx(dx,dy).ravel()
		r_index += 1
outfile = open('phixx_m.npy','wb')
np.save(outfile, phixx_m)
outfile.close()
del phixx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixy_m[r_index,:] = phixy(dx,dy).ravel()
		r_index += 1
outfile = open('phixy_m.npy','wb')
np.save(outfile, phixy_m)
outfile.close()
del phixy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyx_m[r_index,:] = phiyx(dx,dy).ravel()
		r_index += 1
outfile = open('phiyx_m.npy','wb')
np.save(outfile, phiyx_m)
outfile.close()
del phiyx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyy_m[r_index,:] = phiyy(dx,dy).ravel()
		r_index += 1
outfile = open('phiyy_m.npy','wb')
np.save(outfile, phiyy_m)
outfile.close()
del phiyy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixxx_m[r_index,:] = phixxx(dx,dy).ravel()
		r_index += 1
outfile = open('phixxx_m.npy','wb')
np.save(outfile, phixxx_m)
outfile.close()
del phixxx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixxy_m[r_index,:] = phixxy(dx,dy).ravel()
		r_index += 1
outfile = open('phixxy_m.npy','wb')
np.save(outfile, phixxy_m)
outfile.close()
del phixxy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyyx_m[r_index,:] = phiyyx(dx,dy).ravel()
		r_index += 1
outfile = open('phiyyx_m.npy','wb')
np.save(outfile, phiyyx_m)
outfile.close()
del phiyyx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyyy_m[r_index,:] = phiyyy(dx,dy).ravel()
		r_index += 1
outfile = open('phiyyy_m.npy','wb')
np.save(outfile, phiyyy_m)
outfile.close()
del phiyyy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixxxx_m[r_index,:] = phixxxx(dx,dy).ravel()
		r_index += 1
outfile = open('phixxxx_m.npy','wb')
np.save(outfile, phixxxx_m)
outfile.close()
del phixxxx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phixxyy_m[r_index,:] = phixxyy(dx,dy).ravel()
		r_index += 1
outfile = open('phixxyy_m.npy','wb')
np.save(outfile, phixxyy_m)
outfile.close()
del phixxyy_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyyxx_m[r_index,:] = phiyyxx(dx,dy).ravel()
		r_index += 1
outfile = open('phiyyxx_m.npy','wb')
np.save(outfile, phiyyxx_m)
outfile.close()
del phiyyxx_m


r_index = 0
for i in range(N):
	for j in range(N):
		dx,dy = np.meshgrid(Dx[j,:], Dy[i,:], sparse=True)
		phiyyyy_m[r_index,:] = phiyyyy(dx,dy).ravel()
		r_index += 1
outfile = open('phiyyyy_m.npy','wb')
np.save(outfile, phiyyyy_m)
outfile.close()
del phiyyyy_m
