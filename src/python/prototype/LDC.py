#################################################
# Lid Driven Cavity Prototype
#------------------------------------------------
# --Uses Chorin's Projection Method
# --Explicit diffusion (Convection dominant problem)
# --Uniform spacing & time stepping
#################################################
import numpy as np
import scipy.sparse as scysparse
from pdb import set_trace as keyboard
import scipy.sparse.linalg as spysparselinalg  # sparse linear algebra
import scipy.linalg as scylinalg               # non-sparse linear algebra
import pylab as plt
import math
import time
import os
#--------Operators--& Data----------------------#
import operators as op
import ghia_data as gh
#--------LDC Parameters defined----------------#
N = 40
L = 1.0
Re = 1000.0
tf = 30.0

start_time = time.time()
#----------------------------------------------#
def main():
#-----Defining grid----------------------------#	
	dx = L/N
	dx2 = dx*dx

	xh = np.linspace(0.0, L, N+1)
	yh = np.linspace(0.0, L, N+1)

	xh = np.concatenate([xh, [xh[-1]+dx] ])
	yh = np.copy(xh)

	x = np.copy(xh)
	y = np.copy(yh)

	x[1:] = 0.5*( xh[1:]+xh[:-1])
	x[0] = -x[1]
	y = np.copy(x)

	[X,Y] = np.meshgrid(x,y)
	[Xu,Yu]=np.meshgrid(xh,y)
	[Xv,Yv]=np.meshgrid(x,yh)

#-------Initializing fields-----------------------#
	u = np.zeros_like(Xu)
	v = np.zeros_like(Yv)
	p = np.zeros_like(X)

	us = np.copy(u)
	vs = np.copy(v)

	du = np.zeros_like(Xu)
	dv = np.zeros_like(Yv)
	rhs = np.zeros_like(p)

	t = 0.0	
	os.makedirs('Visualize', exist_ok=True)
#------Fixed time problem--------------------------#
	dt = 0.1*min(dx2*Re,dx)
	it = 0

	Nt = N + 2

#-----Creating Sparse Matrix for pressure apriori----#	
	A = scysparse.csr_matrix((Nt*Nt,Nt*Nt),dtype="float64")
	A = op.Pressure_Matrix(Nt)
	vec = np.zeros((Nt*Nt,1))
	ans = np.copy(vec)

#----Time loop---------------------------------------#
	while(t< tf):
		it = it + 1		

#-----Boundary Conditions--------
		v[:, 0] = -v[:,1 ]
		u[:, 0] = 0.0
		v[:,-1] = -v[:,-2]
		u[:,-2] = 0.0
		u[:,-1] = 0.0
		u[0, :] = -u[1, :]
		v[0, :] = 0.0
		u[-1,:] = 2.0 - u[-2, :]
		v[-2,:] = 0.0
		v[-1,:] = 0.0
#-----Convection of U-&-V-----------------------------
		du = (-(op.Dxh(op.Mx(u**2.0))/dx)
		     - (op.Dy( (op.Mxh(v)) * (op.Myh(u)) )/dx))

		dv = (-(op.Dyh(op.My(v**2.0))/dx)
		     - (op.Dx( (op.Myh(u)) * (op.Mxh(v)) )/dx))

#-----Diffusion of U-&-V------------------------------		
		du = du + ((op.Dxh((op.Dx(u))))/(dx2*Re))

		du = du + ((op.Dy( (op.Dyh(u))))/(dx2*Re))

		dv = dv + ((op.Dyh((op.Dy(v))))/(dx2*Re))

		dv = dv + ((op.Dx( (op.Dxh(v))))/(dx2*Re))
#-----Update U-&-V-------------------------------------
		us = du*dt + u

		vs = dv*dt + v
#-----Boundary Conditions------------------------------
		us[:, 0] = 0.0
		us[:,-2] = 0.0
		vs[0, :] = 0.0
		vs[-2,:] = 0.0

#-----Calculate rhs for pressure equation--------------
		rhs = (-op.Dx(us) -op.Dy(vs))*dx/dt
		rhs[:,0] = 0.0
		rhs[:,-1]=0.0
		rhs[0,:] = 0.0
		rhs[-1,:]=0.0
		rhs[1,1] = 0.0####Rooting----------------------

#----------Solving the pressure poisson equation-------
		vec = op.compress_rhs(rhs,Nt)
		ans = spysparselinalg.spsolve(A,vec)
		p = op.expand_p(ans,Nt)

#-----------Updating velocity--------------------------
		u = -(op.Dxh(p)*dt/dx) + us
		v = -(op.Dyh(p)*dt/dx) + vs

		t = t + dt
		if((it % 200)==0):
			print(" t= ",t)

	print("Simulation complete!")
	
	# Get Ghia's benchmark data
	y_ghia, u_ghia, x_ghia, v_ghia = gh.get_ghia_data_re1000()
	
	#-------- Plot 1: U velocity contour --------
	fig = plt.figure(0, figsize=(8,7))
	ax = fig.add_axes([0.15,0.15,0.75,0.75])
	plt.axes(ax)
	levels = np.linspace(np.min(u[:-1,:-1]), np.max(u[:-1,:-1]), 20)
	contour = plt.contourf(Xu[:-1,:-1], Yu[:-1,:-1], u[:-1,:-1], levels=levels, cmap='RdBu_r')
	cbar = plt.colorbar(contour)
	cbar.set_label('U velocity', fontsize=12)
	plt.xlabel('x', fontsize=14)
	plt.ylabel('y', fontsize=14)
	plt.title(f'U Velocity Contour - Lid Driven Cavity (Re={int(Re)})', fontsize=16, fontweight='bold')
	plt.grid(True, alpha=0.3)
	plt.axis('equal')
	plt.savefig("Visualize/U_contour_Re1000.jpg", dpi=300, bbox_inches='tight')	
	plt.close()
	
	#-------- Plot 2: V velocity contour --------
	fig = plt.figure(1, figsize=(8,7))
	ax = fig.add_axes([0.15,0.15,0.75,0.75])
	plt.axes(ax)
	levels = np.linspace(np.min(v[:-1,:-1]), np.max(v[:-1,:-1]), 20)
	contour = plt.contourf(Xv[:-1,:-1], Yv[:-1,:-1], v[:-1,:-1], levels=levels, cmap='RdBu_r')
	cbar = plt.colorbar(contour)
	cbar.set_label('V velocity', fontsize=12)
	plt.xlabel('x', fontsize=14)
	plt.ylabel('y', fontsize=14)
	plt.title(f'V Velocity Contour - Lid Driven Cavity (Re={int(Re)})', fontsize=16, fontweight='bold')
	plt.grid(True, alpha=0.3)
	plt.axis('equal')
	plt.savefig("Visualize/V_contour_Re1000.jpg", dpi=300, bbox_inches='tight')
	plt.close()
	
	#-------- Plot 3: U vs y at centerline (x=0.5) with Ghia's data --------
	# Find the index closest to x=0.5
	centerline_idx = np.argmin(np.abs(xh - 0.5))
	
	fig = plt.figure(2, figsize=(8,7))
	ax = fig.add_axes([0.15,0.15,0.75,0.75])
	
	# Plot computed data
	plt.plot(u[1:-2, centerline_idx], y[1:-2], 'b-', linewidth=2, label='Prototype')
	
	# Plot Ghia's benchmark data
	plt.plot(u_ghia, y_ghia, 'ro', markersize=8, markerfacecolor='none', 
	         markeredgewidth=2, label='Ghia et al. (1982)')
	
	plt.xlabel('U velocity', fontsize=14)
	plt.ylabel('y', fontsize=14)
	plt.title(f'U Velocity along Vertical Centerline (Re={int(Re)})', 
	          fontsize=16, fontweight='bold')
	plt.grid(True, alpha=0.3)
	plt.legend(fontsize=12, loc='best')
	plt.xlim([-0.4, 1.1])
	plt.ylim([0, 1])
	plt.savefig("Visualize/U_vs_y_Re1000.jpg", dpi=300, bbox_inches='tight')
	plt.close()
		
	
	print("All plots saved successfully in folder--Visualize")
	print("  - U_contour.jpg")
	print("  - U_vs_y_comparison.jpg")
	print("  - V_vs_x_comparison.jpg")


main()
end_time = time.time()
print("Runtime = ",end_time-start_time)
