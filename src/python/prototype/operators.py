#######################################################################
# This file contains all the functions for matrix operations
# ----Used in LDC problem
#######################################################################
import numpy as np
import scipy.sparse as scysparse
from pdb import set_trace as keyboard
import scipy.sparse.linalg as spysparselinalg  # sparse linear algebra
import scipy.linalg as scylinalg               # non-sparse linear algebra
import pylab as plt
import math

#--------M_i,j =  (phi_i-1/2,j + phi_i,j)/2--------
def Mxh(phi):
	M = np.zeros_like(phi)
	M[:,:-1] = 0.5*(phi[:,1:] + phi[:,:-1])
	M[:,-1] = M[:,-2]
	return M
#--------M_i,j =  (phi_i,j-1/2 + phi_i,j)/2--------
def Myh(phi):
	M = np.zeros_like(phi)
	M[:-1,:] = 0.5*(phi[1:,:] + phi[:-1,:])
	M[-1,:] = M[-2,:]
	return M
#--------Face-centered x-gradient------------------
def Dxh(phi):
	M = np.zeros_like(phi)
	M[:,:-1] = (phi[:,1:] - phi[:,:-1])
	M[:,-1] = M[:,-2]
	return M
#--------Face-centered y-gradient------------------
def Dyh(phi):
	M = np.zeros_like(phi)
	M[:-1,:] = (phi[1:,:] - phi[:-1,:])
	M[-1,:] = M[-2,:]
	return M
#--------Cell-centered x-neighbors mean------------
def Mx(phi):
	M = np.zeros_like(phi)
	M[:,1:] = 0.5*(phi[:,:-1] + phi[:,1:])
	M[:,0] = M[:,1]
	return M
#--------Cell-centered y-neighbors mean------------
def My(phi):
	M = np.zeros_like(phi)
	M[1:,:] = 0.5*(phi[:-1,:] + phi[1:,:])
	M[0,:] = M[1,:]
	return M
#--------Cell-centered x-gradient------------------
def Dx(phi):
	M = np.zeros_like(phi)
	M[:,1:] = (phi[:,1:] - phi[:,:-1])
	M[:,0] = M[:,1]
	return M
#--------Cell-centered y-gradient------------------
def Dy(phi):
	M = np.zeros_like(phi)
	M[1:,:] = (phi[1:,:] - phi[:-1,:])
	M[0,:] = M[1,:]
	return M
#--------Setting sparse matrix for Pressure--------
def Pressure_Matrix(N):
	A = scysparse.lil_matrix((N*N,N*N),dtype="float64")
	Nt = N*N
	for i in range(0,N):
		for j in range(0,N):
			k = i + j*N
			if(i==0):
				A[k,k] = 1.0
				A[k,k+1] = -1.0
			elif(i==N-1):
				A[k,k] = 1.0
				A[k,k-1] = -1.0
			elif(j==0):
				A[k,k] = 1.0
				A[k,k+N] = -1.0
			elif(j==N-1):
				A[k,k] = 1.0
				A[k,k-N] = -1.0
			elif((i==1)and(j==1)):
				A[k,k] = 1.0
			else:
				A[k,k] = 4.0
				A[k,k-1]=-1.0
				A[k,k+1]=-1.0
				A[k,k-N]=-1.0
				A[k,k+N]=-1.0
	A = A.tocsr()
	return A
#--------2D--->1D flattened--------
def compress_rhs(rhs,N):
	vec = np.zeros((N*N,1))
	for i in range(0,N):
		for j in range(0,N):
			k = i + j*N
			vec[k] = rhs[i,j]
	return vec
#--------1D flattened--->2D--------
def expand_p(vec,N):
	p = np.zeros((N,N))
	for i in range(0,N):
		for j in range(0,N):
			k = i + j*N
			p[i,j] = vec[k]
	return p	
