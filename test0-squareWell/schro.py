# Schroedinger equation :: Variational Method
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adapted from Ian Cooper's work :: University of Sydney
# 	http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Purpose:
# 		To solve or eigenpairs of the 1-D Schroedinger equation 
#		with various potentials and construct a set of basis 
#		functions using the eigenvectors obtained.
# ~~~~~~~~~~~~~
# For later on:
#		Hope to use this code to construct appropriate basis
# 		to solve 2D Schroedinger with particular potential
# 		using a variational technique, i.e. using the eigenbasis
# 		obtained to construct solutions.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Coding Lessons:
# 0. Variables in lists
#		If using a variable list condition to form or fill an 
# 		array in some desired way, e.g.
#			v = [x for x in y],
#		then x is now defined to hold that value from y.
# 		(That is, if you previously defined/filled x, then
#		it is overwritten and now holds that last value from y.)
# 1. np.eigs: 
#		Evals and evecs are NOT ordered necessarily in assumed way. 
#		Be careful to find the indices for which the desired properties 
# 		for the evals are met.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import scipy.sparse as sparse
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Molecular dynamics properties (commented out-- assume unit values)
# #					     # SI units (description)
hbar = 1.055*(1e-34)  # J.s (Planck's constant)
e    = 1.602*1e-19;    # C
me   = 9.109*1e-31;    # kg (electron mass)
mp   = 1.67252*1e-27   # kg (proton mass)
mn   = 1.67482*1e-27   # kg (neutron mass)
eps0 = 8.854*1e-12;    # F/m (vacuum dielectric const)
m    = me;             # kg (mass of particle)
Ese = 1.6*1e-19; # energy scaling factor
Lse = 1e-9;
Cse = -1.0*(hbar**2.0)/(2.0*m)
Cse = Cse/((Lse**2.0)*Ese) # Schro constant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function: square_well(n)
# Input:   n - number of grid points
# Returns: u - array size n: discretized square well potential 
def square_well(n):
	# Single Square well
	xmin = -1.0;
	xmax = 1.0;
	a    = 0.5;
	#n = (n-1)**2+1; # refine more since small domain comp to lattice
	dx = (xmax-xmin)/(n-1)
	v0   = 550;

	# ~~~~~~~~~~~~~~
	# test: comparison with u of syd code:
	xmin = -0.1;
	xmax = 0.1;
	v0=400;
	xc = 0.05; # well width
	n = 1001;

	x = np.linspace(xmin, xmax, n);
	u = [-v0 if abs(y) <= xc else 0 for y in x]

	# ~~~~~~~~~~~~

	#x = np.linspace(xmin, xmax, n);
	#u = [-v0 if abs(y) < a else 0 for y in x]




	return x, u;

	# ax_xmax = xmax*1.1;
	# ax_xmin = xmin*1.1;
	# ax_ymax = 1;
	# ax_ymin = -v0-1;
	# plt.plot(x,u)
	# plt.axis([xmin, xmax, ax_ymin, ax_ymax])
	# plt.show()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Switch statement (in some sense)
	# Call line:
	# 		options[it](n)
	# Purpose: 
	# 	Directs to function associated with itr and supplies argument n
	# 	Just want to easily switch between potentials.
	options = {0 : square_well,
	                1 : lattice,
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	n = 1001;
	x, vpot = options[0](n);
	# ~~~~~~~~~~~~~~~~~~~
	dx = x[1]-x[0];
	print dx

	# Construct centered-difference matrix -Delta
	e = np.ones(n-2)
	A = sparse.spdiags([e, -2.0 * e, e], [-1, 0, 1], n-2,n-2)
	A =A/(dx**2.0)
	
	K = Cse*A;
	U = np.zeros([n-2,n-2])
	np.fill_diagonal(U, vpot[1:n-1])
	# Add on vpot to construct full Hamiltonian
	H = K+ U;
	# # checks: #print H[0:2,0:2]#print H.shape[0], H.shape[1]

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Plot potential:
	xmin=x[0];xmax=x[n-1]; ymin=min(vpot); ymax=max(vpot);
	buffx = 10.0*dx;buffy = 0.05*abs(ymax-ymin)

	plt.plot(x,vpot)
	plt.axis([xmin-buffx, xmax+buffx, ymin-buffy, ymax+buffy])
	plt.xlabel('$x$')
	plt.ylabel(r'$V(x)$');
	plt.title('Plot of Potential $V(x)$')
	plt.show()
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	vals, vecs = np.linalg.eig(H) # Note w[i] has evec v[:,i]

	# Obtain the neg evals
	flag = 0;
	E    = vals[vals<0]
	ind  = np.where(vals < 0)[0]; # get indices for which vals < 0

	
	psi=np.zeros([vecs.shape[0]+2,len(E)])
	for j in range(len(E)):
		psi[1:-1,j] = [k for k in vecs[:,ind[j]]];
		print psi[:,j]
		psi_j2   = np.multiply(psi[:,j],psi[:,j]);
		vol      = integrate.simps(psi_j2,x);
		if vol != 0:
			psi[:,j] = psi[:,j]/np.sqrt(vol)
		else:
			print('Warning! L2 inner product of %d integrated to 0', j)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Plot Eigenfunctions:		
	for j in range(len(E)):
		plt.plot(x, psi[:,j])
		xmin=x[0];xmax=x[n-1]; ymin=min(psi[:,j]); ymax=max(psi[:,j]);
		buffx = 0.0*dx;buffy = 0.05*abs(ymax-ymin)
		plt.axis([xmin-buffx, xmax+buffx, ymin-buffy, ymax+buffy])
		plt.xlabel('$x$')
		plt.ylabel(r'$\psi_{}(x)$'.format(j));
		title_str = 'Eigenfunction $\psi_{}$'.format(j) +' with Eigenvalue $\mu_{} = $'.format(j) + str('%.4f'%E[j]) #'for $x \in [$' + str(xmin) + '$,$' + str(xmax) + '$]$' ;
		plt.title(title_str,fontsize = 12)
		plt.show()
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	np.savetxt("psi-test0.csv", psi, delimiter=","); # save efuncs to file
	



if __name__ == "__main__":
	main();




