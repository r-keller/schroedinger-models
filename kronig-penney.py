# Kronig-Penney Model of 1d Periodic System
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
# Periodic potential well:
#	| v(x) = -v0 for |x| < a/2; 
#	| v(x) =  0  for |x| > a/2
# Periodicity: V(x+v) = V(x)
# V(x) = \sum_n v(x+vn)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Customizable code:
## v0 = np.float(raw_input('Please enter amplitude of pot: v0 >> '));
## a  = np.float(raw_input('Please enter domain of v0 (a in V(x) = -v0 for |x| < a/2): a >> '));
## v  = np.float(raw_input('Please enter  v >= a for V(x) = sum_k[v(x-nv)]: v >> '));
## N  = np.int(raw_input('Please enter number of translates (tiles) N >> '));
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
v0 = 5.0;
a  = 1.0; # non-zero portion of periodic region
v  = 2.0; # full region of periodicity
N = 3*2; # number of translates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct Kronig-Penney Model Potential
# (Not quite necessary for computation of results..but nice nonetheless.)
n = 1000
m = 10; # number of sums in Kronig-Penney sum
b = v; 

dom = N*v;
x = np.linspace(-dom/2, dom/2,n+1)
itr_reg = n/N;

def V(x):
    vpot = np.zeros(len(x))
    for i in range(N):
        vpot[i*(itr_reg):(i+1)*itr_reg] = [-v0 if abs(abs(y-x[0])-i*v) < a else 0 for y in x[i*itr_reg:(i+1)*itr_reg] ]
    return vpot;
vpot = V(x);

ymax = max(vpot); # known: 0
ymin = min(vpot); # known: -v0
xmin = x[0]; xmax = x[n];

plt.plot(x,vpot,'b', x,np.zeros(len(x)),'k--')
plt.axis([xmin, xmax, ymin-1, ymax+1])
plt.show()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Approximate spectrum
e_cutoff = 20; # number of plane waves
k_max    = np.int(np.sqrt(e_cutoff/((2.0*np.pi/v)**2))+0.5);
k_max    = 2*k_max+1;

g = np.zeros([k_max])
e = np.zeros([k_max])
h = np.zeros([k_max,k_max])
#wrk = np.zeros([3*k_max])

g[0] = 0.0;
for i in range(1,k_max,2):
    g[i]   =  (i+1)*np.pi/v;
    g[i+1] = -(i+1)*np.pi/v;

# Loop over [-pi/v, pi/v]:
n = 20;
for i in range(-n, n+1):
    k = np.double(i)*np.pi/(np.double(n)*v);

ij = 0;
for i in range(k_max):
    for j in range(k_max):
        ij=ij+1;
        if i == j:
            h[i][j] = (k + g[i]) * (k+g[i] - v0/v*a);
        else:
            h[i][j] = -v0/v * np.sin( (g[j]-g[i])*a/2.0)/(g[j]-g[i])*2.0;
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Solve using Scipy's lapack algorithm.
vals, vecs, out = la.lapack.dsyev(h)
print vals, vecs, out



