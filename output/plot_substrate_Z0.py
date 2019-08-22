

"""
http://www.mathcancer.org/blog/working-with-physicell-snapshots-in-matlab/

Let’s make an oxygen contour plot through z = 0 μm. First, we find the index corresponding to this z-value:

k = find( MCDS.mesh.Z_coordinates == 0 );

--------------
import scipy.io
fname='output00000000_microenvironment0'
fname='output00000001_microenvironment0'
info_dict = {}
scipy.io.loadmat(fname, info_dict)
M = info_dict['multiscale_microenvironment']
len(M)
M.shape
M.shape[1]
M.shape[1]/(80*80)
M[3,:]
min(M[3,:])
max(M[3,:])
max(M[2,:])
min(M[2,:])
max(M[2,:])
790+790
790+790/80.
(790+790)/80.
M[2,:]
"""

import scipy.io
fname='output00000000_microenvironment0'
info_dict = {}
scipy.io.loadmat(fname, info_dict)
M = info_dict['multiscale_microenvironment']
len(M)
M.shape
M.shape[1]
M.shape[1]/(80*80)
M[3,:]
min(M[3,:])
max(M[3,:])
max(M[2,:])
min(M[2,:])
max(M[2,:])

print('min O2 = ',min(M[4,:]))   # 38.0
print('max O2 = ',max(M[4,:]))
