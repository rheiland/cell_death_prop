# https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html
import pandas as pd
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

#df = pd.read_csv('xy_radial.dat', sep=',',header=None)
df = pd.read_csv('File_20160909_b16f10_aMSH_xy38_simulation_info.csv', sep=',',header=1)
points = df.values
tod = np.arange(len(points))

#plt.scatter(points[:,2], points[:,3], s=10, c=points[:,4])
xp = points[:,2]
yp = points[:,3]
xy_pairs = np.array((xp,yp)).T
tod = points[:,4]

tri = Delaunay(xy_pairs)

plt.triplot(xp, yp, tri.simplices)
#plt.plot(points[:,0], points[:,1], 'o')

for j, p in enumerate(xp):
    # plt.text(xp[j]-0.03, yp[j]+0.03, int(tod[j]), ha='right', size='smaller') # label w TOD
    plt.text(xp[j]-0.03, yp[j]+0.03, j, ha='right', size='smaller') # label w index

#plt.title('Delaunay triangulation, cells labeled with time of death')
plt.title('Delaunay triangulation, cells labeled with index')
plt.xticks([])
plt.yticks([])
plt.show()
