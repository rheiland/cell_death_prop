# https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html
import pandas as pd
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

# 6 cells
#points = np.array([[10,0], [17,3], [17,-3], [20,0], [22,0]  ])
#tod = np.array([15,15,15,20,22])

df = pd.read_csv('xy_radial.dat', sep=',',header=None)
points = df.values
tod = np.arange(len(points))

# 6 cells
#points = np.array([[10,0], [15,0], [17,3], [17,-3], [20,0], [22,0]  ])
#tod = np.array([15,10,15,15,20,22])

# 7 cells
# points = np.array([[10,0], [15,0], [17,3], [17,-3], [16,1], [20,0], [22,0]  ])
# tod = np.array([15,10,15,15,14,20,22])

tri = Delaunay(points)

plt.triplot(points[:,0], points[:,1], tri.simplices)
plt.plot(points[:,0], points[:,1], 'o')

for j, p in enumerate(points):
    plt.text(p[0]-0.03, p[1]+0.03, tod[j], ha='right') # label the points
#for j, s in enumerate(tri.simplices):
#    p = points[s].mean(axis=0)
#    plt.text(p[0], p[1], '#%d' % j, ha='center') # label triangles
#plt.xlim(-0.5, 1.5); plt.ylim(-0.5, 1.5)
plt.title('Delaunay triangulation, cells labeled with time of death')
plt.xticks([])
plt.yticks([])
plt.show()
