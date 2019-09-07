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
x=points[:,0]
y=points[:,1]

tod = np.arange(len(points))

fig = plt.figure()
#print('dir(fig)=',dir(fig))
fig.set_figwidth(7)
fig.set_figheight(7)

#area = (30 * np.random.rand(N))**2  # 0 to 15 point radii
#plt.scatter(x, y, s=area, c=tod, alpha=0.5)
area = 80
plt.scatter(x, y, s=area, c=tod)
#plt.scatter(x, y, c=tod)

plt.xticks([])
plt.yticks([])
plt.axes().set_aspect('equal')
plt.title('Cells colored by time of death')
plt.show()
