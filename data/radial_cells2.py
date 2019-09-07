from math import sin,cos
import random

#print('Cell index,Cell X,Cell Y,Time of Nucleation,Time of death')
cell_idx = 0
xval = yval = 0.
t_nuc = -1
t_death = 0
# print(cell_idx, ',',xval, ',',yval, ',',t_nuc, ',',t_death)
cell_idx = 1

#cell_radius = 10.
cell_radius = 1.
twopi = 6.28
alpha = twopi/8
#num_cells = 
ring_radius = cell_radius * 1.5

for ring in range(11):
    alpha = random.random() * twopi
    # print('alpha=',alpha)
    num_cells_ring = 8 + 2*ring
    alpha_del = twopi / num_cells_ring
    # for cell in range(int(twopi/alpha)):
    for cell in range(num_cells_ring):
        xval = ring_radius * cos(alpha)
        yval = ring_radius * sin(alpha)
        # print(cell_idx, ',',xval, ',',yval, ',',t_nuc, ',',t_death)
        print(xval, ',',yval)
        cell_idx += 1
        alpha += alpha_del
        # print(alpha)
    # alpha = twopi/(8 + ring + 3)
    # print('---> ',alpha)
    ring_radius += cell_radius * 1.2
