from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import csv
import scipy as sp
from scipy import stats
from scipy import interpolate
import cv2
from shapely.geometry import Polygon
from bs4 import BeautifulSoup
import pandas as pd

# configuration file
in_configuration_file = open("configuration.xml", "r")
contents = in_configuration_file.read()
soup = BeautifulSoup(contents, 'xml')
parameters = {'threshold': soup.find_all('death_threshold')[0].get_text(),
              'substrate': soup.find_all('initial_internalized_substrate')[0].get_text(),
              'rate to absorb': soup.find_all('signal_internalization_rate')[0].get_text(),
              'df coefficient': soup.find_all('signal_diffusion_coefficient')[0].get_text()}

# time_comparison
diffs = []
sum_diff = 0
file = open("dead_cells.dat", "r")
for line in file:
    words = line.split()
    diff = (np.abs(int(words[0]) - int(words[1])))  # (|t1-t2|)
    diffs.append([diff, words[2]])

for i in diffs:
    print("Cell index: " + str(i[1]) + ", time difference: " + str(i[0]))
    sum_diff += i[0]
print("Total difference: " + str(sum_diff))

# frac_of_dead_cells
num_of_cells = len(open("dead_cells.dat").readlines())

# simulation
file2 = open("dead_cells.dat", "r")
accumulate_sim = 0
frac_died_until_time_sim = {}
frac_died_until_time_sim_time_scale = {}
sim_time_of_death_per_index = {}
real_time_of_death_per_index = {}
position_per_index = {}
for line in file2:
    words = line.split()
    accumulate_sim += 1
    frac_died_until_time_sim[int(words[0])] = accumulate_sim / num_of_cells
    sim_time_of_death_per_index[int(words[2])] = int(words[0])
    position_per_index[int(words[2])] = words[3]
    real_time_of_death_per_index[int(words[2])] = int(words[1])

# real
time_scale = 0
accumulate_real = 0
frac_died_until_time_real = {}
max_time_of_death = 0
file_name = ""
with open("exp_a.csv", "r") as real_data:
    lines = real_data.readlines()
    last = lines[-1]
    # skip the header row
    for line in lines:
        if line == lines[0]:
            continue
        words = line.split(',')
        accumulate_real += 1
        frac_died_until_time_real[int(words[3])] = accumulate_real / num_of_cells
        if int(words[4]) != 0 and time_scale == 0:
            time_scale = int(int(words[3])/int(words[4]))
        if line is last:
            max_time_of_death = int(words[3])
            file_name = words[6]

#  create fake points every 30
all_times = [x for x in range(0, max_time_of_death+time_scale, time_scale)]
frac_accumulated = 0
last_k = 0
last_x = 0
for x in all_times:
    frac_accumulated = 0
    for k, v in frac_died_until_time_sim.items():
        if k < last_k:
            continue
        if k <= x:
            frac_accumulated = v
        if k > x or k == max_time_of_death:
            if frac_accumulated == 0:
                frac_accumulated = frac_died_until_time_sim_time_scale[last_x]
            frac_died_until_time_sim_time_scale[x] = frac_accumulated
            last_k = k
            last_x = x
            break

# divide x values in max_time_of_death to normalize both sim and real times of death
dataX_sim = list(frac_died_until_time_sim_time_scale.keys())
dataX_sim_normal = []
for x in dataX_sim:
    dataX_sim_normal.append(x / max_time_of_death)
dataY_sim = list(frac_died_until_time_sim_time_scale.values())

dataX_real = list(frac_died_until_time_real.keys())
dataX_real_normal = []
for x in dataX_real:
    dataX_real_normal.append(x / max_time_of_death)
dataY_real = list(frac_died_until_time_real.values())

# calculate area between curves
def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C


def intersection(L1, L2):
    D = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x, y
    else:
        return False

# calculate the area between the sim curve and real curve by dividing the field into polygons and summing their areas
total_area = 0
area = 0
for i in range(len(dataX_real_normal)-1):
    a = (dataX_sim_normal[i], dataY_sim[i])
    b = (dataX_real_normal[i], dataY_real[i])
    c = (dataX_real_normal[i+1], dataY_real[i+1])
    d = (dataX_sim_normal[i+1], dataY_sim[i+1])
    L1 = line(list(a), list(d))
    L2 = line(list(b), list(c))
    # intersection detected
    inter = intersection(L1, L2)
    if inter:
        if (inter[0] >= a[0]) and (inter[0] <= d[0]):
            p1 = Polygon([a, b, inter])
            p2 = Polygon([inter, c, d])
            area = p1.area + p2.area
        else:
            p1 = Polygon([a, b, c, d])
            area = p1.area
    total_area += area

area_below_curve = np.trapz(dataX_real_normal, dataY_real)
areas_ratio = total_area/area_below_curve


line_up, = plt.plot(dataX_real_normal, dataY_real)
line_up.set_label('real experiment')
line_down, = plt.plot(dataX_sim_normal, dataY_sim)
line_down.set_label('simulation')
plt.legend(handles=[line_up, line_down], fontsize='20')
plt.title('Fraction of dead cells over time', fontsize='40')
plt.xlabel('time (min)', color='blue', fontsize='20')
plt.ylabel('fraction', color='blue', fontsize='20')
plt.xlim(0, 1)
plt.xticks(np.arange(0, 1.1, 0.1), fontsize='20')
plt.ylim(0, 1)
plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
plt.text(1, 1, 'Area between curves=\n' + str(total_area), fontsize='20', color='r')

temp_diff = {}
diffs_in_plot = {}
prev = 0
for (time1, frac1) in frac_died_until_time_sim.items():
    if frac1 > 0.75:
        temp_diff[0.75] = prev
        break
    elif (frac1 > 0.5) & (0.5 not in temp_diff):
        temp_diff[0.5] = prev
    elif (frac1 > 0.25) & (0.25 not in temp_diff):
        temp_diff[0.25] = prev
    prev = time1
prev = 0
for (time2, frac2) in frac_died_until_time_real.items():
    if frac2 > 0.75:
        diffs_in_plot[0.75] = temp_diff[0.75] - prev
        break
    elif (frac2 > 0.5) & (0.5 not in diffs_in_plot):
        diffs_in_plot[0.5] = temp_diff[0.5] - prev
    elif (frac2 > 0.25) & (0.25 not in diffs_in_plot):
        diffs_in_plot[0.25] = temp_diff[0.25] - prev
    prev = time2
plt.show()

plt.plot(list(diffs_in_plot.keys()), list(diffs_in_plot.values()))
plt.title('Difference in time of 0.25, 0.5 and 0.75 dead cells', fontsize='40')
plt.xlabel('fraction', fontsize='20', color='r')
plt.ylabel('time difference', fontsize='20', color='r')
plt.xticks(np.array(list(diffs_in_plot.keys())), fontsize='20')
#plt.yticks(np.arange(0, 50, 5), fontsize='20')
plt.show()

# prop neighbors
sim_neighbors_time_diffs = []
real_neighbors_time_diffs = []
normalized_sim_neighbors_time_diffs = []
normalized_real_neighbors_time_diffs = []
diff = 0
distX = 0
distY = 0
dist = 0
with open("exp_b.csv", "r") as real_data_neighbors:
    csv_reader_neighbors = csv.reader(real_data_neighbors, delimiter=',')
    # skip the header
    next(csv_reader_neighbors)
    for line in csv_reader_neighbors:
        cell_index = int(line[0])
        neighbors = line[1].split(',')
        for nei in neighbors:
            diff = np.abs(sim_time_of_death_per_index[cell_index] - sim_time_of_death_per_index[int(nei)])
            sim_neighbors_time_diffs.append(diff)
            distX = float(position_per_index[cell_index].split(',')[0]) - float(
                position_per_index[int(nei)].split(',')[0])
            distY = float(position_per_index[cell_index].split(',')[1]) - float(
                position_per_index[int(nei)].split(',')[1])
            dist = np.sqrt(distX ** 2 + distY ** 2)
            normalized_sim_neighbors_time_diffs.append(diff / dist)
            diff = np.abs(real_time_of_death_per_index[cell_index] - real_time_of_death_per_index[int(nei)])
            real_neighbors_time_diffs.append(diff)
            normalized_real_neighbors_time_diffs.append(diff / dist)

count_sim_diffs = {}
sim_neighbors_time_diffs_set = set(sim_neighbors_time_diffs)
sim_neighbors_time_diffs_list = list(sim_neighbors_time_diffs_set)

count_real_diffs = {}
real_neighbors_time_diffs_set = set(real_neighbors_time_diffs)
real_neighbors_time_diffs_list = list(real_neighbors_time_diffs_set)

# for x in sim_neighbors_time_diffs_list:
#     if x not in real_neighbors_time_diffs_list:
#         real_neighbors_time_diffs_list.append(x)

for x in real_neighbors_time_diffs_list:
    if x not in sim_neighbors_time_diffs_list:
        sim_neighbors_time_diffs_list.append(x)

sim_neighbors_time_diffs_list.sort()
real_neighbors_time_diffs_list.sort()


for x in sim_neighbors_time_diffs_list:
    X = str(x)
    if X not in count_sim_diffs:
        count_sim_diffs[X] = sim_neighbors_time_diffs.count(x)

count_sim_diffs_time_scale = {}
accu_val = 0
for k, v in count_sim_diffs.items():
    accu_val += v
    if int(k) % time_scale == 0:
        count_sim_diffs_time_scale[k] = accu_val
        accu_val = 0

for x in real_neighbors_time_diffs_list:
    X = str(x)
    if X not in count_real_diffs:
        count_real_diffs[X] = real_neighbors_time_diffs.count(x)

for x in count_sim_diffs_time_scale.keys():
    if x not in count_real_diffs.keys():
        count_real_diffs[x] = 0

list(count_real_diffs.keys()).sort()

sim_keys = list(count_sim_diffs_time_scale.keys())
sim_vals = list(count_sim_diffs_time_scale.values())
keys_real = list(count_real_diffs.keys())
vals_real = list(count_real_diffs.values())

y_pos = np.arange(len(sim_keys))
plt.subplot(2, 1, 1)
width = 0.5
sum_all_sim = sum(sim_vals)
sim_vals_fraction = []
for val in sim_vals:
    sim_vals_fraction.append(val/sum_all_sim)
sim_bar = plt.bar(y_pos, sim_vals_fraction, width, color='green')
plt.xticks(y_pos, np.array(sim_keys), fontsize='20')
plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
plt.title('Distribution of time of death differences between neighbors- simulation', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')

plt.subplot(2, 1, 2)
sum_all_real = sum(vals_real)
real_vals_fraction = []
for val in vals_real:
    real_vals_fraction.append(val/sum_all_real)
y_pos_real = np.arange(len(keys_real))
real_bar = plt.bar(y_pos_real, real_vals_fraction, width, color='blue')
plt.title('Distribution of time of death differences between neighbors- real experiment', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(y_pos_real, np.array(keys_real), fontsize='20')
plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
plt.show()

# earth movers distance
emd1 = sum(abs(np.cumsum(sim_vals_fraction)-np.cumsum(real_vals_fraction)))

plt.subplot(2, 1, 1)
cumulative_sim = np.cumsum(sim_vals_fraction)
sim_bar = plt.bar(y_pos, cumulative_sim, width, color='green')
plt.xticks(y_pos, np.array(sim_keys), fontsize='20')
plt.yticks(np.arange(0, max(cumulative_sim) + 0.1, 0.25), fontsize='20')
plt.title('Distribution of time of death differences between neighbors- simulation', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')

plt.subplot(2, 1, 2)
cumulative_real = np.cumsum(real_vals_fraction)
real_bar = plt.bar(y_pos_real, cumulative_real, width, color='blue')
plt.title('Distribution of time of death differences between neighbors- real experiment', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(y_pos_real, np.array(keys_real), fontsize='20')
plt.yticks(np.arange(0, max(cumulative_real) + 0.1, 0.25), fontsize='20')
plt.show()


plt.subplot(2, 1, 1)
plt.hist(normalized_sim_neighbors_time_diffs, density=True, bins=50)
plt.title('Distribution of normalized time of death differences between neighbors- simulation', fontsize='40')
plt.xlabel('time of death differences divided by distance', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(np.arange(0.0, max(normalized_real_neighbors_time_diffs), 0.2), fontsize='20')
plt.yticks(np.arange(0, 10, 2), fontsize='20')

plt.subplot(2, 1, 2)
plt.hist(normalized_real_neighbors_time_diffs, density=True, bins=50)
plt.title('Distribution of normalized time of death differences between neighbors- real experiment', fontsize='40')
plt.xlabel('time of death differences divided by distance', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(np.arange(0.0, max(normalized_real_neighbors_time_diffs), 0.2), fontsize='20')
plt.yticks(np.arange(0, 10, 2), fontsize='20')
plt.show()

# earth movers distance- normalized
emd2 = sum(abs(np.cumsum(normalized_sim_neighbors_time_diffs)-np.cumsum(normalized_real_neighbors_time_diffs)))
print("the EMD of the difference between time of death of each cell and its neighbors- normalized: ", emd2)

# Write into file
measurements_table = pd.DataFrame(
    {'file_name': file_name,
     'Ratio: area between real data curve and simulation curve to area under real data curve': areas_ratio,
     'Earth Movers Distance': [emd1],
     'Parameters': parameters.items()
     })

location = 'Measurements/' + 'measurements for ' + file_name + '.csv'
measurements_table.to_csv(location)

