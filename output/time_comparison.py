import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import csv

# time_comparison
diffs = []
sum = 0
file = open("dead_cells.dat", "r")
for line in file:
    words = line.split()
    diff = (np.abs(int(words[0]) - int(words[1])))   # (|t1-t2|)
    diffs.append([diff, words[2]])

for i in diffs:
    print("Cell index: " + str(i[1]) + ", time difference: " + str(i[0]))
    sum += i[0]
print("Total difference: " + str(sum))

# frac_of_dead_cells
num_of_cells = len(open("dead_cells.dat").readlines())

# simulation
file2 = open("dead_cells.dat", "r")
accumulate_sim = 0
frac_died_until_time_sim = {}
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
accumulate_real = 0
frac_died_until_time_real = {}
with open("test4.csv", "r") as real_data:
    csv_reader = csv.reader(real_data, delimiter=',')
    # skip the header row
    next(csv_reader)
    for line in csv_reader:
        accumulate_real += 1
        frac_died_until_time_real[int(line[3])*15] = accumulate_real / num_of_cells

dataX_sim = list(frac_died_until_time_sim.keys())
dataY_sim = list(frac_died_until_time_sim.values())

dataX_real = list(frac_died_until_time_real.keys())
dataY_real = list(frac_died_until_time_real.values())

line_up, = plt.plot(dataX_real, dataY_real)
line_up.set_label('real experiment')
line_down, = plt.plot(dataX_sim, dataY_sim)
line_down.set_label('simulation')
plt.legend(handles=[line_up, line_down], fontsize='20')
plt.title('Fraction of dead cells over time', fontsize='40')
plt.xlabel('time (min)', color='blue', fontsize='20')
plt.ylabel('fraction', color='blue', fontsize='20')
plt.xlim(0, 225)
plt.xticks(fontsize='20')
plt.ylim(0, 1)
plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')

# calculate field between lines
field_real = np.trapz(dataY_real, dataX_real)
field_sim = np.trapz(dataY_sim, dataX_sim)
field_diff = field_real - field_sim
print(field_real, field_sim)
plt.text(max(dataX_real)-8, max(dataY_real)-0.05, 'Area between curves=\n' + str(field_diff), fontsize='20', color='r')


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
plt.yticks(np.arange(0, 50, 5), fontsize='20')
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
with open("test4b.csv", "r") as real_data_neighbors:
    csv_reader_neighbors = csv.reader(real_data_neighbors, delimiter=',')
    # skip the header
    next(csv_reader_neighbors)
    for line in csv_reader_neighbors:
        cell_index = int(line[0])
        neighbors = line[1].split(',')
        for nei in neighbors:
            diff = np.abs(sim_time_of_death_per_index[cell_index] - sim_time_of_death_per_index[int(nei)])
            sim_neighbors_time_diffs.append(diff)
            distX = float(position_per_index[cell_index].split(',')[0]) - float(position_per_index[int(nei)].split(',')[0])
            distY = float(position_per_index[cell_index].split(',')[1]) - float(position_per_index[int(nei)].split(',')[1])
            dist = np.sqrt(distX**2 + distY**2)
            normalized_sim_neighbors_time_diffs.append(diff/dist)
            diff = np.abs(real_time_of_death_per_index[cell_index] - real_time_of_death_per_index[int(nei)])
            real_neighbors_time_diffs.append(diff)
            normalized_real_neighbors_time_diffs.append(diff/dist)

sim_neighbors_time_diffs.sort()
real_neighbors_time_diffs.sort()

count_sim_diffs = {}
for x in sim_neighbors_time_diffs:
    X = str(x)
    count_sim_diffs[X] = sim_neighbors_time_diffs.count(x)

count_real_diffs = {}
for x in real_neighbors_time_diffs:
    X = str(x)
    count_real_diffs[X] = real_neighbors_time_diffs.count(x)

sim_keys = list(count_sim_diffs.keys())
sim_vals = list(count_sim_diffs.values())
keys_real = list(count_real_diffs.keys())
vals_real = list(count_real_diffs.values())

keys = list(set(sim_keys + keys_real))
int_keys = []
for x in keys:
    int_keys.append(int(x))
int_keys.sort()
keys = []
for x in int_keys:
    keys.append(str(x))

y_pos = np.arange(len(sim_keys))
plt.subplot(2, 1, 1)
width = 0.5
sim_bar = plt.bar(y_pos, sim_vals, width, color='green')
plt.xticks(y_pos, np.array(sim_keys), fontsize='20')
plt.yticks(np.arange(0, max(sim_vals)+100, 100), fontsize='20')
plt.title('Distribution of time of death differences between neighbors- simulation', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')

plt.subplot(2, 1, 2)
y_pos_combined = np.arange(len(keys))
y_pos_real = np.arange(len(keys_real))
real_bar = plt.bar(y_pos_real, vals_real, width, color='blue')
plt.title('Distribution of time of death differences between neighbors- real experiment', fontsize='40')
plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(y_pos_real, np.array(keys_real), fontsize='20')
plt.yticks(np.arange(0, max(sim_vals)+100, 100), fontsize='20')
plt.show()

plt.subplot(2, 1, 1)
plt.hist(normalized_sim_neighbors_time_diffs, density=True, bins=50)
plt.title('Distribution of normalized time of death differences between neighbors- simulation', fontsize='40')
plt.xlabel('time of death differences divided by distance', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(np.arange(0.0, max(normalized_real_neighbors_time_diffs), 0.1), fontsize='20')
plt.yticks(np.arange(0, 45, 5), fontsize='20')

plt.subplot(2, 1, 2)
plt.hist(normalized_real_neighbors_time_diffs, density=True, bins=50)
plt.title('Distribution of normalized time of death differences between neighbors- real experiment', fontsize='40')
plt.xlabel('time of death differences divided by distance', color='red', fontsize='20', horizontalalignment='center')
plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
plt.xticks(np.arange(0.0, max(normalized_real_neighbors_time_diffs), 0.1), fontsize='20')
plt.yticks(np.arange(0, 45, 5), fontsize='20')
plt.show()

