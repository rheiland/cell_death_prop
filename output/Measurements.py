from __future__ import division
import numpy as np
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as p
import csv
import scipy as sp
from scipy import stats
from scipy import interpolate
import cv2
from shapely.geometry import Polygon
from bs4 import BeautifulSoup
import pandas as pd
import os
from datetime import datetime

path = 'C:/Users/cheng/source/repos/cell_death_prop/2/output/Measurements/'


class Measurements:
    def __init__(self, folder_name):
        self.folder_path = path + folder_name
        # create the measurements folder
        if not os.path.exists(self.folder_path):
            os.makedirs(self.folder_path)
        # init variables
        self.parameters = {}
        self.num_of_cells = len(open("dead_cells.dat").readlines())
        self.sim_max_time_of_death = 0
        self.real_max_time_of_death = 0
        self.frac_died_until_time_sim = {}
        self.frac_died_until_time_real = {}
        self.frac_died_until_time_sim_time_scale = {}
        self.sim_time_of_death_per_index = {}
        self.real_time_of_death_per_index = {}
        self.position_per_index = {}
        self.time_scale = 0
        self.file_name = ''
        self.dataX_sim = []
        self.dataX_real = []
        self.dataY_sim = []
        self.dataY_real = []
        self.all_times = []
        self.dataX_sim_normal = []
        self.dataX_real_normal = []
        self.total_area = 0
        self.areas_ratio = 0
        self.area_below_curve = 0
        self.diffs_in_plot = {}
        self.sim_neighbors_time_diffs = []
        self.real_neighbors_time_diffs = []
        self.normalized_sim_neighbors_time_diffs = []
        self.normalized_real_neighbors_time_diffs = []
        self.diff_vs_dist_sim = {}
        self.diff_vs_dist_real = {}
        self.correlation_sim = 0
        self.correlation_real = 0
        self.count_sim_diffs = {}
        self.count_real_diffs = {}
        self.sim_neighbors_time_diffs_list = []
        self.real_neighbors_time_diffs_list = []
        self.count_sim_diffs_time_scale = {}
        self.sim_keys = []
        self.sim_vals = []
        self.keys_real = []
        self.vals_real = []
        self.sim_vals_fraction = []
        self.real_vals_fraction = []
        self.median_sim = 0
        self.median_real = 0
        self.median_normalized_sim = 0
        self.median_normalized_real = 0
        self.emd1 = 0
        self.emd2 = 0
        self.y_pos = []
        self.y_pos_real = []

    def read_files(self):
        self.read_from_config()
        self.read_sim_dead_cells_info()
        self.read_exp_dead_cells_info()
        self.read_neighbors_file()

    # read the fields values from the configuration file
    def read_from_config(self):
        in_configuration_file = open("configuration.xml", "r")
        contents = in_configuration_file.read()
        soup = BeautifulSoup(contents, 'xml')
        self.parameters = {'threshold': soup.find_all('death_threshold')[0].get_text(),
                           'initial substrate': soup.find_all('initial_internalized_substrate')[0].get_text(),
                           'absorption rate': soup.find_all('signal_internalization_rate')[0].get_text(),
                           'df coefficient': soup.find_all('signal_diffusion_coefficient')[0].get_text()}

    # read the time of death for each cell from output simulation file
    def read_sim_dead_cells_info(self):
        accumulate_sim = 0
        with open("dead_cells.dat", "r") as file2:
            lines = file2.readlines()
            last = lines[-1]
            for line in lines:
                words = line.split()
                accumulate_sim += 1
                self.frac_died_until_time_sim[int(words[0])] = accumulate_sim / self.num_of_cells
                self.sim_time_of_death_per_index[int(words[2])] = int(words[0])
                self.position_per_index[int(words[2])] = words[3]
                self.real_time_of_death_per_index[int(words[2])] = int(words[1])
                if line is last:
                    self.sim_max_time_of_death = int(words[0])

    # read the time of death for each cell from output experiment file
    def read_exp_dead_cells_info(self):
        accumulate_real = 0
        with open("exp_a.csv", "r") as real_data:
            lines = real_data.readlines()
            last = lines[-1]
            # skip the header row
            for line in lines:
                if line == lines[0]:
                    continue
                words = line.split(',')
                accumulate_real += 1
                self.frac_died_until_time_real[int(words[3])] = accumulate_real / self.num_of_cells
                if int(words[4]) != 0 and self.time_scale == 0:
                    self.time_scale = int(int(words[3]) / int(words[4]))
                if line is last:
                    self.real_max_time_of_death = int(words[3])
                    self.file_name = words[6]

    # calculate time of death differences between neighbor cells
    def read_neighbors_file(self):
        with open("exp_b.csv", "r") as real_data_neighbors:
            csv_reader_neighbors = csv.reader(real_data_neighbors, delimiter=',')
            # skip the header
            next(csv_reader_neighbors)
            for line in csv_reader_neighbors:
                cell_index = int(line[0])
                neighbors = line[1].split(',')
                for nei in neighbors:
                    if cell_index in self.sim_time_of_death_per_index.keys() and int(
                            nei) in self.sim_time_of_death_per_index.keys():
                        diff = np.abs(self.sim_time_of_death_per_index[cell_index] - self.sim_time_of_death_per_index[int(nei)])
                        self.sim_neighbors_time_diffs.append(diff)
                        distX = np.abs(float(self.position_per_index[cell_index].split(',')[0]) - float(
                            self.position_per_index[int(nei)].split(',')[0]))
                        distY = np.abs(float(self.position_per_index[cell_index].split(',')[1]) - float(
                            self.position_per_index[int(nei)].split(',')[1]))
                        dist = np.sqrt(distX ** 2 + distY ** 2)
                        self.diff_vs_dist_sim[dist] = diff
                        self.normalized_sim_neighbors_time_diffs.append(diff / dist)
                        diff = np.abs(self.real_time_of_death_per_index[cell_index] - self.real_time_of_death_per_index[int(nei)])
                        self.real_neighbors_time_diffs.append(diff)
                        self.normalized_real_neighbors_time_diffs.append(diff / dist)
                        self.diff_vs_dist_real[dist] = diff

    # accumulate all the cells that died between minute i to minute 30i (i runs from 0 to simulation time divided by 30)
    def create_fake_points_every_30(self):
        self.all_times = [x for x in range(0, self.real_max_time_of_death + self.time_scale, self.time_scale)]
        last_k = 0
        last_x = 0
        for x in self.all_times:
            frac_accumulated = 0
            for k, v in self.frac_died_until_time_sim.items():
                if k < last_k:
                    continue
                if k <= x:
                    frac_accumulated = v
                if k > x or k == self.real_max_time_of_death:
                    if frac_accumulated == 0:
                        frac_accumulated = self.frac_died_until_time_sim_time_scale[last_x]
                    self.frac_died_until_time_sim_time_scale[x] = frac_accumulated
                    last_k = k
                    last_x = x
                    break

    # divide x values in max_time_of_death to normalize both sim and real times of death
    def normalize_times_of_death(self):
        # sim
        self.dataX_sim = list(self.frac_died_until_time_sim_time_scale.keys())
        for x in self.dataX_sim:
            self.dataX_sim_normal.append(x / self.real_max_time_of_death)
        self.dataY_sim = list(self.frac_died_until_time_sim_time_scale.values())
        # real
        self.dataX_real = list(self.frac_died_until_time_real.keys())
        for x in self.dataX_real:
            self.dataX_real_normal.append(x / self.real_max_time_of_death)
        self.dataY_real = list(self.frac_died_until_time_real.values())

    def first_graph_area_ratio(self):
        self.calc_area_between_curves()
        line_up, = plt.plot(self.dataX_real_normal, self.dataY_real)
        line_up.set_label('real experiment')
        line_down, = plt.plot(self.dataX_sim_normal, self.dataY_sim)
        line_down.set_label('simulation')
        plt.legend(handles=[line_up, line_down], fontsize='20')
        plt.title('Fraction of dead cells over time', fontsize='40')
        plt.xlabel('time (min)', color='blue', fontsize='20')
        plt.ylabel('fraction', color='blue', fontsize='20')
        plt.xlim(0, 1)
        plt.xticks(np.arange(0, 1.1, 0.1), fontsize='20')
        plt.ylim(0, 1)
        plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
        plt.text(1, 1, 'Ratio of area\nbetween curves=\n' + str(self.areas_ratio), fontsize='20', color='r')
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(29, 27)
        plt.savefig(self.folder_path + '/Fraction of dead cells over time', dpi=50)
        plt.close()

    # receives 2 points, returns the straight line equation
    @staticmethod
    def line(p1, p2):
        a = (p1[1] - p2[1])
        b = (p2[0] - p1[0])
        c = (p1[0] * p2[1] - p2[0] * p1[1])
        return a, b, -c

    # returns the intersection point of 2 points and False if they don't intersect
    @staticmethod
    def intersection(l1, l2):
        d = l1[0] * l2[1] - l1[1] * l2[0]
        dx = l1[2] * l2[1] - l1[1] * l2[2]
        dy = l1[0] * l2[2] - l1[2] * l2[0]
        if d != 0:
            x = dx / d
            y = dy / d
            return x, y
        else:
            return False

    # Dividing the area between curves into polygons according to intersections, calculating each polygon's area and summing all areas
    def calc_area_between_curves(self):
        area = 0
        for i in range(len(self.dataX_sim_normal) - 1):
            a = (self.dataX_sim_normal[i], self.dataY_sim[i])
            b = (self.dataX_real_normal[i], self.dataY_real[i])
            c = (self.dataX_real_normal[i + 1], self.dataY_real[i + 1])
            d = (self.dataX_sim_normal[i + 1], self.dataY_sim[i + 1])
            L1 = self.line(list(a), list(d))
            L2 = self.line(list(b), list(c))
             # intersection detected
            inter = self.intersection(L1, L2)
            if inter:
                if (inter[0] >= a[0]) and (inter[0] <= d[0]):
                    p1 = Polygon([a, b, inter])
                    p2 = Polygon([inter, c, d])
                    area = p1.area + p2.area
                else:
                    p1 = Polygon([a, b, c, d])
                    area = p1.area
            self.total_area += area

        self.area_below_curve = np.trapz(self.dataX_real_normal, self.dataY_real)
        self.areas_ratio = self.total_area / self.area_below_curve

    def create_dead_cells_by_fractions(self):
        temp_diff = {}
        prev = 0
        for (time1, frac1) in self.frac_died_until_time_sim.items():
            if frac1 > 0.75:
                temp_diff[0.75] = prev
                break
            elif (frac1 > 0.5) & (0.5 not in temp_diff):
                temp_diff[0.5] = prev
            elif (frac1 > 0.25) & (0.25 not in temp_diff):
                temp_diff[0.25] = prev
            prev = time1
        prev = 0
        for (time2, frac2) in self.frac_died_until_time_real.items():
            if frac2 > 0.75:
                self.diffs_in_plot[0.75] = temp_diff[0.75] - prev
                break
            elif (frac2 > 0.5) & (0.5 not in self.diffs_in_plot):
                self.diffs_in_plot[0.5] = temp_diff[0.5] - prev
            elif (frac2 > 0.25) & (0.25 not in self.diffs_in_plot):
                self.diffs_in_plot[0.25] = temp_diff[0.25] - prev
            prev = time2

    def second_graph_dead_cells_by_fractions(self):
        self.create_dead_cells_by_fractions()
        plt.plot(list(self.diffs_in_plot.keys()), list(self.diffs_in_plot.values()))
        plt.title('Difference in time of 0.25, 0.5 and 0.75 dead cells', fontsize='40')
        plt.xlabel('fraction', fontsize='20', color='r')
        plt.ylabel('time difference', fontsize='20', color='r')
        plt.xticks(np.array(list(self.diffs_in_plot.keys())), fontsize='20')
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(29, 27)
        plt.savefig(self.folder_path + '/Difference in time of 025, 05 and 075 dead cells', dpi=50)
        plt.close()

    # create a matrix of values for sim heat map
    # go over diff vs dist sim
    def create_matrix_for_heat_map(self, skips_x, skips_y, diff_vs_dist, max_time_of_death):
        max_dist = int(max(diff_vs_dist.keys()) + 1)
        y_bins = np.array(range(0, max_time_of_death + skips_x, skips_x))
        x_bins = np.array(range(0, max_dist + skips_y, skips_y))
        matrix_rows = len(y_bins)
        matrix_cols = len(x_bins)
        matrix = [[0 for x in range(matrix_cols - 1)] for y in range(matrix_rows - 1)]
        data_frame = pd.DataFrame(columns=["x", "y"])
        data_frame["x"] = diff_vs_dist.keys()
        data_frame["y"] = diff_vs_dist.values()
        tuples_x = []
        tuples_y = []
        x_labels = []
        y_labels = []
        for i in range(matrix_rows - 1):
            tuples_y.append((y_bins[i], y_bins[i + 1]))
            y_labels.append(str((y_bins[i], y_bins[i + 1])))
        for j in range(matrix_cols - 1):
            tuples_x.append((x_bins[j], x_bins[j + 1]))
            x_labels.append(str((x_bins[j], x_bins[j + 1])))

        tuples_combined = [[x, y] for x in tuples_x for y in tuples_y]

        # go over data frame and count points that fall within each square
        for ranges in tuples_combined:
            size = data_frame[((ranges[0][0] <= data_frame["x"]) & (data_frame["x"] <= ranges[0][1])) &
                              ((ranges[1][0] <= data_frame["y"]) & (data_frame["y"] <= ranges[1][1]))].shape[0]
            matrix[(ranges[1][0]) // skips_x][ranges[0][0] // skips_y] = size
        return matrix, x_labels, y_labels

    # remove information of cells that are not actually neighbors, trim most right cols and most up rows
    def reduce_matrix_size(self, matrix, threshold_row, threshold_col):
        first_matrix = []
        for row in matrix:
            first = row[:(len(matrix[0]) // 2)]
            first_matrix.append(first)

        matrix_trim_rows = []
        found_row = False
        # trim less than 2 sum rows from above
        for row in reversed(first_matrix):
            if sum(row) > threshold_row:
                found_row = True
            if found_row:
                matrix_trim_rows.append(row)

        matrix_trim_rows.reverse()
        first_matrix = matrix_trim_rows

        # trim less than 2 sum cols from the left
        found_col = False
        matrix_trim_cols = []
        for row in reversed(np.transpose(first_matrix)):
            if sum(row) > threshold_col:
                found_col = True
            if found_col:
                matrix_trim_cols.append(row)

        matrix_trim_cols.reverse()
        first_matrix = np.transpose(matrix_trim_cols)
        first_matrix_rows = len(first_matrix)
        first_matrix_cols = len(first_matrix[0])
        return first_matrix, first_matrix_rows, first_matrix_cols

    # create a heat map show casting the correlation between:
    # neighbor cells time of death differences (propagation time) and distance between neighbors
    def create_heat_map(self, x_labels, y_labels, matrix, matrix_rows, matrix_cols, sim_or_real):
        x_labels_reduced = x_labels[:(len(x_labels) + 1) // 2]
        fig, ax = plt.subplots()
        matsum = np.sum(matrix)
        reds = cm.get_cmap('Reds', 19)
        cmap = colors.ListedColormap(reds(range(19)))
        bounds = [-10, 0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
        bounds = [x / matsum for x in bounds]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        array_norm = [x / matsum for x in matrix]
        im = ax.imshow(np.array(array_norm), interpolation="none", cmap=cmap, norm=norm)
        ax.set_xlabel('Distance', fontsize='20', color='red')
        ax.set_ylabel('Time of death differences', fontsize='20', color='red')
        ax.set_xticks(np.arange(matrix_cols))
        ax.set_yticks(np.arange(matrix_rows))
        ax.set_xticklabels(x_labels_reduced)
        ax.set_yticklabels(y_labels)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        for i in range(matrix_rows):
            for j in range(matrix_cols):
                text = ax.text(j, i, round((matrix[i][j] / matsum), 3), ha="center", va="center", color="black", fontsize='15')

        ax.set_title(sim_or_real + ' - Neighbor cells time of death difference vs. distance', fontsize='30')
        plt.gca().invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.5)
        fig.set_size_inches(29, 27, forward=True)
        plt.colorbar(im, cax=cax)
        fig.savefig(self.folder_path + '/' + sim_or_real + ' - Heatmap of cells time of death difference vs distance', dpi=100, pad_inches=1.0)
        plt.close()

    # call all relevant functions for creating the heat maps
    def create_all_heat_maps(self):
        matrix_sim_5, sim_5_x_labels, sim_5_y_labels = self.create_matrix_for_heat_map(5, 20, self.diff_vs_dist_sim, self.sim_max_time_of_death)
        matrix_sim_30, sim_30_x_labels, sim_30_y_labels = self.create_matrix_for_heat_map(30, 20, self.diff_vs_dist_sim, self.sim_max_time_of_death)
        matrix_real, real_x_labels, real_y_labels = self.create_matrix_for_heat_map(30, 20, self.diff_vs_dist_real, self.real_max_time_of_death)
        reduced_matrix_sim_5, reduced_matrix_sim_5_rows, reduced_matrix_sim_5_cols = self.reduce_matrix_size(matrix_sim_5, 9, 9)
        reduced_matrix_sim_30, reduced_matrix_sim_30_rows, reduced_matrix_sim_30_cols = self.reduce_matrix_size(matrix_sim_30, 4, 4)
        reduced_matrix_real, reduced_matrix_real_rows, reduced_matrix_real_cols = self.reduce_matrix_size(matrix_real, 11, 11)
        self.create_heat_map(sim_5_x_labels, sim_5_y_labels, reduced_matrix_sim_5, reduced_matrix_sim_5_rows,
                        reduced_matrix_sim_5_cols, 'Simulation')
        self.create_heat_map(sim_30_x_labels, sim_30_y_labels, reduced_matrix_sim_30, reduced_matrix_sim_30_rows,
                        reduced_matrix_sim_30_cols, 'Simulation')
        self.create_heat_map(real_x_labels, real_y_labels, reduced_matrix_real, reduced_matrix_real_rows,
                        reduced_matrix_real_cols, 'Experiment')


    def count_differences(self):
        # get rid of duplicates
        sim_neighbors_time_diffs_set = set(self.sim_neighbors_time_diffs)
        self.sim_neighbors_time_diffs_list = list(sim_neighbors_time_diffs_set)
        real_neighbors_time_diffs_set = set(self.real_neighbors_time_diffs)
        self.real_neighbors_time_diffs_list = list(real_neighbors_time_diffs_set)

        for nei in self.real_neighbors_time_diffs_list:
            if nei not in self.sim_neighbors_time_diffs_list:
                self.sim_neighbors_time_diffs_list.append(nei)

        self.sim_neighbors_time_diffs_list.sort()
        self.real_neighbors_time_diffs_list.sort()

        for nei in self.sim_neighbors_time_diffs_list:
            str_nei = str(nei)
            if str_nei not in self.count_sim_diffs:
                self.count_sim_diffs[str_nei] = self.sim_neighbors_time_diffs.count(nei)

        # accumulate
        accu_val = 0
        for k, v in self.count_sim_diffs.items():
            accu_val += v
            if int(k) % self.time_scale == 0:
                self.count_sim_diffs_time_scale[k] = accu_val
                accu_val = 0

        for nei in self.real_neighbors_time_diffs_list:
            str_nei = str(nei)
            if str_nei not in self.count_real_diffs:
                self.count_real_diffs[str_nei] = self.real_neighbors_time_diffs.count(nei)

        for nei in self.count_sim_diffs_time_scale.keys():
            if nei not in self.count_real_diffs.keys():
                self.count_real_diffs[nei] = 0

        list(self.count_real_diffs.keys()).sort()
        self.sim_keys = list(self.count_sim_diffs_time_scale.keys())
        self.sim_vals = list(self.count_sim_diffs_time_scale.values())
        self.keys_real = list(self.count_real_diffs.keys())
        self.vals_real = list(self.count_real_diffs.values())

    def vals_fraction(self, vals):
        sum_all = sum(vals)
        return [(val / sum_all) for val in vals]

    def third_graph(self):
        self.count_differences()
        self.y_pos = np.arange(len(self.sim_keys))
        plt.subplot(2, 1, 1)
        width = 0.5
        self.sim_vals_fraction = self.vals_fraction(self.sim_vals)
        self.median_sim = np.median(self.sim_vals_fraction)
        sim_bar = plt.bar(self.y_pos, self.sim_vals_fraction, width, color='green')
        plt.xticks(self.y_pos, np.array(self.sim_keys), fontsize='20')
        plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
        plt.title('Distribution of time of death differences between neighbors- simulation', fontsize='40')
        plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
        plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')

        plt.subplot(2, 1, 2)
        self.real_vals_fraction = self.vals_fraction(self.vals_real)
        self.median_real = np.median(self.real_vals_fraction)
        self.y_pos_real = np.arange(len(self.keys_real))
        real_bar = plt.bar(self.y_pos_real, self.real_vals_fraction, width, color='blue')
        plt.title('Distribution of time of death differences between neighbors- real experiment', fontsize='40')
        plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
        plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')
        plt.xticks(self.y_pos_real, np.array(self.keys_real), fontsize='20')
        plt.yticks(np.arange(0, 1.25, 0.25), fontsize='20')
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(29, 27)
        plt.savefig(self.folder_path + '/Distribution of time of death differences between neighbors', dpi=50)
        plt.close()
        # emd1
        self.emd1 = sum(abs(np.cumsum(self.sim_vals_fraction) - np.cumsum(self.real_vals_fraction)))

    def cumulative(self, vals_fraction, keys, sim_or_real, y_pos):
        cumulative = np.cumsum(vals_fraction)
        width = 0.5
        bar = plt.bar(y_pos, cumulative, width, color='green')
        plt.xticks(y_pos, np.array(keys), fontsize='20')
        plt.yticks(np.arange(0, max(cumulative) + 0.1, 0.25), fontsize='20')
        plt.title('Distribution of time of death differences between neighbors- ' + sim_or_real, fontsize='40')
        plt.xlabel('time of death differences', color='red', fontsize='20', horizontalalignment='center')
        plt.ylabel('frequency', color='red', fontsize='20', horizontalalignment='center')
        return bar

    def forth_graph(self):
        plt.subplot(2, 1, 1)
        sim_bar = self.cumulative(self.sim_vals_fraction, self.sim_keys, 'Simulation', self.y_pos)
        plt.subplot(2, 1, 2)
        real_bar = self.cumulative(self.real_vals_fraction, self.keys_real, 'Experiment', self.y_pos_real)
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(29, 27)
        plt.savefig(path + 'Distribution of time of death differences between neighbors- cummulative', dpi=50)
        plt.close()

        plt.subplot(2, 1, 1)
        plt.hist(self.normalized_sim_neighbors_time_diffs, density=True, bins=50)
        plt.title('Distribution of normalized time of death differences between neighbors- simulation', fontsize='40')
        plt.xlabel('time of death differences divided by distance', color='red', fontsize='20',
                   horizontalalignment='center')
        plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
        plt.xticks(np.arange(0.0, max(self.normalized_real_neighbors_time_diffs), 0.2), fontsize='20')
        plt.yticks(np.arange(0, 10, 2), fontsize='20')
        self.median_normalized_sim = np.median(self.normalized_sim_neighbors_time_diffs)
        self.median_normalized_real = np.median(self.normalized_real_neighbors_time_diffs)

        plt.subplot(2, 1, 2)
        plt.hist(self.normalized_real_neighbors_time_diffs, density=True, bins=50)
        plt.title('Distribution of normalized time of death differences between neighbors- real experiment',
                  fontsize='40')
        plt.xlabel('time of death differences divided by distance', color='red', fontsize='20',
                   horizontalalignment='center')
        plt.ylabel('fraction', color='red', fontsize='20', horizontalalignment='center')
        plt.xticks(np.arange(0.0, max(self.normalized_real_neighbors_time_diffs), 0.2), fontsize='20')
        plt.yticks(np.arange(0, 10, 2), fontsize='20')
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(29, 27)
        plt.savefig(self.folder_path + '/Distribution of normalized time of death differences between neighbors', dpi=50)
        plt.close()
        # earth movers distance- normalized
        self.emd2 = sum(abs(np.cumsum(self.normalized_sim_neighbors_time_diffs) - np.cumsum(self.normalized_real_neighbors_time_diffs)))
        print("the EMD of the difference between time of death of each cell and its neighbors- normalized: ", self.emd2)

    def write_to_file(self):
        # xml parameters
        first_param = list(self.parameters.items())[0]
        second_param = list(self.parameters.items())[1]
        third_param = list(self.parameters.items())[2]
        forth_param = list(self.parameters.items())[3]
        median = np.abs(self.median_sim.round(5) - self.median_real.round(5)) / self.median_real.round(5)
        median_normalized = np.abs(self.median_normalized_sim.round(5) - self.median_normalized_real.round(5)) / \
                            self.median_normalized_real.round(5)
        # Write into file
        file_name = self.file_name[0:-1]
        measurements_table = pd.DataFrame(
            {'File_name': file_name,
             'Areas Ratio': self.areas_ratio,
             'Earth Movers Distance': [self.emd1],
             'Median of time diff': median,
             'Median of normalized time diff': median_normalized,
             'Correlation - distance and time of death difference - sim': self.correlation_sim[0],
             'Correlation - distance and time of death difference - real': self.correlation_real[0],
             first_param[0]: first_param[1],
             second_param[0]: second_param[1],
             third_param[0]: third_param[1],
             forth_param[0]: forth_param[1],
             })
        date_and_time = datetime.now()
        dt_string = date_and_time.strftime("%d-%m-%Y-%H%M")
        location = self.folder_path + 'measurements for ' + file_name + ' ' + dt_string + '.csv'
        measurements_table.to_csv(location)

    # pearson correlation for both simulation and experiment - btw time difference and distance
    def correlation(self):
        self.correlation_sim = sp.stats.pearsonr(list(self.diff_vs_dist_sim.keys()),
                                                 list(self.diff_vs_dist_sim.values()))
        self.correlation_real = sp.stats.pearsonr(list(self.diff_vs_dist_real.keys()),
                                                  list(self.diff_vs_dist_real.values()))
        return self.correlation_sim, self.correlation_real


# change according to experiment
exp = 11
set_num = 4

measurements_obj = Measurements(exp + '/Set' + set_num + '/')
measurements_obj.read_files()
measurements_obj.create_fake_points_every_30()
measurements_obj.normalize_times_of_death()
measurements_obj.first_graph_area_ratio()
measurements_obj.second_graph_dead_cells_by_fractions()
measurements_obj.create_all_heat_maps()
correlation_sim, correlation_real = measurements_obj.correlation()
print("correlation_sim: ", correlation_sim)
print("correlation_real: ", correlation_real)
measurements_obj.third_graph()
measurements_obj.forth_graph()
measurements_obj.write_to_file()
