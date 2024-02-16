import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import math


def fill_axes_for_subplot(axes, row, col, x_values, y_values, plot_title):
    """function filling plots of a subplot depending on number of rows and columns of subplot"""
    def fill_axes_for_subplot_1d(row_or_col):
        """ function filling plots of a subplot if row number == 1 or column number == 1"""
        axes[row_or_col].plot(x_values, y_values)
        axes[row_or_col].set_title(plot_title)
        axes[row_or_col].set_xlabel("Time (ms)")
        axes[row_or_col].set_ylabel("ADC (AU)")

    if range_row == 1:
        fill_axes_for_subplot_1d(col)
    elif range_row == 1:
        fill_axes_for_subplot_1d(row)
    else:
        axes[row][col].plot(x_values, y_values)
        axes[row][col].set_title(plot_title)
        axes[row][col].set_xlabel("Time (ms)")
        axes[row][col].set_ylabel("ADC (AU)")

# path where file to open is stored
os.chdir("/home/milewski/these/mesures/")
# name of the file to open
data_file = "M8-1_G1_GUIDose22.2Gy.bin"

dt = np.dtype('uint16')

f = open(data_file, "rb")
data = f.read()
z_data = np.frombuffer(data, dt)

CPT = []
CPT1 = []
CPT2 = []
CPT3 = []
CPT4 = []
CPT5 = []
CPT6 = []
CPT7 = []
CPT8 = []

nb_event = np.size(z_data) // 27

# retrieves values of each strip and store them in lists
for num_event in range(nb_event):
    CPT1.append(np.uint32(((z_data[3 + 0 * 2 + num_event * 27]) << 16) + (z_data[4 + 0 * 2 + num_event * 27])))
    CPT2.append(np.uint32(((z_data[3 + 1 * 2 + num_event * 27]) << 16) + (z_data[4 + 1 * 2 + num_event * 27])))
    CPT3.append(np.uint32(((z_data[3 + 2 * 2 + num_event * 27]) << 16) + (z_data[4 + 2 * 2 + num_event * 27])))
    CPT4.append(np.uint32(((z_data[3 + 3 * 2 + num_event * 27]) << 16) + (z_data[4 + 3 * 2 + num_event * 27])))
    CPT5.append(np.uint32(((z_data[3 + 4 * 2 + num_event * 27]) << 16) + (z_data[4 + 4 * 2 + num_event * 27])))
    CPT6.append(np.uint32(((z_data[3 + 5 * 2 + num_event * 27]) << 16) + (z_data[4 + 5 * 2 + num_event * 27])))
    CPT7.append(np.uint32(((z_data[3 + 6 * 2 + num_event * 27]) << 16) + (z_data[4 + 6 * 2 + num_event * 27])))
    CPT8.append(np.uint32(((z_data[3 + 7 * 2 + num_event * 27]) << 16) + (z_data[4 + 7 * 2 + num_event * 27])))

# Subplots of 8 strips
# to plot less than 8 strips remove the strips not wanted in dic_strips_to_plot
dic_strips_to_plot = {'Strip 1 ': CPT1, 'Strip 2': CPT2, 'Strip 3': CPT3, 'Strip 4': CPT4, 'Strip 5': CPT5,
                      'Strip 6': CPT6, 'Strip 7': CPT7, 'Strip 8': CPT8}
range_row = math.ceil(len(dic_strips_to_plot) / 2)
range_col = 2

if range_row == 1:
    fig_height_inches = 5
else:
    fig_height_inches = 12
fig, axes = plt.subplots(range_row, range_col, figsize=(16, fig_height_inches))
fig.suptitle('ADC response')
index_strip = 0

# loop filling each space of the subplot with a strip response
for row in range(0, range_row):
    for col in range(0, range_col):
        if len(dic_strips_to_plot) == index_strip:
            break
        x_values = [i for i in range(nb_event)]
        y_values = list(dic_strips_to_plot.values())[index_strip]
        plot_title = list(dic_strips_to_plot.keys())[index_strip]
        fill_axes_for_subplot(axes, row, col, x_values, y_values, plot_title)
        plt.subplots_adjust(top=0.9, wspace = 0.4, hspace= 0.6)
        index_strip += 1

# suppresses last plot in subplot if number of plots odd (as number of columns = 2)
if len(dic_strips_to_plot) % 2 != 0:
    if range_row == 1:
        fig.delaxes(axes[range_col - 1])
    elif range_col == 1:
        fig.delaxes(axes[range_row - 1])
    else:
        fig.delaxes(axes[range_row - 1][range_col - 1])

# # 8 strips on same plot
# to plot less than 8 strips remove the strips not wanted in dic_list_strips_to_plot
dic_strips_to_plot = {'Strip 2': CPT2, 'Strip 3': CPT3, 'Strip 5': CPT5, 'Strip 6': CPT6,
                      'Strip 7': CPT7,'Strip 8': CPT8}
plt.figure(figsize=(14, 7))
for strip in dic_strips_to_plot.keys():
    plt.title('ADC response')
    figure = plt.plot([i for i in range(nb_event)], dic_strips_to_plot[strip], label = strip)
    # plt.yticks(np.arange(0, 350, 10))
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    plt.xlabel("Time (ms)")
    plt.ylabel("AU")


# plot the mean 8 strips responses
# calculates the mean values of each event for the 8 strips and retrieves the mean values calculated in mean_strip_lists
# to remove 1 strip of the mean calculation : remove it from strips_array definition below
strips_array = np.array([CPT1, CPT2, CPT3, CPT4, CPT5, CPT6, CPT7, CPT8])
mean_strips_list = []
for event in range(0, nb_event):
    mean_strips_list.append(np.mean(strips_array[:, event]))

plt.figure(figsize=(14, 7))
plt.title('Mean ADC response')
figure = plt.plot([i for i in range(nb_event)], mean_strips_list)
plt.xlabel("Time (ms)")
plt.ylabel("AU")

plt.show()
