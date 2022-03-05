# -*- coding: utf-8 -*-
from pda import PSDAnalyser
import numpy
import math
# import os
import matplotlib.pyplot as plt
import numpy as np



pd = PSDAnalyser()

file = "C:\\Users\\hughr\\MATLAB Drive\\PSD-Analyser\\MS2000 Export.xlsx"
# xldata = pd.load_excel(file, header = None)
xldata = pd.load_excel(file)

''' Check through keys '''
stub = 'Result Between'
x_ind, header = [(i,k) for i,k in enumerate(xldata.keys()) if stub in str(k)][0]
y_ind = [0]

''' Get list of columns containing bin sizes '''
bin_occurrences = [y+1 for y,el in enumerate(xldata[header]) if el == header]
if bin_occurrences:
    y_ind.extend(bin_occurrences)

''' Grab column indices of bin sizes '''
bins_row = xldata.iloc[0]
bins_start = x_ind
bins_end = [i for i,el in enumerate(bins_row) if not numpy.isnan(bins_row[i]) and numpy.isnan(bins_row[i+1])][0]

''' Get all sets of bin sizes '''
bin_sets = []
for y in y_ind:
    bins = [el for el in xldata.iloc[y][bins_start:bins_end+1]]
    bin_sets.append(bins)

''' Grab all rows with data '''
names_raw = {i:el for i,el in enumerate(xldata['Sample Name']) if str(el) != 'Sample Name'}

''' Reject any numpy "nan", i.e. empty cells in original Excel file '''
names = {}
for k,v in names_raw.items():
    i = k
    name = v
    if type(name) == str:
        names[k] = v
    elif (type(name) == float and not numpy.isnan(name)) or type(name) == int:
        names[k] = str(v)

''' Build all data sets '''
y_ind_plus = [el for el in y_ind]
y_ind_plus.append(len(xldata['Sample Name'])+2)
blocks = {}
for i,yi in enumerate(y_ind_plus[:-1]):
    # print(yi)
    # within_block[yi] = [(j,el) for j,el in names if j>yi and j<y_ind_plus[i+1]]
    blocks[yi] = [j for j,el in names.items() if j>yi and j<y_ind_plus[i+1]]
    # print(within_block)
    # print(len(within_block))

''' Grab each PSD '''

psds = {}
for i,(k,v) in enumerate(blocks.items()):
    for el in v:
        bins = bin_sets[i]
        psd_data = [el for el in xldata.iloc[el][bins_start+1:bins_end+1]]
        psds[el] = (bins, psd_data)


yval = 182
all_bins = psds[yval][0]
all_y = psds[yval][1]

all_widths = [all_bins[i+1]-el for i,el in enumerate(all_bins[:-1])]
all_x = [np.mean([all_bins[i+1],el]) for i,el in enumerate(all_bins[:-1])]

i1,i2 = 20,60
x = all_x[i1:i2]
y = all_y[i1:i2]
widths = all_widths[i1:i2]

fig, ax = plt.subplots()
# x = np.arange(len(data))
# for i in range(len(data[0])):
#     y = [d[i] for d in data]
#     b = ax.bar(x + i * dimw, y, dimw, bottom=0.001)

# ax.set_xscale('log')

# ax.set_xlabel('x')
# ax.set_ylabel('y')

# widths = [0.13*el for el in bins[1:]]
# widths = [bins]
ax.bar(x, height = y, width = widths, facecolor = 'white', edgecolor = 'black')
ax.scatter(x,y,marker = '*')
plt.show()
