# -*- coding: utf-8 -*-
from pda import PSDAnalyser
import numpy
# import os

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
bins_end   = [i for i,el in enumerate(bins_row) if not numpy.isnan(bins_row[i]) and numpy.isnan(bins_row[i+1])][0]

''' Get all sets of bin sizes '''
bin_sets = []
for y in y_ind:
    bins = [el for el in xldata.iloc[y][bins_start:bins_end+1]]
    bin_sets.append(bins)

''' Grab all rows with data '''
name_rows = [(i,el) for i,el in enumerate(xldata['Sample Name']) if type(el) == str]

''' Build all data sets '''
y_ind_plus = [el for el in y_ind]
y_ind_plus.append(len(xldata['Sample Name'])+1)
for yi in y_ind_plus[:-1]:
    print(yi)
    