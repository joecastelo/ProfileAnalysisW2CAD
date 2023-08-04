from tkinter import filedialog,Tk
from pylinac import FieldAnalysis
import matplotlib.pyplot  as plt
from scipy.interpolate import interp1d
import numpy as np

def openfile_dialog():
    root = Tk()
    root.withdraw()
    selected = filedialog.askopenfile()
    return selected

dcm = openfile_dialog()

content = dcm.readlines()


import pandas as pd 
import itertools

def get_data_between_header(linesneeded, header_elements):
    dc = {}
    for index,e in header_elements:
        try:
            indxs = [i for i,f in header_elements]
            between_header_elements_data_range = index, indxs[indxs.index(index)+1]
            bddr = between_header_elements_data_range ## just for simplification 
            data_for = linesneeded[bddr[0] + 1 : bddr[1]]
            ## test data2cad
            if not e in dc.keys():
                dc[e] = data_for
            else:
                dc[e].extend(data_for)
        except:
            print('all elements are analysed')
    return dc
def cleandata(content):
    linesneeded = [l for l in content if ('<' in l) or ('%DPTH'in l ) or ("%FLSZ" in l)]
    field_sizes = [(i,f) for i,f in enumerate(linesneeded) if '%FLSZ' in f ]
    fs_dict = get_data_between_header(linesneeded, field_sizes)
    dc1 = {'fieldsize': [], 'data': []}
    ## clean the data for each field
    
    for i,fs in enumerate(fs_dict.keys()):
        linesneeded = fs_dict[fs]

        depths = [d for d in linesneeded if '%DPTH' in d]
        # if i == 0:
        #     print(linesneeded)
        dc={}
        print(i,fs)
        for d in depths:
            
            try:
                inner_dc= {'d_from_cax': [], 'dose' : []}
                next_s = depths[depths.index(d) + 1]
                
                between_depths_data_range = linesneeded.index(d), linesneeded.index(next_s)
                bddr = between_depths_data_range ## just for simplification 
                data_for_depth = linesneeded[bddr[0] + 1 : bddr[1]]
                ## test data2cad
                data = [parser_forw2cad(l) for l in data_for_depth]
                inner_dc['d_from_cax'].append([d[0] for d in data])
                inner_dc['dose'].append([d[1] for d in data])
                dc[d] = inner_dc
            except:
                print('all elements are analysed')
        dc1['fieldsize'].append(fs)
        dc1['data'].append(dc)
    return pd.DataFrame(dc1)
def parser_forw2cad(l):
    cleanse = l.strip('<').strip('>\n').split(' ')
    return (float(cleanse[0]),float(cleanse[-1]))
df = cleandata(content)
df.head(10)


depth = '013' ## dmax
p_depth  = f'%DPTH {depth}\n'
x = list(df[df['fieldsize'] == '%FLSZ 100*100\n']['data'].iloc[0][p_depth]['d_from_cax'][0])
xr = np.arange(min(x), max(x), .25)

y = list(df[df['fieldsize'] == '%FLSZ 100*100\n']['data'].iloc[0][p_depth]['dose'][0])
lin = interp1d(x,y)
yr = list(lin(xr))
plt.scatter(xr,yr)



def findnearest_index(arr, value):
    # calculate the difference array
    difference_array = np.absolute(np.array(arr)-value)

    # find the index of minimum element from the array
    index = difference_array.argmin()
    return index

xr = list(xr)
# element to which nearest value is to be found
half = 50 
x_positive = [s for s in xr if s >= 0]
indexes_positive = xr.index(x_positive[0]), xr.index(x_positive[-1])
y_positive = yr[indexes_positive[0]:indexes_positive[-1] + 1]
x_negative = [s for s in xr if s <= 0]
indexes_negative = xr.index(x_negative[0]), xr.index(x_negative[-1])
y_negative = yr[indexes_negative[0]:indexes_negative[-1] +1]
plt.scatter(x_positive,y_positive)
plt.scatter(x_negative,y_negative)

#### Symm Calculation

half_neg = findnearest_index(y_negative, 50)
half_pos = findnearest_index(y_positive, 50)
print("Y = 50% at positive and negative : ", max(x_positive[0:half_pos]), min(x_negative[half_neg::]))
area_pos = np.trapz(y = y_positive[0:half_pos], x = x_positive[0:half_pos])
area_neg = np.trapz(y = y_negative[half_neg::], x = x_negative[half_neg::])
sym = 100*(area_neg - area_pos)/(area_pos + area_neg)
print(f"Beam symmetry = {sym} %")
print(f"Right Area = {area_pos}")
print(f"Left Area = {area_neg}")


plt.scatter(xr,yr)
plt.axvline(x = x_positive[half_pos])
plt.axvline(x = x_negative[half_neg])
plt.axvline(x = 0)



#### Flatness

depth = '100' ## 10 cm
p_depth  = f'%DPTH {depth}\n'
field_size = 300
x = list(df[df['fieldsize'] == f'%FLSZ {field_size}*{field_size}\n']['data'].iloc[0][p_depth]['d_from_cax'][0])
xr = np.arange(min(x), max(x), .25)

y = list(df[df['fieldsize'] == f'%FLSZ {field_size}*{field_size}\n']['data'].iloc[0][p_depth]['dose'][0])
lin = interp1d(x,y)
yr = list(lin(xr))
plt.scatter(xr,yr)

xr = list(xr)
# element to which nearest value is to be found
half = 50 
x_positive = [s for s in xr if s >= 0]
indexes_positive = xr.index(x_positive[0]), xr.index(x_positive[-1])
y_positive = yr[indexes_positive[0]:indexes_positive[-1] + 1]
x_negative = [s for s in xr if s <= 0]
indexes_negative = xr.index(x_negative[0]), xr.index(x_negative[-1])
y_negative = yr[indexes_negative[0]:indexes_negative[-1] +1]
plt.scatter(x_positive,y_positive)
plt.scatter(x_negative,y_negative)

half_neg = findnearest_index(x_negative, -1*40*field_size/100)
half_pos = findnearest_index(x_positive, 40*field_size/100)
print( half_neg , half_pos)
profile_80 = y_positive[0:half_pos]
profile_80_n = y_negative[half_neg::]
maxes = max(profile_80), max(profile_80_n)
mins = min(profile_80), min(profile_80_n)
print("X = 80% at positive and negative : ", max(maxes), min(mins))

flat = 100*(max(maxes)- min(mins))/(max(maxes)+ min(mins))
print(f"Beam flatness = {flat} %")


plt.scatter(xr,yr)
plt.scatter(x_positive[0:half_pos], profile_80)
plt.scatter(x_negative[half_neg::], profile_80_n)

plt.axvline(x = x_positive[half_pos])
plt.axvline(x = x_negative[half_neg])
