import os
import sys
from textwrap import fill 
from matplotlib import colorbar
from matplotlib import colors
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


# param = "param-" + angle + "-" + J + "-" + DMI +"-"+ J_twist + "-" + J_prist + "-" + J_intra + "-" + DMI_sub
path = os.getcwd() +   "/cells-coords.cfg"

try:
    open(path)
except:
    print("cannot find " + path)
    

print(path)
x_pos_bottom = []
y_pos_bottom = []
mc_1_avg_bottom = []
mc_1_avg_bottom_x = []
mc_1_avg_bottom_y = []

x_width = 60 # np.max(x_pos_bottom)
y_width = 60 #np.max(y_pos_bottom)
x_min = 0
y_min = 0
pos_data = np.genfromtxt(path, dtype = float, unpack = True, delimiter='\t', skip_header = 13)
# print(mag_data[2])
normalise = 1.0
count = 0
for length in pos_data[5]:
    if(length != 230):
        
        x_pos_bottom.append(0.1*pos_data[3][count])
        y_pos_bottom.append(0.1*pos_data[4][count])
        
    count += 1


path = os.getcwd() +   "/cells-00000000.cfg"
try:
    open(path)
except:
    print("cannot find " + path)
    
mag_data = np.genfromtxt(path, dtype = float, unpack = True, delimiter='\t', skip_header = 13)
# print(mag_data[2])
normalise = 1.0
count = 0
for length in mag_data[2]:
    if(count % 2 == 1):
        
        mc_1_avg_bottom.append(mag_data[7][count])
        mc_1_avg_bottom_x.append(mag_data[5][count])
        mc_1_avg_bottom_y.append(mag_data[6][count])   
    count += 1

grid_x, grid_y = np.mgrid[x_min:x_width:50j, y_min:y_width:50j]

grid_mx_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
grid_my_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
grid_m_avg  = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")

mag_min = (np.min(grid_m_avg))
mag_max = (np.max(grid_m_avg))
if(mag_min < 0.0):
    range = mag_max - mag_min
    if(abs(mag_min)/ range > 0.1):
        nodes = [0, abs(mag_min)/ range, abs(mag_min)/ range + 0.333*mag_max/range, abs(mag_min)/ range + 0.667*mag_max/range, 1.0]
        cmap1 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #blue-white-purple
        norm1 = colors.Normalize(vmin= mag_min, vmax= mag_max)
    else:
        nodes = [0, 0.333, 0.667, 1.0]
        cmap1 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
        norm1 = colors.Normalize(vmin= 0.0, vmax= mag_max)
else:
    nodes = [0, 0.333, 0.667, 1.0]
    cmap1 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
    norm1 = colors.Normalize(vmin= mag_min, vmax= mag_max)
        
plt.subplot(1,3, 1)
plt.xticks([])
plt.yticks([])
#grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
plt.imshow(grid_mx_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )

plt.subplot(1,3, 2)
plt.xticks([])
plt.yticks([])
#grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
plt.imshow(grid_my_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )

plt.subplot(1,3, 3)
plt.xticks([])
plt.yticks([])
#grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
plt.imshow(grid_m_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
#plt.colorbar( cmap = cmap1, shrink = 0.5, aspect= 10, norm=norm1 )

#plt.savefig(param + "-" + K + "-" + index +"-z-layers-avg-ori.png", bbox_inches = 'tight')
plt.savefig("FGT-FC-510mT-dipole.png", bbox_inches = 'tight')
plt.clf()







