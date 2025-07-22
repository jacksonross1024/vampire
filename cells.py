

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

#inter



for K in [""]:
    for i in range(10,12,1):
    
        index = f"{i:08d}"
        
            
        # param = "param-" + angle + "-" + J + "-" + DMI +"-"+ J_twist + "-" + J_prist + "-" + J_intra + "-" + DMI_sub
        path = os.getcwd() +   "/cells-" + index + ".txt"
        try:
            open(path)
        except:
            print("cannot find " + path)
            continue
        
        print(path)
        x_pos_bottom = []
        y_pos_bottom = []
        mc_1_avg_bottom = []
        mc_2_avg_bottom = []
        mc_3_avg_bottom = []
        mc_4_avg_bottom = []
        mc_1_avg_bottom_x = []
        mc_2_avg_bottom_x = []
        mc_3_avg_bottom_x = []
        mc_4_avg_bottom_x = []
        mc_1_avg_bottom_y = []
        mc_2_avg_bottom_y = []
        mc_3_avg_bottom_y = []
        mc_4_avg_bottom_y = []
        count = 0
        x_width = 60 # np.max(x_pos_bottom)
        y_width = 60 #np.max(y_pos_bottom)
        x_min = 0
        y_min = 0
        mag_data = np.genfromtxt(path, dtype = float, unpack = True, delimiter='\t', skip_header = 1)
        # print(mag_data[2])
        normalise = 1.0
        for length in mag_data[2]:
            if(length == 0):
                if(0.1*mag_data[0][count] > float(x_width)):
                    count += 1
                    continue
                if(0.1*mag_data[1][count] > float(y_width)):
                    count += 1
                    continue
                if(0.1*mag_data[0][count] < float(x_min)):
                    count += 1
                    continue
                if(0.1*mag_data[1][count] < float(y_min)):
                    count += 1
                    continue
                if(mag_data[2][count] > float(0.0)):
                    count += 1
                    continue
                x_pos_bottom.append(0.1*mag_data[0][count])
                y_pos_bottom.append(0.1*mag_data[1][count])
                # mc_1_avg_bottom.append(3.14*mag_data[5][count]*mag_data[6][count]*mag_data[7][count]/normalise)
                
                # mc_1_avg_bottom_x.append(3.14*mag_data[3][count]*mag_data[6][count]*mag_data[7][count]/normalise)
            
                # mc_1_avg_bottom_y.append(3.14*mag_data[4][count]*mag_data[6][count]*mag_data[7][count]/normalise)

                mc_1_avg_bottom.append(mag_data[20][count]*mag_data[21][count])
                
                mc_1_avg_bottom_x.append(mag_data[18][count]*mag_data[21][count])
            
                mc_1_avg_bottom_y.append(mag_data[19][count]*mag_data[21][count])
                
            count += 1
        

        # #J_inter_weight_bottom = np.add(np.multiply(mc_1_avg_bottom,0.183), np.add(np.multiply(mc_2_avg_bottom,0.186), np.add(np.multiply(mc_3_avg_bottom,0.269), np.multiply(mc_4_avg_bottom, 0.43))))
        # J_inter_avg_bottom = np.add(np.multiply(mc_1_avg_bottom,1.0), np.add(np.multiply(mc_2_avg_bottom,1.0), np.add(np.multiply(mc_3_avg_bottom,1.0), np.multiply(mc_4_avg_bottom, 1.0))))
        # #J_inter_weight_bottom_x = np.add(np.multiply(mc_1_avg_bottom_x,0.183), np.add(np.multiply(mc_2_avg_bottom_x,0.186), np.add(np.multiply(mc_3_avg_bottom_x,0.269), np.multiply(mc_4_avg_bottom_x, 0.43))))
        # J_inter_avg_bottom_x = np.add(np.multiply(mc_1_avg_bottom_x,1.0), np.add(np.multiply(mc_2_avg_bottom_x,1.0), np.add(np.multiply(mc_3_avg_bottom_x,1.0), np.multiply(mc_4_avg_bottom_x, 1.0))))
        # #J_inter_weight_bottom_y = np.add(np.multiply(mc_1_avg_bottom_y,0.183), np.add(np.multiply(mc_2_avg_bottom_y,0.186), np.add(np.multiply(mc_3_avg_bottom_y,0.269), np.multiply(mc_4_avg_bottom_y, 0.43))))
        # J_inter_avg_bottom_y = np.add(np.multiply(mc_1_avg_bottom_y,1.0), np.add(np.multiply(mc_2_avg_bottom_y,1.0), np.add(np.multiply(mc_3_avg_bottom_y,1.0), np.multiply(mc_4_avg_bottom_y, 1.0))))
        
        
        grid_x, grid_y = np.mgrid[x_min:x_width:50j, y_min:y_width:50j]
        #J_inter_weight_bottom = mc_1_avg_bottom*0.183 + mc_2_avg_bottom*0.186 + mc_3_avg_bottom*0.269 + mc_4_avg_bottom*0.43
        #J_inter_avg_bottom = mc_1_avg_bottom*0.25 + mc_2_avg_bottom*0.25 + mc_3_avg_bottom*0.25 + mc_4_avg_bottom*0.25

        #x_pos_bottom = x_pos[::2]
        #y_pos_bottom = y_pos[::2]
        #print(mag_data_bottom[2])

        #J_inter_weight_bottom = J_inter_weight[::2]
        #J_inter_avg_bottom = J_inter_avg[::2]
        #mc_1_avg_bottom = mc_1_avg[::2]
        #mc_2_avg_bottom = mc_2_avg[::2]
        #mc_3_avg_bottom = mc_3_avg[::2]
        #mc_4_avg_bottom = mc_4_avg[::2]
        
        grid_mx_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
        grid_my_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
        grid_m_avg  = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
       # grid_alloy  = griddata((0.1*alloy_data[0],0.1*alloy_data[1]), alloy_data[2], (grid_x, grid_y), fill_value=1.0, method="cubic")
    
        #print(np.max(grid_J1,1))
        #nodes = [0, 0.5, 0.55, 0.65, 1.0]
        #cmap = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"])))

        
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
        plt.title(K + "; " + index )
        #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
        plt.imshow(grid_my_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
        
        plt.subplot(1,3, 3)
        plt.xticks([])
        plt.yticks([])
        #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
        plt.imshow(grid_m_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
        #plt.colorbar( cmap = cmap1, shrink = 0.5, aspect= 10, norm=norm1 )

   
        #plt.savefig(param + "-" + K + "-" + index +"-z-layers-avg-ori.png", bbox_inches = 'tight')
        plt.savefig("FGT-FC-510mT-" + index + ".png", bbox_inches = 'tight')
        plt.clf()







