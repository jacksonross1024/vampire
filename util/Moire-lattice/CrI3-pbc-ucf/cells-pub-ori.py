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



count = 0
#for angle in ["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9"]:
for angle in ["0.5"]:
#        for J_twist in ["1.0"]:
#            for index in ["0", "1", "2", "3"]:

    for J in ["0.04"]:
        for DMI in ["1"]:
            for J_twist in ["2.44"]:
                for J_intra in ["1"]:
                    for J_prist in ["6"]:
                        for K in ["0.5K-1nm"]:
                            for index in ["00000000", "00000001"]:
                    #for index in ["5K-AF-hr"]: #, "00000000", "00000006", "00000016"]:
                            
                    #for index in ["5K-AF-hr", "5K-FC-hrs"]:
                        #index = "%08.f" % i
                                param = index #"param-" + angle + "-" + J + "-" + DMI + "-" + J_twist + "-" + J_prist + "-" + J_intra
                                path = os.getcwd() + "/cells-" + index + ".txt"
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
                                x_width = 300# np.max(x_pos_bottom)
                                y_width = 300 #np.max(y_pos_bottom)
                                normalise = 25.0

                                mag_data = np.genfromtxt(path, dtype = float, unpack = True, delimiter='\t', skip_header = 1)
                                # print(mag_data[2])
                                for length in mag_data[2]:
                                    if(length == 0):
                                        if(0.1*mag_data[0][count] > float(x_width)):
                                            count += 1
                                            continue
                                        if(0.1*mag_data[1][count] > float(y_width)):
                                            count += 1
                                            continue
                                        x_pos_bottom.append(0.1*mag_data[0][count])
                                        y_pos_bottom.append(0.1*mag_data[1][count])
                                        mc_1_avg_bottom.append(2.98*mag_data[5][count]*mag_data[6][count]*mag_data[7][count]/normalise)
                                        mc_2_avg_bottom.append(2.98*mag_data[10][count]*mag_data[11][count]*mag_data[12][count]/normalise)
                                        mc_3_avg_bottom.append(2.98*mag_data[15][count]*mag_data[16][count]*mag_data[17][count]/normalise)
                                        mc_4_avg_bottom.append(2.98*mag_data[20][count]*mag_data[21][count]*mag_data[22][count]/normalise)
                                        mc_1_avg_bottom_x.append(2.98*mag_data[3][count]*mag_data[6][count]*mag_data[7][count]/normalise)
                                        mc_2_avg_bottom_x.append(2.98*mag_data[8][count]*mag_data[11][count]*mag_data[12][count]/normalise)
                                        mc_3_avg_bottom_x.append(2.98*mag_data[13][count]*mag_data[16][count]*mag_data[17][count]/normalise)
                                        mc_4_avg_bottom_x.append(2.98*mag_data[18][count]*mag_data[21][count]*mag_data[22][count]/normalise)
                                        mc_1_avg_bottom_y.append(2.98*mag_data[4][count]*mag_data[6][count]*mag_data[7][count]/normalise)
                                        mc_2_avg_bottom_y.append(2.98*mag_data[9][count]*mag_data[11][count]*mag_data[12][count]/normalise)
                                        mc_3_avg_bottom_y.append(2.98*mag_data[14][count]*mag_data[16][count]*mag_data[17][count]/normalise)
                                        mc_4_avg_bottom_y.append(2.98*mag_data[19][count]*mag_data[21][count]*mag_data[22][count]/normalise)
                                    count += 1
                                

                                #J_inter_weight_bottom = np.add(np.multiply(mc_1_avg_bottom,0.183), np.add(np.multiply(mc_2_avg_bottom,0.186), np.add(np.multiply(mc_3_avg_bottom,0.269), np.multiply(mc_4_avg_bottom, 0.43))))
                                J_inter_avg_bottom = np.add(np.multiply(mc_1_avg_bottom,1.0), np.add(np.multiply(mc_2_avg_bottom,1.0), np.add(np.multiply(mc_3_avg_bottom,1.0), np.multiply(mc_4_avg_bottom, 1.0))))
                                #J_inter_weight_bottom_x = np.add(np.multiply(mc_1_avg_bottom_x,0.183), np.add(np.multiply(mc_2_avg_bottom_x,0.186), np.add(np.multiply(mc_3_avg_bottom_x,0.269), np.multiply(mc_4_avg_bottom_x, 0.43))))
                                J_inter_avg_bottom_x = np.add(np.multiply(mc_1_avg_bottom_x,1.0), np.add(np.multiply(mc_2_avg_bottom_x,1.0), np.add(np.multiply(mc_3_avg_bottom_x,1.0), np.multiply(mc_4_avg_bottom_x, 1.0))))
                                #J_inter_weight_bottom_y = np.add(np.multiply(mc_1_avg_bottom_y,0.183), np.add(np.multiply(mc_2_avg_bottom_y,0.186), np.add(np.multiply(mc_3_avg_bottom_y,0.269), np.multiply(mc_4_avg_bottom_y, 0.43))))
                                J_inter_avg_bottom_y = np.add(np.multiply(mc_1_avg_bottom_y,1.0), np.add(np.multiply(mc_2_avg_bottom_y,1.0), np.add(np.multiply(mc_3_avg_bottom_y,1.0), np.multiply(mc_4_avg_bottom_y, 1.0))))
                                
                                grid_x, grid_y = np.mgrid[0:x_width:75j, 0:y_width:75j]
                                
                                grid_mx_avg = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_my_avg = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_m_avg  = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_1x_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_1y_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_1z_avg = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_2x_avg = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_2y_avg = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_2z_avg = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_3x_avg = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_3y_avg = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_3z_avg = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_4x_avg = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_4y_avg = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                grid_4z_avg = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                
                                mag_min = -25 #np.min(grid_m_avg)
                                mag_max = 25 #abs(np.max(grid_m_avg))
                                
                                if(mag_min < 0.0):
                                    range = mag_max - mag_min
                                    if(abs(mag_min)/ range > 0.1):
                                        nodes = [0, abs(mag_min)/ range, abs(mag_min)/ range + 0.333*mag_max/range, abs(mag_min)/ range + 0.667*mag_max/range, 1.0]
                                        cmap = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #blue-white-purple
                                        norm = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                    else:
                                        nodes = [0, 0.333, 0.667, 1.0]
                                        cmap = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                        norm = colors.Normalize(vmin= 0.0, vmax= mag_max)
                                else:
                                    nodes = [0, 0.333, 0.667, 1.0]
                                    cmap = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                    norm = colors.Normalize(vmin= mag_min, vmax= mag_max)

                                plt.xticks([])
                                plt.yticks([])

                                plt.plot(111)
                                plt.gcf().set_size_inches(6, 6)
                                plt.title(param + "; " + K + "; " + index)
                                plt.imshow(grid_m_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap,norm=norm, interpolation='gaussian' )
                         
                                plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                
                                #plt.savefig(param + "/" + param + "-" + K + "-" + index +"-z-avg-1nm.png", bbox_inches = 'tight')
                                plt.savefig(param + ".png", bbox_inches = 'tight')
                                plt.clf()
                                
    
                                # np.savetxt("Cr1_J_inter_map.txt", (grid_J1.T ), delimiter='\t', fmt='%.4f')
                                plt.gcf().set_size_inches(6, 8)
                                plt.subplot(5,3,9)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_m_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap,norm=norm, interpolation='gaussian' )
                                plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                      
                                
                                mag_min = np.min(grid_1z_avg)
                                mag_max = abs(np.max(grid_1z_avg))
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
                                        
                                plt.subplot(5,3, 1)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_1x_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
                                
                                plt.subplot(5,3, 2)
                                plt.xticks([])
                                plt.yticks([])
                                plt.title(param + "; " + K + "; " + index)
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_1y_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
                                
                                plt.subplot(5,3, 3)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_1_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_1z_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap1,norm=norm1, interpolation='gaussian' )
                                plt.colorbar( cmap = cmap1, shrink = 0.5, aspect= 10, norm=norm1 )
                                
                                mag_min = -15 #np.min(grid_2z_avg)
                                mag_max = 15 #abs(np.max(grid_2z_avg))
                                if(mag_min < 0.0):
                                    range = mag_max - mag_min
                                    if(abs(mag_min)/ range > 0.1):
                                        nodes = [0, abs(mag_min)/ range, abs(mag_min)/ range + 0.333*mag_max/range, abs(mag_min)/ range + 0.667*mag_max/range, 1.0]
                                        cmap2 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #blue-white-purple
                                        norm2 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                    else:
                                        nodes = [0, 0.333, 0.667, 1.0]
                                        cmap2 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                        norm2 = colors.Normalize(vmin= 0.0, vmax= mag_max)
                                else:
                                    nodes = [0, 0.333, 0.667, 1.0]
                                    cmap2 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                    norm2 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                        
                                plt.subplot(5,3, 4)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_2x_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap2,norm=norm2, interpolation='gaussian' )
                                
                                plt.subplot(5,3, 5)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_2y_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap2,norm=norm2, interpolation='gaussian' )
                                
                                plt.subplot(5,3, 6)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_2_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_2z_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap2,norm=norm2, interpolation='gaussian' )
                                plt.colorbar( cmap = cmap2, shrink = 0.5, aspect= 10, norm=norm2 )
                                
                                mag_min = -15 #np.min(grid_3z_avg)
                                mag_max = 15 # abs(np.max(grid_3z_avg))
                                if(mag_min < 0.0):
                                    range = mag_max - mag_min
                                    if(abs(mag_min)/ range > 0.1):
                                        nodes = [0, abs(mag_min)/ range, abs(mag_min)/ range + 0.333*mag_max/range, abs(mag_min)/ range + 0.667*mag_max/range, 1.0]
                                        cmap3 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #blue-white-purple
                                        norm3 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                    else:
                                        nodes = [0, 0.333, 0.667, 1.0]
                                        cmap3 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                        norm3 = colors.Normalize(vmin= 0.0, vmax= mag_max)
                                else:
                                    nodes = [0, 0.333, 0.667, 1.0]
                                    cmap3 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                    norm3 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                        
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                
                                plt.subplot(5,3, 7)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_mx_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap,norm=norm, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                
                                
                                plt.subplot(5,3, 8)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), J_inter_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_my_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap,norm=norm, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                
                                
                                mag_min = np.min(grid_4z_avg)
                                mag_max = abs(np.max(grid_4z_avg))
                                if(mag_min < 0.0):
                                    range = mag_max - mag_min
                                    if(abs(mag_min)/ range > 0.1):
                                        nodes = [0, abs(mag_min)/ range, abs(mag_min)/ range + 0.333*mag_max/range, abs(mag_min)/range + 0.667*mag_max/range, 1.0]
                                        cmap4 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #blue-white-purple
                                        norm4 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                    else:
                                        nodes = [0, 0.333, 0.667, 1.0]
                                        cmap4 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                        norm4 = colors.Normalize(vmin= 0.0, vmax= mag_max)
                                else:
                                    nodes = [0, 0.333, 0.667, 1.0]
                                    cmap4 = LinearSegmentedColormap.from_list("mycolors", list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"]))) #white-purple
                                    norm4 = colors.Normalize(vmin= mag_min, vmax= mag_max)
                                        
                                
                                plt.subplot(5,3, 10)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_3x_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap3,norm=norm3, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                
                                
                                plt.subplot(5,3, 11)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_3y_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap3,norm=norm3, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                               
                                
                                plt.subplot(5,3, 12)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_3_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_3z_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap3, norm=norm3, interpolation='gaussian' )
                                plt.colorbar( cmap = cmap3, shrink = 0.5, aspect= 10, norm=norm3 )
                                
                                plt.subplot(5,3, 13)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom_x, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_4x_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap4,norm=norm4, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                #plt.savefig(param + "-" + index +"-z-layers-avg.png", bbox_inches = 'tight')
                                #plt.clf()
                                
                                plt.subplot(5,3, 14)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom_y, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_4y_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap4,norm=norm4, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10, norm=norm  )
                                #plt.savefig(param + "-" + index +"-z-layers-avg.png", bbox_inches = 'tight')
                               # plt.clf()
                                
                                plt.subplot(5,3, 15)
                                plt.xticks([])
                                plt.yticks([])
                                #grid_J1 = griddata((x_pos_bottom,y_pos_bottom), mc_4_avg_bottom, (grid_x, grid_y), fill_value=0.0, method="cubic")
                                plt.imshow(grid_4z_avg.T, extent=(0,x_width, 0,x_width), origin='lower',cmap=cmap4,norm=norm4, interpolation='gaussian' )
                                #plt.colorbar( cmap = cmap, shrink = 0.5, aspect= 10  )
                                plt.colorbar( cmap = cmap4, shrink = 0.5, aspect= 10, norm=norm4  )
                                #plt.savefig(param + "-" + K + "-" + index +"-z-layers-avg-ori.png", bbox_inches = 'tight')
                                #plt.savefig(param + "/" + param + "-" + K + "-" + index +"-z-layers-avg-1nm.png", bbox_inches = 'tight')
                                plt.savefig(param + "-2.0degree-510mT-0K.png", bbox_inches = 'tight')
                                plt.clf()
                        
                        

    
