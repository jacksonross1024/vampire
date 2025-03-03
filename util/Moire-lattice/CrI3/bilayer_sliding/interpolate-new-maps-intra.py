
import os 
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt
from scipy.ndimage import spline_filter

grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]
grid_J = np.array((40,40))
grid_Dx = np.array((40,40))
grid_Dy = np.array((40,40))
grid_Dz = np.array((40,40))
plt.plot()
plt.gcf().set_size_inches(6, 8)
#print(x_pos[...])
J_avg = [[[]]]
Dx_avg = [[[]]]
Dy_avg = [[[]]]
Dz_avg = [[[]]]

file1 = 'Cr1_intra'
file2 = 'Cr2_intra'
file3 = 'Cr3_intra'
file4 = 'Cr4_intra'
#for file in ["Cr1_intra", "Cr2_intra", "Cr3_intra", "Cr4_intra"]:
path1 = os.getcwd() + "/" + file1 + ".txt"
path2 = os.getcwd() + "/" + file2 + ".txt"
path3 = os.getcwd() + "/" + file3 + ".txt"
path4 = os.getcwd() + "/" + file4 + ".txt"


#with open(file + '_data.txt', 'wb') as f:
f1 = open(file1 + '_data.txt', 'wb')
f2 = open(file2 + '_data.txt', 'wb')
f3 = open(file3 + '_data.txt', 'wb')
f4 = open(file4 + '_data.txt', 'wb')
intra_data1 = np.genfromtxt(path1, dtype=float, delimiter='\t', unpack='False')
intra_data2 = np.genfromtxt(path2, dtype=float, delimiter='\t', unpack='False')
intra_data3 = np.genfromtxt(path3, dtype=float, delimiter='\t', unpack='False')
intra_data4 = np.genfromtxt(path4, dtype=float, delimiter='\t', unpack='False')
for i in range(12):
    grid_Ji_1  = griddata((intra_data1[8][i::12], intra_data1[9][i::12]), intra_data1[4][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dx1_1 = griddata((intra_data1[8][i::12], intra_data1[9][i::12]),intra_data1[5][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dy1_1 = griddata((intra_data1[8][i::12], intra_data1[9][i::12]),intra_data1[6][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dz1_1 = griddata((intra_data1[8][i::12], intra_data1[9][i::12]),intra_data1[7][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)

    grid_Ji_2  = griddata((intra_data2[8][i::12], intra_data2[9][i::12]), intra_data2[4][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dx1_2 = griddata((intra_data2[8][i::12], intra_data2[9][i::12]),intra_data2[5][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dy1_2 = griddata((intra_data2[8][i::12], intra_data2[9][i::12]),intra_data2[6][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dz1_2 = griddata((intra_data2[8][i::12], intra_data2[9][i::12]),intra_data2[7][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)

    grid_Ji_3  = griddata((intra_data3[8][i::12], intra_data3[9][i::12]), intra_data3[4][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dx1_3 = griddata((intra_data3[8][i::12], intra_data3[9][i::12]),intra_data3[5][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dy1_3 = griddata((intra_data3[8][i::12], intra_data3[9][i::12]),intra_data3[6][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dz1_3 = griddata((intra_data3[8][i::12], intra_data3[9][i::12]),intra_data3[7][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)

    grid_Ji_4  = griddata((intra_data4[8][i::12], intra_data4[9][i::12]), intra_data4[4][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dx1_4 = griddata((intra_data4[8][i::12], intra_data4[9][i::12]),intra_data4[5][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dy1_4 = griddata((intra_data4[8][i::12], intra_data4[9][i::12]),intra_data4[6][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    grid_Dz1_4 = griddata((intra_data4[8][i::12], intra_data4[9][i::12]),intra_data4[7][i::12], (grid_x, grid_y), method="cubic", fill_value=0.0)
    
    grid_J_avg = np.add(np.add(np.add(grid_Ji_1, grid_Ji_2), grid_Ji_3), grid_Ji_4)/4.0
    if(i < 3 or i > 8): 
        grid_Dx_avg = np.add(np.add(np.add(grid_Dx1_1, -1*grid_Dx1_2), -1*grid_Dx1_3), grid_Dx1_1)/4.0
        grid_Dy_avg = np.add(np.add(np.add(grid_Dy1_1, -1*grid_Dy1_2), -1*grid_Dy1_3), grid_Dy1_4)/4.0
        grid_Dz_avg = np.add(np.add(np.add(grid_Dz1_1, -1*grid_Dz1_2), -1*grid_Dz1_3), grid_Dz1_4)/4.0

        x_vec1 = np.full((np.size(grid_J_avg),1), intra_data1[1][i])
        y_vec1 = np.full((np.size(grid_J_avg),1), intra_data1[2][i])
        x_vec2 = np.full((np.size(grid_J_avg),1), intra_data2[1][i])
        y_vec2 = np.full((np.size(grid_J_avg),1), intra_data2[2][i])
        x_vec3 = np.full((np.size(grid_J_avg),1), intra_data3[1][i])
        y_vec3 = np.full((np.size(grid_J_avg),1), intra_data3[2][i])
        x_vec4 = np.full((np.size(grid_J_avg),1), intra_data4[1][i])
        y_vec4 = np.full((np.size(grid_J_avg),1), intra_data4[2][i])

        x_shift = np.reshape(grid_x, (np.size(grid_J_avg),1))
        y_shift = np.reshape(grid_y, (np.size(grid_J_avg),1))

        J_data1 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data1 = np.reshape(grid_Dx_avg.T, (np.size(grid_Dx_avg),1))
        Dy_data1 = np.reshape(grid_Dy_avg.T, (np.size(grid_Dy_avg),1))
        Dz_data1 = np.reshape(grid_Dz_avg.T, (np.size(grid_Dz_avg),1))

        J_data2 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data2 = np.reshape(grid_Dx_avg.T, (np.size(grid_Dx_avg),1))
        Dy_data2 = np.reshape(grid_Dy_avg.T, (np.size(grid_Dy_avg),1))
        Dz_data2 = np.reshape(grid_Dz_avg.T, (np.size(grid_Dz_avg),1))

        J_data3 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data3 = np.reshape(grid_Dx_avg.T, (np.size(grid_Dx_avg),1))
        Dy_data3 = np.reshape(grid_Dy_avg.T, (np.size(grid_Dy_avg),1))
        Dz_data3 = np.reshape(grid_Dz_avg.T, (np.size(grid_Dz_avg),1))

        J_data4 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data4 = np.reshape(grid_Dx_avg.T, (np.size(grid_Dx_avg),1))
        Dy_data4 = np.reshape(grid_Dy_avg.T, (np.size(grid_Dy_avg),1))
        Dz_data4 = np.reshape(grid_Dz_avg.T, (np.size(grid_Dz_avg),1))

        np.savetxt(file1 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data1, Dx_data1, Dy_data1, Dz_data1), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file2 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data2, Dx_data2, Dy_data2, Dz_data2), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file3 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data3, Dx_data3, Dy_data3, Dz_data3), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file4 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data4, Dx_data4, Dy_data4, Dz_data4), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f1, (np.concatenate((x_vec1, y_vec1, J_data1, Dx_data1, Dy_data1, Dz_data1), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f2, (np.concatenate((x_vec2, y_vec2, J_data2, -1*Dx_data2, -1*Dy_data2, -1*Dz_data2), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f3, (np.concatenate((x_vec3, y_vec3, J_data3, -1*Dx_data3, -1*Dy_data3, -1*Dz_data3), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f4, (np.concatenate((x_vec4, y_vec4, J_data4, Dx_data4, Dy_data4, Dz_data4), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')

    elif(i < 9):
        # grid_Dx_avg = grid_Dx1_1 # np.add(np.add(np.add(grid_Dx1_1, -1*grid_Dx1_2), -1*grid_Dx1_3), grid_Dx1_1)/4.0
        # grid_Dy_avg = grid_Dy1_1 #np.add(np.add(np.add(grid_Dy1_1, -1*grid_Dy1_2), -1*grid_Dy1_3), grid_Dy1_4)/4.0
        # grid_Dz_avg = np.add(np.add(np.add(grid_Dz1_1, -1*grid_Dz1_2), -1*grid_Dz1_3), grid_Dz1_4)/4.0

        x_vec1 = np.full((np.size(grid_J_avg),1), intra_data1[1][i])
        y_vec1 = np.full((np.size(grid_J_avg),1), intra_data1[2][i])
        x_vec2 = np.full((np.size(grid_J_avg),1), intra_data2[1][i])
        y_vec2 = np.full((np.size(grid_J_avg),1), intra_data2[2][i])
        x_vec3 = np.full((np.size(grid_J_avg),1), intra_data3[1][i])
        y_vec3 = np.full((np.size(grid_J_avg),1), intra_data3[2][i])
        x_vec4 = np.full((np.size(grid_J_avg),1), intra_data4[1][i])
        y_vec4 = np.full((np.size(grid_J_avg),1), intra_data4[2][i])

        x_shift = np.reshape(grid_x, (np.size(grid_J_avg),1))
        y_shift = np.reshape(grid_y, (np.size(grid_J_avg),1))

        J_data1 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data1 = np.reshape(grid_Dx1_1.T, (np.size(grid_Dx1_1),1))
        Dy_data1 = np.reshape(grid_Dy1_1.T, (np.size(grid_Dy1_1),1))
        Dz_data1 = np.reshape(grid_Dz1_1.T, (np.size(grid_Dz1_1),1))

        J_data2 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data2 = np.reshape(grid_Dx1_2.T, (np.size(grid_Dx1_2),1))
        Dy_data2 = np.reshape(grid_Dy1_2.T, (np.size(grid_Dy1_2),1))
        Dz_data2 = np.reshape(grid_Dz1_2.T, (np.size(grid_Dz1_2),1))

        J_data3 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data3 = np.reshape(grid_Dx1_1.T, (np.size(grid_Dx1_1),1))
        Dy_data3 = np.reshape(grid_Dy1_1.T, (np.size(grid_Dy1_1),1))
        Dz_data3 = np.reshape(grid_Dz1_1.T, (np.size(grid_Dz1_1),1))

        J_data4 = np.reshape(grid_J_avg.T, (np.size(grid_J_avg),1))
        Dx_data4 = np.reshape(grid_Dx1_2.T, (np.size(grid_Dx1_2),1))
        Dy_data4 = np.reshape(grid_Dy1_2.T, (np.size(grid_Dy1_2),1))
        Dz_data4 = np.reshape(grid_Dz1_2.T, (np.size(grid_Dz1_2),1))

        np.savetxt(file1 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data1, Dx_data1, Dy_data1, Dz_data1), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file2 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data2, Dx_data2, Dy_data2, Dz_data2), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file3 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data3, -1*Dx_data3, -1*Dy_data3, -1*Dz_data3), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(file4 + "_" + str(i) + "avg.txt", (np.concatenate((x_shift, y_shift, J_data4, -1*Dx_data4, -1*Dy_data4, -1*Dz_data4), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f1, (np.concatenate((x_vec1, y_vec1, J_data1, Dx_data1, Dy_data1, Dz_data1), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f2, (np.concatenate((x_vec2, y_vec2, J_data2, Dx_data2, Dy_data2, Dz_data2), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f3, (np.concatenate((x_vec3, y_vec3, J_data3, -1*Dx_data3, -1*Dy_data3, -1*Dz_data3), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')
        np.savetxt(f4, (np.concatenate((x_vec4, y_vec4, J_data4, -1*Dx_data4, -1*Dy_data4, -1*Dz_data4), axis=1)) , delimiter='\t', newline='\n', fmt='%.4f')





