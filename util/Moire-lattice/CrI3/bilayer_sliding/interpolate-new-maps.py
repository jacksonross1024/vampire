
import os 
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt





#intralayer 
# grid_x, grid_y = np.mgrid[0:0.9:100j, 0:1.0:100j]
# for i in range(1,13):
#     path = os.getcwd() + "/Cr1_intra_"+str(i)+".txt"
#     print(path)
#     x_pos = np.genfromtxt(path, dtype=None, usecols = 8)
#     y_pos = np.genfromtxt(path, dtype=None, usecols = 9)
#     J_inter = np.genfromtxt(path, dtype=None, usecols = 4)

#     #print(x_pos[...])
#     grid_z = griddata((x_pos,y_pos),J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)

#     np.savetxt("Cr1_intra_"+str(i)+"_map.txt", (grid_z.T), delimiter='\t', fmt='%.4f')

#     interp = CloughTocher2DInterpolator(list(zip(x_pos, y_pos)), J_inter)
#     Z = interp(grid_x, grid_y)

#     plt.plot()
#     plt.imshow(grid_z.T, extent=(0,0.9,0,1.0), origin='lower')
#     plt.title("Cr1_intra_"+str(i))
#     plt.gcf().set_size_inches(6, 6)

#     plt.savefig("Cr1_intra_"+str(i)+".png", bbox_inches = 'tight')

#inter


path = os.getcwd() + "/Cr1_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)

grid_x, grid_y = np.mgrid[-7.3:7.3:200j, -7.3:7.3:200j]
#print(x_pos[...])
grid_J1 = griddata((x_pos,y_pos),J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr1_J_inter_map.txt", (grid_J1.T ), delimiter='\t', fmt='%.4f')
grid_D1x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr1_Dx_inter_map.txt", (grid_D1x.T ), delimiter='\t', fmt='%.4f')
grid_D1y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr1_Dy_inter_map.txt", (grid_D1y.T ), delimiter='\t', fmt='%.4f')
grid_D1z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr1_Dz_inter_map.txt", (grid_D1z.T ), delimiter='\t', fmt='%.4f')

path = os.getcwd() + "/Cr2_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)
y_pos = y_pos*-1

grid_x, grid_y = np.mgrid[-7.3:7.3:200j, -7.3:7.3:200j]
#print(x_pos[...])
grid_J2 = griddata((x_pos,y_pos),J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr2_J_inter_map.txt", (-1*(grid_J2.T) ), delimiter='\t', fmt='%.4f')
grid_D2x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr2_Dx_inter_map.txt", -1*(grid_D2x.T ), delimiter='\t', fmt='%.4f')
y_pos = y_pos*-1
grid_D2y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr2_Dy_inter_map.txt", -1*(grid_D2y.T ), delimiter='\t', fmt='%.4f')
grid_D2z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr2_Dz_inter_map.txt", -1*(grid_D2z.T ), delimiter='\t', fmt='%.4f')

grid_J = np.add(grid_J1, grid_J2,  dtype=np.float64)
grid_Dx = np.add(grid_D1x, -1*grid_D2x,  dtype=np.float64)
grid_Dy = np.add(grid_D1y, -1*grid_D2y,  dtype=np.float64)
grid_Dz = np.add(grid_D1z, -1*grid_D2z,  dtype=np.float64)

path = os.getcwd() + "/Cr3_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)
x_pos = x_pos * -1

grid_x, grid_y = np.mgrid[-7.3:7.3:200j, -7.3:7.3:200j]
#print(x_pos[...])
grid_J3 = griddata((x_pos,y_pos), J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr3_J_inter_map.txt", (grid_J3.T ), delimiter='\t', fmt='%.4f')
grid_D3x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr3_Dx_inter_map.txt", -1*(grid_D3x.T ), delimiter='\t', fmt='%.4f')
grid_D3y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr3_Dy_inter_map.txt", (grid_D3y.T ), delimiter='\t', fmt='%.4f')
grid_D3z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr3_Dz_inter_map.txt", -1*(grid_D3z.T ), delimiter='\t', fmt='%.4f')

grid_J = np.add(grid_J, grid_J3,  dtype=np.float64)
grid_Dx = np.add(grid_Dx, -1*grid_D3x,  dtype=np.float64)
grid_Dy = np.add(grid_Dy, grid_D3y,  dtype=np.float64)
grid_Dz = np.add(grid_Dz, -1*grid_D3z,  dtype=np.float64)

path = os.getcwd() + "/Cr4_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)

x_pos = x_pos*-1
y_pos = y_pos*-1

grid_x, grid_y = np.mgrid[-7.3:7.3:200j, -7.3:7.3:200j]
#print(x_pos[...])
grid_J4 = griddata((x_pos,y_pos), J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr4_J_inter_map.txt", (grid_J4.T ), delimiter='\t', fmt='%.4f')
grid_D4x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr4_Dx_inter_map.txt", (grid_D4x.T ), delimiter='\t', fmt='%.4f')
grid_D4y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr4_Dy_inter_map.txt", (grid_D4y.T ), delimiter='\t', fmt='%.4f')
grid_D4z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
np.savetxt("Cr4_Dz_inter_map.txt", (grid_D4z.T ), delimiter='\t', fmt='%.4f')

grid_J = np.add(grid_J, grid_J4,  dtype=np.float64)/4.0
grid_Dx = np.add(grid_Dx, grid_D4x,  dtype=np.float64)/4.0
grid_Dy = np.add(grid_Dy, grid_D4y,  dtype=np.float64)/4.0
grid_Dz = np.add(grid_Dz, grid_D4z,  dtype=np.float64)/4.0

np.savetxt("J_total_inter_map.txt", (grid_J.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dx_total_inter_map.txt", (grid_Dx.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dy_total_inter_map.txt", (grid_Dy.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dz_total_inter_map.txt", (grid_Dz.T), delimiter='\t', fmt='%.4f')

np.savetxt("Cr1_inter_map.txt", np.ravel(grid_J.T), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dx_inter_map_avg.txt", np.ravel(grid_Dx.T ), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dy_inter_map_avg.txt", np.ravel(grid_Dy.T ), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dz_inter_map_avg.txt", np.ravel(grid_Dz.T ), delimiter='\t', fmt='%.4f')

np.savetxt("Cr2_inter_map.txt", np.ravel(np.flip(grid_J.T,0)), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dx_inter_map_avg.txt", np.ravel(-1*(np.flip(grid_Dx.T, 0) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dy_inter_map_avg.txt", np.ravel(-1*(grid_Dy.T )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dz_inter_map_avg.txt", np.ravel(-1*(grid_Dz.T )), delimiter='\t', fmt='%.4f')

np.savetxt("Cr3_inter_map.txt", np.ravel(np.flip(grid_J.T,1)), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dx_inter_map_avg.txt", np.ravel(-1*np.flip(grid_Dx.T,1 )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dy_inter_map_avg.txt", np.ravel(np.flip(grid_Dy.T,1 )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dz_inter_map_avg.txt", np.ravel(-1*np.flip(grid_Dz.T,1 )), delimiter='\t', fmt='%.4f')

np.savetxt("Cr4_inter_map.txt", np.ravel(np.flip(grid_J.T, (0,1))), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dx_inter_map_avg.txt", np.ravel(np.flip(grid_Dx.T, (0,1) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dy_inter_map_avg.txt", np.ravel(np.flip(grid_Dy.T, (0,1) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dz_inter_map_avg.txt", np.ravel(np.flip(grid_Dz.T, (0,1) )), delimiter='\t', fmt='%.4f')


# plt.subplot(221)
# plt.imshow(grid_z1.T, extent=(-7.3,7.3,-7.3,7.3), origin='lower')
# plt.title("Cr1_inter")
# plt.subplot(222)
# plt.imshow(grid_z2.T, extent=(-7.3,7.3,-7.3,7.3), origin='lower')
# plt.title("Cr2_inter")
# plt.subplot(223)
# plt.imshow(grid_z3.T, extent=(-7.3,7.3,-7.3,7.3), origin='lower')
# plt.title("Cr3_inter")
# plt.subplot(224)
# plt.imshow(grid_z4.T, extent=(-7.3,7.3,-7.3,7.3), origin='lower')
# plt.title("Cr4_inter")

# plt.gcf().set_size_inches(6, 6)

# plt.savefig("Cr1_inter.png", bbox_inches = 'tight')
