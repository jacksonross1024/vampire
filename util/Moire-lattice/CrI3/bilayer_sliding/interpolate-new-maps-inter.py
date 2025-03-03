
import os 
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt
from scipy.signal import cspline2d

path = os.getcwd() + "/Cr1_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)

grid_x, grid_y = np.mgrid[-7.276:7.276:100j, -7.402:7.402:100j]
grid_x2, grid_y2 = np.mgrid[-7.276:7.276:200j, -7.402:7.402:200j]
#print(x_pos[...])
grid_J1 = griddata((x_pos,y_pos),J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_J1 = cspline2d(grid_J1, 100.0, -0.00001)
np.savetxt("Cr1_J_inter_map_2.txt", (grid_J1.T ), delimiter='\t', fmt='%.4f')
grid_D1x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D1x = cspline2d(grid_D1x, 100.0, -0.00001)
np.savetxt("Cr1_Dx_inter_map_2.txt", (grid_D1x.T), delimiter='\t', fmt='%.4f')
grid_D1y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D1y = cspline2d(grid_D1y, 100.0, -0.00001)
np.savetxt("Cr1_Dy_inter_map_2.txt", (grid_D1y.T ), delimiter='\t', fmt='%.4f')
grid_D1z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D1z = cspline2d(grid_D1z, 100.0, -0.00001)
np.savetxt("Cr1_Dz_inter_map_2.txt", (grid_D1z.T ), delimiter='\t', fmt='%.4f')

path = os.getcwd() + "/Cr2_inter.txt"
print(path)
x_pos = np.genfromtxt(path, dtype=None, usecols = 1)
y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
J_inter = np.genfromtxt(path, dtype=None, usecols = 4)
Dx_inter = np.genfromtxt(path, dtype=None, usecols = 5)
Dy_inter = np.genfromtxt(path, dtype=None, usecols = 6)
Dz_inter = np.genfromtxt(path, dtype=None, usecols = 7)
y_pos = y_pos*-1

#grid_x, grid_y = np.mgrid[np.min(x_pos):np.max(x_pos):200j, np.min(y_pos):np.max(y_pos):200j]
#print(x_pos[...])
grid_J2 = griddata((x_pos,y_pos),J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_J2 = cspline2d(grid_J2, 100.0, -0.00001)
np.savetxt("Cr2_J_inter_map_2.txt", (-1*(grid_J2.T) ), delimiter='\t', fmt='%.4f')
grid_D2x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D2x = cspline2d(grid_D2x, 100.0, -0.00001)
np.savetxt("Cr2_Dx_inter_map_2.txt", -1*(grid_D2x.T ), delimiter='\t', fmt='%.4f')
y_pos = y_pos*-1
grid_D2y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D2y = cspline2d(grid_D2y, 100.0, -0.00001)
np.savetxt("Cr2_Dy_inter_map_2.txt", -1*(grid_D2y.T ), delimiter='\t', fmt='%.4f')
grid_D2z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D2z = cspline2d(grid_D2z, 100.0, -0.00001)
np.savetxt("Cr2_Dz_inter_map_2.txt", -1*(grid_D2z.T ), delimiter='\t', fmt='%.4f')

grid_J  = np.add(grid_J1, grid_J2,  dtype=np.float64)
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

#grid_x, grid_y = np.mgrid[np.min(x_pos):np.max(x_pos):200j, np.min(y_pos):np.max(y_pos):200j]
#print(x_pos[...])
grid_J3 = griddata((x_pos,y_pos), J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_J3 = cspline2d(grid_J3, 100.0, -0.00001)
np.savetxt("Cr3_J_inter_map_2.txt", (grid_J3.T ), delimiter='\t', fmt='%.4f')
grid_D3x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D3x = cspline2d(grid_D3x, 100.0, -0.00001)
np.savetxt("Cr3_Dx_inter_map_2.txt", -1*(grid_D3x.T ), delimiter='\t', fmt='%.4f')
grid_D3y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D3y = cspline2d(grid_D3y, 100.0, -0.00001)
np.savetxt("Cr3_Dy_inter_map_2.txt", (grid_D3y.T ), delimiter='\t', fmt='%.4f')
grid_D3z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D3z = cspline2d(grid_D3z, 100.0, -0.00001)
np.savetxt("Cr3_Dz_inter_map_2.txt", -1*(grid_D3z.T ), delimiter='\t', fmt='%.4f')

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

#grid_x, grid_y = np.mgrid[np.min(x_pos):np.max(x_pos):200j, np.min(y_pos):np.max(y_pos):200j]
#print(x_pos[...])
grid_J4 = griddata((x_pos,y_pos), J_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_J4 = cspline2d(grid_J4, 100.0, -0.00001)
np.savetxt("Cr4_J_inter_map_2.txt", (grid_J4.T ), delimiter='\t', fmt='%.4f')
grid_D4x = griddata((x_pos,y_pos), Dx_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D4x = cspline2d(grid_D4x, 100.0, -0.00001)
np.savetxt("Cr4_Dx_inter_map_2.txt", (grid_D4x.T ), delimiter='\t', fmt='%.4f')
grid_D4y = griddata((x_pos,y_pos), Dy_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D4y = cspline2d(grid_D4y, 100.0, -0.00001)
np.savetxt("Cr4_Dy_inter_map_2.txt", (grid_D4y.T ), delimiter='\t', fmt='%.4f')
grid_D4z = griddata((x_pos,y_pos), Dz_inter, (grid_x, grid_y), method="cubic", fill_value=0.0)
grid_D4z = cspline2d(grid_D4z, 100.0, -0.00001)
np.savetxt("Cr4_Dz_inter_map_2.txt", (grid_D4z.T ), delimiter='\t', fmt='%.4f')

grid_J = np.add(grid_J, grid_J4,  dtype=np.float64)/4.0
grid_Dx = np.add(grid_Dx, grid_D4x,  dtype=np.float64)/4.0
grid_Dy = np.add(grid_Dy, grid_D4y,  dtype=np.float64)/4.0
grid_Dz = np.add(grid_Dz, grid_D4z,  dtype=np.float64)/4.0

np.savetxt("J_total_inter_map_2.txt", (grid_J.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dx_total_inter_map_2.txt", (grid_Dx.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dy_total_inter_map_2.txt", (grid_Dy.T), delimiter='\t', fmt='%.4f')
np.savetxt("Dz_total_inter_map_2.txt", (grid_Dz.T), delimiter='\t', fmt='%.4f')

np.savetxt("Cr1_inter_map_2.txt", np.ravel(grid_J.T), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dx_inter_map_2_avg.txt", np.ravel(grid_Dx.T ), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dy_inter_map_2_avg.txt", np.ravel(grid_Dy.T ), delimiter='\t', fmt='%.4f')
np.savetxt("Cr1_Dz_inter_map_2_avg.txt", np.ravel(grid_Dz.T ), delimiter='\t', fmt='%.4f')

np.savetxt("Cr2_inter_map_2.txt", np.ravel(np.flip(grid_J.T,0)), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dx_inter_map_2_avg.txt", np.ravel(-1*(np.flip(grid_Dx.T, 0) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dy_inter_map_2_avg.txt", np.ravel(-1*(grid_Dy.T )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr2_Dz_inter_map_2_avg.txt", np.ravel(-1*(grid_Dz.T )), delimiter='\t', fmt='%.4f')

np.savetxt("Cr3_inter_map_2.txt", np.ravel(np.flip(grid_J.T,1)), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dx_inter_map_2_avg.txt", np.ravel(-1*np.flip(grid_Dx.T,1 )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dy_inter_map_2_avg.txt", np.ravel(np.flip(grid_Dy.T,1 )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr3_Dz_inter_map_2_avg.txt", np.ravel(-1*np.flip(grid_Dz.T,1 )), delimiter='\t', fmt='%.4f')

np.savetxt("Cr4_inter_map_2.txt", np.ravel(np.flip(grid_J.T, (0,1))), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dx_inter_map_2_avg.txt", np.ravel(np.flip(grid_Dx.T, (0,1) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dy_inter_map_2_avg.txt", np.ravel(np.flip(grid_Dy.T, (0,1) )), delimiter='\t', fmt='%.4f')
np.savetxt("Cr4_Dz_inter_map_2_avg.txt", np.ravel(np.flip(grid_Dz.T, (0,1) )), delimiter='\t', fmt='%.4f')



plt.subplot(422)
plt.imshow(grid_J.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(424)
plt.imshow(np.flip(grid_J.T,0), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(426)
plt.imshow(np.flip(grid_J.T,1), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(428)
plt.imshow(np.flip(grid_J.T, (0,1)), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_inter")

plt.subplot(421)
plt.imshow(grid_J1.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(423)
plt.imshow(np.flip(grid_J2.T,0), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(425)
plt.imshow(np.flip(grid_J3.T,1), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(427)
plt.imshow(np.flip(grid_J4.T, (0,1)), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_inter")

plt.gcf().set_size_inches(6, 8)

plt.savefig("Cr_J_inter.png", bbox_inches = 'tight')

plt.subplot(422)
plt.imshow(grid_Dx.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(424)
plt.imshow(-1*(np.flip(grid_Dx.T, 0) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(426)
plt.imshow(-1*np.flip(grid_Dx.T,1 ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(428)
plt.imshow(np.flip(grid_Dx.T, (0,1) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_Dy_inter")

plt.subplot(421)
plt.imshow(grid_D1x.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(423)
plt.imshow(-1*(np.flip(grid_D2x.T, 0) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(425)
plt.imshow(-1*np.flip(grid_D3x.T,1 ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(427)
plt.imshow(np.flip(grid_D4x.T, (0,1) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_Dy_inter")

plt.gcf().set_size_inches(6, 6)

plt.savefig("Cr_Dx_inter.png", bbox_inches = 'tight')

plt.subplot(221)
plt.imshow(grid_Dy.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(222)
plt.imshow(-1*(grid_Dy.T ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(223)
plt.imshow(np.flip(grid_Dy.T,1 ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(224)
plt.imshow(np.flip(grid_Dy.T, (0,1) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_Dy_inter")

plt.gcf().set_size_inches(6, 6)

plt.savefig("Cr_Dy_inter.png", bbox_inches = 'tight')

plt.subplot(221)
plt.imshow(grid_Dz.T, extent=(-7.276,7.276,-7.402,7.402), origin='lower', interpolation='none')
plt.title("Cr1_inter")
plt.subplot(222)
plt.imshow(-1*(grid_Dz.T ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr2_inter")
plt.subplot(223)
plt.imshow(-1*np.flip(grid_Dz.T,1 ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr3_inter")
plt.subplot(224)
plt.imshow(np.flip(grid_Dz.T, (0,1) ), extent=(-7.276,7.276,-7.402,7.402), origin='lower')
plt.title("Cr4_Dz_inter")

plt.gcf().set_size_inches(6, 6)

plt.savefig("Cr_Dz_inter.png", bbox_inches = 'tight')
