
from math import sqrt
import numpy as np
import os 
import matplotlib.pyplot as plt
import sys


for angle in ["0.5", "1.1", "2.0"]:
    for J in ["0.2", "0.25", "0.3"]:
        for DMI in ["4", "7", "10"]:
            for J_twist in ["1.0", "0.9"]:
                for J_intra in ["1.0", "1.1"]:
                    for index in ["1", "3"]:
                        param = "param-" + angle + "-" + J + "-" + DMI + "-" + J_twist + "-" + J_intra
                        path = os.getcwd() + "/" + param + "/cells-0000000" + index + ".txt"
                        path = os.getcwd() + "/dipole-cells-test.data"
                        print(path)

                        cell_data = np.genfromtxt(path, dtype=float, usecols = (3,4,5,19,0), unpack=True, delimiter = '\t')
                        # y_pos = np.genfromtxt(path, dtype=None, usecols = 2)
                        # m_z = np.genfromtxt(path, dtype=None, usecols = 26)
                        # M = np.genfromtxt(path, dtype=None, 

                        print(np.size(cell_data[0]))
                        M_max = 1
                        M_corr = [[[0,0]]]
                        print(cell_data[1])

                        x_cells = int(100/5)
                        y_cells = int(100/5)
                        print(x_cells*y_cells*5)
                        cell_size = 5.0
                        for i in range(x_cells):
                            M_corr.append([[0,0]])
                            for j in range(y_cells):
                                M_corr[i].append([0,0])
                                for k in range(int(x_cells*y_cells*5)):
                                    if(cell_data[2][k] == 19.62):
                                        dx = i*5 + 2.5 - 0.1*cell_data[0][k]
                                        dy = j*5 + 2.5 - 0.1*cell_data[1][k]
                                        dr = sqrt(dx*dx+dy*dy)
                                        #print(str(dx) + ", " + str(dy))
                                        if(dr > 10 ):
                                            w_i = 1/pow(dr,3.0)
                                        else:
                                            w_i = 1
                                            continue
                                        # print(dr, w_i, cell_data[2][k], cell_data[3][k])

                                        M_corr[i][j][0] += cell_data[3][k]*w_i
                                        M_corr[i][j][1] += w_i 
                                        print(w_i)

                        M_max = 0
                        counter = 0
                        corr_max = 0.0
                        orig_out = sys.stdout
                        f = open('corr-dipole.txt', 'w')
                        sys.stdout = f
                        for i in range(x_cells):
                            for j in range(y_cells):#
                                if(M_corr[i][j][1] > 0):
                                    M_val = M_corr[i][j][0] / M_corr[i][j][1]
                                    M_corr[i][j][0] = M_val
                                    if(M_max < M_val):
                                        M_max = M_val
                                    #cell_data[4][counter] = M_val
                                    counter+=1
                                    print(i, j, M_val)

                        sys.stdout = orig_out
                        print("Corr max: " + str(M_max) + " w " + str(counter))
                        f.close()
