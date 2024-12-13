#!/bin/bash

for J_twist in "1.0" 
do 
    for J_intra in "1.0" "1.2"
    do 

for angle in "0.5" "1.1" "2.0"
do
    for J in "0.2" "0.25" "0.3" "0.35"
    do
        for DMI in "2" "4" "8"
        do

                    ./main $angle 9.9 $J $DMI $J_twist $J_intra
                   # mv config_energy_atomic.txt config-energy-atomic-$angle-$J-$DMI-$J_twist-$J_intra-Jinter.txt
                    mv config_energy_cells.txt config-energy-cells-$angle-$J-$DMI-$J_twist-$J_intra-Jintra.txt
                    #cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

                    # cat input-mc > input
                    # mpirun -np 16 ../../../vampire-parallel
                    # ../../vdc/vdc --cells --cell-size 10,10,30
                    # mv cells-00000000.txt cells-$angle-5K-$J-$DMI-$J_twist-$J_intra-serial-FC.txt
                    # mv cells-00000001.txt cells-$angle-5K-$J-$DMI-$J_twist-$J_intra-serial-AF.txt
                    # ../../vdc/vdc --cells --cell-size 100,100,30
                    # mv cells-00000000.txt cells-lowres-$angle-5K-$J-$DMI-$J_twist-$J_intra-serial-FC.txt
                    # mv cells-00000001.txt cells-lowres-$angle-5K-$J-$DMI-$J_twist-$J_intra-serial-AF.txt
                    # cat input-dipole > input
                    # mpirun -np 16 ../../../vampire-parallel
                    # paste cells-coords.cfg cells-00000000.cfg > dipole-cells-$angle-5K-$J-$DMI-$J_twist-$J_intra-serial-AF.data
                done
            done
        done
    done
done
