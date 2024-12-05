#!/bin/bash


for angle in "0.5" "1.1" "2.0"
do
    for J in "0.2" "0.25" "0.3"
    do
        for DMI in "7"
        do
            for J_twist in "0.975" "0.95" "0.9"
            do 
            ./main $angle 9.9 --dmi $J $DMI $J_twist 
            mv config_energy_atomic.txt config-energy-atomic-$angle-$J-$DMI-$J_twist.txt
            mv config_energy_cells.txt config-energy-cells-$angle-$J-$DMI-$J_twist.txt
            cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

            cat input-mc > input
            mpirun -np 16 ../../../vampire-parallel
            ../../vdc/vdc --cells --cell-size 20,20,30
            mv cells-00000000.txt cells-$angle-5K-$J-$DMI-$J_twist-FC.txt
            mv cells-00000001.txt cells-$angle-5K-$J-$DMI-$J_twist-AF.txt
            ../../vdc/vdc --cells --cell-size 100,100,30
            mv cells-00000000.txt cells-lowres-$angle-5K-$J-$DMI-$J_twist-FC.txt
            mv cells-00000001.txt cells-lowres-$angle-5K-$J-$DMI-$J_twist-AF.txt
            cat input-dipole > input
            mpirun -np 16 ../../../vampire-parallel
            paste cells-coords.cfg cells-00000000.cfg > dipole-cells-$angle-5K-$J-$DMI-$J_twist-AF.data
            done
        done
    done
done
