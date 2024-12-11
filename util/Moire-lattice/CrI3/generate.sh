#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=CrI3-setup
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --account=e89-ed_p
#SBATCH --partition=serial
#SBATCH --qos=serial
module load PrgEnv-gnu
export OMP_NUM_THREADS=1

rm -r param*
for angle in "1.1"
do
    for J in "0.2" "0.25" "0.3" "0.35"
    do
        for DMI in  "4" "7" "10"
        do
            for J_twist in "1.0" "0.9"
            do 
                for J_intra in "1.0" "1.2"
                do
                mkdir param-$angle-$J-$DMI-$J_twist-$J_intra
                ../main $angle 9.9 --dmi $J $DMI $J_twist $J_intra
                mv config_energy_atomic.txt  param-$angle-$J-$DMI-$J_twist-$J_intra/config-energy-atomic-$angle-$J-$DMI-$J_twist_$J_intra.txt
                mv config_energy_cells.txt  param-$angle-$J-$DMI-$J_twist-$J_intra/config-energy-cells-$angle-$J-$DMI-$J_twist_$J_intra.txt
                mv header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf  param-$angle-$J-$DMI-$J_twist-$J_intra/

                cp input-mc input-dipole CrI3.mat  param-$angle-$J-$DMI-$J_twist-$J_intra/
                cd param-$angle-$J-$DMI-$J_twist-$J_intra/
                cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf
                
                #sbatch vmpr.sh
                cd ../

                done
            done
        done
    done
done
