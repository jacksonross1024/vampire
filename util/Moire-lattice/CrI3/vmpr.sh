#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=CrI3-dipole
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1

#SBATCH --account=e89-ed_p
#SBATCH --partition=standard
#SBATCH --qos=taskfarm

module load PrgEnv-gnu
export OMP_NUM_THREADS=1

input-mc-20K > input
srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel
../../../../vdc/vdc --cells --cell-size 50,50,30
mv cells-00000000.txt cells-20K-AF.txt
mv cells-00000001.txt cells-20K-FC.txt
../../../../vdc/vdc --cells --cell-size 10,10,30
mv cells-00000000.txt cells-20K-AF-hr.txt
mv cells-00000001.txt cells-20K-FC-hr.txt

cat input-dipole-20K > input
srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel

paste cells-coords.cfg cells-00000000.cfg > dipole-cells-20K-AF.data

input-mc-5K > input
srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel
../../../../vdc/vdc --cells --cell-size 50,50,30
mv cells-00000000.txt cells-5K-AF.txt
mv cells-00000001.txt cells-5K-FC.txt
../../../../vdc/vdc --cells --cell-size 10,10,30
mv cells-00000000.txt cells-5K-AF-hr.txt
mv cells-00000001.txt cells-5K-FC-hr.txt

cat input-dipole-5K > input
srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel

paste cells-coords.cfg cells-00000000.cfg > dipole-cells-5K-AF.data

