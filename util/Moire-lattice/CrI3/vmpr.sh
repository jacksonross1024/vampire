
#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=CrI3-dipole
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1

#SBATCH --account=e89-ed_p
#SBATCH --partition=standard
#SBATCH --qos=standard
module load PrgEnv-gnu
export OMP_NUM_THREADS=1

mc > input

srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel
../../vdc/vdc --cells --cell-size 10,10,30

cat input-dipole > input
srun --distribution=block:cyclic:block /work/e89/e89/jross71/vampire/moire/vampire-parallel

paste cells-coords.cfg cells-00000000.cfg > dipole-cells-AF.data
rm CrI3.ucf

