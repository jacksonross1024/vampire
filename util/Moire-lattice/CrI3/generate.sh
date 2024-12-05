#!/bin/bash
# make
./main 1.1 9.9 --dmi 0.25 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.25-10-AF.txt

../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.25-10-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.25-10-AF.data


./main 1.1 9.9 --dmi 0.25 7.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.25-7-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.25-7-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.25-7-AF.data


./main 1.1 9.9 --dmi 0.25 4.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.25-4-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.25-4-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.25-4-AF.data


./main 1.1 9.9 --dmi 0.30 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.3-10-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.3-10-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.3-10-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.3-10-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.3-10-AF.data


./main 1.1 9.9 --dmi 0.3 7.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.3-7-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.3-7-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.3-7-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.3-7-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.3-7-AF.data

./main 1.1 9.9 --dmi 0.3 4.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.3-4-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.3-4-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.3-4-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.3-4-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.3-4-AF.data


./main 1.1 9.9 --dmi 0.2 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.2-10-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.2-10-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-1.1-lowres-5K-0.2-10-FC.txt
mv cells-00000001.txt cells-1.1-lowres-5K-0.2-10-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.2-10-AF.data


./main 1.1 9.9 --dmi 0.2 7.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.2-7-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.2-7-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.2-7-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.2-7-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.2-7-AF.data


./main 1.1 9.9 --dmi 0.2 4.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-1.1-5K-0.2-4-FC.txt
mv cells-00000001.txt cells-1.1-5K-0.2-4-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-1.1-5K-0.2-4-FC.txt
mv cells-00000001.txt cells-lowres-1.1-5K-0.2-4-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-1.1-5K-0.2-4-AF.data


./main 2.0 9.9 --dmi 0.25 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-2.0-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-2.0-5K-0.25-10-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-2.0-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-lowres-2.0-5K-0.25-10-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-2.0-5K-0.25-10-AF.data


./main 2.0 9.9 --dmi 0.25 7.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-2.0-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-2.0-5K-0.25-7-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-2.0-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-lowres-2.0-5K-0.25-7-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-2.0-5K-0.25-7-AF.data


./main 2.0 9.9 --dmi 0.25 4.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-2.0-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-2.0-5K-0.25-4-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-2.0-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-lowres-2.0-5K-0.25-4-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-2.0-5K-0.25-4-AF.data



./main 0.5 9.9 --dmi 0.25 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-0.5-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-0.5-5K-0.25-10-AF.txt

../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-0.5-5K-0.25-10-FC.txt
mv cells-00000001.txt cells-lowres-0.5-5K-0.25-10-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-0.5-5K-0.25-10-AF.data


./main 0.5 9.9 --dmi 0.25 7.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-0.5-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-0.5-5K-0.25-7-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-0.5-5K-0.25-7-FC.txt
mv cells-00000001.txt cells-lowres-0.5-5K-0.25-7-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-0.5-5K-0.25-7-AF.data


./main 0.5 9.9 --dmi 0.25 4.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
cat input-mc > input
mpirun -np 16 ../../../vampire-parallel
cat input-dipole > input
mpirun -np 16 ../../../vampire-parallel

../../vdc/vdc --cells --cell-size = 10,10,30
mv cells-00000000.txt cells-0.5-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-0.5-5K-0.25-4-AF.txt
../../vdc/vdc --cells --cell-size = 100,100,30
mv cells-00000000.txt cells-lowres-0.5-5K-0.25-4-FC.txt
mv cells-00000001.txt cells-lowres-0.5-5K-0.25-4-AF.txt
paste cells-coords.cfg cells-00000000.cfg > dipole-cells-0.5-5K-0.25-4-AF.data


