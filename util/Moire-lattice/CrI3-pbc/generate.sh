#!/bin/bash
make
./main 0.0 9.9 --dmi 0.2 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.2J

./main 0.0 9.9 --dmi 0.225 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.225J

./main 0.0 9.9 --dmi 0.25 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.25J

./main 0.0 9.9 --dmi 0.275 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.275J

./main 0.5 9.9 --dmi 0.2 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.2J-0.5

./main 0.5 9.9 --dmi 0.225 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.225J-0.5

./main 0.5 9.9 --dmi 0.25 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.25J-0.5

./main 0.5 9.9 --dmi 0.275 10.0
cat header.ucf atom_positions.xyz header_interactions.ucf  interactions.ucf > CrI3.ucf
mpirun -np 8 ../../../vampire-parallel
mv output output-0.275J-0.5