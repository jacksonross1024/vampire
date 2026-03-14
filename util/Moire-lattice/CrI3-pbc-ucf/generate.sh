#!/bin/bash



./main 1.1 9.9 --dmi 0.085 1 1 5.5 1.45 1 1
cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

mpirun -np 16 ../../../vampire-parallel
../../vdc/vdc --cells --cell-size = 10,10,30
mkdir xcross-1-AFM
mv cells* xcross-1-AFM

./main 1.1 9.9 --dmi 0.085 1 1 5.5 1.45 3 1
cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

mpirun -np 16 ../../../vampire-parallel
../../vdc/vdc --cells --cell-size = 10,10,30
mkdir xcross-3-AFM
mv cells* xcross-3-AFM


./main 1.1 9.9 --dmi 0.085 0 1 5.5 1.45 1 1
cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

mpirun -np 16 ../../../vampire-parallel
../../vdc/vdc --cells --cell-size = 10,10,30
mkdir xcross-1-AFM-0
mv cells* xcross-1-AFM-0

./main 1.1 9.9 --dmi 0.085 0 1 5.5 1.45 3 1
cat header.ucf atom_positions.xyz header_interactions.ucf interactions.ucf > CrI3.ucf

mpirun -np 16 ../../../vampire-parallel
../../vdc/vdc --cells --cell-size = 10,10,30
mkdir xcross-3-AFM-0
mv cells* xcross-3-AFM-0