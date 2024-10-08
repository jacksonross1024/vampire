#!/bin/bash
make
./main 0.0 9.9 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.0-DMI-3NN/CrI3.ucf

./main 0.0 7.7 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.0-DMI-2NN/CrI3.ucf

./main 0.0 9.9 
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.0-noDMI-3NN/CrI3.ucf

./main 0.0 7.7
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.0-noDMI-2NN/CrI3.ucf

./main 0.5 9.9 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.5-DMI-3NN/CrI3.ucf

./main 0.5 7.7 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.5-DMI-2NN/CrI3.ucf

./main 0.5 9.9 
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.5-noDMI-3NN/CrI3.ucf

./main 0.5 7.7
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-0.5-noDMI-2NN/CrI3.ucf

./main 1.1 9.9 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.1-DMI-3NN/CrI3.ucf

./main 1.1 7.7 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.1-DMI-2NN/CrI3.ucf

./main 1.1 9.9 
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.1-noDMI-3NN/CrI3.ucf

./main 1.1 7.7
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.1-noDMI-2NN/CrI3.ucf

./main 1.41 9.9 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.41-DMI-3NN/CrI3.ucf

./main 1.41 7.7 --dmi
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.41-DMI-2NN/CrI3.ucf

./main 1.41 9.9 
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.41-noDMI-3NN/CrI3.ucf

./main 1.41 7.7
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > A-1.41-noDMI-2NN/CrI3.ucf
