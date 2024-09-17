#!/bin/bash
make
./main
cat header.ucf atom_positions.ucf header_interactions.ucf  interactions.ucf > CrI3.ucf
