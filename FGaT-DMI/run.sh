
set -e 


cat input-FC > input 
echo "sim:applied-field-strength = 0.3" >> input
mpirun -np 8 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,50

mv *chk cells*.txt 300mT/


cat input-FC > input 
echo "sim:applied-field-strength = 0.35" >> input
mpirun -np 8 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,50

mv *chk cells*.txt 350mT/



