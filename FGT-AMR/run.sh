

mkdir -p spin-resistance 

cat input-FC > input
mpirun -np 16 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,21
cat cells-00000000.txt > cells-FC.txt

cat input-SOT > input
mpirun -np 16 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,21

cat input-amr > input
mpirun -np 16 ../vampire-parallel
