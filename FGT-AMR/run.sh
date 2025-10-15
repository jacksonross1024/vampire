

mkdir -p spin-resistance 
mkdir FC-chk

##field cooling
cat input-FC > input
mpirun -np 16 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,21
cat cells-00000000.txt > cells-FC.txt
cp *chk FC-chk/

##amr
cat input-amr > input
mpirun -np 16 ../vampire-parallel

##sot
cat input-SOT > input
mpirun -np 16 ../vampire-parallel
../util/vdc/vdc --cells --cell-size = 10,10,21

##amr
#cat input-amr > input
#mpirun -np 16 ../vampire-parallel
