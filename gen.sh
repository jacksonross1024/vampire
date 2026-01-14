
for integrator in  "llg-heun" "llg-midpoint"
do
for damping in "1" "0.1" "0.01" "0.001"
do 
    for temp in "1"  "40"  "80" "120" "160" "200"
    do
cat FGaT-mara.mat > FGaT-mara-dampvar.mat
echo "material[1]:damping-constant=" $damping >> FGaT-mara-dampvar.mat
echo "material[2]:damping-constant=" $damping >> FGaT-mara-dampvar.mat

cat input-timestep > input
echo "sim:equilibration-temperature=" $temp >> input
echo "sim:temperature=" $temp >> input

echo "sim:integrator = " $integrator >> input
mpirun -np 8 ./vampire-parallel
cat output > output-$temp-K-$damping-alpha-$integrator-drift

done
done
done