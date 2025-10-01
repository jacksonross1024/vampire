
set -e


##Al2O3
cat input-10ps > input
echo "material:file= FGaT-mara-Al2O3.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-10ps-Al2O3.mat


cat input-100ps > input
echo "material:file= FGaT-mara-Al2O3.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-100ps-Al2O3.mat


cat input-1000ps > input
echo "material:file= FGaT-mara-Al2O3.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-1000ps-Al2O3.mat


##SiO2
cat input-10ps > input
echo "material:file= FGaT-mara-SiO2.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-10ps-SiO2.mat


cat input-100ps > input
echo "material:file= FGaT-mara-SiO2.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-100ps-SiO2.mat


cat input-1000ps > input
echo "material:file= FGaT-mara-SiO2.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-1000ps-SiO2.mat


##hBN
cat input-10ps > input
echo "material:file= FGaT-mara-hBN.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-10ps-hBN.mat


cat input-100ps > input
echo "material:file= FGaT-mara-hBN.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-100ps-hBN.mat


cat input-1000ps > input
echo "material:file= FGaT-mara-hBN.mat" >> input
../vampire-serial
cat vertical_temperature_profile.dat > vertical_temperature_profile-1000ps-hBN.mat