


set ylabel "y position (nm)"
set xlabel "x position (nm)"

set ytics 50 out 
set xtics 50 out 
set mytics 5 
set mxtics 5

chck(S,s, L, l, x) = (S == s && L == l) ? (x):(-100.0)

set colorbox
set cblabel "meV/mu_B^2"
set style fill solid noborder 
normalise = 2.98*2.98


set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.6'#614199')

file = sprintf("config_energy_atomic")
print(file.".txt")

set xrange [190:350]
set yrange [475:590]

set term pngcairo size 1600,1600
set output sprintf("%s-J-intra-twist1NN-3NN.png", file)
set multiplot layout 4,4


print("plotting loaded data")

layer = 1
species = 1

set title sprintf("1NN/interaction")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($7)/normalise/$6) w boxxy palette notitle 

layer = 1
species = 1

set title sprintf("1NN, net")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($7)/normalise) w boxxy palette notitle 


set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
species = 12
set title sprintf("<1NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($14+$24)/normalise/3/($13+$23)) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
species = 12
set title sprintf("<1NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($14+$24)/normalise/($13+$23)) w boxxy palette notitle 


layer = 1
species = 1
set title sprintf("2NN/interaction")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($12)/normalise/$11) w boxxy palette notitle 

layer = 1
species = 1
set title sprintf("2NN, net")
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer+1)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($12)/normalise) w boxxy palette notitle 

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
species = 12
set title sprintf("<2NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($15+$25)/normalise/6/($13+$23)) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
species = 12
set title sprintf("<2NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($15+$25)/normalise/($13+$23)) w boxxy palette notitle 


llayer = 1
species = 1
set title sprintf("3NN/interaction")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise/$16) w boxxy palette notitle 

layer = 1
species = 1
set title sprintf("3NN, net")
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer+1)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise) w boxxy palette notitle 

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
species = 12
set title sprintf("<3NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($16+$26)/normalise/($13+$23)/3) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
species = 12
set title sprintf("<3NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($16+$26)/normalise/($13+$23)) w boxxy palette notitle 



layer = 1
species = 1
set title sprintf("<J_{intra}>/interaction")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17+$12+$7)/normalise/($16+$11+$6)) w boxxy palette notitle 

layer = 1
species = 2
set title sprintf("<J_{intra}>")
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17+$12+$7)/normalise) w boxxy palette notitle 

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
species = 12
set title sprintf("<J_{intra}>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($16+$26+$14+$15+$25+$24)/normalise/3/($13+$23)/6/3) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
species = 12
set title sprintf("<J_{intra}>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($16+$26+$14+$15+$25+$24)/normalise/($13+$23)) w boxxy palette notitle 


unset multiplot 

set term pngcairo size 800,1600
set output sprintf("%s-J-inter-twist-1NN-3NN-cells.png", file)
set multiplot layout 4,2

set palette defined (-0.04 '#81B1CB', 0.0 '#f5f0f7', 0.06 '#7A5FA9',  0.12 '#614199')
species = 12
set title sprintf("<1NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($18+$28)/normalise/($17+$27)) w boxxy palette notitle 

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.2 '#614199')
species = 12
set title sprintf("<1NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($18+$28)/normalise/($13+$23)) w boxxy palette notitle 


set palette defined (-0.08 '#81B1CB', 0.0 '#f5f0f7', 0.05 '#7A5FA9',  0.1 '#614199')
species = 12
set title sprintf("<2NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($20+$30)/normalise/($19+$29)) w boxxy palette notitle 

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.2 '#614199')
species = 12
set title sprintf("<2NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($20+$30)/normalise/($13+$23)) w boxxy palette notitle 



set palette defined (-0.05 '#81B1CB', 0.0 '#f5f0f7')
species = 12
set title sprintf("<3NN>/<interactions>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($22+$32)/normalise/($21+$31)) w boxxy palette notitle 

set palette defined (-0.3 '#81B1CB', 0.0 '#f5f0f7')
species = 12
set title sprintf("<3NN>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($22+$32)/normalise/($13+$23)) w boxxy palette notitle 


set palette defined (-0.02 '#81B1CB', 0.0 '#f5f0f7', 0.005 '#7A5FA9',  0.015 '#614199')
species = 12
set title sprintf("<J_{inter}>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($22+$32+$20+$18+$30+$28)/normalise/($21+$31+$19+$17+$29+$27)) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.075 '#7A5FA9',  0.15 '#614199')
species = 12
set title sprintf("<J_{inter}>, species: %.f", species)
plot "config_energy_cells.txt" u ($1):($2):(0.693):(0.6002):(($22+$32+$20+$18+$30+$28)/normalise/($13+$23)) w boxxy palette notitle 


unset multiplot 


set term pngcairo size 1600,1600
set output sprintf("%s-J-inter-twist-1NN-3NN-atomistic.png", file)
set multiplot layout 4,4

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
layer = 1
species = 1
set title sprintf("<1NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($37/normalise/$36) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
layer = 1
species = 1
set title sprintf("<1NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($37/normalise) w boxxy palette notitle 

set palette defined (-0.1 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.25 '#614199')
layer = 1
species = 2
set title sprintf("<1NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($37/normalise/$36) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.4 '#614199')
layer = 1
species = 2
set title sprintf("<1NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($37/normalise) w boxxy palette notitle 


set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.2 '#614199')
layer = 1
species = 1
set title sprintf("<2NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($42/normalise/$41) w boxxy palette notitle 

set palette defined (-0.3 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.6 '#614199')
layer = 1
species = 1
set title sprintf("<2NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($42/normalise) w boxxy palette notitle 

set palette defined (-0.2 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.2 '#614199')
layer = 1
species = 2
set title sprintf("<2NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($42/normalise/$41) w boxxy palette notitle 

set palette defined (-0.3 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.6 '#614199')
layer = 1
species = 2
set title sprintf("<2NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($42/normalise) w boxxy palette notitle 

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.15 '#614199')
layer = 1
species = 1
set title sprintf("<3NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($47)/normalise/$46) w boxxy palette notitle 

set palette defined (-0.7 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.3 '#614199')
layer = 1
species = 1
set title sprintf("<3NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($47)/normalise) w boxxy palette notitle 

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.15 '#614199')
layer = 1
species = 2
set title sprintf("<3NN>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($47)/normalise/$46) w boxxy palette notitle 

set palette defined (-0.7 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.3 '#614199')
layer = 1
species = 2
set title sprintf("<3NN>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($47)/normalise) w boxxy palette notitle 


set palette defined (-0.06 '#81B1CB', 0.0 '#f5f0f7', 0.05 '#7A5FA9',  0.1 '#614199')
layer = 1
species = 1
set title sprintf("<J_{inter}>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise/($36+$41+$46)) w boxxy palette notitle 

set palette defined (-0.5 '#81B1CB', 0.0 '#f5f0f7', 0.2 '#7A5FA9',  0.6 '#614199')
species = 1
set title sprintf("<J_{inter}>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

set palette defined (-0.06 '#81B1CB', 0.0 '#f5f0f7', 0.05 '#7A5FA9',  0.1 '#614199')
species = 2
set title sprintf("<J_{inter}>/<interactions>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise/($36+$41+$46)) w boxxy palette notitle 

set palette defined (-0.5 '#81B1CB', 0.0 '#f5f0f7', 0.2 '#7A5FA9',  0.6 '#614199')
species = 2
set title sprintf("<J_{inter}>, species: %.f", species)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

unset multiplot