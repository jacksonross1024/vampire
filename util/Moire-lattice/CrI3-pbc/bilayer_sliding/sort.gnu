set palette defined ( 0 '#04233B',\
                      1 '#293E4D',\
                      2 '#20547C',\
                      3 '#77B1DE',\
                      4 '#A9C7DE',\
                      5 '#EFD5B2',\
                      6 '#EFBD79',\
                      7 '#C07F28',\
                      8 '#785F3D',\
                      9 '#5B3605')

set term pngcairo size 900,900

set palette defined (  0 'blue', 1 'white',  2 'red')

a0x = 6.93
a1y = 6.002
a1x = -3.465

set colorbox 

set xlabel "x (A)"
set ylabel "y (A)"
set ytics out nomirror
set xtics out nomirror


a0 = 9.9
set xrange [-a0:a0]
set yrange [-a0:a0]

set output "Jinter.png"
set cbrange [-0.8:0.8]

set style fill solid noborder

set multiplot layout 2,2

plot "Cr1_inter.txt" us 2:3:(0.2):(0.4):($5) with boxxy palette title "Cr1"

plot "Cr2_inter.txt" us 2:3:(0.2):(0.4):($5) with boxxy palette title "Cr2"

plot "Cr3_inter.txt" us ($2):($3):(0.2):(0.4):($5) with boxxy palette title "Cr3"

plot "Cr4_inter.txt" us ($2):($3):(0.2):(0.4):($5) with boxxy palette title "Cr4"

unset multiplot

set output "Dinter.png"
set cbrange [-0.1:0.1]

set style fill solid noborder

set multiplot layout 4,3
file1 = "Cr1_inter.txt"

plot file1 us 2:3:(0.2):(0.4):($6) with boxxy palette title "Cr1_x"
plot file1 us 2:3:(0.2):(0.4):($7) with boxxy palette title "Cr1_y"
plot file1 us 2:3:(0.2):(0.4):($8) with boxxy palette title "Cr1_z"

file2 = "Cr2_inter.txt"
plot file2 us 2:3:(0.2):(0.4):($6) with boxxy palette title "Cr2_x"
plot file2 us 2:3:(0.2):(0.4):($7) with boxxy palette title "Cr2_y"
plot file2 us 2:3:(0.2):(0.4):($8) with boxxy palette title "Cr2_z"

file3 = "Cr3_inter.txt"
plot file3 us 2:3:(0.2):(0.4):($6) with boxxy palette title "Cr3_x"
plot file3 us 2:3:(0.2):(0.4):($7) with boxxy palette title "Cr3_y"
plot file3 us 2:3:(0.2):(0.4):($8) with boxxy palette title "Cr3_z"

file4 = "Cr4_inter.txt"
plot file4 us 2:3:(0.2):(0.4):($6) with boxxy palette title "Cr4_x"
plot file4 us 2:3:(0.2):(0.4):($7) with boxxy palette title "Cr4_y"
plot file4 us 2:3:(0.2):(0.4):($8) with boxxy palette title "Cr4_z"

unset multiplot

set output "Jintra_Cr1.png"
set cbrange [0.6:3]

set style fill solid noborder

set multiplot layout 3,3

plot "Cr1_intra_1.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_2.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_3.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [0.3:0.9]

plot "Cr1_intra_4.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_6.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_8.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [-0.2:0.2]

plot "Cr1_intra_10.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr1_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

unset multiplot


set output "Jintra_Cr2.png"
set cbrange [0.6:3]

set style fill solid noborder

set multiplot layout 3,3

plot "Cr2_intra_1.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_2.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_3.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [0.3:0.9]

plot "Cr2_intra_4.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_6.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_8.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [-0.2:0.2]

plot "Cr2_intra_10.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr2_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

unset multiplot


set output "Jintra_Cr3.png"
set cbrange [0.6:3]

set style fill solid noborder

set multiplot layout 3,3

plot "Cr3_intra_1.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_2.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_3.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [0.3:0.9]

plot "Cr3_intra_4.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_6.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_8.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [-0.2:0.2]

plot "Cr3_intra_10.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr3_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

unset multiplot


set output "Jintra_Cr4.png"
set cbrange [0.6:3]

set style fill solid noborder

set multiplot layout 3,3

plot "Cr4_intra_1.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_2.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_3.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "1NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [0.3:0.9]

plot "Cr4_intra_4.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_6.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_8.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "2NN+/-_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

set cbrange [-0.2:0.2]

plot "Cr4_intra_10.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{90}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{210}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

plot "Cr4_intra_11.txt" us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette title "3NN_{330}",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($5) with boxxy palette notitle "Cr1",\

unset multiplot


set output "Dintra_Cr1_1NN.png"
set cbrange [-0.1:0.1]

set style fill solid noborder

set multiplot layout 3,3

file1 = "Cr1_intra_1.txt"

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file2 = "Cr1_intra_2.txt"
plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file3 = "Cr1_intra_3.txt"
plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

unset multiplot


set output "Dintra_Cr2_1NN.png"
set cbrange [-0.1:0.1]

set style fill solid noborder

set multiplot layout 3,3

file1 = "Cr2_intra_1.txt"

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file2 = "Cr2_intra_2.txt"
plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file3 = "Cr2_intra_3.txt"
plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

unset multiplot


set output "Dintra_Cr3_1NN.png"
set cbrange [-0.1:0.1]

set style fill solid noborder

set multiplot layout 3,3

file1 = "Cr3_intra_1.txt"

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file2 = "Cr3_intra_2.txt"
plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file3 = "Cr3_intra_3.txt"
plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

unset multiplot


set output "Dintra_Cr4_1NN.png"
set cbrange [-0.1:0.1]

set style fill solid noborder

set multiplot layout 3,3

file1 = "Cr4_intra_1.txt"

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file1 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file2 = "Cr4_intra_2.txt"
plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file2 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

file3 = "Cr4_intra_3.txt"
plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette title "Cr1_x",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($6) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette title "Cr1_y",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($7) with boxxy palette notitle "Cr1_x",\

plot file3 us ($9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette title "Cr1_z",\
"" us (-$9*a0x+a1x*$10):($10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us (-$9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\
"" us ($9*a0x-a1x*$10):(-$10*a1y):(0.35):(0.3):($8) with boxxy palette notitle "Cr1_x",\

unset multiplot


