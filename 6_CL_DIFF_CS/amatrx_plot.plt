#set xlabel " strain"
#set ylabel " stress"
m="./amatrx_plot.out"
set terminal x11 0
#set terminal png
#set output 'passive_vs_active_shear.png'
#set terminal postscript eps enhanced color font 'Helvetica,12'
#set output 'passive_vs_active_shear.eps'
#set nokey
set key bottom right
#set grid
set pm3d map
set palette model HSV rgbformulae 0,5,15
#set pal gray
#set log y
#set log x
#plot m using 1:3  with lines title 'active'
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
splot m matrix
