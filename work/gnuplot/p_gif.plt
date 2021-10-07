set cblabel ""
#set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.8333 1 0 0, 1 1 1 1 )
#set term gif animate
set palette gray
 #set terminal postscript eps size 600,480 dashed enhanced color font 'Roman,10'
set terminal pngcairo font "Times,14" enhanced size 800,450


set pm3d 
set pm3d map
set xlabel "{nx}"
set ylabel "{ny"
#set autoscale z


time =  sprintf('./aimp_all.png')
set output time
#set title ttl
#set output time
#set xr[0:40]
#set yr[0:15]
#set cbrange[0.9:1.2]
set size ratio -1
 #set pm3d interpolate 10,10
unset key
str = sprintf("./Lenna_grayscale_pm3d.txt")
splot str with pm3d






