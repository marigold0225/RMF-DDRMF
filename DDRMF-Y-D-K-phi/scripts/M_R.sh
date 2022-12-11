#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red")
echo "set terminal png truecolor
      set output \"./gnuplot/M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{U_k=-120}\"
      set xrange [ 8.0 : 16.0 ]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./L=47.2/M_R.txt\" u 4:3 t \"L = 47.2\" w l lw $lw lc \"${color[0]}\", \\
           \"./L=50/M_R.txt\" u 4:3 t \"L = 50\" w l lw $lw lc \"${color[1]}\", \\
           \"./L=60/M_R.txt\" u 4:3 t \"L = 60\" w l lw $lw lc \"${color[2]}\", \\
           \"./L=70/M_R.txt\" u 4:3 t \"L = 70\" w l lw $lw lc \"${color[3]}\", \\
           \"./L=80/M_R.txt\" u 4:3 t \"L = 80\" w l lw $lw lc \"${color[4]}\", \\
           \"./L=90/M_R.txt\" u 4:3 t \"L = 90\" w l lw $lw lc \"${color[5]}\",\\
           \"./L=100/M_R.txt\" u 4:3 t \"L = 100\" w l lw $lw lc \"${color[6]}\",\\
           \"./L=110/M_R.txt\" u 4:3 t \"L = 110\" w l lw $lw lc \"${color[7]}\"" | gnuplot
