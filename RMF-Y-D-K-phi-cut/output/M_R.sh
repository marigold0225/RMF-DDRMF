#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{M-R}\"
      set xrange [ 8.0 : 16.0 ]
	  set yrange [0:2.2]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./FSUGOLD/Uk=-120/c_s=0/M_R.txt\" u 4:3 t \"c_s=0\" w l lw $lw lc \"${color[1]}\", \\
           \"./FSUGOLD/Uk=-140/c_s=0/M_R.txt\" u 4:3 t \"c_s=0\" w l lw $lw lc \"${color[2]}\", \\
           \"./FSUGOLD/Uk=-160/c_s=0/M_R.txt\" u 4:3 t \"c_s=0\" w l lw $lw lc \"${color[3]}\", \\
           \"./FSUGOLD/Uk=-120/c_s=0.15/M_R.txt\" u 4:3 t \"c_s=0.15\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"./FSUGOLD/Uk=-140/c_s=0.15/M_R.txt\" u 4:3 t \"c_s=0.15\" w l lw $lw dt 2 lc \"${color[2]}\", \\
           \"./FSUGOLD/Uk=-160/c_s=0.15/M_R.txt\" u 4:3 t \"c_s=0.15\" w l lw $lw dt 2 lc \"${color[3]}\"" | gnuplot
