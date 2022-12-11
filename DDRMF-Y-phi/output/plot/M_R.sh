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
      set xrange [ 10 : 13.0 ]
	  set yrange [0.25:2.2]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./PKDD/1/M_R.txt\" u 4:3 t \"PKDD\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"./DDME1/1/M_R.txt\" u 4:3 t \"DDEM1\" w l lw $lw dt 2 lc \"${color[2]}\", \\
		   \"./DDME2/1/M_R.txt\" u 4:3 t \"DDEM2\" w l lw $lw dt 2 lc \"${color[3]}\", \\
           \"./TW99/1/M_R.txt\" u 4:3 t \"TW99\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
