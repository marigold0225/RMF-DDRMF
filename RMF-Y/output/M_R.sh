#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red")
echo "set terminal png truecolor
      set output \"./M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{U_k=-120}\"
      set xrange [ 10 : 12.0 ]
	  set yrange [0.2:1.4]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./FSUgold/M_R.dat\" u 2:3 t \"FSUgold\" w l lw $lw lc \"${color[0]}\", \\
           \"./IUFSU/M_R.txt\" u 5:4 t \"IUFSU\" w l lw $lw lc \"${color[1]}\", \\
           \"./TM1/M_R.txt\" u 5:4 t \"TM1\" w l lw $lw lc \"${color[2]}\"" | gnuplot
