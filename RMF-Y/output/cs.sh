#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red")
echo "set terminal png truecolor
      set output \"./cs.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{CS^2}\"
      set xrange [ 0 : 2 ]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./FSUgold/EOS.txt\" u 1:4 t \"FSUgold\" w l lw $lw lc \"${color[0]}\", \\
           \"./IUFSU/EOS.txt\" u 1:4 t \"IUFSU\" w l lw $lw lc \"${color[1]}\", \\
           \"./TM1/EOS.txt\" u 1:4 t \"TM1\" w l lw $lw lc \"${color[2]}\"" | gnuplot
