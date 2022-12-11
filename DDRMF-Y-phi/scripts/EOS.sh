#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red")
echo "set terminal png truecolor
      set output \"./gnuplot/EOS.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{U_k=-120}\"
      set xrange [ 0 : 1500 ]
      set xlabel \"Energy\"
      set ylabel \"Pressure\"
      plot \"./L=47.2/EOS.txt\" u 2:3 t \"L = 47.2\" w l lw $lw lc \"${color[0]}\", \\
           \"./L=50/EOS.txt\" u 2:3 t \"L = 50\" w l lw $lw lc \"${color[1]}\", \\
           \"./L=60/EOS.txt\" u 2:3 t \"L = 60\" w l lw $lw lc \"${color[2]}\", \\
           \"./L=70/EOS.txt\" u 2:3 t \"L = 70\" w l lw $lw lc \"${color[3]}\", \\
           \"./L=80/EOS.txt\" u 2:3 t \"L = 80\" w l lw $lw lc \"${color[4]}\", \\
           \"./L=90/EOS.txt\" u 2:3 t \"L = 90\" w l lw $lw lc \"${color[5]}\",\\
           \"./L=100/EOS.txt\" u 2:3 t \"L = 100\" w l lw $lw lc \"${color[6]}\",\\
           \"./L=110/EOS.txt\" u 2:3 t \"L = 110\" w l lw $lw lc \"${color[7]}\"" | gnuplot
