#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red")
echo "set terminal png truecolor
      set output \"./EOS.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{EOS}\"
      set xrange [ 0 : 1500 ]
      set xlabel \"Energy\"
      set ylabel \"Pressure\"
	  plot \"./M_delta=1112/pa1&gr=0&s=0.15/EOS.txt\" u 2:3 t \"M=1112\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"./M_delta=1112/pa1&gr=gr&s=0.15/EOS.txt\" u 2:3 t \"M=1112\" w l lw $lw lc \"${color[1]}\", \\
           \"./M_delta=1232/pa1&gr=0&s=0.15/EOS.txt\" u 2:3 t \"M=1232\" w l lw $lw dt 2 lc \"${color[2]}\", \\
           \"./M_delta=1232/pa1&gr=gr&s=0.15/EOS.txt\" u 2:3 t \"M=1232\" w l lw $lw lc \"${color[2]}\", \\
           \"./M_delta=1352/pa1&gr=0&s=0.15/EOS.txt\" u 2:3 t \"M=1352\" w l lw $lw dt 2 lc \"${color[6]}\", \\
           \"./M_delta=1352/pa1&gr=gr&s=0.15/EOS.txt\" u 2:3 t \"M=1352\" w l lw $lw lc \"${color[6]}\",\\
           \"./M_delta=1112/pa1&gr=0&s=0/EOS.txt\" u 2:3 t \"M=1112\" w l lw $lw dt 2 lc \"${color[1]}\",\\
           \"./M_delta=1112/pa1&gr=gr&s=0/EOS.txt\" u 2:3 t \"M=1112\" w l lw $lw lc \"${color[1]}\",\\
           \"./M_delta=1232/pa1&gr=0&s=0/EOS.txt\" u 2:3 t \"M=1232\" w l lw $lw dt 2 lc \"${color[2]}\",\\
           \"./M_delta=1232/pa1&gr=gr&s=0/EOS.txt\" u 2:3 t \"M=1232\" w l lw $lw lc \"${color[2]}\",\\
           \"./M_delta=1352/pa1&gr=0&s=0/EOS.txt\" u 2:3 t \"M=1352\" w l lw $lw dt 2 lc \"${color[6]}\",\\
           \"./M_delta=1352/pa1&gr=gr&s=0/EOS.txt\" u 2:3 t \"M=1352\" w l lw $lw lc \"${color[6]}\"" | gnuplot
