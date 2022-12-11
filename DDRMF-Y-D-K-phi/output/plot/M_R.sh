#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{M-R-150}\"
      set xrange [ 10 : 15.0 ]
	  set yrange [0:2.5]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"../DDME1/npeuY/M_R.txt\" u 4:3 t \"DDME1-npeuY\" w l lw $lw  lc \"${color[1]}\", \\
           \"../DDME1/delta=-150/npeuYD/M_R.txt\" u 4:3 t \"DDME1-npeuYD\" w l lw $lw dt '...' lc \"${color[1]}\", \\
           \"../DDME1/delta=-150/Uk=-120/M_R.txt\" u 4:3 t \"DDME1-120\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../DDME1/delta=-150/Uk=-140/M_R.txt\" u 4:3 t \"DDME1-140\" w l lw $lw dt '..-' lc \"${color[1]}\", \\
           \"../DDME1/delta=-150/Uk=-160/M_R.txt\" u 4:3 t \"DDME1-160\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"../DDME2/npeuY/M_R.txt\" u 4:3 t \"DDME2-npeuY\" w l lw $lw lc \"${color[2]}\", \\
           \"../DDME2/delta=-150/npeuYD/M_R.txt\" u 4:3 t \"DDME2-npeuYD\" w l lw $lw dt '...' lc \"${color[2]}\", \\
           \"../DDME2/delta=-150/Uk=-120/M_R.txt\" u 4:3 t \"DDME2--120\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../DDME2/delta=-150/Uk=-140/M_R.txt\" u 4:3 t \"DDME2--140\" w l lw $lw dt '..-' lc \"${color[2]}\", \\
           \"../DDME2/delta=-150/Uk=-160/M_R.txt\" u 4:3 t \"DDME2--160\" w l lw $lw dt 2 lc \"${color[2]}\",\\
           \"../DD2/npeuY/M_R.txt\" u 4:3 t \"DD2-npeuY\" w l lw $lw lc \"${color[3]}\",\\
           \"../DD2/delta=-150/npeuYD/M_R.txt\" u 4:3 t \"DD2-npeuYD\" w l lw $lw dt '...' lc \"${color[3]}\",\\
           \"../DD2/delta=-150/Uk=-120/M_R.txt\" u 4:3 t \"DD2--120\" w l lw $lw dt '--' lc \"${color[3]}\",\\
           \"../DD2/delta=-150/Uk=-140/M_R.txt\" u 4:3 t \"DD2--140\" w l lw $lw dt '..-' lc \"${color[3]}\",\\
           \"../DD2/delta=-150/Uk=-160/M_R.txt\" u 4:3 t \"DD2--160\" w l lw $lw dt 2 lc \"${color[3]}\",\\
           \"../DDMEX/npeuY/M_R.txt\" u 4:3 t \"DDMEX-npeuY\" w l lw $lw lc \"${color[6]}\",\\
           \"../DDMEX/delta=-150/npeuYD/M_R.txt\" u 4:3 t \"DDMEX-npeuYD\" w l lw $lw dt '...' lc \"${color[6]}\",\\
           \"../DDMEX/delta=-150/Uk=-120/M_R.txt\" u 4:3 t \"DDMEX--120\" w l lw $lw dt '--' lc \"${color[6]}\",\\
           \"../DDMEX/delta=-150/Uk=-140/M_R.txt\" u 4:3 t \"DDMEX--140\" w l lw $lw dt '..-' lc \"${color[6]}\",\\
           \"../DDMEX/delta=-150/Uk=-160/M_R.txt\" u 4:3 t \"DDMEX--160\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
