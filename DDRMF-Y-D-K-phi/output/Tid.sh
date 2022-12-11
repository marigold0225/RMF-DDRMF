#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./Tid.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{Tid}\"
      set xrange [1 : 3]
	  set logscale y
	  set yrange [1:2000]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"./TD/DDME1/npeuY/output_tov.dat\" u 3:4 t \"DDME1-npeuY\" w l lw $lw  lc \"${color[1]}\", \\
           \"./TD/DDME1/Ud=-100/npeuYD/output_tov.dat\" u 3:4 t \"DDME1-npeuYD\" w l lw $lw dt '...' lc \"${color[1]}\", \\
           \"./TD/DDME1/Ud=-100/Uk=-120/output_tov.dat\" u 3:4 t \"DDME1-120\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"./TD/DDME1/Ud=-100/Uk=-140/output_tov.dat\" u 3:4 t \"DDME1-140\" w l lw $lw dt '..-' lc \"${color[1]}\", \\
           \"./TD/DDME2/npeuY/output_tov.dat\" u 3:4 t \"DDME2-npeuY\" w l lw $lw lc \"${color[2]}\", \\
           \"./TD/DDME2/Ud=-100/npeuYD/output_tov.dat\" u 3:4 t \"DDME2-npeuYD\" w l lw $lw dt '...' lc \"${color[2]}\", \\
           \"./TD/DDME2/Ud=-100/Uk=-120/output_tov.dat\" u 3:4 t \"DDME2--120\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"./TD/DDME2/Ud=-100/Uk=-140/output_tov.dat\" u 3:4 t \"DDME2--140\" w l lw $lw dt '..-' lc \"${color[2]}\", \\
           \"./TD/DDME2/Ud=-100/Uk=-160/output_tov.dat\" u 3:4 t \"DDME2--160\" w l lw $lw dt 2 lc \"${color[2]}\",\\
           \"./TD/DDMEX/npeuY/output_tov.dat\" u 3:4 t \"DDMEX-npeuY\" w l lw $lw lc \"${color[6]}\",\\
           \"./TD/DDMEX/Ud=-100/npeuYD/output_tov.dat\" u 3:4 t \"DDMEX-npeuYD\" w l lw $lw dt '...' lc \"${color[6]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-120/output_tov.dat\" u 3:4 t \"DDMEX--120\" w l lw $lw dt '--' lc \"${color[6]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-140/output_tov.dat\" u 3:4 t \"DDMEX--140\" w l lw $lw dt '..-' lc \"${color[6]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-160/output_tov.dat\" u 3:4 t \"DDMEX--160\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
