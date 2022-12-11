#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./DDME1-Tid_M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{Tid}\"
      set xrange [10 : 14]
	  set yrange [0:2.5]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"../DDME1/npeuY/output_tov.dat\" u 2:3 t \"DDME1-npeuY\" w l lw $lw  lc \"${color[0]}\", \\
           \"../DDME1/X_sig_d=1.05/npeuYD/output_tov.dat\" u 2:3 t \"DDME1-npeuYD\" w l lw $lw dt '...' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-120/output_tov.dat\" u 2:3 t \"DDME1-120\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-140/output_tov.dat\" u 2:3 t \"DDME1-140\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.1/npeuYD/output_tov.dat\" u 2:3 t \"DDME1-npeuYD\" w l lw $lw dt '...' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-120/output_tov.dat\" u 2:3 t \"DDME1-120\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-140/output_tov.dat\" u 2:3 t \"DDME1-140\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
		   \"../DDME1/X_sig_d=1.15/npeuYD/output_tov.dat\" u 2:3 t \"DDME1-npeuYD\" w l lw $lw dt '...' lc \"${color[6]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-120/output_tov.dat\" u 2:3 t \"DDME1--120\" w l lw $lw dt '--' lc \"${color[6]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-140/output_tov.dat\" u 2:3 t \"DDME1--140\" w l lw $lw dt '.-' lc \"${color[6]}\"" | gnuplot
