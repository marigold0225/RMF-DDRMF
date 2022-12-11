#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./1112-DDME1-M_R.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{M-R}\"
      set xrange [ 10 : 14.0 ]
	  set yrange [0:2.5]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"../1112-DDME1/npeuY/M_R.txt\" u 4:3 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"../1112-DDME1/X_sig_d=1.05/npeuYD/M_R.txt\" u 4:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../1112-DDME1/X_sig_d=1.05/Uk=-120/M_R.txt\" u 4:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../1112-DDME1/X_sig_d=1.05/Uk=-140/M_R.txt\" u 4:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[1]}\", \\
           \"../1112-DDME1/X_sig_d=1.10/npeuYD/M_R.txt\" u 4:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../1112-DDME1/X_sig_d=1.10/Uk=-120/M_R.txt\" u 4:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../1112-DDME1/X_sig_d=1.10/Uk=-140/M_R.txt\" u 4:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[2]}\", \\
           \"../1112-DDME1/X_sig_d=1.15/npeuYD/M_R.txt\" u 4:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[6]}\", \\
           \"../1112-DDME1/X_sig_d=1.15/Uk=-120/M_R.txt\" u 4:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../1112-DDME1/X_sig_d=1.15/Uk=-140/M_R.txt\" u 4:3 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
