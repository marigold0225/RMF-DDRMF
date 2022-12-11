#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./DDME1-parameter.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{effect_mass}\"
      set xrange [ 0 : 7.5 ]
	  set yrange [0:0.5]
      set xlabel \"R\"
      set ylabel \"M\"
      plot \"../DDME1/npeuY/parameter.txt\" u 1:2 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"../DDME1/X_sig_d=1.05/npeuYD/parameter.txt\" u 1:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-120/parameter.txt\" u 1:8 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-140/parameter.txt\" u 1:8 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.1/npeuYD/parameter.txt\" u 1:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-120/parameter.txt\" u 1:8 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-140/parameter.txt\" u 1:8 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.15/npeuYD/parameter.txt\" u 1:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[5]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-120/parameter.txt\" u 1:8 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[5]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-140/parameter.txt\" u 1:8 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[5]}\"" | gnuplot
