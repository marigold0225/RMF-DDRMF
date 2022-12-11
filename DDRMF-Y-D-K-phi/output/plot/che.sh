#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./che.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{chemical potential}\"
      set xrange [ 0 : 7 ]
	  set yrange [0:500]
      set xlabel \"{/Symbol r_B}\"
      set ylabel \"chemical potential\"
      plot \"../DDME1/X_sig_d=1.05/Uk=-120/parameter.txt\" u 1:5 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
      	   \"../DDME1/X_sig_d=1.05/Uk=-120/parameter.txt\" u 1:6 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-140/parameter.txt\" u 1:5 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.05/Uk=-140/parameter.txt\" u 1:6 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-120/parameter.txt\" u 1:5 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-120/parameter.txt\" u 1:6 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-140/parameter.txt\" u 1:5 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.1/Uk=-140/parameter.txt\" u 1:6 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[2]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-120/parameter.txt\" u 1:5 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-120/parameter.txt\" u 1:6 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-140/parameter.txt\" u 1:5 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[6]}\", \\
           \"../DDME1/X_sig_d=1.15/Uk=-140/parameter.txt\" u 1:6 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
