#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./cut-IUFSU-output_tov.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{M-R}\"
      set xrange [ 8 : 15.0 ]
	  set yrange [0:2.5]
      set xlabel \"R\"
      set ylabel \"M\"
      plot  \"../IUFSU/npeuY/output_tov.dat\" u 2:3 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"../IUFSU/X_sig_d=1.05/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-160/output_tov.dat\" u 2:3 t \"Uk=-160\" w l lw $lw dt 2 lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.10/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-160/output_tov.dat\" u 2:3 t \"Uk=-160\" w l lw $lw dt 2 lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.15/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-160/output_tov.dat\" u 2:3 t \"Uk=-160\" w l lw $lw dt 2 lc \"${color[6]}\", \\
	  	   \"../cut-IUFSU/npeuY/output_tov.dat\" u 2:3 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"../cut-IUFSU/X_sig_d=1.05/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../cut-IUFSU/X_sig_d=1.05/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../cut-IUFSU/X_sig_d=1.05/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[1]}\", \\
           \"../cut-IUFSU/X_sig_d=1.10/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../cut-IUFSU/X_sig_d=1.10/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../cut-IUFSU/X_sig_d=1.10/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[2]}\", \\
           \"../cut-IUFSU/X_sig_d=1.15/npeuYD/output_tov.dat\" u 2:3 t \"npeuYD\" w l lw $lw dt '--' lc \"${color[6]}\", \\
           \"../cut-IUFSU/X_sig_d=1.15/Uk=-120/output_tov.dat\" u 2:3 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../cut-IUFSU/X_sig_d=1.15/Uk=-140/output_tov.dat\" u 2:3 t \"Uk=-140\" w l lw $lw dt 2 lc \"${color[6]}\"" | gnuplot
