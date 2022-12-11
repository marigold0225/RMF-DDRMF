#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "light-blue" "dark-blue" "dark-green" "web-green" "dark-cyan" "dark-yellow ")
echo "set terminal png truecolor
      set output \"./IUFSU_EOS.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{DDRMFEOS}\"
      set xrange [ 0 : 1200 ]
      set xlabel \"Energy\"
      set ylabel \"Pressure\"
     plot \"../IUFSU/npeuY/EOS.dat\" u 3:2 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"../IUFSU/X_sig_d=1.05/npeuYD/EOS.dat\" u 3:2 t \"npeuYD\" w l lw $lw lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-120/EOS.dat\" u 3:2 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-140/EOS.dat\" u 3:2 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.05/Uk=-160/EOS.dat\" u 3:2 t \"Uk=-140\" w l lw $lw dt '--' lc \"${color[1]}\", \\
           \"../IUFSU/X_sig_d=1.10/npeuYD/EOS.dat\" u 3:2 t \"npeuYD\" w l lw $lw lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-120/EOS.dat\" u 3:2 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-140/EOS.dat\" u 3:2 t \"Uk=-140\" w l lw $lw dt '..' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.10/Uk=-160/EOS.dat\" u 3:2 t \"Uk=-140\" w l lw $lw dt '--' lc \"${color[2]}\", \\
           \"../IUFSU/X_sig_d=1.15/npeuYD/EOS.dat\" u 3:2 t \"npeuYD\" w l lw $lw lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-120/EOS.dat\" u 3:2 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-140/EOS.dat\" u 3:2 t \"Uk=-120\" w l lw $lw dt '.-' lc \"${color[6]}\", \\
           \"../IUFSU/X_sig_d=1.15/Uk=-160/EOS.dat\" u 3:2 t \"Uk=-140\" w l lw $lw dt '--' lc \"${color[6]}\"" | gnuplot
