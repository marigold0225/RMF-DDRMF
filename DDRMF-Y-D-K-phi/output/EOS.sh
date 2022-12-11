#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "dark-green")
echo "set terminal png truecolor
      set output \"./DDMEX-EOS.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"{DDMEX-EOS}\"
      set xrange [ 0 : 1500 ]
      set xlabel \"Energy\"
      set ylabel \"Pressure\"
      plot \"./TD/DDMEX/npeuY/EOS.dat\" u 3:2 t \"npeuY\" w l lw $lw lc \"${color[0]}\", \\
           \"./TD/DDMEX/Ud=-50/npeuYD/EOS.dat\" u 3:2 t \"-50-npeuYD\" w l lw $lw lc \"${color[1]}\", \\
           \"./TD/DDMEX/Ud=-50/Uk=-120/EOS.dat\" u 3:2 t \"-50-120\" w l lw $lw lc \"${color[2]}\", \\
           \"./TD/DDMEX/Ud=-50/Uk=-140/EOS.dat\" u 3:2 t \"-50-140\" w l lw $lw lc \"${color[3]}\", \\
           \"./TD/DDMEX/Ud=-50/Uk=-160/EOS.dat\" u 3:2 t \"-50-160\" w l lw $lw lc \"${color[4]}\", \\
           \"./TD/DDMEX/Ud=-100/npeuYD/EOS.dat\" u 3:2 t \"-100-npeuYD\" w l lw $lw dt 2 lc \"${color[5]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-120/EOS.dat\" u 3:2 t \"-100-120\" w l lw $lw dt 2 lc \"${color[6]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-140/EOS.dat\" u 3:2 t \"-100-140\" w l lw $lw dt 2 lc \"${color[7]}\",\\
           \"./TD/DDMEX/Ud=-100/Uk=-160/EOS.dat\" u 3:2 t \"-100-160\" w l lw $lw dt 2 lc \"${color[8]}\"" | gnuplot
