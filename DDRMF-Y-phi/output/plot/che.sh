#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "dark-yellow" "dark-blue")
echo "set terminal png truecolor
      set output \"./che.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"chemical\"
	  set xtics 0.2
	  set mxtics 1
      set xrange [ 0 : 1.2 ]
	  set yrange [800:1600]
	  set xlabel \"{/Symbol r_b}\"
      set ylabel \"{chemical potential}\"
      plot \"./che.txt\" u 1:2 t \"n\" w l lw $lw lc \"${color[0]}\", \\
           \"./che.txt\" u 1:3 t \"p\" w l lw $lw lc \"${color[1]}\", \\
           \"./che.txt\" u 1:4 t \"e\" w l lw $lw lc \"${color[2]}\", \\
           \"./che.txt\" u 1:5 t \"u\" w l lw $lw lc \"${color[3]}\", \\
           \"./che.txt\" u 1:6 t \"{/Symbol L}\" w l lw $lw lc \"${color[4]}\", \\
           \"./che.txt\" u 1:7 t \"{/Symbol S^+}\" w l lw $lw lc \"${color[5]}\",\\
           \"./che.txt\" u 1:8 t \"{/Symbol S^0}\" w l lw $lw lc \"${color[6]}\",\\
           \"./che.txt\" u 1:9 t \"{/Symbol S^-}\" w l lw $lw lc \"${color[7]}\",\\
           \"./che.txt\" u 1:10 t \"{/Symbol X^0}\" w l lw $lw lc \"${color[8]}\",\\
           \"./che.txt\" u 1:11 t \"{/Symbol X^-}\" w l lw $lw lc \"${color[9]}\"" | gnuplot
