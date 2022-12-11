#!/bin/bash
lw=2
pt=(5 7 9 11 13 15)
ps=1.0
color=("black" "blue" "red" "orange" "magenta" "cyan" "green" "dark-red" "dark-yellow" "dark-blue")
echo "set terminal png truecolor
      set output \"./population.png\"
      set term pngcairo size 1024,768
      set grid
	  set title \"density\"
	  set xtics 0.2
	  set mxtics 1
      set xrange [ 0 : 1.2 ]
	  set logscale y
	  set yrange [0.001:1]
	  set xlabel \"{/Symbol r_b}\"
      set ylabel \"{population}\"
      plot \"./density.txt\" u 1:2 t \"n\" w l lw $lw lc \"${color[0]}\", \\
           \"./density.txt\" u 1:3 t \"p\" w l lw $lw lc \"${color[1]}\", \\
           \"./density.txt\" u 1:4 t \"e\" w l lw $lw lc \"${color[2]}\", \\
           \"./density.txt\" u 1:5 t \"u\" w l lw $lw lc \"${color[3]}\", \\
           \"./density.txt\" u 1:6 t \"{/Symbol L}\" w l lw $lw lc \"${color[4]}\", \\
           \"./density.txt\" u 1:7 t \"{/Symbol S^+}\" w l lw $lw lc \"${color[5]}\",\\
           \"./density.txt\" u 1:8 t \"{/Symbol S^0}\" w l lw $lw lc \"${color[6]}\",\\
           \"./density.txt\" u 1:9 t \"{/Symbol S^-}\" w l lw $lw lc \"${color[7]}\",\\
           \"./density.txt\" u 1:10 t \"{/Symbol X^0}\" w l lw $lw lc \"${color[8]}\",\\
           \"./density.txt\" u 1:11 t \"{/Symbol X^-}\" w l lw $lw lc \"${color[9]}\"" | gnuplot
