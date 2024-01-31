# for plotting rmsd in gnuplot
set yr [0:2.6]
set xr [0:660]
plot "run0/rmsd.xvg" u ($0/10):2 w l t "rep0","run1/rmsd.xvg" u ($0/10):2 w l t "rep1","run2/rmsd.xvg" u ($0/10):2 w l t "rep2","run3/rmsd.xvg" u ($0/10):2 w l t "rep3"
set xlabel "time (ns)"
set ylabel "rmsd (nm)"
repl
