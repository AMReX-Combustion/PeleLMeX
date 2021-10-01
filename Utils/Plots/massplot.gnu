set terminal qt size 700,500
set xlabel font "Helvetica,12"
set tics font ",12"
set tmargin 2
set multiplot layout 2,2 rowsfirst
set grid

set label 1 'mass [kg]' at graph 0.45,1.05 font ',12'
set xlabel 'time [ms]'
plot 'temporals/tempMass' u ($2*1000):3 w lp notitle

set label 1 'dmdt [kg/s]' at graph 0.45,1.05 font ',12'
set xlabel 'time [ms]'
plot 'temporals/tempMass' u ($2*1000):4 w lp notitle
replot 'temporals/tempMass' u ($2*1000):5 w lp notitle

set label 1 'balance [kg/s]' at graph 0.45,1.05 font ',12'
set xlabel 'time [ms]'
plot 'temporals/tempMass' u ($2*1000):6 w lp notitle

pause mouse close
unset multiplot
