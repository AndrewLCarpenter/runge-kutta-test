# Gnuplot script file for plotting data from "fort.80" file
# run at command line by > load "Boldrighini_Franceschini.p"
set key outside;
set xrange[+85:100] ;
set yrange[-6:6] ;
plot  "fort.80a" using 1:2 with lines, \
      "" using 1:3 with lines, \
      "" using 1:4 with lines, \
      "" using 1:5 with lines, \
      "" using 1:6 with lines, \
      "fort.80b" using 1:2 with lines, \
      "" using 1:3 with lines, \
      "" using 1:4 with lines, \
      "" using 1:5 with lines, \
      "" using 1:6 with lines, \
      "fort.80c" using 1:2 with lines, \
      "" using 1:3 with lines, \
      "" using 1:4 with lines, \
      "" using 1:5 with lines, \
      "" using 1:6 with lines;
