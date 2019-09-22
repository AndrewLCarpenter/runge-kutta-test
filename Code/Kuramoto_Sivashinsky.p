# Gnuplot script file for plotting data from "fort.80" file
# run at command line by > load "Kuramoto_Sivashinsky.p"
set autoscale;
set key outside;
plot 'fort.80' using 1:02 with lines, \
     'fort.80' using 1:17 with lines, \
     'fort.80' using 1:33 with lines, \
     'fort.80' using 1:49 with lines, \
     'fort.80' using 1:65 with lines#, \
    #'fort.80' using 1:81 with lines, \
    #'fort.80' using 1:97 with lines, \
    #'fort.80' using 1:113 with lines, \
    #'fort.80' using 1:129 with lines ;
