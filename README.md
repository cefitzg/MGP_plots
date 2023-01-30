Code for analysis and recreating figures in "Temporal Dynamics of Faculty Hiring in Mathematics" by Cody FitzGerald, Yitong Huang, Katelyn Plaiser-Leisman, and Chad M. Topaz. 

log_reg.R runs the statistical analysis on logisticdata.Rdata (Section 3.3)

MGP_plots.m is a MATLAB script (written and tested in 2021b) that does the analysis and 
plots the figures 1-9 found in the manuscript. MGP_plots.m reads MGPdata.csv and uses the helper script school_list_compute.m. 

fig10.m plots Figure 10 and uses target6.csv, a non-anonymized version of MGPdata.csv.

Things to check: 

Did we plot Fig. 4 right based on Pepper's description? Is this clear? 


