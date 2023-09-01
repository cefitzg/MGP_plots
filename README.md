# Code and Data for "Temporal Dynamics of Faculty Hiring in Mathematics" by Cody FitzGerald, Yitong Huang, Katelyn Plaiser Leisman, and Chad M. Topaz. 

## Data

MGPdata.csv is an anonmyized data set that was scrapped from the Mathematics Genealogy Project (https://www.genealogy.math.ndsu.nodak.edu/). logisticdata.Rdata is a version of MGPdata.csv used for the logistic regression analysis. 

## Code 
log_reg.R runs the statistical analysis on logisticdata.Rdata (Section 3.3) 

school_list_compute.m computes which university math departments we consider "elite" and which departments we consider "well-placing."

MGP_plots.m is a MATLAB script (written and tested in 2021b) that does the analysis and 
plots the figures 1-10 found in the manuscript. MGP_plots.m reads MGPdata.csv and uses the helper script school_list_compute.m.





