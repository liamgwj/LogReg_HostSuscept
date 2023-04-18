Liam GW Johnson 2022-02-02
This directory contains the R scripts needed to reproduce the analysis from the manuscript "Logistic regression is ineffective for predicting host susceptibility based on phylogenetic distance".

'01_plr-simulation.r' runs a 36000-iteration simulation that takes about 20 minutes on my machine (4 cores, 8 threads, 2.7GHz; 30GB RAM).

'00_plr-functions.r' contains a number of custom functions needed to run the simulation; it's called from source by the simulation script.

'02_plr-stats-plots.r' reads in the simulation output and generates the figure included in the manuscript. It should be run from a fresh R session after completing the simulation.

The scripts are written assuming that they all sit together in the same directory, from which R is called and to which all output is written.
