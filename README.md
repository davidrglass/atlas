Code related to An integrated multi-omic single cell atlas of human B cell identity: link pending 

These scripts allow processing of fcs files and generation of the main figures from the manuscript. Instructions for each script are written in script headers.

screen_preprocess.R processes fcs files from the surface screen. Files can be found here: https://flowrepository.org/id/FR-FCM-Z2MA
Requires:
fcs files
screen_metadata.scv

other_preprocess.R process fcs files from the other datasets in the manuscript. Files can be found here: https://flowrepository.org/id/FR-FCM-Z2MC
Requires: 
fcs files
figure_4_metadata.csv
figure_5_biosynthesis_metadata.csv
figure_5_metabolism_metadata.csv
figure_5_signaling_metadata.csv
figure_6_metadata.csv

screen_figures.R generates subfigures from figures 1-3.
Requires:
go_parsed.csv

figure_4.R figure_5.R figure_6.R generate subfigures from figures 4-6.
