Files ending in .Rmd are R markdown files that run on free R software (see e.g., https://rmarkdown.rstudio.com/ ).

Files that start with "COP" create analyses that match with the specified section of the Statistical Analysis Plan.

"paper figures.Rmd" creates the figures in the paper.

reformat_paragraph.tex is LaTex code called by CoP_SAP_6.4_Appendix_Final.Rmd

All data are in one file in the data directory. Results were created using R version 4 (4.0.0-4.6.1).

R should be run from the same directory as the files and the data file should be in the /data subdirectory of that directory so that the data file can be found. The folder intermediate_outcomes is intentionally empty, but will be filled with intermediate results if some .Rmd files are run.

