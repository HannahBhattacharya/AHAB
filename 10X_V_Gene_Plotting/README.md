WHAT IT DOES:
This series of R scripts are designed to manipulate BCRseq data from the 
10X Genomics Chromium platform. It starts with single chains export file 
from 10X Genomics Loupe VDJ Browser as input and pulls VH or VL data from
it (see filenames). VH family and gene plots are rendered using ggplot2
(no ggplot knowledge required if you can live with how they're rendered).
This way you can "gate" on subsets of interest and only look at corresponding
BCRs.

REQUISITE NEWBIE APOLOGY:
I'm a total novice in R, feel free to make suggestions for streamlined code.
I frequently sacrificed "pretty" for "functional" because ain't nobody got
time for perfect.
