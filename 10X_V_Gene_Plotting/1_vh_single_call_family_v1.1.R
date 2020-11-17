#This script takes single chains export file from 
#10X Genomics Loupe VDJ Browser as input and pulls 
#only the VH data from it.
#CSV output has filename, chain, clonotype_id, contig_id,
#single VH call (single chain), proportion, and vh_family

#Things you may need to revise are marked:
######ALL CAPS#########

library(tidyverse)
library(magrittr)
library(readr)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

######DEFINE DIR IN########
#this is where files come from
dir.in <- "~/desktop/10X_BCR_pipelines/input"
######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*.csv",
                        full.names = TRUE)
file.list

make.single.vh.call.csv <- function(file, the.dir) {
  df <- data.frame(read.csv(file, header = TRUE))
  #keep these calls here, it fails on select without them
  library(dplyr)
  library(tidyr)
  #vh retains multiple calls
  #vh_single_call only keeps first VH call
  vh <- select(filter(df, chain == "IGH"), 
               c(chain,
                 clonotype_ids,
                 contig_ids,
                 v_genes,
                 proportion))
                                             
  vh.single.call <- vh %>% separate(v_genes, c("single_call"),
                                    sep = "[;]", extra = "drop", 
                                    fill = "left")
  
  #this adds column for VH family
  vh.fam <- vh.single.call %>% 
    separate(single_call, c("vh.family"), 
             sep = "[-]", extra = "drop", 
             fill = "left")
  vh.family <- as.factor(vh.fam$vh.family)
  vh.family.df <- data.frame(vh.single.call, vh.family)
 
  #this generates filename for csv
  filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                     "usage_")) [2,1])
  
  #filename_final <- as.character((data.frame(strsplit(filename2, "usage_")) [2,1]
  vh.single.call.freq.filename <- paste(dir.out,
                                        "vh_single_call_freq_",
                                        filename.final, sep="")
  #could add this back in to see file that still has multiple vh calls listed
  #vh_freq_filename <- paste(dir_out,filename_final,"_vh_freq",".csv",sep="")
  #write.csv(vh, file = vh_freq_filename)
  library(tibble)
  #this adds sample ID column
  file.id <- as.character(data.frame(strsplit(filename.final, 
                                              ".csv")) [1,1])
  filename <- c(rep_len(file.id, 
                        length(as.character(vh.family.df$vh.family))))
  vh.fam.id <- data.frame(filename, vh.family.df)
  write.csv(vh.fam.id, file = vh.single.call.freq.filename)
}

#ignore the error that pops up after running lapply, this happens
#when the subset isn't present in a given sample
lapply(file.list,
       FUN = make.single.vh.call.csv,
       the.dir = dir.out)

############### STOP HERE #################
#--------------------------------------------
#test code to see how this is working on single file
library(tidyverse)
library(magrittr)
library(readr)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

#this is where files come from
dir.in <- "~/desktop/10X_BCR_pipelines/input"
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file <- "~/desktop/10X_BCR_pipelines/input/VDJ_gene-usage_4025-RB-1_cd21low_1.csv"
file.list <- read.csv(file)

df <- data.frame(file.list, header = TRUE)
#vh retains multiple calls
#vh_single_call only keeps first VH call
vh <- select(filter(df, chain == "IGH"), c(chain,clonotype_ids,v_genes,proportion))
vh.single.call <- vh %>% separate(v_genes, c("single_call"),
                                  sep = "[;]", extra = "drop", 
                                  fill = "left")

#this adds column for VH family
vh.fam <- vh.single.call %>% 
  separate(single_call, c("vh.family"), 
            sep = "[-]", extra = "drop", 
            fill = "left")
vh.family <- as.factor(vh.fam$vh.family)
vh.family.df <- data.frame(vh.single.call, vh.family)

#this generates filename for csv
filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                   "usage_")) [2,1])
#filename1 <- as.character((data.frame(strsplit(file, "/")) [7,1]))
#filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
#filename_final <- as.character((data.frame(strsplit(filename2, "usage_")) [2,1]))
vh.single.call.freq.filename <- paste(dir.out,
                                      "vh_single_call_freq_",
                                      filename.final, sep="")
#could add this back in to see file that still has multiple vh calls listed
#vh.freq.filename <- paste(dir.out,
 #                         "vh_freq_",
#                          filename.final, sep="")
  #write.csv(vh, file = vh.freq.filename)
library(tibble)
#this adds sample ID column
file.id <- as.character(data.frame(strsplit(filename.final, 
                                            ".csv")) [1,1])
filename <- c(rep_len(file.id, 
                      length(as.character(vh.family.df$vh.family))))
vh.fam.id <- data.frame(filename, vh.family.df)
write.csv(vh.fam.id, file = vh.single.call.freq.filename)
