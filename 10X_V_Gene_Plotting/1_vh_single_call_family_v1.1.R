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
  vh <- select(filter(df, chain == "IGH"), c(chain,
                                             clonotype_ids,
                                             contig_ids,
                                             v_genes,
                                             proportion))
  vh.single.call <- vh %>% separate(v_genes, c("single_call"),
                                    sep = "[;]", extra = "drop", 
                                    fill = "left")
  
  #this adds column for VH family
  vh.fam <- vh.single.call %>% 
    separate(single_call, c("vh_family"), 
             sep = "[-]", extra = "drop", 
             fill = "left")
  vh_family <- as.factor(vh.fam$vh_family)
  vh.family.df <- data.frame(vh.single.call, vh_family)
 
  #this generates filename for csv
  filename.final <- basename(file)
  vh.single.call.freq.filename <- paste(dir.out,filename.final,
                                        "_vh_single_call_freq",".csv",sep="")
  #could add this back in to see file that still has multiple vh calls listed
  #vh_freq_filename <- paste(dir_out,filename_final,"_vh_freq",".csv",sep="")
  #write.csv(vh, file = vh_freq_filename)
  library(tibble)
  #this adds sample ID column
  filename <- c(rep_len(filename.final, 
                        length(as.character(vh.family.df$vh_family))))
  vh.fam.id <- data.frame(filename, vh.family.df)
  write.csv(vh.fam.id, file = vh.single.call.freq.filename)
}

#ignore the error that pops up after running lapply, this happens
#when the subset isn't present in a given sample
lapply(file.list,
       FUN = make.single.vh.call.csv,
       the.dir = dir.out)

#--------------------------------------------
#test code to see how this is working on single file
library(tidyverse)
library(magrittr)
library(readr)

check_create_dir <- function(the_dir) {
  if (!dir.exists(the_dir)) {
    dir.create(the_dir, recursive = TRUE) }
}

#this is where files come from
dir_in <- "~/desktop/10X_BCR_pipelines/input"
#this is where output files go
dir_out <- "~/desktop/10X_BCR_pipelines/output/vh/"
directories <- c(dir_in, dir_out)
lapply(directories,
       FUN = check_create_dir)

library(dplyr)
library(tidyr)

file <- read.csv("~/desktop/10X_BCR_pipelines/input/VDJ_gene-usage_4025-RB-1_cd21low_1.csv")

df <- data.frame(file, header = TRUE)
#vh retains multiple calls
#vh_single_call only keeps first VH call
vh <- select(filter(df, chain == "IGH"), c(chain,clonotype_ids,v_genes,proportion))
vh_single_call <- vh %>% separate(v_genes, c("single_call"),
                                  sep = "[;]", extra = "drop", 
                                  fill = "left")

#this adds column for VH family
vh_fam <- vh_single_call %>% 
  separate(single_call, c("vh_family"), 
            sep = "[-]", extra = "drop", 
            fill = "left")
vh_family <- as.factor(vh_fam$vh_family)
vh_family_df <- data.frame(vh_single_call, vh_family)

?basename 
#this generates filename for csv
#ignore that the naming code below doesn't work for an individual file
filename <- basename(as.character(file))
filename1 <- as.character((data.frame(strsplit(file, "/")) [7,1]))
filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
filename_final <- as.character((data.frame(strsplit(filename2, "usage_")) [2,1]))
vh_single_call_freq_filename <- paste(dir_out,filename_final,"_vh_single_call_freq",".csv",sep="")
#could add this back in to see file that still has multiple vh calls listed
  #vh_freq_filename <- paste(dir_out,filename_final,"_vh_freq",".csv",sep="")
  #write.csv(vh, file = vh_freq_filename)
library(tibble)
#this adds sample ID column
filename_final <- "bullshit"
filename <- c(rep_len(filename_final, length(as.character(vh_family_df$vh_family))))
vh_fam_id <- data.frame(filename, vh_family_df)
write.csv(vh_fam_id, file = vh_single_call_freq_filename)
