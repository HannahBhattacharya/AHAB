#This script takes single chains export file from 
#10X Genomics Loupe VDJ Browser as input and pulls 
#only the VL data from it. VL includes both kappa and lambda.
#CSV output has filename, chain, clonotype_id, contig_id,
#single VL call (single chain), proportion, and vl_family

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
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*.csv",
                        full.names = TRUE)
file.list

make.single.vlk.call.csv <- function(file, the.dir) {
  df <- data.frame(read.csv(file, header = TRUE))
  #keep these calls here, it fails on select without them
  library(dplyr)
  library(tidyr)
  #vlk retains multiple calls
  #vlk_single_call only keeps first VL or VK call
  vlk <- select(filter(df, chain == "IGL" | chain == "IGK"), 
                c(chain,
                  clonotype_ids,
                  contig_ids,
                  v_genes,
                  proportion))
                                                               
  vlk.single.call <- vlk %>% separate(v_genes, c("single_call"),
                                    sep = "[;]", extra = "drop", 
                                    fill = "left")
  
  #this adds column for VL family
  vlk.fam <- vlk.single.call %>% 
    separate(single_call, c("vlk.family"), 
             sep = "[-]", extra = "drop", 
             fill = "left")
  vlk.family <- as.factor(vlk.fam$vlk.family)
  vlk.family.df <- data.frame(vlk.single.call, vlk.family)
 
  #this generates filename for csv
  filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                     "usage_")) [2,1])
  vlk.single.call.freq.filename <- paste(dir.out,
                                         "vlk_single_call_freq_",
                                         filename.final, sep="")
  #could add this back in to see file that still has multiple vl or vk calls listed
  #vlk_freq_filename <- paste(dir_out,filename_final,"_vlk_freq",".csv",sep="")
  #write.csv(vlk, file = vlk_freq_filename)
  library(tibble)
  #this adds sample ID column
  file.id <- as.character(data.frame(strsplit(filename.final, 
                                              ".csv")) [1,1])
  filename <- c(rep_len(file.id, 
                        length(as.character(vlk.family.df$vlk.family))))
  vlk.fam.id <- data.frame(filename, vlk.family.df)
  write.csv(vlk.fam.id, file = vlk.single.call.freq.filename)
}

#ignore error that pops up after running lapply, this happens when that subset
#wasn't identified in a given sample
lapply(file.list,
       FUN = make.single.vlk.call.csv,
       the.dir = dir.out)

########## STOP HERE ##################
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
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file <- "~/desktop/10X_BCR_pipelines/input/VDJ_gene-usage_4025-RB-1_cd21low_1.csv"
file.list <- read.csv(file)

df <- data.frame(file.list, header = TRUE)
#vlk retains multiple calls
#vlk_single_call only keeps first VL call
vlk <- select(filter(df, chain == "IGL" | chain == "IGK"), 
              c(chain,clonotype_ids,v_genes,proportion))
vlk.single.call <- vlk %>% separate(v_genes, c("single_call"),
                                  sep = "[;]", extra = "drop", 
                                  fill = "left")

#this adds column for VH family
vlk.fam <- vlk.single.call %>% 
  separate(single_call, c("vlk.family"), 
            sep = "[-]", extra = "drop", 
            fill = "left")
vlk.family <- as.factor(vlk.fam$vlk.family)
vlk.family.df <- data.frame(vlk.single.call, vlk.family)
  
#this generates filename for csv
#position call different for this line than for using lapply function above
filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                   "usage_")) [2,1])
#filename1 <- as.character((data.frame(strsplit(file, "/")) [7,1]))
#filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
#filename_final <- as.character((data.frame(strsplit(filename2, "usage_")) [2,1]))
vlk.single.call.freq.filename <- paste(dir.out,
                                       "vlk.single.call.freq_",
                                       filename.final, sep="")
#could add this back in to see file that still has multiple vh calls listed
#vlk_freq_filename <- paste(dir_out,
#                           "vlk_freq",
#                           filename_final, sep="")
  #write.csv(vlk, file = vlk.freq.filename)
library(tibble)
#this adds sample ID column
file.id <- as.character(data.frame(strsplit(filename.final, 
                                            ".csv")) [1,1])
filename <- c(rep_len(file.id, 
                      length(as.character(vlk.family.df$vlk.family))))
vlk.fam.id <- data.frame(filename, vlk.family.df)
write.csv(vlk.fam.id, file = vlk.single.call.freq.filename)
