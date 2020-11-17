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
  vlk <- select(filter(df, chain == "IGL" | chain == "IGK"), c(chain,
                                                               clonotype_ids,
                                                               contig_ids,
                                                               v_genes,
                                                               proportion))
  vlk.single.call <- vlk %>% separate(v_genes, c("single_call"),
                                    sep = "[;]", extra = "drop", 
                                    fill = "left")
  
  #this adds column for VL family
  vlk.fam <- vlk.single.call %>% 
    separate(single_call, c("vlk_family"), 
             sep = "[-]", extra = "drop", 
             fill = "left")
  vlk_family <- as.factor(vlk.fam$vlk_family)
  vlk.family.df <- data.frame(vlk.single.call, vlk_family)
 
  #this generates filename for csv
  filename.final <- basename(file)
  vlk.single.call.freq.filename <- paste(dir.out,filename.final,
                                         "_vlk_single_call_freq",".csv",sep="")
  #could add this back in to see file that still has multiple vl or vk calls listed
  #vlk_freq_filename <- paste(dir_out,filename_final,"_vlk_freq",".csv",sep="")
  #write.csv(vlk, file = vlk_freq_filename)
  library(tibble)
  #this adds sample ID column
  filename <- c(rep_len(filename.final, 
                        length(as.character(vlk.family.df$vlk_family))))
  vlk.fam.id <- data.frame(filename, vlk.family.df)
  write.csv(vlk.fam.id, file = vlk.single.call.freq.filename)
}

#ignore error that pops up after running lapply, this happens when that subset
#wasn't identified in a given sample
lapply(file.list,
       FUN = make.single.vlk.call.csv,
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
dir_out <- "~/desktop/10X_BCR_pipelines/output/vlk/"
directories <- c(dir_in, dir_out)
lapply(directories,
       FUN = check_create_dir)

library(dplyr)
library(tidyr)

file <- read.csv("~/desktop/10X_BCR_pipelines/input/VDJ_gene-usage_RB4025_1_cd21low_1.csv")
df <- data.frame(file, header = TRUE)
file <- "~/desktop/10X_BCR_pipelines/input/VDJ_gene-usage_RB4025_1_cd21low_1.csv"
#vlk retains multiple calls
#vlk_single_call only keeps first VL call
vlk <- select(filter(df, chain == "IGL" | chain == "IGK"), c(chain,clonotype_ids,v_genes,proportion))
vlk_single_call <- vlk %>% separate(v_genes, c("single_call"),
                                  sep = "[;]", extra = "drop", 
                                  fill = "left")

#this adds column for VH family
vlk_fam <- vlk_single_call %>% 
  separate(single_call, c("vlk_family"), 
            sep = "[-]", extra = "drop", 
            fill = "left")
vlk_family <- as.factor(vlk_fam$vlk_family)
vlk_family_df <- data.frame(vlk_single_call, vlk_family)
  
#this generates filename for csv
#position call different for this line than for using lapply function above
#ignore that this part doesn't work
filename1 <- as.character((data.frame(strsplit(file, "/")) [7,1]))
filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
filename_final <- as.character((data.frame(strsplit(filename2, "usage_")) [2,1]))
vlk_single_call_freq_filename <- paste(dir_out,filename_final,"_vlk_single_call_freq",".csv",sep="")
#could add this back in to see file that still has multiple vh calls listed
  #vlk_freq_filename <- paste(dir_out,filename_final,"_vlk_freq",".csv",sep="")
  #write.csv(vl, file = vlk_freq_filename)
library(tibble)
#this adds sample ID column
filename <- c(rep_len(filename_final, length(as.character(vlk_family_df$vlk_family))))
vlk_fam_id <- data.frame(filename, vlk_family_df)
write.csv(vlk_fam_id, file = vlk_single_call_freq_filename)
