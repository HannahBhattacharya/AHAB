#This code takes output files from 1_vlk_single_call_family Rscript
#this code collapses data for each sample to show counts and
#proportions of VH family calls and VH gene ID calls

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
dir.in <- "~/desktop/10X_BCR_pipelines/output/vh"

######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/freq_summary/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*.csv",
                        full.names = TRUE)
file.list
vh.fam.list <- c("IGHV1", "IGHV2", "IGHV3", "IGHV3/OR16", 
                 "IGHV4", "IGHV5", "IGHV6", "IGHV7")
#need to check the list below for completeness!
vh.gene.list <- c("IGHV1-18", "IGHV1-2", "IGHV1-24", "IGHV1-3", 
                "IGHV1-45", "IGHV1-46", "IGHV1-58", "IGHV1-69", 
                "IGHV1-69-2", "IGHV1-69D", "IGHV1-8", "IGHV2-26",
                "IGHV2-5", "IGHV2-70", "IGHV2-70D", "IGHV3-11", 
                "IGHV3-13", "IGHV3-15", "IGHV3-20", "IGHV3-21", 
                "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-43",
                "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64", 
                "IGHV3-64D", "IGHV3-66", "IGHV3-7", "IGHV3-72", 
                "IGHV3-73", "IGHV3-74", "IGHV3/OR16-12", 
                "IGHV4-28", "IGHV4-30-2", "IGHV4-31", "IGHV4-34", 
                "IGHV4-39", "IGHV4-4", "IGHV4-59",  "IGHV4-61", 
                "IGHV5-10-1", "IGHV5-51", "IGHV6-1", "IGHV7-4-1", 
                "IGHV7-81")

#This function calculates proportion of each VH gene
#and each VH gene family and writes these out as CSV files
sum.vh.freq.make.csv <- function(file, the.dir) {
  df <- data.frame(read.csv(file, header = TRUE))
  df2 <- data.frame(read.csv(file, header = TRUE))
  #keep these calls here, it fails on select without them
  library(dplyr)
  library(tidyr)
  #may not need to call tibble
  library(tibble)
  
  #this makes summary file with proportions by VH family
  sum.file.vh.fam.df <- df %>%
    group_by(vh.family) %>%
    summarise(
      observations = n(),
      vh.fam.freq.sum = sum(proportion)
    )
  sum.file.vh.fam.df

  #this makes summary file with proportions by VH gene id
  sum.file.vh.gene.df <- df %>%
    group_by(single_call) %>%
    summarise(
      observations = n(),
      vh.gene.freq.sum = sum(proportion)
    )
  sum.file.vh.gene.df
  
  #this generates filename for csv files
  filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                     "freq_")) [2,1])
  vh.fam.freq.sum.filename <- paste(dir.out,
                                    "vh_fam_freq_sum_",
                                    filename.final, sep="")
  vh.gene.freq.sum.filename <- paste(dir.out,
                                   "vh_gene_freq_sum_",
                                   filename.final, sep="")

  #this adds sample ID column to both summary files
  file.id <- as.character(data.frame(strsplit(filename.final, 
                                              ".csv")) [1,1])
  
  file.gene <- c(rep_len(file.id, 
                         length(as.character(sum.file.vh.gene.df$single_call))))
  
  file.fam <- c(rep_len(file.id,
                        length(as.character(sum.file.vh.fam.df$vh.family))))
  sum.file.vh.fam.df <- data.frame(file.fam, sum.file.vh.fam.df)
  sum.file.vh.gene.df <- data.frame(file.gene, sum.file.vh.gene.df)
  
  #this adds full possible VH family list to the file
  #so that frequencies calculate properly when
  #averaged across donors
  vh.fam.list.df <- data.frame(vh.fam.list)
  increase.id.len <- c(rep_len(file.id, 
                               length(as.character(vh.fam.list.df$vh.fam.list))))
  add.id.full.vh.fam.df <- data.frame(vh.fam.list.df, increase.id.len)
  full.df <- merge(add.id.full.vh.fam.df, sum.file.vh.fam.df,
                   by.x = "vh.fam.list",
                   by.y = "vh.family",
                   all.x = TRUE)
  drops <- c("file.fam")
  full.df <- full.df[ , !(names(full.df) %in% drops)]
  #changing column order
  full.df <- full.df[,c(2,1,3,4)]
  #changing column headers
  colnames(full.df) [1:4] <- c("file.id", "vh.family", "count", "frequency")
  
  #turning NA to 0 so averages calculate properly
  vh.fam.freq.df <- full.df %>%
    complete(file.id, nesting(vh.family),
             fill = list(count = 0,
                         frequency = 0))
  
  #this adds full possible VH ID list to the file
  #so that frequencies calculate properly when
  #averaged across donors
  vh.gene.list.df <- data.frame(vh.gene.list)
  increase.id.len <- c(rep_len(file.id, 
                               length(as.character
                                      (vh.gene.list.df$vh.gene.list))))
  add.id.full.vh.gene.df <- data.frame(vh.gene.list.df, increase.id.len)
  full.df <- merge(add.id.full.vh.gene.df, sum.file.vh.gene.df,
                   by.x = "vh.gene.list",
                   by.y = "single_call",
                   all.x = TRUE)
  drops <- c("file.id")
  full.df <- full.df[ , !(names(full.df) %in% drops)]
  #changing column order
  full.df <- full.df[,c(2,1,3,4)]
  #changing column header names
  colnames(full.df) [1:4] <- c("file.id", "vh.gene", "count", "frequency")
  
  #turning NA to 0 so averages calculate properly
  vh.gene.freq.df <- full.df %>%
    complete(file.id, nesting(vh.gene),
             fill = list(count = 0,
                         frequency = 0))
  
  #write csv files out
  write.csv(vh.fam.freq.df, file = vh.fam.freq.sum.filename)
  write.csv(vh.gene.freq.df, file = vh.gene.freq.sum.filename)
}

#if it throws an error here, shut down RStudio completely and reopen
#this RScript and try running again
lapply(file.list,
       FUN = sum.vh.freq.make.csv,
       the.dir = dir.out)

########STOP HERE########################
#------------------------------------------------
#this code is test code using single file to check for issues
#with code above

library(tidyverse)
library(magrittr)
library(readr)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

#this is where files come from
dir.in <- "~/desktop/10X_BCR_pipelines/output/vh"
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/freq_summary/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file <- "~/desktop/10X_BCR_pipelines/output/vh/vh_single_call_freq_4025-RB-1_all_1.csv"
file.list <- read.csv(file)
df <- data.frame(file.list, header = TRUE)

vh.fam.list <- c("IGHV1", "IGHV2", "IGHV3", "IGHV3/OR16", 
                 "IGHV4", "IGHV5", "IGHV6", "IGHV7")
#need to check the list below for completeness!
vh.gene.list <- c("IGHV1-18", "IGHV1-2", "IGHV1-24", "IGHV1-3", 
                "IGHV1-45", "IGHV1-46", "IGHV1-58", "IGHV1-69", 
                "IGHV1-69-2", "IGHV1-69D", "IGHV1-8", "IGHV2-26",
                "IGHV2-5", "IGHV2-70", "IGHV2-70D", "IGHV3-11", 
                "IGHV3-13", "IGHV3-15", "IGHV3-20", "IGHV3-21", 
                "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-43",
                "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64", 
                "IGHV3-64D", "IGHV3-66", "IGHV3-7", "IGHV3-72", 
                "IGHV3-73", "IGHV3-74", "IGHV3/OR16-12", 
                "IGHV4-28", "IGHV4-30-2", "IGHV4-31", "IGHV4-34", 
                "IGHV4-39", "IGHV4-4", "IGHV4-59",  "IGHV4-61", 
                "IGHV5-10-1", "IGHV5-51", "IGHV6-1", "IGHV7-4-1", 
                "IGHV7-81")

#keep these calls here, it fails on select without them
library(dplyr)
library(tidyr)
#may not need to call tibble
library(tibble)

#this makes summary file with proportions by VH family
#if it throws an error here, shut down RStudio completely and reopen
#this RScript and try running again
sum.file.vh.fam.df <- df %>%
  group_by(vh.family) %>%
  summarise(
    observations = n(),
    vh.fam.freq.sum = sum(proportion)
  )
sum.file.vh.fam.df

#this makes summary file with proportions by VH gene id
sum.file.vh.gene.df <- df %>%
  group_by(single_call) %>%
  summarise(
    observations = n(),
    vh.id.freq.sum = sum(proportion)
  )
sum.file.vh.gene.df

#this generates filename for csv files
#the indices change below based on whether it's using
#a list of files in lapply or not
filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                   "freq_")) [2,1])
#filename1 <- as.character((data.frame(strsplit(file, "/")) [6,1]))
#filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
#filename_final <- as.character((data.frame(strsplit(filename2, "vh_")) [1,1]))
vh.fam.freq.sum.filename <- paste(dir.out,
                                  "vh_fam_freq_sum_",
                                  filename.final, sep="")
vh.gene.freq.sum.filename <- paste(dir.out,
                                 "vh_gene_freq_sum_",
                                 filename.final, sep="")

#this adds sample ID column to both summary files
file.id <- as.character(data.frame(strsplit(filename.final, 
                                            ".csv")) [1,1])

file.gene <- c(rep_len(file.id, 
                       length(as.character(sum.file.vh.gene.df$single_call))))

file.fam <- c(rep_len(file.id,
                      length(as.character(sum.file.vh.fam.df$vh.family))))
sum.file.vh.fam.df <- data.frame(file.fam, sum.file.vh.fam.df)
sum.file.vh.gene.df <- data.frame(file.gene, sum.file.vh.gene.df)

#this adds full possible VH family list to the file
#so that frequencies calculate properly when
#averaged across donors
vh.fam.list.df <- data.frame(vh.fam.list)
increase.id.len <- c(rep_len(filename.final, 
                             length(as.character(vh.fam.list.df$vh.fam.list))))
add.id.full.vh.fam.df <- data.frame(vh.fam.list.df, increase.id.len)
full.df <- merge(add.id.full.vh.fam.df, sum.file.vh.fam.df,
                 by.x = "vh.fam.list",
                 by.y = "vh.family",
                 all.x = TRUE)
drops <- c("file_fam")
full.df <- full.df[ , !(names(full.df) %in% drops)]
full.df <- full.df[,c(2,1,3,4)]
colnames(full.df) [1:4] <- c("file.id", "vh.family", "count", "frequency")

vh.fam.freq.df <- full.df %>%
  complete(file.id, nesting(vh.family),
           fill = list(count = 0,
                       frequency = 0))

#this adds full possible VH ID list to the file
#so that frequencies calculate properly when
#averaged across donors
vh.gene.list.df <- data.frame(vh.gene.list)
increase.id.len <- c(rep_len(filename.final, 
                             length(as.character(vh.gene.list.df$vh.gene.list))))
add.id.full.vh.gene.df <- data.frame(vh.gene.list.df, increase.id.len)
full.df <- merge(add.id.full.vh.gene.df, sum.file.vh.gene.df,
                 by.x = "vh.gene.list",
                 by.y = "single_call",
                 all.x = TRUE)
drops <- c("file.gene")
full.df <- full.df[ , !(names(full.df) %in% drops)]
full.df <- full.df[,c(2,1,3,4)]
colnames(full.df) [1:4] <- c("file.id", "vh.gene", "count", "frequency")

vh.gene.freq.df <- full.df %>%
  complete(file.id, nesting(vh.gene),
           fill = list(count = 0,
                       frequency = 0))

#write csv files out
write.csv(vh.fam.freq.df, file = vh.fam.freq.sum.filename)
write.csv(vh.gene.freq.df, file = vh.gene.freq.sum.filename)
