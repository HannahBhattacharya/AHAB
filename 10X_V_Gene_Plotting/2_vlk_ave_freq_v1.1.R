#Run the 1_vlk... Rscript first to generate the input files needed here.
#This code takes output files from 1_vh_single_call_family Rscript
#and collapses data for each sample to show counts and
#proportions of VL or VK family calls and VL or VK gene ID calls

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
dir.in <- "~/desktop/10X_BCR_pipelines/output/vlk"

######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/freq_summary/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*.csv",
                        full.names = TRUE)
file.list
vlk.fam.list <- c("IGKV1", "IGKV1D", "IGKV2", "IGKV2D", 
                  "IGKV3", "IGKV3/OR2", "IGKV3D", "IGKV4",
                  "IGKV5", "IGKV6", "IGKV6D", "IGLV1", 
                  "IGLV2", "IGLV3", "IGLV4", "IGLV5", "IGLV6",
                  "IGLV7", "IGLV8", "IGLV9", "IGLV10")
#need to check the list below for completeness!

vlk.gene.list <- c("IGKV1-5", "IGKV1-6", "IGKV1-8", "IGKV1D-8", 
                 "IGKV1-9", "IGKV1-12", "IGKV1D-12", "IGKV1D-13",
                 "IGKV1-16", "IGKV1D-16", "IGKV1-17", "IGKV1D-17", 
                 "IGKV1-27", "IGKV1-33", "IGKV1D-33", "IGKV1-37", 
                 "IGKV1D-37", "IGKV1-39", "IGKV1-42", "IGKV1-43",
                 "IGKV2-18", "IGKV2-24", "IGKV2D-24", "IGKV2-26", 
                 "IGKV2D-26", "IGKV2-28", "IGKV2D-28", "IGKV2-29", 
                 "IGKV2D-29", "IGKV2-30", "IGKV2D-30", "IGKV2-40", 
                 "IGKV2D-40",  "IGKV3D-7", "IGKV3-11", "IGKV3D-11", 
                 "IGKV3-15", "IGKV3D-15", "IGKV3-20", "IGKV3D-20",
                 "IGKV3/OR2-268", "IGKV4-1", "IGKV5-2", "IGKV6-21", 
                 "IGKV6D-21", "IGKV6D-41", 
                 "IGLV1-36", "IGLV1-40", "IGLV1-41", "IGLV1-44", 
                 "IGLV1-47", "IGLV1-50", "IGLV1-51", "IGLV2-8", 
                 "IGLV2-11", "IGLV2-14", "IGLV2-18", "IGLV2-23", 
                 "IGLV2-33", "IGLV3-1", "IGLV3-9", "IGLV3-10", 
                 "IGLV3-12", "IGLV3-16", "IGLV3-19", "IGLV3-21", 
                 "IGLV3-22", "IGLV3-25", "IGLV3-27", "IGLV3-32", 
                 "IGLV4-3", "IGLV4-60", "IGLV4-69", "IGLV5-37",
                 "IGLV5-39", "IGLV5-45", "IGLV5-48", "IGLV5-52", 
                 "IGLV6-57", "IGLV7-43",  "IGLV7-46", "IGLV8-61", 
                 "IGLV9-49", "IGLV10-54", "IGLV11-55")

sum.vlk.freq.make.csv <- function(file, the.dir) {
  df <- data.frame(read.csv(file, header = TRUE))
  df2 <- data.frame(read.csv(file, header = TRUE))
  #keep these calls here, it fails on select without them
  library(dplyr)
  library(tidyr)
  #may not need to call tibble
  library(tibble)
  
  #this makes summary file with proportions by VH family
  sum.file.vlk.fam.df <- df %>%
    group_by(vlk.family) %>%
    summarise(
      observations = n(),
      vlk.fam.freq.sum = sum(proportion)
    )
  sum.file.vlk.fam.df

  #this makes summary file with proportions by VL or VK gene id
  sum.file.vlk.gene.df <- df %>%
    group_by(single_call) %>%
    summarise(
      observations = n(),
      vlk.gene.freq.sum = sum(proportion)
    )
  sum.file.vlk.gene.df
  
  #this generates filename for csv files
  filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                     "freq_")) [2,1])
  vlk.fam.freq.sum.filename <- paste(dir.out,
                                     "vlk_fam_freq_sum_",
                                     filename.final, sep="")
  vlk.gene.freq.sum.filename <- paste(dir.out,
                                    "vlk_gene_freq_sum_",
                                    filename.final, sep="")

  #this adds sample ID column to both summary files
  file.id <- as.character(data.frame(strsplit(filename.final, 
                                              ".csv")) [1,1])
  file.gene <- c(rep_len(file.id, 
                         length(as.character
                                (sum.file.vlk.gene.df$single_call))))
  file.fam <- c(rep_len(file.id,
                        length(as.character
                               (sum.file.vlk.fam.df$vlk.family))))
  sum.file.vlk.fam.df <- data.frame(file.fam, sum.file.vlk.fam.df)
  sum.file.vlk.gene.df <- data.frame(file.gene, sum.file.vlk.gene.df)
 
  #this adds full possible VH family list to the file
  #so that frequencies calculate properly when
  #averaged across donors
  vlk.fam.list.df <- data.frame(vlk.fam.list)
  increase.id.len <- c(rep_len(file.id, 
                               length(as.character
                                      (vlk.fam.list.df$vlk.fam.list))))
  add.id.full.vlk.fam.df <- data.frame(vlk.fam.list.df, increase.id.len)
  full.df <- merge(add.id.full.vlk.fam.df, sum.file.vlk.fam.df,
                   by.x = "vlk.fam.list",
                   by.y = "vlk.family",
                   all.x = TRUE)
  drops <- c("file.fam")
  full.df <- full.df[ , !(names(full.df) %in% drops)]
  #changing column order
  full.df <- full.df[,c(2,1,3,4)]
  #changing column headers
  colnames(full.df) [1:4] <- c("file.id", "vlk.family", 
                               "count", "frequency")
  
  #turning NA to 0 so averages calculate properly
  vlk.fam.freq.df <- full.df %>%
    complete(file.id, nesting(vlk.family),
             fill = list(count = 0,
                         frequency = 0))
  
  #this adds full possible VL or VK ID list to the file
  #so that frequencies calculate properly when
  #averaged across donors
  vlk.gene.list.df <- data.frame(vlk.gene.list)
  increase.id.len <- c(rep_len(file.id, 
                               length(as.character
                                      (vlk.gene.list.df$vlk.gene.list))))
  add.id.full.vlk.gene.df <- data.frame(vlk.gene.list.df, increase.id.len)
  full.df <- merge(add.id.full.vlk.gene.df, sum.file.vlk.gene.df,
                   by.x = "vlk.gene.list",
                   by.y = "single_call",
                   all.x = TRUE)
  drops <- c("file.gene")
  full.df <- full.df[ , !(names(full.df) %in% drops)]
  #changin column order
  full.df <- full.df[,c(2,1,3,4)]
  #changing column header names
  colnames(full.df) [1:4] <- c("file.id", "vlk.gene", 
                               "count", "frequency")
  
  #turning NA to 0 so averages calculate properly
  vlk.gene.freq.df <- full.df %>%
    complete(file.id, nesting(vlk.gene),
             fill = list(count = 0,
                         frequency = 0))
  
  #write csv files out
  write.csv(vlk.fam.freq.df, file = vlk.fam.freq.sum.filename)
  write.csv(vlk.gene.freq.df, file = vlk.gene.freq.sum.filename)
}

#if it throws an error here, shut down RStudio completely and reopen
#this RScript and try running again
lapply(file.list,
       FUN = sum.vlk.freq.make.csv,
       the.dir = dir.out)

###########STOP HERE##################
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
dir.in <- "~/desktop/10X_BCR_pipelines/output/vlk"
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/freq_summary/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file <- "~/desktop/10X_BCR_pipelines/output/vlk/vlk_single_call_freq_4025-RB-1_all_1.csv"
file.list <- read.csv(file)
df <- data.frame(file.list, header = TRUE)

#re-run this code with summary file that includes all files
#after making them to be sure I didn't miss any
levels(df$vlk.family)
levels(df$single_call)
vlk.fam.list <- c("IGKV1", "IGKV1D", "IGKV2", "IGKV2D", "IGKV3", 
                  "IGKV3/OR2", "IGKV3D", "IGKV4", "IGKV5", "IGKV6",
                  "IGKV6D", 
                  "IGLV1", "IGLV2", "IGLV3", "IGLV4", "IGLV5", 
                  "IGLV6", "IGLV7", "IGLV8", "IGLV9", "IGLV10",
                  "test")

#need to check the list below for completeness!
vlk.gene.list <- c("IGKV1-5", "IGKV1-6", "IGKV1-8", "IGKV1D-8", 
                 "IGKV1-9", "IGKV1-12", "IGKV1D-12", "IGKV1D-13",
                 "IGKV1-16", "IGKV1D-16", "IGKV1-17", "IGKV1D-17", 
                 "IGKV1-27", "IGKV1-33", "IGKV1D-33", "IGKV1-37", 
                 "IGKV1D-37", "IGKV1-39", "IGKV1-42", "IGKV1-43",
                 "IGKV2-18", "IGKV2-24", "IGKV2D-24", "IGKV2-26", 
                 "IGKV2D-26", "IGKV2-28", "IGKV2D-28", "IGKV2-29", 
                 "IGKV2D-29", "IGKV2-30", "IGKV2D-30", "IGKV2-40", 
                 "IGKV2D-40",  "IGKV3D-7", "IGKV3-11", "IGKV3D-11", 
                 "IGKV3-15", "IGKV3D-15", "IGKV3-20", "IGKV3D-20",
                 "IGKV3/OR2-268", "IGKV4-1", "IGKV5-2", "IGKV6-21", 
                 "IGKV6D-21", "IGKV6D-41", 
                 "IGLV1-36", "IGLV1-40", "IGLV1-41", "IGLV1-44", 
                 "IGLV1-47", "IGLV1-50", "IGLV1-51", "IGLV2-8", 
                 "IGLV2-11", "IGLV2-14", "IGLV2-18", "IGLV2-23", 
                 "IGLV2-33", "IGLV3-1", "IGLV3-9", "IGLV3-10", 
                 "IGLV3-12", "IGLV3-16", "IGLV3-19", "IGLV3-21", 
                 "IGLV3-22", "IGLV3-25", "IGLV3-27", "IGLV3-32", 
                 "IGLV4-3", "IGLV4-60", "IGLV4-69", "IGLV5-37",
                 "IGLV5-39", "IGLV5-45", "IGLV5-48", "IGLV5-52", 
                 "IGLV6-57", "IGLV7-43",  "IGLV7-46", "IGLV8-61", 
                 "IGLV9-49", "IGLV10-54", "IGLV11-55")

#keep these calls here, it fails on select without them
library(dplyr)
library(tidyr)
#may not need to call tibble
library(tibble)

#this makes summary file with proportions by VL or VK family
#if it throws an error here, shut down RStudio completely and reopen
#this RScript and try running again
sum.file.vlk.fam.df <- df %>%
  group_by(vlk.family) %>%
  summarise(
    observations = n(),
    vlk.fam.freq.sum = sum(proportion)
  )
sum.file.vlk.fam.df

#this makes summary file with proportions by VL or VK gene id
sum.file.vlk.gene.df <- df %>%
  group_by(single_call) %>%
  summarise(
    observations = n(),
    vlk.gene.freq.sum = sum(proportion)
  )
sum.file.vlk.gene.df

#this generates filename for csv files
#the indices change below based on whether it's using
#a list of files in lapply or not
filename.final <- as.character(data.frame(strsplit(basename(file), 
                                                   "freq_")) [2,1])
#filename1 <- as.character((data.frame(strsplit(file, "/")) [6,1]))
#filename2 <- as.character((data.frame(strsplit(filename1, "[.]")) [1,1]))
#filename_final <- as.character((data.frame(strsplit(filename2, "vlk_")) [1,1]))
vlk.fam.freq.sum.filename <- paste(dir.out,
                                   "vlk_fam_freq_sum_",
                                   filename.final, sep="")
vlk.gene.freq.sum.filename <- paste(dir.out,
                                  "vlk_id_freq_sum_",
                                  filename.final, sep="")

#this adds sample ID column to both summary files
file.id <- as.character(data.frame(strsplit(filename.final, 
                                            ".csv")) [1,1])
file.gene <- c(rep_len(file.id, 
                       length(as.character(sum.file.vlk.gene.df$single_call))))
file.fam <- c(rep_len(file.id,
                      length(as.character(sum.file.vlk.fam.df$vlk.family))))
sum.file.vlk.fam.df <- data.frame(file.fam, sum.file.vlk.fam.df)
sum.file.vlk.gene.df <- data.frame(file.gene, sum.file.vlk.gene.df)

#this adds full possible VL or VK family list to the file
#so that frequencies calculate properly when
#averaged across donors
vlk.fam.list.df <- data.frame(vlk.fam.list)
increase.id.len <- c(rep_len(file.id, 
                             length(as.character(vlk.fam.list.df$vlk.fam.list))))
add.id.full.vlk.fam.df <- data.frame(vlk.fam.list.df, increase.id.len)
full.df <- merge(add.id.full.vlk.fam.df, sum.file.vlk.fam.df,
                 by.x = "vlk.fam.list",
                 by.y = "vlk.family",
                 all.x = TRUE)
drops <- c("file.fam")
full.df <- full.df[ , !(names(full.df) %in% drops)]
full.df <- full.df[,c(2,1,3,4)]
colnames(full.df) [1:4] <- c("file.id", "vlk.family", "count", "frequency")

vlk.fam.freq.df <- full.df %>%
  complete(file.id, nesting(vlk.family),
           fill = list(count = 0,
                       frequency = 0))

#this adds full possible VH ID list to the file
#so that frequencies calculate properly when
#averaged across donors
vlk.gene.list.df <- data.frame(vlk.gene.list)
increase.id.len <- c(rep_len(file.id, 
                             length(as.character(vlk.gene.list.df$vlk.gene.list))))
add.id.full.vlk.gene.df <- data.frame(vlk.gene.list.df, increase.id.len)
full.df <- merge(add.id.full.vlk.gene.df, sum.file.vlk.gene.df,
                 by.x = "vlk.gene.list",
                 by.y = "single_call",
                 all.x = TRUE)
drops <- c("file.gene")
full.df <- full.df[ , !(names(full.df) %in% drops)]
full.df <- full.df[,c(2,1,3,4)]
colnames(full.df) [1:4] <- c("file.id", "vlk.gene", "count", "frequency")

vlk.gene.freq.df <- full.df %>%
  complete(file.id, nesting(vlk.gene),
           fill = list(count = 0,
                       frequency = 0))

#write csv files out
write.csv(vlk.fam.freq.df, file = vlk.fam.freq.sum.filename)
write.csv(vlk.gene.freq.df, file = vlk.gene.freq.sum.filename)
