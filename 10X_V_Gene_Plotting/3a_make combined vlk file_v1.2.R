#This script combines all files into single summary file
#that includes VL/VK family proportion and counts per sample.
#It also adds columns to use as grouping variables
#in ggplot and removes groups by sample for which n<10 sequences
#to avoid inappropriate skewing of the data.

#A summary table is output as a CSV file to show the count
#and frequency of all Vgenes combined per sample (file.id)
#and subset.dx. This allows you to determine if one sample
#is skewing the data set.

#Things you may need to revise are marked:
######ALL CAPS#########

library(tidyverse)
library(magrittr)
library(readr)

#error here if you try to run as one block of code with above
check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

######DEFINE DIR IN########
#this is where files come from
dir.in <- "~/desktop/10X_BCR_pipelines/output/vlk/freq_summary/"

######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/graphs"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*vlk_fam_freq_sum",
                        full.names = TRUE)

#this puts all the data into one larger dataframe containing all samples
#this puts all the data into one larger dataframe containing all samples
data.list <- lapply(file.list, read.csv)
smaller.df <- subset(do.call(rbind, data.list), select = -c(X))

##########NEED TO REPLACE SUBSET NAMES FOR YOUR SUBSET NAMES#############
########E.G. REPLACE "all" AND "cd21low" ACCORDINGLY
#this adds subset group id column
smaller.df$subset.group <- ifelse(grepl("all", smaller.df$file.id), 
                                  "all", "cd21low")
#if additional subset groups are present, could do the following:
#smaller.df$subset.group <- ifelse(grepl("all", 
#                                        smaller.df$file.id), "all",
#                                  ifelse(grepl("cd21low", smaller.df$file.id), 
#                                         "cd21low", "other"))

########### EDIT SAMPLE METADATA BELOW #######################
#adding T1D vs. CTL and B cell subset grouping columns
library(stringr)
combined.vlk.df <- smaller.df %>%
  mutate(dx_group = case_when(
    #all
    smaller.df$file.id == "4025-RB-1_all_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-1_all_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-1_all_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-2_all_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-2_all_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-2_all_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-3_all_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-3_all_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-3_all_5" ~ "FDR",
    smaller.df$file.id == "4025-RB-4_all_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-4_all_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-4_all_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-5_all_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-5_all_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-5_all_5" ~ "FDR",
    smaller.df$file.id == "4025-RB-6_all_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-6_all_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-6_all_5" ~ "FDR",
    #cd21low
    smaller.df$file.id == "4025-RB-1_cd21low_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-1_cd21low_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-1_cd21low_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-2_cd21low_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-2_cd21low_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-2_cd21low_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-3_cd21low_1" ~ "FDR",
    smaller.df$file.id == "4025-RB-3_cd21low_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-3_cd21low_5" ~ "FDR",
    smaller.df$file.id == "4025-RB-4_cd21low_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-4_cd21low_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-4_cd21low_5" ~ "T1D",
    smaller.df$file.id == "4025-RB-5_cd21low_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-5_cd21low_2" ~ "FDR",
    smaller.df$file.id == "4025-RB-5_cd21low_5" ~ "FDR",
    smaller.df$file.id == "4025-RB-6_cd21low_1" ~ "T1D",
    smaller.df$file.id == "4025-RB-6_cd21low_2" ~ "T1D",
    smaller.df$file.id == "4025-RB-6_cd21low_5" ~ "FDR"
    )
    )

#this creates a column with subset_dx as added column
combined.vlk.df <- cbind(combined.vlk.df, 
                         data.frame(paste(combined.vlk.df$subset.group, 
                                          combined.vlk.df$dx_group, 
                                          sep = "_")))
colnames(combined.vlk.df) [6:7] <- c("dx.group", "subset.dx")
write.csv(combined.vlk.df, 
          file = paste(dir.out, "/combined_vlk_fam_freq.csv", 
                       sep = ""))

#combine VLK.fam and subset.dx into single column
vlk.fam.subset.dx <- (paste(combined.vlk.df$vlk.family,
                           combined.vlk.df$subset.dx, sep = "_"))
combined.vlk.fam.subset.dx.df <- cbind(vlk.fam.subset.dx, 
                                       combined.vlk.df)
#reorganize columns in df
combined.vlk.fam.subset.dx.df <- combined.vlk.fam.subset.dx.df[,c(2:8,1)]

smaller.df <- data.frame(combined.vlk.fam.subset.dx.df[, 1:4],
                         combined.vlk.fam.subset.dx.df[, 7])
#rename column in df
colnames(smaller.df) [5] <- c("subset.dx")

#eliminating all or cd21low from file.id
smaller.df <- smaller.df %>%
  #split by "_" and give arbitrary column names to each
  #of 3 new columns that were added
  separate(file.id, c("A", "B", "C"), "_") %>%
  unite(file.id, A, C, sep = "_") %>%
  #ditching column B
  select(!(B)) %>%
  mutate(file.id = as.factor(file.id))

vl.df <- data.frame(filter(smaller.df, grepl(
  "IGLV", smaller.df$vlk.family)))

vk.df <- data.frame(filter(smaller.df, grepl(
  "IGKV", smaller.df$vlk.family)))

write.csv(vl.df, file = paste(dir.out, 
                              "/filtered_combined_vl_fam_freq.csv", 
                              sep = ""))
write.csv(vk.df, file = paste(dir.out, 
                              "/filtered_combined_vk_fam_freq.csv", 
                              sep = ""))

#summary table output of n per group by sample 
BCR.count.df<- vl.df %>%
  group_by(subset.dx, file.id) %>%
  summarise(counts = sum(count, na.rm = TRUE)) %>%
  mutate(freq = counts / sum(counts))
write.csv(BCR.count.df, file = paste(dir.out, 
                                     "/VL_count_table_by_sample_subsetdx.csv", 
                                     sep = ""))

BCR.count.df<- vk.df %>%
  group_by(subset.dx, file.id) %>%
  summarise(counts = sum(count, na.rm = TRUE)) %>%
  mutate(freq = counts / sum(counts))
write.csv(BCR.count.df, file = paste(dir.out, 
                                     "/VK_count_table_by_sample_subsetdx.csv", 
                                     sep = "")) 

#reorganizing structure of df to make group comparisons easier
#by adding each as a separate factor in df
vlk.fam.wide.df <- smaller.df %>%
  pivot_wider(
    names_from = subset.dx,
    values_from = c(frequency, count)
  ) %>%
  #include this call or the output is a tbl, dataframe, and something else
  data.frame()

#generating mean values by vlk.gene per group (B cell subset)
#make sure it generates these correctly, if not, restart RStudio
#this has been a common glitch for me.
all.FDR.mean.freq <- vlk.fam.wide.df %>%
  group_by(vlk.family) %>%
  summarize(mean(frequency_all_FDR, na.rm = TRUE))

all.T1D.mean.freq <- vlk.fam.wide.df %>%
  group_by(vlk.family) %>%
  summarize(mean(frequency_all_T1D, na.rm = TRUE))

cd21low.FDR.mean.freq <- vlk.fam.wide.df %>%
  group_by(vlk.family) %>%
  summarize(mean(frequency_cd21low_FDR, na.rm = TRUE))

cd21low.T1D.mean.freq <- vlk.fam.wide.df %>%
  group_by(vlk.family) %>%
  summarize(mean(frequency_cd21low_T1D, na.rm = TRUE))

mean.Vlk.fam.freq <- merge(merge(merge(all.FDR.mean.freq, 
                                      all.T1D.mean.freq, 
                                      by = "vlk.family", all = TRUE),
                                cd21low.FDR.mean.freq, 
                                by = "vlk.family", all = TRUE),
                          cd21low.T1D.mean.freq, 
                          by = "vlk.family", all = TRUE)

colnames(mean.Vlk.fam.freq) [2:5] <- c("all.FDR.mean", 
                                       "all.T1D.mean",
                                       "cd21low.FDR.mean",
                                       "cd21low.T1D.mean")

#this pulls VH gene id's for which cd21low T1D > cd21low.FDR
#frequency is TRUE, and writes csv for use with ggplot2
cd21low.T1D.vlk.plot.df <- subset(
  mean.Vlk.fam.freq, cd21low.T1D.mean > cd21low.FDR.mean)

cd21low.T1D.vl.plot.df <- data.frame(filter(cd21low.T1D.vlk.plot.df, grepl(
  "IGLV", cd21low.T1D.vlk.plot.df$vlk.family)))

cd21low.T1D.vk.plot.df <- data.frame(filter(cd21low.T1D.vlk.plot.df, grepl(
  "IGKV", cd21low.T1D.vlk.plot.df$vlk.family)))

write.csv(cd21low.T1D.vl.plot.df, 
          file = paste(dir.out, 
                       "/vl_fam_plot_when_cd21low_T1D_increased.csv", 
                       sep = ""))
write.csv(cd21low.T1D.vk.plot.df, 
          file = paste(dir.out, 
                       "/vk_fam_plot_when_cd21low_T1D_increased.csv", 
                       sep = ""))

#----------------------------
#STOP HERE
#not sure i need to do below?
#code below adds zero rows for missing VL/Vk family rows per sample

new.df <- combined.vlk.df %>% 
  complete(file.id, nesting(subset.group, vlk.family), 
           fill = list(frequency = 0,
                       count = 0))
#data %>% complete(vlk.family, nesting(item_id, item_name), fill = list(vlk.fam.freq.sum = 0, observations = 0))

write.csv(combined.vlk.df, file = paste(dir.out, "/combined_Vlk_fam_freq.csv", sep = ""))
#data_agg <- aggregate(value ~ index, smaller.df, mean)