#Run the 2_vlk... Rscript first to generate the input files needed here.
#This script combines all files into single summary file
#that includes VH gene ID proportion and counts per sample.
#It also adds columns to use as grouping variables
#in ggplot and removes groups by sample for which n<10 sequences
#to avoid inappropriate skewing of the data.

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

file.list <- list.files(dir.in, pattern = "*vlk_gene_freq_sum",
                        full.names = TRUE)
file.list

#this puts all the data into one larger dataframe containing all samples
data.list <- lapply(file.list, read.csv)
smaller.df <- subset(do.call(rbind, data.list), select = -c(X))


########## NEED TO REPLACE SUBSET NAMES FOR YOUR SUBSET NAMES #############
######## E.G. REPLACE "all" AND "cd21low" ACCORDINGLY #####################
#this adds subset group id column
smaller.df$subset.group <- ifelse(grepl("all", smaller.df$file.id), "all", "cd21low")
#if additional subset groups are present, could do the following:
#smaller.df$subset.group <- ifelse(grepl("all", smaller.df$file.id), "all", 
 #                                 ifelse(grepl("cd21low", smaller.df$file.id), "cd21low", "other"))

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

#this creates a column with subset_dx as added column and writes csv
combined.vlk.df <- cbind(combined.vlk.df, 
                         data.frame(paste(combined.vlk.df$subset.group, 
                                          combined.vlk.df$dx_group, 
                                          sep = "_")))
colnames(combined.vlk.df) [6:7] <- c("dx.group", "subset.dx")
write.csv(combined.vlk.df, file = paste(dir.out, "/combined_vlk_gene_freq.csv", sep = ""))
#data_agg <- aggregate(value ~ index, smaller.df, mean)

#combine Vlk.id and subset.dx into single column
vlk.gene.subset.dx <- (paste(combined.vlk.df$vh.gene,
                            combined.vlk.df$subset.dx, sep = "_"))
combined.vlk.gene.subset.dx.df <- cbind(vlk.gene.subset.dx, combined.vlk.df)
#reorganize columns in df
combined.vlk.gene.subset.dx.df <- combined.vlk.gene.subset.dx.df[,c(2:8,1)]

smaller.df <- data.frame(combined.vlk.gene.subset.dx.df[, 1:4],
                         combined.vlk.gene.subset.dx.df[, 7])
#rename column in df
colnames(smaller.df) [5] <- c("subset.dx")

#eliminate groups in any sample for which n<10  
smaller.df <- smaller.df %>%
  group_by(file.id) %>%
  filter(sum(count) >= 10)

#eliminating all or cd21low from file.id
smaller.df <- smaller.df %>%
  #split by "_" and give arbitrary column names to each
  #of 3 new columns that were added
  separate(file.id, c("A", "B", "C"), "_") %>%
  unite(file.id, A, C, sep = "_") %>%
  #ditching column B
  select(!(B)) %>%
  mutate(file.id = as.factor(file.id))

#need this to ensure it's only type data frame, otherwise
#str lists tbl and other stuff too.
smaller.df <- data.frame(smaller.df)

vl.df <- data.frame(filter(smaller.df, grepl(
  "IGLV", smaller.df$vlk.gene)))

vk.df <- data.frame(filter(smaller.df, grepl(
  "IGKV", smaller.df$vlk.gene)))

write.csv(vl.df, file = paste(dir.out, 
                              "/filtered_combined_vl_gene_freq.csv", 
                              sep = ""))
write.csv(vk.df, file = paste(dir.out, 
                              "/filtered_combined_vk_gene_freq.csv", 
                              sep = ""))

#reorganizing structure of df to make group comparisons easier
#by adding each as a separate factor in df
vlk.gene.wide.df <- smaller.df %>%
  pivot_wider(
    names_from = subset.dx,
    values_from = c(frequency, count)
  ) %>%
  #include this call or the output is a tbl, dataframe, and something else
  data.frame()

######### RENAME SUBSET NAMES BELOW ACCORDINGLY FOR YOUR DATA ###########
#I left my subset names because i thought it was better than 1 vs. 2
#for you to follow what is happening and is also nice to have that 
#defined in your output CSV files. You'll need to edit through the
#end of this code.

#Generating mean values by vlk.gene per group (B cell subset)
#keep dplyr call here, it was screwing up without it
#watch to be sure these dataframes are made correctly
#this part is glitchy for some reason
library(dplyr)
all.FDR.mean.freq <- vlk.gene.wide.df %>%
  group_by(vlk.gene) %>%
  summarize(mean(frequency_all_FDR, na.rm = TRUE))

all.T1D.mean.freq <- vlk.gene.wide.df %>%
  group_by(vlk.gene) %>%
  summarize(mean(frequency_all_T1D, na.rm = TRUE))

cd21low.FDR.mean.freq <- vlk.gene.wide.df %>%
  group_by(vlk.gene) %>%
  summarize(mean(frequency_cd21low_FDR, na.rm = TRUE))

cd21low.T1D.mean.freq <- vlk.gene.wide.df %>%
  group_by(vlk.gene) %>%
  summarize(mean(frequency_cd21low_T1D, na.rm = TRUE))

mean.Vlk.gene.freq <- merge(merge(merge(all.FDR.mean.freq, 
                                       all.T1D.mean.freq, 
                                       by = "vlk.gene", all = TRUE),
                                 cd21low.FDR.mean.freq, 
                                 by = "vlk.gene", all = TRUE),
                           cd21low.T1D.mean.freq, 
                           by = "vlk.gene", all = TRUE)

colnames(mean.Vlk.gene.freq) [2:5] <- c("all.FDR.mean", 
                                        "all.T1D.mean",
                                        "cd21low.FDR.mean", 
                                        "cd21low.T1D.mean")

#this pulls VH gene id's for which cd21low T1D > cd21low.FDR
#frequency is TRUE and writes csv for use with ggplot2
cd21low.T1D.vlk.plot.df <- subset(
  mean.Vlk.gene.freq, cd21low.T1D.mean > cd21low.FDR.mean)

Vlk.plot.list <- data.frame(cd21low.T1D.vlk.plot.df$vlk.gene)
colnames(Vlk.plot.list) [1] <- c("cd21low.T1D.vlk")

Vl.plot.list <- data.frame(filter(Vlk.plot.list, grepl(
  "IGLV", Vlk.plot.list$cd21low.T1D.vlk)))
colnames(Vl.plot.list) [1] <- c("cd21low.T1D.vl")

Vk.plot.list <- data.frame(filter(Vlk.plot.list, grepl(
  "IGKV", Vlk.plot.list$cd21low.T1D.vlk)))
colnames(Vk.plot.list) [1] <- c("cd21low.T1D.vk")

#write.csv(mean.VH.gene.freq, file = paste(dir.out, "/mean_Vlk_gene_freq_by_subset.csv", sep = ""))
#write.csv(Vlk.plot.list, file = paste(dir.out, "/Vlk_plot_when_cd21low_T1D_increased.csv", sep = ""))
write.csv(Vl.plot.list, 
          file = paste(dir.out, 
                       "/vl_gene_plot_when_cd21low_T1D_increased.csv", 
                       sep = ""))
write.csv(Vk.plot.list, 
          file = paste(dir.out, 
                       "/vk_gene_plot_when_cd21low_T1D_increased.csv", 
                       sep = ""))
