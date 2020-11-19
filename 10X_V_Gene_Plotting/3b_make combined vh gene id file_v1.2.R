#Run the 2_vh... Rscript first to generate the input files needed here.
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
dir.in <- "~/desktop/10X_BCR_pipelines/output/vh/freq_summary/"
######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/graphs"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)

file.list <- list.files(dir.in, pattern = "*vh_gene_freq_sum",
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
         #                         ifelse(grepl("cd21low", smaller.df$file.id), "cd21low", "other"))

########### EDIT SAMPLE METADATA BELOW #######################
#adding T1D vs. CTL and B cell subset grouping columns
library(stringr)
combined.vh.df <- smaller.df %>%
  mutate(dx.group = case_when(
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
combined.vh.df <- cbind(combined.vh.df, data.frame(paste(combined.vh.df$subset.group, combined.vh.df$dx.group, sep = "_")))
colnames(combined.vh.df) [7] <- c("subset.dx")
write.csv(combined.vh.df, file = paste(dir.out, "/combined_vh_gene_freq.csv", sep = ""))
#data_agg <- aggregate(value ~ index, smaller.df, mean)

#combine VH.id and subset.dx into single column
vh.gene.subset.dx <- (paste(combined.vh.df$vh.gene,
                            combined.vh.df$subset.dx, sep = "_"))
combined.vh.gene.subset.dx.df <- cbind(vh.gene.subset.dx, combined.vh.df)
#reorganize columns in df
combined.vh.gene.subset.dx.df <- combined.vh.gene.subset.dx.df[,c(2:8,1)]

smaller.df <- data.frame(combined.vh.gene.subset.dx.df[, 1:4],
                         combined.vh.gene.subset.dx.df[, 7])
#rename column in df
colnames(smaller.df) [5] <- c("subset.dx")

#eliminate groups in any sample for which n<10  
smaller.df <- smaller.df %>%
  group_by(file.id) %>%
  filter(sum(count) >= 10)

#eliminating all or T1D from file.id
smaller.df <- smaller.df %>%
  #split by "_" and give arbitrary column names to each
  #of 3 new columns that were added
  separate(file.id, c("A", "B", "C"), "_") %>%
  unite(file.id, A, C, sep = "_") %>%
  #ditching column B
  select(!(B)) %>%
  mutate(file.id = as.factor(file.id))

write.csv(smaller.df, file = paste(dir.out, "/filtered_combined_vh_gene_freq.csv", sep = ""))
#consider writing summary table output of n per group by sample here 

#reorganizing structure of df to make group comparisons easier
#by adding each as a separate factor in df
vh.gene.wide.df <- smaller.df %>%
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

#Generating mean values by vh.gene per group (B cell subset)
#keep dplyr call here, it was screwing up without it
#watch to see that these dataframes are generated properly,
#this part is glitchy for some reason
library(dplyr)
all.FDR.mean.freq <- vh.gene.wide.df %>%
  group_by(vh.gene) %>%
  summarise(mean(frequency_all_FDR, na.rm = TRUE))

all.T1D.mean.freq <- vh.gene.wide.df %>%
  group_by(vh.gene) %>%
  summarize(mean(frequency_all_T1D, na.rm = TRUE))

cd21low.FDR.mean.freq <- vh.gene.wide.df %>%
  group_by(vh.gene) %>%
  summarize(mean(frequency_cd21low_FDR, na.rm = TRUE))

cd21low.T1D.mean.freq <- vh.gene.wide.df %>%
  group_by(vh.gene) %>%
  summarize(mean(frequency_cd21low_T1D, na.rm = TRUE))

mean.VH.gene.freq <- merge(merge(merge(all.FDR.mean.freq, 
                                       all.T1D.mean.freq, 
                                       by = "vh.gene", all = TRUE),
                                 cd21low.FDR.mean.freq, 
                                 by = "vh.gene", all = TRUE),
                           cd21low.T1D.mean.freq, 
                           by = "vh.gene", all = TRUE)

colnames(mean.VH.gene.freq) [2:5] <- c("all.FDR.mean", 
                                       "all.T1D.mean",
                                       "cd21low.FDR.mean", 
                                       "cd21low.T1D.mean")

#this pulls VH gene id's for which cd21low T1D > cd21low.FDR
#frequency is TRUE, and writes csv for use with ggplot2
cd21low.T1D.vh.plot.df <- subset(
  mean.VH.gene.freq, cd21low.T1D.mean > cd21low.FDR.mean)

VH.plot.list <- data.frame(cd21low.T1D.vh.plot.df$vh.gene)
#write.csv(mean.VH.gene.freq, file = paste(dir.out, "/mean_VH_gene_freq_by_subset.csv", sep = ""))
write.csv(VH.plot.list, file = paste(dir.out, "/vh_gene_plot_when_cd21low_T1D_increased.csv", sep = ""))
