#This script takes single chains export file 
#from 10X Genomics Loupe VDJ Browser as input
#and pulls only the VH data from it. CSV output has 
#filename, chain, clonotype_id, 
#single Vk call (single chain), proportion, 
#and vk_family

########################
#SEE ALL CAPS FOR THINGS YOU'RE LIKELY TO NEED
#TO EDIT BASED ON YOUR SPECIFIC NEEDS

library(tidyverse)
library(magrittr)
library(readr)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

#########################
#DEFINE FILEPATH WHERE LCB (LOUPE CELL BROWSER)
#FILES ARE COMING FROM
dir.in <- "~/desktop/10X_BCR_pipelines/input"

####################
#DEFINE FILEPATH OUT
dir.out <- "~/desktop/10X_BCR_pipelines/REDCap/out/"

directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)
library(stringr)

#####################################
#THIS NEEDS TO POINT TO LCB FILE PATH
#THIS IS HOW YOU FILTER DOWN TO B CELL
#SUBSET, E.G. CD21LOW
#"cd21low.+csv" will find any file with the cd21low
# subset name in the csv file name
subset.name <- "all"
file.list <- list.files(dir.in, pattern = "all.+csv",
                        full.names = TRUE)
file.list

#this puts all the data into one list of data frames,
#one df for each samples
data.list <- lapply(file.list, read.csv)

#pull file id without leading path, subset or .csv
####################
#MAY NEED TO ALTER BASED ON FILE PATH HIERARCHY/NAME
make.names <- function(file) {
  filename <- as.character((data.frame
                            (strsplit(file, "/")) [7,1]))
  filename <- as.character((data.frame
                            (strsplit(filename, "usage_")) [2,1]))
  filename <- as.character((data.frame
                            (strsplit(filename, ".csv")) [1,1]))
  filename <- gsub(subset.name, "", filename)
  filename <- gsub("__", "_", filename)
}
filenames <- lapply(file.list, make.names)
filenames
#insert filename information into each
#data frame in list
data.list <- Map(cbind, data.list, file_id = filenames)

##########################################
#ALTER THIS TO INCLUDE WHATEVER METADATA YOU NEED
#FOR DOWNSTREAM PLOTTING; WHAT GROUPS DO YOU WANT
#TO COMPARE? KEEP TWO GROUPS NAMED AS FDR OR T1D
#OR EDIT ALL DOWNSTREAM VARIABLE NAMES IN CODE
#define dx group based on order of samples in list
disease <- c("FDR","T1D","T1D",
             "FDR", "FDR","T1D",
             "FDR","T1D","FDR",
             "T1D","FDR","T1D",
             "T1D","FDR","FDR",
             "T1D","T1D","FDR")

subset.list <- c(subset.name,subset.name,subset.name,
                 subset.name,subset.name,subset.name,
                 subset.name,subset.name,subset.name,
                 subset.name,subset.name,subset.name,
                 subset.name,subset.name,subset.name,
                 subset.name,subset.name,subset.name)

data.list <- Map(cbind, data.list, dx_group = disease)
data.list <- Map(cbind, data.list, subset_name = subset.list)
#this creates a new list of data frames that have
#dx_group == FDR; returns zeros for all values in 
#data frames do not have dx_group == FDR (e.g. T1D)
#I can't figure out how to avoid this, just live with them
#they will get tossed in a minute by rbind
FDR.data.list <- lapply(data.list, function(x) subset(x, dx_group == "FDR"))
T1D.data.list <- lapply(data.list, function(x) subset(x, dx_group == "T1D"))

#put FDR samples into one dataframe and T1D
#samples into a 2nd dataframe
#omit X column that we don't need
#eliminating proportion and frequency so they aren't
#inadvertantly used later
FDR.df <- subset(do.call(rbind, FDR.data.list), select = -c(X, proportion, frequency))
T1D.df <- subset(do.call(rbind, T1D.data.list), select = -c(X, proportion, frequency))
#removing things we don't need anymore
rm(data.list, FDR.data.list, T1D.data.list,
   disease, subset.list)

#need to put separate row for each barcode and repeat
#other data. Otherwise the proportion calls later are
#totally messed up
T1D.df <- separate_rows(T1D.df, barcodes, sep = ";")
FDR.df <- separate_rows(FDR.df, barcodes, sep = ";")

#filtering based on which heavy or light chain 
#each cell expresses; make sure you change 1:14
#to reflect the total # of columns you have 
#1:16 if you have 16 columns
T1D.vh <- select(filter(T1D.df, chain == "IGH"), c(1:14))
T1D.vk <- select(filter(T1D.df, chain == "IGK"), c(1:14))
T1D.vl <- select(filter(T1D.df, chain == "IGL"), c(1:14))
FDR.vh <- select(filter(FDR.df, chain == "IGH"), c(1:14))
FDR.vk <- select(filter(FDR.df, chain == "IGK"), c(1:14))
FDR.vl <- select(filter(FDR.df, chain == "IGL"), c(1:14))

#combine all together into one list of dataframes
df.list <- list(T1D.vh, T1D.vk, T1D.vl,
             FDR.vh, FDR.vk, FDR.vl)

#pull file id without leading path, subset or .csv
make.single.calls <- function(x) {
  vh.single.call <- x %>% 
    separate(v_genes,
             c("v_genes"),
             sep = "[;]", extra = "drop", 
             fill = "left")
  return(vh.single.call)
}
df.list <- lapply(df.list, make.single.calls)

#separating elements back out of list of 
#dataframes again
T1D.vh <- data.frame(df.list[1])
T1D.vk <- data.frame(df.list[2])
T1D.vl <- data.frame(df.list[3])
FDR.vh <- data.frame(df.list[4])
FDR.vk <- data.frame(df.list[5])
FDR.vl <- data.frame(df.list[6])

#write out T1D csv files
file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_T1D_vh.csv")
write_csv(T1D.vh, path = file.name, 
          col_names = TRUE)

file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_T1D_vl.csv")
write_csv(T1D.vl, path = file.name, 
          col_names = TRUE)
file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_T1D_vk.csv")
write_csv(T1D.vk, path = file.name, 
          col_names = TRUE)

#Write out FDR files
file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_FDR_vh.csv")
write_csv(FDR.vh, path = file.name, 
          col_names = TRUE)

file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_FDR_vl.csv")
write_csv(FDR.vl, path = file.name, 
          col_names = TRUE)
file.name <- paste(dir.out, sep = "",
                   as.character(subset.name), 
                   "_FDR_vk.csv")
write_csv(FDR.vk, path = file.name, 
          col_names = TRUE)
#remove things we no longer need
rm(df.list, filenames, FDR.df, T1D.df, subset.name,
   make.names, make.single.calls)
###########################################
#DEFINE FILE PATH OF ANTIGEN-SPECIFIC BCR CSV FILE
#THAT COMES FROM REDCAP BCR DATABASE EXPORT
#Make anti-insulin BCRs data frame
#this is where csv files come from
dir.in <- "~/desktop/10X_BCR_pipelines/REDCap"

##############################
#DEFINE FILE PATH OUT
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/REDCap/graphs/"

directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

file.list <- list.files(dir.in, pattern = "*.csv",
                        full.names = TRUE)
file.list
ins.df <- read.csv(file.list, header = TRUE)

#################################
#REVISE THIS TO INCLUDE WHATEVER COLUMNS OF CSV FILE
#YOU WANT TO BRING INTO THIS ANALYSIS, BUT SHOULDN'T
#NEED TO CHANGE
ins.df <- data.frame(ins.df$Monoclonal.hybridoma.line.ID,
                     ins.df$VL..kappa.or.lambda..gene.id..scroll.down.for.lambda.,
                     ins.df$VH.gene.ID,
                     ins.df$Laboratory.that.generated.this.mAb,
                     ins.df$Autoantigen.reactivity.)

colnames(ins.df)[1:5] <- c("mAb_line", "vk_vl", "vh",
                           "lab_source", "reactivity")

#keep these so NA added and blank cells not
#included as a factor downstream
ins.df$vh <- ifelse(grepl
                ("IGH", ins.df$vh),
                ins.df$vh, NA)
ins.df$vk <- ifelse(grepl
                 ("IGK", ins.df$vk_vl),
                 ins.df$vk_vl, NA)
ins.df$vl <- ifelse(grepl
                    ("IGL", ins.df$vk_vl),
                    ins.df$vk_vl, NA)

###########################
#EDIT VARIABLES/NAMES BELOW ACCORDINGLY
#Make data frames for Chi square test later

#vh
FDR.vh.small <- data.frame(subset
                           (FDR.vh, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(FDR.vh.small)[1:3] <- c("cell_id", "v_genes", "group")

T1D.vh.small <- data.frame(subset
                           (T1D.vh, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(T1D.vh.small)[1:3] <- c("cell_id", "v_genes", "group")

ins.vh.small <- data.frame(subset
                           (ins.df, 
                             select = (c(mAb_line,
                                         vh,
                                         reactivity))))
colnames(ins.vh.small)[1:3] <- c("cell_id", "v_genes", "group")

vh.df <- rbind(FDR.vh.small, T1D.vh.small, ins.vh.small)

#Vkappa
FDR.vk.small <- data.frame(subset
                           (FDR.vk, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(FDR.vk.small)[1:3] <- c("cell_id", "v_genes", "group")

T1D.vk.small <- data.frame(subset
                           (T1D.vk, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(T1D.vk.small)[1:3] <- c("cell_id", "v_genes", "group")

ins.vk.small <- data.frame(subset
                           (ins.df, 
                             select = (c(mAb_line,
                                         vk,
                                         reactivity))))
colnames(ins.vk.small)[1:3] <- c("cell_id", "v_genes", "group")

vk.df <- rbind(FDR.vk.small, T1D.vk.small, ins.vk.small)

#Vlambda
FDR.vl.small <- data.frame(subset
                           (FDR.vl, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(FDR.vl.small)[1:3] <- c("cell_id", "v_genes", "group")

T1D.vl.small <- data.frame(subset
                           (T1D.vl, 
                             select = (c(barcodes,
                                         v_genes,
                                         dx_group))))
colnames(T1D.vl.small)[1:3] <- c("cell_id", "v_genes", "group")

ins.vl.small <- data.frame(subset
                           (ins.df, 
                             select = (c(mAb_line,
                                         vl,
                                         reactivity))))
colnames(ins.vl.small)[1:3] <- c("cell_id", "v_genes", "group")

vl.df <- rbind(FDR.vl.small, T1D.vl.small, ins.vl.small)

rm(FDR.vh.small, FDR.vk.small, FDR.vl.small, 
   T1D.vh.small, T1D.vk.small, T1D.vl.small,
   ins.vh.small, ins.vk.small, ins.vl.small)

###########################
#Calculate proportions for anti-insulin VH, VL, and VK
#Need to ensure NA removed, ensure sum(df$proportion = 100)
df <- ins.df

df <- df %>%
  filter(!is.na(vh)) %>%
  group_by(vh) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$vh))))
sum(df$proportion)
ins.vh.prop <- df
ins.vh.prop$group <- "ins.vh"

df <- ins.df
df <- df %>%
  filter(!is.na(vk)) %>%
  group_by(vk) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$vk))))
sum(df$proportion)
ins.vk.prop <- df
ins.vk.prop$group <- "ins.vk"

df <- ins.df
df <- df %>%
  filter(!is.na(vl)) %>%
  group_by(vl) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$vl))))
sum(df$proportion)
ins.vl.prop <- df
ins.vl.prop$group <- "ins.vl"
##############################
#Calculate proportion for FDR
#vh, vk, and vl
#some day when i have more time i'll figure out
#how to write this as a function to use with lapply

#heavy chain
df <- FDR.vh
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
FDR.vh.prop <- df
FDR.vh.prop$group <- "FDR.vh"

#kappa
df <- FDR.vk
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
FDR.vk.prop <- df
FDR.vk.prop$group <- "FDR.vk"

#lambda
df <- FDR.vl
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
FDR.vl.prop <- df
FDR.vl.prop$group <- "FDR.vl"

########
#T1D
#heavy chain
df <- T1D.vh
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
T1D.vh.prop <- df
T1D.vh.prop$group <- "T1D.vh"

#kappa
df <- T1D.vk
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
T1D.vk.prop <- df
T1D.vk.prop$group <- "T1D.vk"

#lambda
df <- T1D.vl
df <- data.frame(df[, c("v_genes")],
                 stringsAsFactors = TRUE)
colnames(df)[1] <- "v_genes"
df <- df %>%
  filter(!is.na(v_genes)) %>%
  group_by(v_genes) %>%
  summarise(proportion = 100 * n()/(nrow(df) - sum(is.na(df$v_genes))))
sum(df$proportion)
T1D.vl.prop <- df
T1D.vl.prop$group <- "T1D.vl"

#############################
#Make combined data frame with VH proportion for
#each VH gene for FDR, T1D, and insulin BCRs
colnames(FDR.vh.prop)[1] <- "vh"
colnames(FDR.vk.prop)[1] <- "vk"
colnames(FDR.vl.prop)[1] <- "vl"
colnames(T1D.vh.prop)[1] <- "vh"
colnames(T1D.vk.prop)[1] <- "vk"
colnames(T1D.vl.prop)[1] <- "vl"

vh.df.prop <- rbind(FDR.vh.prop, T1D.vh.prop, ins.vh.prop)
vk.df.prop <- rbind(FDR.vk.prop, T1D.vk.prop, ins.vk.prop)
vl.df.prop <- rbind(FDR.vl.prop, T1D.vl.prop, ins.vl.prop)

rm(FDR.vh.prop, FDR.vk.prop, FDR.vl.prop, ins.vh.prop,
  ins.vk.prop, ins.vl.prop, T1D.vh.prop, T1D.vk.prop, 
  T1D.vl.prop)

rm(FDR.vh, FDR.vk, FDR.vl,
   T1D.vh, T1D.vk, T1D.vl,
   ins.df, df)
##############################################
#Make plots with ggplot2
#this plots the average frequency of each VH gene 
#per each subset_dx group
library(ggplot2)

#vh
data <- vh.df.prop
ggplot(data, aes(x=factor(vh), y=proportion, fill=group)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values = c("lightcoral", 
                               "darkorchid1", 
                               "cyan3")) +
  theme(axis.text.x=element_text(size = 8, angle=90, hjust=1))
ggsave("vh_gene.pdf", plot = last_plot(), 
       width = 8,
       height = 4,
       path = dir.out)

#vk
data <- vk.df.prop
ggplot(data, aes(x=factor(vk), y=proportion, fill=group)) +
  geom_bar(stat='identity', position='dodge') + 
  scale_fill_manual(values = c("lightcoral", 
                               "darkorchid1", 
                               "cyan3")) +
  theme(axis.text.x=element_text(size = 8, angle=90, hjust=1))
ggsave("vk_gene.pdf", plot = last_plot(), 
       width = 8,
       height = 4,
       path = dir.out)

#vl
data <- vl.df.prop
ggplot(data, aes(x=factor(vl), y=proportion, fill=group)) +
  geom_bar(stat='identity', position='dodge') + 
  scale_fill_manual(values = c("lightcoral", 
                               "darkorchid1", 
                               "cyan3")) +
  theme(axis.text.x=element_text(size = 8, angle=90, hjust=1))
ggsave("vl_gene.pdf", plot = last_plot(), 
       width = 8,
       height = 4,
       path = dir.out)

##################################
#Chi-square test of independence
#FDR vs T1D
library(stats)
chisq.test.stats <- function(x) {
  data <- x %>% filter(group != "Insulin")
  chi <- chisq.test(table(data$group, data$v_genes))
  return(chi)
}
data.list <-list(vh.df, vk.df, vl.df)

FDRvsT1D.chisq <- lapply(data.list, 
                         FUN =chisq.test.stats)
names(FDRvsT1D.chisq) <- c("vh", "vk", "vl")

#Chi-square test of independence
#FDR vs Ins
chisq.test.stats <- function(x) {
  data <- x %>% filter(group != "T1D")
  chi <- chisq.test(table(data$group, data$v_genes))
  return(chi)
}
data.list <-list(vh.df, vk.df, vl.df)

FDRvsIns.chisq <- lapply(data.list, 
                         FUN =chisq.test.stats)
names(FDRvsIns.chisq) <- c("vh", "vk", "vl")

#Chi-square test of independence
#T1D vs Ins
chisq.test.stats <- function(x) {
  data <- x %>% filter(group != "FDR")
  chi <- chisq.test(table(data$group, data$v_genes))
  return(chi)
}
data.list <-list(vh.df, vk.df, vl.df)

T1DvsIns.chisq <- lapply(data.list, 
                         FUN =chisq.test.stats)
names(T1DvsIns.chisq) <- c("vh", "vk", "vl")

##############
#Summarise Chi square test p values into data frame

#VH
vh.chisq.df<- data.frame(FDRvsIns = FDRvsIns.chisq$vh$p.value)
vh.chisq.df<- cbind(vh.chisq.df, 
                    data.frame(T1DvsIns = T1DvsIns.chisq$vh$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.chisq$vh$p.value))
vh.chisq.df$v_gene <- "vh"

#Vk
vk.chisq.df<- data.frame(FDRvsIns = FDRvsIns.chisq$vk$p.value)
vk.chisq.df<- cbind(vk.chisq.df, 
                    data.frame(T1DvsIns = T1DvsIns.chisq$vk$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.chisq$vk$p.value))
vk.chisq.df$v_gene <- "vk"

#Vl
vl.chisq.df<- data.frame(FDRvsIns = FDRvsIns.chisq$vl$p.value)
vl.chisq.df<- cbind(vl.chisq.df, 
                    data.frame(T1DvsIns = T1DvsIns.chisq$vl$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.chisq$vl$p.value))
vl.chisq.df$v_gene <- "vl"

#single data frame with all v_genes
v.chisq.df <- rbind(vh.chisq.df, vk.chisq.df, vl.chisq.df) 
rm(vh.chisq.df, vk.chisq.df, vl.chisq.df)
rm(FDRvsIns.chisq, FDRvsT1D.chisq, T1DvsIns.chisq)
#########################
#Fisher's exact test of independence

#library(stats)

#FDR vs T1D
fisher.test.stats <- function(x) {
  data <- x %>% filter(group != "Insulin")
  fisher <- fisher.test(table(data$group, data$v_genes),
                        simulate.p.value = TRUE,
                        B = 10000)
  return(fisher)
}
data.list <- list(vh.df, vk.df, vl.df)

FDRvsT1D.fisher <- lapply(data.list, 
                         FUN =fisher.test.stats)
names(FDRvsT1D.fisher) <- c("vh", "vk", "vl")

#FDR vs. Ins
fisher.test.stats <- function(x) {
  data <- x %>% filter(group != "T1D")
  fisher <- fisher.test(table(data$group, data$v_genes),
                        simulate.p.value = TRUE,
                        B = 10000)
  return(fisher)
}
data.list <- list(vh.df, vk.df, vl.df)

FDRvsIns.fisher <- lapply(data.list, 
                          FUN =fisher.test.stats)
names(FDRvsIns.fisher) <- c("vh", "vk", "vl")

#T1D vs. Ins
fisher.test.stats <- function(x) {
  data <- x %>% filter(group != "FDR")
  fisher <- fisher.test(table(data$group, data$v_genes),
                        simulate.p.value = TRUE,
                        B = 10000)
  return(fisher)
}
data.list <- list(vh.df, vk.df, vl.df)

T1DvsIns.fisher <- lapply(data.list, 
                          FUN =fisher.test.stats)
names(T1DvsIns.fisher) <- c("vh", "vk", "vl")

##############
#Summarise Fisher's exact test p values into data frame

#VH
vh.fisher.df<- data.frame(FDRvsIns = FDRvsIns.fisher$vh$p.value)
vh.fisher.df<- cbind(vh.fisher.df, 
                    data.frame(T1DvsIns = T1DvsIns.fisher$vh$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.fisher$vh$p.value))
vh.fisher.df$v_gene <- "vh"

#Vk
vk.fisher.df<- data.frame(FDRvsIns = FDRvsIns.fisher$vk$p.value)
vk.fisher.df<- cbind(vk.fisher.df, 
                    data.frame(T1DvsIns = T1DvsIns.fisher$vk$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.fisher$vk$p.value))
vk.fisher.df$v_gene <- "vk"

#Vl
vl.fisher.df<- data.frame(FDRvsIns = FDRvsIns.fisher$vl$p.value)
vl.fisher.df<- cbind(vl.fisher.df, 
                    data.frame(T1DvsIns = T1DvsIns.fisher$vl$p.value),
                    data.frame(FDRvsT1D = FDRvsT1D.fisher$vl$p.value))
vl.fisher.df$v_gene <- "vl"

#single data frame with all v_genes
v.fisher.df <- rbind(vh.fisher.df, vk.fisher.df, vl.fisher.df) 
rm(vh.fisher.df, vk.fisher.df, vl.fisher.df)
rm(FDRvsT1D.fisher, FDRvsIns.fisher, T1DvsIns.fisher)

#write tables out
write.csv(v.chisq.df, 
            file = paste(dir.out, "v_chisq.csv", sep = ""))
write.csv(v.fisher.df, 
            file = paste(dir.out, "v_fisher.csv", sep = ""))

##################################
#STOP HERE
#STOP HERE
#STOP HERE
#STOP HERE
#Keep this code in case it's useful later

##################################
#If you want to graph in excel, use this code
#note: this data structure does not work with ggplot2
vh_join <- full_join(x=FDR_vh_prop, y=ins_vh_prop,
                     by = c("v_genes" = "vh"))
colnames(vh_join)[1] <- "v_gene"
colnames(vh_join)[2] <- "FDR_all"
colnames(vh_join)[3] <- "ins_BCRs"

#Make combined data frame with Vl proportion for
#each Vl gene for FDR and insulin BCRs
vl_join <- full_join(x=FDR_vl_prop, y=ins_vl_prop,
                     by = c("v_genes" = "vl"))
colnames(vl_join)[1] <- "v_gene"
colnames(vl_join)[2] <- "FDR_all"
colnames(vl_join)[3] <- "ins_BCRs"

#Make combined data frame with Vk proportion for
#each Vl gene for FDR and insulin BCRs
vk_join <- full_join(FDR_vk_prop, ins_vk_prop,
                     by = c("v_genes" = "vk"))
colnames(vk_join)[1] <- "v_gene"
colnames(vk_join)[2] <- "FDR_all"
colnames(vk_join)[3] <- "ins_BCRs"

#write out csv files
vh_name <- paste(dir_out, sep = "", "vh_prop.csv")
write_csv(vh_join, path = vh_name, 
          col_names = TRUE)
vl_name <- paste(dir_out, sep = "", "vl_prop.csv")
write_csv(vl_join, path = vl_name, 
          col_names = TRUE)
vk_name <- paste(dir_out, sep = "", "vk_prop.csv")
write_csv(vh_join, path = vk_name, 
          col_names = TRUE)


dx_group <- c("4025-RB-1_1" ~ "FDR", 
              "4025-RB-1_2" ~ "T1D",
              "4025-RB-1_5" ~ "T1D",
              "4025-RB-2_1" ~ "FDR",
              "4025-RB-2_2" ~ "FDR",
              "4025-RB-2_5" ~ "T1D",
              "4025-RB-3_1" ~ "FDR",
              "4025-RB-3_2" ~ "T1D",
              "4025-RB-3_5" ~ "FDR",
              "4025-RB-4_1" ~ "T1D",
              "4025-RB-4_2" ~ "FDR",
              "4025-RB-4_5" ~ "T1D",
              "4025-RB-5_1" ~ "T1D",
              "4025-RB-5_2" ~ "FDR",
              "4025-RB-5_5" ~ "FDR",
              "4025-RB-6_1" ~ "T1D",
              "4025-RB-6_2" ~ "T1D",
              "4025-RB-6_5" ~ "FDR")

#this will change the valued in dx_group that are 
#file IDs for the disease group they belong to.
data_list <- Map(cbind, data_list, dx_group = filenames)
data_list1 <- lapply(seq_along(data_list), 
                     function(i, x){
                       x[[i]]$dx_group <- disease[i]
                       return (x[[i]])
                     }, data_list
)

# Now the names work - but I pass in datalist instead 
data_list2 <- lapply(names(data_list), 
                     function(n, x){
                       x[[n]]$dx_group <- n
                       return (x[[n]])
                     }, data_list
)
data_list <- data_list1
rm(data_list1, data_list2)
