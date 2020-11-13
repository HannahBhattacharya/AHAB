#This script takes 6 columns of CDR data output from
#mAb REDCap database, removes () that define anchor residues,
#adds three flanking residues: AAAsequenceIII 
#to ensure e.g. CDR1 vs. CDR2 alignment occurs 
#independently of length variance, and concatenates 
#all 3 VH CDRs into one string.
#This will be used in hierarchical clustering to
#identify groups of BCRs that have similar CDR aa motifs.

library(tidyverse)
library(magrittr)
library(readr)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

#this is where csv files come from
dir.in <- "~/desktop/10X_BCR_pipelines/CDR"

#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/CDR/out/"

directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

library(dplyr)
library(tidyr)
library(stringr)

file.list <- list.files(dir.in, pattern = "*.csv",
                            full.names = TRUE)
file.list
df <- read.csv(file.list, header = TRUE)

####OPTIONAL: 
#remove delimiters but retain anchor residues:
#remove_delimiters <- function(data){
#  data <- str_remove_all(data, "[(]")
#  data <- str_remove_all(data, "[)]")
#}
#new_df <- data.frame(df[1], lapply(df[2:7], remove_delimiters) )

###OPTIONAL:
#Alternatively, remove anchor residues
#This requires both ends of CDR to be defined with 
# ( )
remove.anchor.residues <- function(data){
  data2 <- sub(".*?)", "", data)
  data <- sub("\\(.*","", data2)
}
df <- data.frame(df[1], lapply(df[2:7], remove.anchor.residues) )

#concatenate all six CDR together into one column
df$CDR.concat <- paste(df$VH.CDR1.sequence,
                       df$VH.CDR2.sequence,
                       df$VH.CDR3.sequence,
                       df$VL.CDR1.sequence,
                       df$VL.CDR2.sequence,
                       df$VL.CDR3.sequence)

data <- data.frame(df$Monoclonal.hybridoma.line.ID,
                          df$CDR.concat)

#remove ( ) and spaces if they remain
test<- str_remove_all(data$df.CDR.concat, "[(]")
test<- str_remove_all(test, "[)]")
test<- str_remove_all(test, "[ ]")
data$df.CDR.concat <- test

colnames(data)[1] <- "name"
colnames(data)[2] <- "seq"

filename<- paste(dir.out,"VH_VL_CDR123_concat",".fasta",sep="")
#function to write out fasta file
write.fasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

write.fasta(data, filename)

###################################################
#Use this to add "fake" anchor residues to help MegaX
#do a better job with alignment and tree calculation

df <- read.csv(file.list, header = TRUE)

# Add AAA and III as anchor residues at beginning and end of
#CDR to help with alignment
#This requires both ends of CDR to be defined with 
# ( )
add.fake.anchor.residues <- function(data){
  data2 <- sub(".*?)", "AAA", data)
  data <- sub("\\(.*","III", data2)
}
df <- data.frame(df[1], lapply(df[2:7], add.fake.anchor.residues) )

#concatenate all six CDR together into one column
df$CDR.concat <- paste(df$VH.CDR1.sequence,
                       df$VH.CDR2.sequence,
                       df$VH.CDR3.sequence,
                       df$VL.CDR1.sequence,
                       df$VL.CDR2.sequence,
                       df$VL.CDR3.sequence)

data <- data.frame(df$Monoclonal.hybridoma.line.ID,
                   df$CDR.concat)

#remove ( ) and spaces if they remain
test<- str_remove_all(data$df.CDR.concat, "[(]")
test<- str_remove_all(test, "[)]")
test<- str_remove_all(test, "[ ]")
data$df.CDR.concat <- test

colnames(data)[1] <- "name"
colnames(data)[2] <- "seq"

filename<- paste(dir.out,"VH_VL_CDR123_concat_add_anchor",".fasta",sep="")

#write out fasta file
write.fasta(data, filename)
