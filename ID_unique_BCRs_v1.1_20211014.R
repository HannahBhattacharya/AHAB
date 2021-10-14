#Test_Sequences2_Metadata_V3

#This code will identify duplicate CDR3 amino acid sequences between
#two lists of BCR sequences. These BCR sequence lists were created 
#by running nucleotide sequences through IMGT/HighV-QUEST to generate
#output CSV files which contain columns including CDR3 AA sequences.

#revise the companion metadata R script to refer to directories and file paths
#and "Source" (button at top right) this R script prior to running this pipeline

library(dplyr)

#database
db.name <- paste0(dir.in, db.filename, sep = "")
AAseqs <- read.csv(db.name)

#subsets the Sequence ID and CDR3 of sequences in existing database
db1.CDR3 <- subset(AAseqs, 
                   select = c(Sequence.ID, V.GENE.and.allele, CDR3.IMGT))

#differentiates heavy and light chains
#chain <- grep("IGKV-", db1.CDR3$V.GENE.and.allele, value = TRUE)
#chain <- grep(pattern = "Homsap IGKV-", db1.CDR3$V.GENE.and.allele, value = TRUE)

#removes sequences without CDR3s
db2.CDR3 <- db1.CDR3 %>%
  filter(CDR3.IMGT != "") 

#query
query.name <- paste0(dir.in, query.filename, sep = "")
AAseqs <- read.csv(query.name)

#subsets the Sequence ID and CDR3 of sequences in existing database
query1.CDR3 <- subset(AAseqs, 
                   select = c(Sequence.ID, V.GENE.and.allele, CDR3.IMGT))

#removes sequences without CDR3s
query2.CDR3 <- query1.CDR3 %>%
  filter(CDR3.IMGT != "")

#if you want to test the database against itself, use this code
#sets the query sequence as the database to find matches within the existing sequences
#query2.CDR3 <- db2.CDR3

#creates the initial data frame that will display all matches
IgH.df <- data.frame(HYB.line = "", CDR3.AA = "", query.ID = "")

#this loop runs the query sequence against each known sequence to identify matches
i <- 1
which(query2.CDR3$CDR3.IMGT[i] %in% db2.CDR3$CDR3.IMGT)
#runs loop against the number of sequences in query
while (i <= length(query2.CDR3$CDR3.IMGT)){
  #identifies the position numbers of all matches
  bcr.row <- which(db2.CDR3$CDR3.IMGT %in% query2.CDR3$CDR3.IMGT[i])
  #calls out the CDR3 of the match sequences
  CDR3.AA <- db2.CDR3$CDR3.IMGT[bcr.row]
  #calls out the name of the match sequences
  HYB.line <- db2.CDR3$Sequence.ID[bcr.row]
  #inserts the names and CDR3s of the matches in the data frame
  df.2 <- data.frame(HYB.line, CDR3.AA)
  #inserts the name of the query sequence in the data frame
  df.2$query.ID <- query2.CDR3$Sequence.ID[i]
  #binds the initial data frame with the data frame of all matches
  IgH.df <- rbind(IgH.df, df.2)
  #informs the loop to iterate through all sequences in the database
  i <- i+1
}

#exports the dataframe with all matches as a CSV file
write.csv(IgH.df, paste0(dir.out, "victory.csv", sep = ""))

#if you want to generate a toy data set, use this
#define test BCR CDR3
#query2.CDR3 <- subset(AAseqs,
#subset = Sequence.number[] == c(1,6), 
#select = c(Sequence.ID, CDR3.IMGT))
