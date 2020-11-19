#Run the 3a_vlk... and 3a_vlk... Rscripts first 
#to generate the input files needed here.

#This will generate Vgene usage plots using ggplot2. 
#This will likely require some level of R/ggplot2 experience.

#Things you may need to revise are marked:
######ALL CAPS#########

library(pacman)
# We can then used pacman to load the other packages we need for this project. These packages include,
# * **tidyverse**: this package contains many data editing and plotting tools, including ggplot2
# * **scales**: allows users to manipulate the plot axes more easily
# * **RColorBrewer**: allows you to select from a number of color palettes for plotting within ggplot2
# * **geofacet**: create plots based on a geographical location (will make more sense after some examples)
# * **ggpubr**: easily export ggplot results into PDF, PNG, and other picture formats
# * **ggalt**: some extra goodies to use in ggplot
# * **ggplotAssist**: a plot-&-click RStudio addin that simplifies the use of ggplot
pacman::p_load(tidyverse,ggplot2,scales,RColorBrewer,geofacet,
               ggpubr,ggalt,ggplotAssist)

library(plyr)
library(repr) #to set plot dimensions
options(repr.plot.width=10, repr.plot.height=6)

check.create.dir <- function(the.dir) {
  if (!dir.exists(the.dir)) {
    dir.create(the.dir, recursive = TRUE) }
}

######DEFINE DIR IN########
#this is where files come from
dir.in <- "~/desktop/10X_BCR_pipelines/output/vlk/graphs"
######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vlk/graphs/"

directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

#kappa
file.fam <- paste0(dir.in, "/filtered_combined_vk_fam_freq.csv")
file.gene <- paste0(dir.in,"/filtered_combined_vk_gene_freq.csv")

#this plots the average frequency of each VLK family per each subset.dx group
data <- data.frame(read.csv(file.fam, header = TRUE))
aesthetic <- aes(x=factor(vlk.family), y=frequency, 
                 group=subset.dx, fill=subset.dx)

data.summary <- function(data, varname, groupnames) {
  require(plyr)
  summary.func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data.sum<-ddply(data, groupnames, .fun=summary.func, varname)
  data.sum <- rename(data.sum, c("mean" = varname))
  return(data.sum)
}

df <- data.summary(data, varname="frequency", 
                   groupnames=c("subset.dx", "vlk.family"))

geometry <- stat_summary(fun="mean", geom="bar", position="dodge")
geometry.error <- geom_errorbar(aes(ymin=frequency-sd, 
                                    ymax=frequency+sd), width=.1, 
                                position=position_dodge(.9))

#note: am using df (not data) as the dataframe for the plot now to include error
ggplot(df, aesthetic) + geometry + geometry.error + 
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("vk_family_ave.pdf", plot = last_plot(), path = dir.out)

#--------------------------------------------
#this plots the average frequency of each VLK gene per each subset.dx group
data <- data.frame(read.csv(file.gene, header = TRUE))
aesthetic <- aes(x=factor(vlk.gene), y=frequency, 
                 group=subset.dx, fill=subset.dx)

df <- data.summary(data, varname="frequency", 
                   groupnames=c("subset.dx", "vlk.gene"))

geometry <- stat_summary(fun="mean", geom="bar", position="dodge")
geometry.error <- geom_errorbar(aes(ymin=frequency-sd, 
                                    ymax=frequency+sd), width=.1, 
                                position=position_dodge(.9))

#to plot without error bars:
ggplot(data, aesthetic) + geometry + 
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("vk_gene_ave.pdf", plot = last_plot(), path = dir.out)

#to plot with error bars
#ggplot(df, aesthetic) + geometry + geometry_error + theme(axis.text.x=element_text(angle=90, hjust=1))

#--------------------------------------
#to plot geom boxplot based on cd21low.T1D > cd21low.FDR for Vk and Vl separately
file.enriched.subset <- paste0(dir.out, "/vk_gene_plot_when_cd21low_T1D_increased.csv")
data.enriched.subset <- data.frame(read.csv(file.enriched.subset, header = TRUE))
colnames(data.enriched.subset) [2] <- c("v.gene")
data.enriched.subset.char <- as.character(data.enriched.subset$v.gene)

#make sure to call this again here
data <- data.frame(read.csv(file.gene, header = TRUE))
colnames(data) [3] <- c("v.gene")

library(dplyr)
data <- filter(data, v.gene %in% data.enriched.subset.char)

aesthetic <- aes(x=factor(v.gene), y=frequency, fill=subset.dx)
geometry <- geom_boxplot(position="dodge")

ggplot(data, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))
ggsave("vk_gene_plot_when_cd21low_T1D_increased.pdf", 
       plot = last_plot(), path = dir.out)
