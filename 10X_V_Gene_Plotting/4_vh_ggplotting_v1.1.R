#Run the 3a_vlk... and 3a_vlk... Rscripts first 
#to generate the input files needed here.

#This will generate Vgene usage plots using ggplot2.
#This will likely require some level of R/ggplot2 experience

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
dir.in <- "~/desktop/10X_BCR_pipelines/output/vh/graphs"

######DEFINE DIR OUT########
#this is where output files go
dir.out <- "~/desktop/10X_BCR_pipelines/output/vh/graphs/"
directories <- c(dir.in, dir.out)
lapply(directories,
       FUN = check.create.dir)

file.fam <- paste0(dir.in, "/filtered_combined_vh_fam_freq.csv")
file.gene <- paste0(dir.in,"/filtered_combined_vh_gene_freq.csv")

library(dplyr)
library(tidyr)

#this plots the average frequency of each VH family per each subset.dx group
data <- data.frame(read.csv(file.fam, header = TRUE))
aesthetic <- aes(x=factor(vh.family), y=frequency, 
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
                    groupnames=c("subset.dx", "vh.family"))

geometry <- stat_summary(fun="mean", geom="bar", position="dodge")
geometry.error <- geom_errorbar(aes(ymin=frequency-sd, ymax=frequency+sd), 
                                width=.1, position=position_dodge(.9))

#note: am using df (not data) as the dataframe for the plot now to include error
ggplot(df, aesthetic) + geometry + geometry.error + 
  theme(text = element_text(size=24),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave("vh_family_ave.pdf", plot = last_plot(), path = dir.out)
#--------------------------------------------
#this plots the average frequency of each VH gene per each subset.dx group
data <- data.frame(read.csv(file.gene, header = TRUE))
aesthetic <- aes(x=factor(vh.gene), y=frequency, group=subset.dx, fill=subset.dx)

df <- data.summary(data, varname="frequency", 
                   groupnames=c("subset.dx", "vh.gene"))

geometry <- stat_summary(fun="mean", geom="bar", position="dodge")
geometry.error <- geom_errorbar(aes(ymin=frequency-sd, ymax=frequency+sd), 
                                width=.1, position=position_dodge(.9))

#to plot without error bars:
ggplot(data, aesthetic) + geometry + 
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("vh_gene_ave.pdf", plot = last_plot(), path = dir.out)

#to plot with error bars
ggplot(df, aesthetic) + geometry + geometry.error +
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))

#--------------------------------
#to plot geom boxplot based on cd21low.T1D > cd21low.FDR
data.vh.gene <- data.frame(read.csv(file.gene, header = TRUE))

file.cd21low.T1D <- paste0(dir.in, "/vh_gene_plot_when_cd21low_T1D_increased.csv")
cd21low.T1D.increased.df <- data.frame(read.csv(file.cd21low.T1D,
                                                header = TRUE))
colnames(cd21low.T1D.increased.df) [2] <- c("vh.gene")
cd21low.T1D.increased.char <- as.character(cd21low.T1D.increased.df$vh.gene)
str(cd21low.T1D.increased.char)

library(dplyr)
data <- filter(data.vh.gene, vh.gene %in% cd21low.T1D.increased.char)

aesthetic <- aes(x=factor(vh.gene), y=frequency, fill=subset.dx)
geometry <- geom_boxplot(position="dodge")

ggplot(data, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))
ggsave("vh_gene_increased_T1D_cd21low.pdf", plot = last_plot(), path = dir.out)
#------------------------------------------------
#rest of graphs below probably not so helpful

#this plots the average frequency of each VH gene per each subset.dx group, 
#but I've removed the factor levels that don't show such dramatic 
#differences on the graph above for visual simplicity
data.subset <- data %>%
  subset(vh.gene == "IGHV1-18" | vh.gene == "IGHV1-69" | 
           vh.gene == "IGHV3-11" | vh.gene == "IGHV3-15" | 
           vh.gene == "IGHV3-20" | vh.gene == "IGHV3-23" |
           vh.gene == "IGHV3-30" | vh.gene == "IGHV3-43" | 
           vh.gene == "IGHV3-7" | vh.gene == "IGHV3-72" | 
           vh.gene == "IGHV3-73" | vh.gene == "IGHV4-30-2" |
           vh.gene == "IGHV4-4" | vh.gene == "IGHV4-59" | 
           vh.gene == "IGHV5-10-1" | vh.gene == "IGHV7-4-1")
aesthetic <- aes(x=factor(vh.gene), y=frequency, 
                 group=subset.dx, fill=subset.dx)

df <- data.summary(data.subset, varname="frequency", 
                   groupnames=c("subset.dx", "vh.gene"))

geometry <- stat_summary(fun="mean", geom="bar", position="dodge")
geometry.error <- geom_errorbar(aes(ymin=frequency-sd, ymax=frequency+sd), 
                                width=.1, position=position_dodge(.9))

#to plot bar graph without error bars:
ggplot(df, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))

#to plot bar graph with error bars:
ggplot(df, aesthetic) + geometry + 
  geometry_error + theme(text = element_text(size=24), 
                         axis.text.x=element_text(angle=90, hjust=1))

#--------------------------------------
#to plot dotplot
data.subset.all <- data %>%
  subset(vh.gene == "IGHV1-18" | vh.gene == "IGHV1-69" |
           vh.gene == "IGHV3-11" | vh.gene == "IGHV3-15" | 
           vh.gene == "IGHV3-20" | vh.gene == "IGHV3-23" |
           vh.gene == "IGHV3-30" | vh.gene == "IGHV3-43" | 
           vh.gene == "IGHV3-7" | vh.gene == "IGHV3-72" | 
           vh.gene == "IGHV3-73" | vh.gene == "IGHV4-30-2" |
           vh.gene == "IGHV4-4" | vh.gene == "IGHV4-59" | 
           vh.gene == "IGHV5-10-1" | vh.gene == "IGHV7-4-1")

aesthetic <- aes(x=factor(vh.gene), y=frequency, fill=subset.dx)
geometry <- geom_dotplot(binaxis="y", binwidth=.01, position="dodge")

ggplot(data.subset, aesthetic) + geometry + 
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggplot(data.subset.vh3, aesthetic) + geometry + 
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggplot(data.subset.notvh3, aesthetic) + geometry + 
  theme(axis.text.x=element_text(angle=90, hjust=1))

#--------------------------------------
#to plot geom boxplot
data.subset.all <- data %>%
  subset(vh.gene == "IGHV1-18" | vh.gene == "IGHV1-69" | 
           vh.gene == "IGHV3-11" | vh.gene == "IGHV3-15" | 
           vh.gene == "IGHV3-20" | vh.gene == "IGHV3-23" |
           vh.gene == "IGHV3-30" | vh.gene == "IGHV3-43" | 
           vh.gene == "IGHV3-7" | vh.gene == "IGHV3-72" | 
           vh.gene == "IGHV3-73" | vh.gene == "IGHV4-30-2" |
           vh.gene == "IGHV4-4" | vh.gene == "IGHV4-59" | 
           vh.gene == "IGHV5-10-1" | vh.gene == "IGHV7-4-1")

aesthetic <- aes(x=factor(vh.gene), y=frequency, fill=subset.dx)
geometry <- geom_boxplot(position="dodge")

ggplot(data.subset, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))
ggplot(data.subset.vh3, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))
ggplot(data.subset.notvh3, aesthetic) + geometry + 
  theme(text = element_text(size=24), 
        axis.text.x=element_text(angle=90, hjust=1))