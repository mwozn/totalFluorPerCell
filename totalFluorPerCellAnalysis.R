# total mitochondrial fluorescence per cell - Michael Wozny 20230320
# Make plots from the outputs of totalFluorPerCell.jim
#

# MODIFY zSteps ACCORDING TO DATASET (MAKE MORE GENERAL/ROBUST LATER)
zSteps <- 19

library(dunn.test)

setwd("~/data/ERMES/revisions3/mdm10_Fluorescence/images/")

# dir names within working dir with subdir mitoMeasurements containing .csv files 
strainNames = c("wdt","qua")

# Search directories for .csv files with total mitochondrial fluorescence measurements
for (i in 1:length(strainNames)){
  pathToMitoCSVFiles <- paste("./", strainNames[i], "/mitoMeasurements", sep = "")
  mitoCSVfiles <- c()
  mitoCSVfiles <<- rbind(mitoCSVfiles,list.files(pathToMitoCSVFiles, pattern = 'resultsMeasure*', full.names = TRUE))
  if (i == 1){
    mitoCSVdata_all <- c()
  }
  for (j in 1:length(mitoCSVfiles[1,])){
    
    if (j == 1){
      cellIdx <- 1 
    }
    
    # Import data
    mitoCSVdata <- read.csv(file=mitoCSVfiles[j], header=TRUE, stringsAsFactors = FALSE)
    
    # Add the date from the .csv input to column 5
    mitoCSVdata[,9] <- substr(mitoCSVfiles[j],3,5)
    colnames(mitoCSVdata)[9] <- "strain"
    
    # Keep track of .csv files for each cell
    mitoCSVdata[,10] <- cellIdx
    colnames(mitoCSVdata)[10] <- "cellIdx"
    cellIdx <- cellIdx+1
    
    # Collect data into mitoCSVdata_all
    mitoCSVdata_all <<- rbind(mitoCSVdata_all,mitoCSVdata)
  }
}

# Search directories for .csv files with cell fluorescence (i.e. background excluding mitochondrial fluorescence) measurements
for (i in 1:length(strainNames)){
  pathToBkgrCSVFiles <- paste("./", strainNames[i], "/mitoMeasurements", sep = "")
  backgroundCSVfiles <- c()
  backgroundCSVfiles <<- rbind(backgroundCSVfiles,list.files(pathToBkgrCSVFiles, pattern = 'resultsBackground*', full.names = TRUE))
  if (i == 1){
    backgroundCSVdata_all <- c()
  }
  for (j in 1:length(backgroundCSVfiles[1,])){
    
    if (j == 1){
      cellIdx <- 1 
    }
    
    # Import data
    backgroundCSVdata <- read.csv(file=backgroundCSVfiles[j], header=TRUE, stringsAsFactors = FALSE)
    
    # Add the date from the .csv input to column 5
    backgroundCSVdata[,9] <- substr(backgroundCSVfiles[j],3,5)
    colnames(backgroundCSVdata)[9] <- "strain"
    
    # Keep track of .csv files for each cell
    backgroundCSVdata[,10] <- cellIdx
    colnames(backgroundCSVdata)[10] <- "cellIdx"
    cellIdx <- cellIdx+1
    
    # Collect data into backgroundCSVdata_all
    backgroundCSVdata_all <<- rbind(backgroundCSVdata_all,backgroundCSVdata)
  }
}

for (i in 1:length(mitoCSVdata_all[,1])){
  # CTCF = Integrated Density – (Area of selected cell X Mean fluorescence of background readings)
  # https://theolb.readthedocs.io/en/latest/imaging/measuring-cell-fluorescence-using-imagej.html
  mitoCSVdata_all[,11] <- mitoCSVdata_all$IntDen - (mitoCSVdata_all$Area * backgroundCSVdata_all$Mean)
}

colnames(mitoCSVdata_all)[11] <- "CorrectedIntDen"

# Sum CorrectedIntDen 
idx <- 0
for (i in seq(from = 1, to = nrow(mitoCSVdata_all), by = zSteps)){
  if (i == 1){
    intDenCorrMitoTotal <- c()
    intDenCorrMitoTotal_all <- c()
  }
  intDenCorrMitoTotal <- sum(mitoCSVdata_all$CorrectedIntDen[i:(i+(zSteps-1))])
  intDenCorrMitoTotal <- cbind(intDenCorrMitoTotal,mitoCSVdata_all$strain[i])
  intDenCorrMitoTotal <- cbind(intDenCorrMitoTotal,mitoCSVdata_all$cellIdx[i])
  colnames(intDenCorrMitoTotal)[1] <- "mitoTotalIntDen"
  colnames(intDenCorrMitoTotal)[2] <- "strain"
  colnames(intDenCorrMitoTotal)[3] <- "cellIdx"
  intDenCorrMitoTotal_all <<- rbind(intDenCorrMitoTotal_all,intDenCorrMitoTotal)
  idx <- idx + 1
}

# Little prince elephant plots from Michah Allen: https://micahallen.org/2018/03/15/introducing-raincloud-plots/

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

library(plyr)
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# Upper and lower SD bounds around the mean
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)
lmad <- function(x) median(x) - mad(x)
umad <- function(x) median(x) + mad(x)

intDenCorrMitoTotal_all <- data.frame(intDenCorrMitoTotal_all)
intDenCorrMitoTotal_all$mitoTotalIntDen <- as.numeric(intDenCorrMitoTotal_all$mitoTotalIntDen)


# Summary stats for mean of all data
sumld <- ddply(intDenCorrMitoTotal_all, ~strain, summarise, mean = mean(mitoTotalIntDen), median = median(mitoTotalIntDen), lower = lb(mitoTotalIntDen), upper = ub(mitoTotalIntDen), n = length(mitoTotalIntDen))

# Summary stats for medians for each experiment
sumExp1 <- ddply(intDenCorrMitoTotal_all, ~strain, summarise, mean = mean(mitoTotalIntDen), median = median(mitoTotalIntDen), lower = lb(mitoTotalIntDen), upper = ub(mitoTotalIntDen), n = length(mitoTotalIntDen))

# Summary stats for Mean of Medians
sumAll <- rbind(sumExp1)#,sumExp2,sumExp3)
sumAllMedians <- ddply(sumAll, ~strain, summarise, meanOfMedians = mean(median), lower = lb(median), upper = ub(median))

# Order data for x-axis
intDenCorrMitoTotal_all$strain <- factor(intDenCorrMitoTotal_all$strain, levels = c("wdt"
                                                                                    #, "tri"
                                                                                    , "qua"))

write.csv(intDenCorrMitoTotal_all , "Path to export the DataFrame\\File Name.csv", row.names=FALSE)

g <- ggplot(data = intDenCorrMitoTotal_all, aes(y = mitoTotalIntDen, x = strain, #fill = strain
)) +
  
  # Half-violin of pooled replicates for each strain
  geom_flat_violin(aes(y = mitoTotalIntDen, 
                       #fill = strain, 
                       color = strain
  ), position = position_nudge(x = .3, y = 0) ) +
  
  # Jitter-points of data points coloured by ExperimentDate
  geom_point(aes(y = mitoTotalIntDen, 
                 color = strain
  ), position = position_jitterdodge(dodge.width = 0.5), size = 1, alpha = 0.5, shape =16) +
  
  # For black outlines of median points
  geom_point(data = sumExp1, aes(x = strain, y = median), color = "black", position = position_nudge(x = -0.17), size = 3.2) +
  
  # Medians and MAD error bars for experimental replicates
  geom_errorbar(data = sumExp1, aes(ymin = lower, ymax = upper, y = median), position = position_nudge(x = -0.17), width = 0) +
  geom_point(data = sumExp1, aes(x = strain, y = median), color = (values = c("#7570b3","#1b9e77")), position = position_nudge(x = -0.17), size = 2.5) +
  
  # Label number of data-points
  #geom_text(data = sumExp1, aes(x = strain, label = n, y = 5.15),position = position_nudge(x = -0.17),vjust = 0, color = (values = c("#1b9e77"))) +
  
  expand_limits(x = 5.25) +
  guides(fill = "none") +
  guides(color = "none") +
  scale_color_manual(values = c("#1b9e77","#7570b3")) +
  scale_fill_brewer(palette = "PRGn") +
  theme_bw() +
  raincloud_theme +
  # Remove x-label 'strain'
  labs(x = "", y = "Total Mitochondrial Fluorescence per Cell") 
 
# Save as 300 dpi .tff
tiff("wdt_qua.tiff", units="in", width=5, height=4.5, res=300)
g
dev.off()

write.table(sumAll, paste("summary_stats_wtd_qua.csv", sep = ""), sep = ",")

library(dunn.test)

# Pairwise comparison of the treatment distributions
table_dunn_Exp1 <- dunn.test(intDenCorrMitoTotal_all$mitoTotalIntDen, intDenCorrMitoTotal_all$strain, method="bonferroni")
write.table(table_dunn_Exp1, paste("table_dunn_wtd_qua.csv", sep = ""), sep = ",")
