setwd("~/Documents/University/Manuscripts/Brochet_strain-interactions/Datasets/Invivo_pairwise_transfer")
# Packages --------

library(ggplot2)
#library(dplyr)
library(data.table)
library(patchwork)
library(reshape)
library(rcompanion)
library(magrittr)
library(tidyr)
library(forcats)
library(tidyverse)
library(circular)
library(bpnreg)
library(RColorBrewer)
library(ggpmisc)
#library(plyr)
#library(CircStats)

# Check treatments for contamination ----

data <- read.table("IVPT2_3P.txt", h=T, check.names = FALSE) #table with raw counts results

# Re-shape dataset and calculate proportions

data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
prop <- melt(data, id.vars = "Strain")

start <- 0
end <- (unlist(gregexpr('_', prop$variable))[4] - 1)
prop$variable <- substr(prop$variable, start, end)
prop$variable <- sub("_S1_L$","",prop$variable)
prop$variable <- sub("_S1$","",prop$variable)
prop$variable <- sub("_S1_$","",prop$variable)


# Add replicate column

prop <-  prop %>% 
  extract(variable, into = c("variable", "replicate"), "(.*)_([^_]+)$")
#prop$variable <- sub("[P]", "", prop$variable)
#prop$variable <- as.numeric(prop$variable)
#prop <- prop[order((prop$variable)) , ]
prop$variable <- as.factor(prop$variable)

# Make stacked barplots (check contaminations qualitatively)

colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000")


pdf("stacked_barplots.pdf")

# Iterate over samples (=variable)
prop$variable <- as.factor(prop$variable)

for (variableName in levels(prop$variable)) {
  
  # Get sample specific data
  
  variableData <- prop %>% 
    filter(variable == variableName) %>% 
    mutate(plot_strain= as.factor(1:n()))
  
  p <- ggplot(variableData, aes(x = replicate, y = value, fill = factor(Strain), color = Strain)) + 
    geom_bar(stat = "identity") +
    ggtitle(variableName) +
    xlab("Replicates") + 
    ylab("Percentage") +
    theme_minimal() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    #scale_x_discrete(labels=variableData$strain_replicate) +
    theme(
      legend.position = "",
      #axis.ticks.x = element_blank(),
      #axis.text.x = element_blank()
    )
  
  print(p)
  
}

dev.off()


# Create a dataframe in which you only have the non-focal strains
# to quantify contaminations (check contaminations qualitatively)

prop3 <- data.frame(Strain=prop$Strain, variable=prop$variable, value=prop$value, replicate=prop$replicate)

# Add variable only column

prop3$variable_num <- prop3$variable

prop3$variable <- sub("[P]0_", "", prop3$variable)
prop3$variable <- sub("[P]1_", "", prop3$variable)
prop3$variable <- sub("[P]2_", "", prop3$variable)


pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE) # table containing information about which strain is included in which treatment (pair)

pairs_info$variable <- as.factor(pairs_info$variable)
prop3$variable <- as.factor(prop3$variable)

pairwise_pairs1 <- head(prop3, 0) %>% mutate(col_type = factor(0,levels = c("mono","co")))

for (variableName in levels(prop3$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  
  pairwise_pairs1 <- bind_rows(
    pairwise_pairs1,
    prop3 %>% filter(variable == variableName, !(Strain %in% variableData)) %>%
      mutate(col_type = ifelse(length(variableData) == 1, "mono", "co"))
  )
}

write.csv(pairwise_pairs1, "contamination_quant.csv")

# Sum per sample

pairwise_pairs1_sum <- data.frame(aggregate(pairwise_pairs1$value, by=list(variable=pairwise_pairs1$variable_num, replicate=pairwise_pairs1$replicate), FUN=sum))

write.csv(pairwise_pairs1_sum, "contamination_quant_sum.csv")

# We can use this file to evaluate contamination levels in each sample and 
# set a threshold to exclude contaminated samples. 

# Analysis CFUs data ----

testDataMelt <- read.csv("CFUs.csv")

# Adapt dataframe, transform in log10

start <- 0
end <- (unlist(gregexpr('_', testDataMelt$variable))[2] - 0)
testDataMelt$variable <- substr(testDataMelt$variable, start, end)
testDataMelt <- na.omit(testDataMelt)
testDataMelt$CFUs <- log10(testDataMelt$CFUs)

# Calculate mean and sd across replicates

setDT(testDataMelt)
testDataMelt_agg <- testDataMelt[, .(mean = mean(CFUs), 
                                     se = sd(CFUs)/.N), 
                                 by = .(variable)]

# Plot CFUs data

pdf("CFUs.pdf")
ggplot(testDataMelt_agg, 
       aes(x = variable, y = mean)) + 
  geom_bar(stat = "identity") +
  # Add error bars (here +/- 1.96 SE)
  geom_errorbar(aes(ymax = mean + se, 
                    ymin = mean - se)) +
  xlab("Treatment") + 
  ylab("CFU/mL") +
  coord_cartesian(ylim=c(6,10)) +
  labs(fill = "") + 
  scale_fill_grey() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=3))
dev.off()

# Normalise using CFUs ----

CFUs <- read.csv("CFUs.csv")

data <- read.table("IVPT2_3P.txt", h=T, check.names = FALSE) #table with raw counts results

# Re-shape dataset and calculate proportions

data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
prop <- melt(data, id.vars = "Strain")

start <- 0
end <- (unlist(gregexpr('_', prop$variable))[4] - 1)
prop$variable <- substr(prop$variable, start, end)
prop$variable <- sub("_S1_L$","",prop$variable)
prop$variable <- sub("_S1$","",prop$variable)
prop$variable <- sub("_S1_$","",prop$variable)



# Add CFUs column to dataset

for (i in 1:nrow(prop)){    
  for (j in 1:nrow(CFUs)){ 
    if(prop[i,2]==as.character(CFUs[j,1])){  # if variable match
      prop[i,4]<-CFUs[j,2]
      print(i)
      break
    }
    
  }
}

colnames(prop)[4] <- "CFUs"

# Multiply counts*CFUs to obtain normalised counts
prop$Normalised_counts <- prop$value * prop$CFUs

# Dynamics over time ----
# Create a dataframe (norm_pairs) including only the focal strains for each treatment

prop1 <-  prop %>% 
  extract(variable, into = c("variable", "replicate"), "(.*)_([^_]+)$")


prop1 <-  prop1 %>% 
  extract(variable, into = c("transfer", "variable"), "(.*)_([^_]+)$")

prop1 <- na.omit(prop1)

#prop$variable <- sub("[P]", "", prop$variable)
prop1$variable <- as.numeric(prop1$variable)
prop1 <- prop1[order((prop1$variable)) , ]
prop1$variable <- as.factor(prop1$variable)

pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE)

pairs_info$variable <- as.factor(pairs_info$variable)
prop1$variable <- as.factor(prop1$variable)

norm_pairs <- head(prop1, 0) 

for (variableName in levels(prop1$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  norm_pairs <- bind_rows(
    norm_pairs,
    prop1 %>% filter(variable == variableName, Strain %in% variableData) 
  )
}

write.csv(norm_pairs, "norm_pairs.csv")

norm_pairs$Strain <- as.factor(norm_pairs$Strain)
norm_pairs$variable <- as.factor(norm_pairs$variable)

# Lineplots ----

pdf("lineplots.pdf")

colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000")

norm_pairs$variable <- as.factor(norm_pairs$variable)
norm_pairs$transfer <- as.factor(norm_pairs$transfer)


# Iterate over samples (=variable)

for (variableName in levels(norm_pairs$variable)) {
  
  # Get sample specific data
  
  variableData <- norm_pairs %>% 
    filter(variable == variableName) 
  
  p <- qplot(transfer, Normalised_counts, data=variableData, color = Strain) + 
    scale_y_log10() +
    coord_cartesian(ylim=c(1e+04,1e+10)) +
    #ggtitle(variableName) +
    xlab("Transfer") + 
    ylab("Normalised counts") +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      #axis.ticks.x = element_blank(),
      #axis.text.x = element_blank()
      axis.text=element_text(size=14),
      axis.title=element_text(size=20)
    ) +
    stat_summary(fun.y=mean, geom="line", aes(group = Strain, color = Strain))
  
  print(p)
  
}

dev.off()

# LOD calc ----

# Create LOD dataset

# Dataframe with number of reads of focal strains

data1 <- read.table("IVPT2_3P.txt", h=T, check.names = FALSE) #table with raw counts results
data1 <- melt(data1, id.vars = "Strain")
start <- 0
end <- (unlist(gregexpr('_', data1$variable))[4] - 1)
data1$variable <- substr(data1$variable, start, end)
data1$variable <- sub("_S1_L$","",data1$variable)
data1$variable <- sub("_S1$","",data1$variable)
data1$variable <- sub("_S1_$","",data1$variable)

data1$variable_num <- data1$variable

data1$variable <- as.factor(data1$variable)
data1$variable_num <- data1$variable

data1$variable <- sub("[P]0_", "", data1$variable)
data1$variable <- sub("[P]1_", "", data1$variable)
data1$variable <- sub("[P]2_", "", data1$variable)
data1$variable <- sub("([^_]+)$", "", data1$variable)
data1$variable <- sub("_", "", data1$variable)

pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE)

pairs_info$variable <- as.factor(pairs_info$variable)
data1$variable <- as.factor(data1$variable)

norm_pairs1 <- head(data1, 0) 

for (variableName in levels(data1$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  norm_pairs1 <- bind_rows(
    norm_pairs1,
    data1 %>% filter(variable == variableName, Strain %in% variableData) 
  )
}

# Sum by sample

norm_pairs1_sum <- data.frame(aggregate(norm_pairs1$value, by=list(variable=norm_pairs1$variable_num), FUN=sum))
colnames(norm_pairs1_sum)[2] <- "total_reads"

CFUs <- read.csv("CFUs.csv")
norm_pairs1_sum$total_reads <- as.factor(norm_pairs1_sum$total_reads)


for (i in 1:nrow(norm_pairs1_sum)){    
  for (j in 1:nrow(CFUs)){ 
    if(norm_pairs1_sum[i,1]==as.character(CFUs[j,1])){  # if variable match
      norm_pairs1_sum[i,3]<-CFUs[j,2]
      print(i)
      break
    }
    
  }
}

colnames(norm_pairs1_sum)[3] <- "CFUs"

norm_pairs1_sum <- na.omit(norm_pairs1_sum)
norm_pairs1_sum$total_reads <- as.character(norm_pairs1_sum$total_reads)

norm_pairs1_sum$total_reads <- as.numeric(norm_pairs1_sum$total_reads)

norm_pairs1_sum$LOD <- (norm_pairs1_sum$CFUs / norm_pairs1_sum$total_reads)


# Line plots with LOD

# Add LOD information to lineplots dataset

norm_pairs$variable_name <- paste(norm_pairs$transfer,norm_pairs$variable,norm_pairs$replicate,sep="_")

for (i in 1:nrow(norm_pairs)){    
  for (j in 1:nrow(norm_pairs1_sum)){ 
    if(norm_pairs[i,8]==as.character(norm_pairs1_sum[j,1])){  # if variable match
      norm_pairs[i,9]<-norm_pairs1_sum[j,4]
      print(i)
      break
    }
    
  }
}

colnames(norm_pairs)[9] <- "LOD"

pdf("lineplots_LOD.pdf")

colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000")

norm_pairs$variable <- as.factor(norm_pairs$variable)
norm_pairs$transfer <- as.factor(norm_pairs$transfer)


# Iterate over samples (=variable)

for (variableName in levels(norm_pairs$variable)) {
  
  
  # Get sample specific data
  
  variableData <- norm_pairs %>% 
    filter(variable == variableName) 
  
  library(dplyr)
  
  variableData$transfer <- sub("P", "", variableData$transfer)
  
  variableData$transfer <- as.character(variableData$transfer)
  variableData$transfer <- as.numeric(variableData$transfer)
  
  test <- variableData %>% group_by(transfer) %>% summarise(minLOD = min(LOD), maxLOD = max(LOD))
  y_coordinates <- c(test$maxLOD, rev(test$minLOD))
  x_coordinates <- as.numeric(c(test$transfer, rev(test$transfer)))
  
  positions <- data.frame(x = x_coordinates,
                          y = y_coordinates)
  
  z_coordinates <- c(test$minLOD, rep(0,length(test$minLOD)))

  positions1 <- data.frame(x = x_coordinates,
                           y =z_coordinates)
  

  p1 = ggplot(variableData, aes(transfer, Normalised_counts)) + 
    geom_polygon(aes(x=x,y=y), data = positions, fill='grey', alpha =0.5)+
    geom_polygon(aes(x=x,y=y), data = positions1, fill='grey', alpha =0.8)+
    geom_point(aes(color=Strain))+
    stat_summary(fun=mean, geom="line", aes(group = Strain, color = Strain)) +
    scale_y_log10() +
    coord_cartesian(ylim=c(1e+2,1e+10)) +
    #ggtitle(variableName) +
    xlab("Transfer") + 
    ylab("Normalised counts") +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      #axis.ticks.x = element_blank(),
      #axis.text.x = element_blank()
      axis.text=element_text(size=14),
      axis.title=element_text(size=20)
    ) +
    scale_x_continuous(breaks=c(0,1,2,3),labels=c("P0", "P1", "P2", "P3"))

    print(p1)
}


dev.off()





