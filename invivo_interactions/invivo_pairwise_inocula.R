setwd("~/Documents/University/Manuscripts/Brochet_strain-interactions/Datasets/invivo_interactions/inocula_analysis")
# Packages --------

library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)
library(reshape)
library(rcompanion)
library(magrittr)
library(tidyr)
library(forcats)
library(tidyverse)



# Examination CFUs data ----

testDataMelt <- read.csv("CFUs_inocula.csv")
start <- 0
end <- (unlist(gregexpr('_', testDataMelt$variable))[1] - 0)

testDataMelt$variable <- as.factor(testDataMelt$variable)

testDataMelt <- testDataMelt[order((testDataMelt$variable)) , ]


testDataMelt$CFUs <- log10(testDataMelt$CFUs)


pdf("CFUs_inocula.pdf")
ggplot(testDataMelt, 
       aes(x = variable, y = CFUs)) + 
  geom_bar(stat = "identity") +
  xlab("Treatment") + 
  ylab("CFU/mL") +
  coord_cartesian(ylim=c(6,10.5)) +
  labs(fill = "") + 
  scale_fill_grey() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=3))
dev.off()



# Investigate contamination qualitatively ----

data <- read.table("invivo_pairwise_inocula.txt", h=T, check.names = FALSE)

# Re-shape dataset and calculate proportions
data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
prop <- melt(data, id.vars = "Strain")

# Change samples names
start <- 0
end <- (unlist(gregexpr('_', prop$variable))[1] - 1)
prop$variable <- substr(prop$variable, start, end)
prop$variable <- sub("_$","",prop$variable)
prop$variable <- sub("inoc", "", prop$variable)

prop$variable <- factor(prop$variable, 1:78)

prop <- prop[order((prop$variable)) , ]

subset <- subset(prop, prop$variable=="67" | prop$variable== "68" | prop$variable=="69" | prop$variable=="70" | prop$variable=="71" | prop$variable=="72" | prop$variable=="73" | prop$variable=="74" | prop$variable=="75" | prop$variable=="76" | prop$variable=="77")

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

pdf("stacked_barplots_inocula1.pdf")

ggplot(prop, aes(x = variable, y = value, fill = factor(Strain), color = Strain)) + 
  geom_bar(stat = "identity", width = 1) +
  #ggtitle(variableName) +
  xlab("Strain") + 
  ylab("Percentage") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #scale_x_discrete(labels=variableData$strain_replicate) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    legend.position = ""
  )

dev.off()


# Create a dataframe in which you only have the non-focal strains
# to better quantify contaminations
prop <- na.omit(prop)

prop3 <- data.frame(Strain=prop$Strain, variable=prop$variable, value=prop$value)

pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE)

pairs_info$variable <- as.factor(pairs_info$variable)
prop3$variable <- as.factor(prop3$variable)

pairwise_pairs1 <- head(prop3, 0) %>% mutate(col_type = factor(0,levels = c("mono","co")))

for (variableName in levels(prop$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  
  pairwise_pairs1 <- bind_rows(
    pairwise_pairs1,
    prop %>% filter(variable == variableName, !(Strain %in% variableData)) %>%
      mutate(col_type = ifelse(length(variableData) == 1, "mono", "co"))
  )
}

write.csv(pairwise_pairs1, "contamination_quant_inocula.csv")

# Create dataframe with only focal strains
pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE)

pairs_info$variable <- as.factor(pairs_info$variable)
prop$variable <- as.factor(prop$variable)

norm_pairs <- head(prop, 0) %>% mutate(col_type = factor(0,levels = c("mono","co")))
#counts_pairs <- head(data1, 0) %>% mutate(col_type = factor(0,levels = c("mono","co")))

for (variableName in levels(prop$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  norm_pairs <- bind_rows(
    norm_pairs,
    prop %>% filter(variable == variableName, Strain %in% variableData) %>%
      mutate(col_type = ifelse(length(variableData) == 1, "mono", "co"))
  )
}

write.csv(norm_pairs,"dataframe_inocula.csv")

# Correlate inocula proportions with pairs proportions (after 10 days)

norm_pairs_10days <- read.csv("norm_pairs.csv") #file from 10days experiment

norm_pairs_10days$X <- NULL

# add info inocula

for (i in 1:nrow(norm_pairs_10days)){    
  for (j in 1:nrow(norm_pairs)){ 
    if(norm_pairs_10days[i,1]==as.character(norm_pairs[j,1])&(norm_pairs_10days[i,2]==as.character(norm_pairs[j,2]))){  # if variable match
      norm_pairs_10days[i,8]<-norm_pairs[j,3]
      print(i)
      break
    }
    
  }
}

colnames(norm_pairs_10days)[8] <- "inocula_prop"

norm_pairs_10days <- norm_pairs_10days[norm_pairs_10days$col_type != "mono", ]

test <- norm_pairs_10days[cumsum(rle(norm_pairs_10days$replicate)$lengths),]

#write.csv(norm_pairs_10days,"norm_pairs_10days.csv")

pdf("prop_corr-inocula-end1.pdf")
library(ggpmisc)
ggplot(test, aes(x = inocula_prop, y = value)) +
  geom_point() +
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  ylab("Proportions inocula") +
  xlab("Proportions after 10 days") +
  theme_bw()
dev.off()



