setwd("~/Documents/University/Manuscripts/Brochet_strain-interactions/Datasets/invivo_interactions")
# Packages --------

library(ggplot2)
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
library(multcomp)

# Analysis CFUs data ----

testDataMelt <- read.csv("CFUs.csv")

# Adapt dataframe, transform in log10

start <- 0
end <- (unlist(gregexpr('_', testDataMelt$variable))[1] - 0)
testDataMelt$variable <- substr(testDataMelt$variable, start, end)
testDataMelt$variable <- sub("_$","",testDataMelt$variable)
testDataMelt <- na.omit(testDataMelt)
testDataMelt$variable <- sub("[P]", "", testDataMelt$variable)
testDataMelt$variable <- as.numeric(testDataMelt$variable)
testDataMelt <- testDataMelt[order((testDataMelt$variable)) , ]
testDataMelt$variable <- as.factor(testDataMelt$variable)
testDataMelt$CFUs <- log10(testDataMelt$CFUs)

# Calculate mean and sd across replicates

setDT(testDataMelt)
testDataMelt_agg <- testDataMelt[, .(mean = mean(CFUs), 
                  se = sd(CFUs)/.N), 
              by = .(variable)]

# Plot CFUs data

library(multcompView)
library(ggpubr)

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
  theme(axis.text.x = element_text(angle = 45, size=3))
dev.off()


# Check treatments for contamination ----

data <- read.table("invivo_pairwise.txt", h=T, check.names = FALSE) #table with raw counts results

# Re-shape dataset and calculate proportions

data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
prop <- melt(data, id.vars = "Strain")

start <- 0
end <- (unlist(gregexpr('_', prop$variable))[2] - 1)
prop$variable <- substr(prop$variable, start, end)
prop$variable <- sub("_$","",prop$variable)

# Add replicate column

prop <-  prop %>% 
  extract(variable, into = c("variable", "replicate"), "(.*)_([^_]+)$")
prop$variable <- sub("[P]", "", prop$variable)
prop$variable <- as.numeric(prop$variable)
prop <- prop[order((prop$variable)) , ]
prop$variable <- as.factor(prop$variable)

# Make stacked barplots (check contaminations qualitatively - visualize them in stacked barplots)

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
    xlab("Strain") + 
    ylab("Percentage") +
    theme_minimal() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels=variableData$strain_replicate) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = ""
    )
  
  print(p)
  
}

dev.off()

# Create a dataframe in which you only have the non-focal strains
# to quantify contaminations (check contaminations qualitatively)

prop3 <- data.frame(Strain=prop$Strain, variable=prop$variable, value=prop$value, replicate=prop$replicate)

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


pairwise_pairs1_sum <- data.frame(aggregate(pairwise_pairs1$value, by=list(variable=pairwise_pairs1$variable, replicate=pairwise_pairs1$replicate), FUN=sum))


write.csv(pairwise_pairs1_sum, "contamination_quant_sum.csv")


# We can use this file to evaluate contamination levels in each sample and 
# set a threshold to exclude contaminated samples. For now I set this 
# threshold to =<10% contamination to define a sample as contaminated and exclude it.
# Now work on new file for analysis = invivo_pairwise_up.txt

# Normalise using CFUs ----

data <- read.table("invivo_pairwise_up.txt", h=T, check.names = FALSE) # updated dataset without contaminated samples

# Re-shape dataset

data[, -1] <- lapply( data[ , -1], function(x) x/sum(x, na.rm=TRUE) )
prop <- melt(data, id.vars = "Strain")
start <- 0
end <- (unlist(gregexpr('_', prop$variable))[2] - 1)
prop$variable <- substr(prop$variable, start, end)
prop$variable <- sub("_$","",prop$variable)

CFUs <- read.csv("CFUs.csv")

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

# Strains interactions ----
# I will use the method suggested by German = calculate the effect of A on B
# taking the log2 of the normalised counts of A (by replicate) in presence of B
# divided by the mean of the normalised counts of A in mono-colonisation.

# Create a dataframe (norm_pairs) including only the focal strains for each treatment

prop <-  prop %>% 
  extract(variable, into = c("variable", "replicate"), "(.*)_([^_]+)$")
prop$variable <- sub("[P]", "", prop$variable)
prop$variable <- as.numeric(prop$variable)
prop <- prop[order((prop$variable)) , ]
prop$variable <- as.factor(prop$variable)

pairs_info <- read.table("pairs_info.txt", header = TRUE, check.names = FALSE)

pairs_info$variable <- as.factor(pairs_info$variable)
prop$variable <- as.factor(prop$variable)

norm_pairs <- head(prop, 0) %>% mutate(col_type = factor(0,levels = c("mono","co")))

for (variableName in levels(prop$variable)) {
  variableData <- pairs_info %>% filter(variable == variableName) %>% pull('Strain')
  norm_pairs <- bind_rows(
    norm_pairs,
    prop %>% filter(variable == variableName, Strain %in% variableData) %>%
      mutate(col_type = ifelse(length(variableData) == 1, "mono", "co"))
  )
}

write.csv(norm_pairs, "norm_pairs.csv")
norm_pairs <- read.csv("norm_pairs.csv")

norm_pairs$Strain <- as.factor(norm_pairs$Strain)
norm_pairs$variable <- as.factor(norm_pairs$variable)

# Calculate effects of A on B and vice versa

german_df <- data.frame(effect_of = character(0), effect_on = character(0), log2FC = numeric(0), variable = character(0))

for (strain in levels(norm_pairs$Strain)) {
  focal_strain <- norm_pairs %>% filter(Strain == strain)
  co_variables <- focal_strain %>% filter(col_type == "co") %>% pull('variable') %>% unique()
  co_variables1 <- focal_strain %>% filter(col_type =="co")
  mono_variables <- focal_strain %>% filter(col_type == "mono")

  for (i in co_variables) {
    pair_name <- pairs_info %>% filter(variable == i, Strain != strain) %>% pull('Strain')
    data_frame <- bind_rows(mono_variables, focal_strain %>% filter(variable == i))
    
    average_mono <- mean(mono_variables$Normalised_counts)
    
    average_mono <- rep(average_mono, length(data_frame$Normalised_counts[data_frame$col_type=="co"]))
    
    log2FC <- log2((data_frame$Normalised_counts[data_frame$col_type=="co"])/(average_mono))
    
    #log2FC2 <- log2((data_frame$Normalised_counts[data_frame$col_type=="mono"])/(average_mono))
    
    subet <- subset(co_variables1, co_variables1$variable==i)
    
    german_df <- bind_rows(
      german_df,
      data.frame(
        effect_of = pair_name,
        effect_on = rep(mono_variables$Strain,length(data_frame$Normalised_counts[data_frame$col_type=="co"])),
        log2FC = log2FC,
        variable = subet$variable))
      deduped.data <- unique(german_df[ , 1:4 ])
    
    }
}

# Do the same for mono-colonisations

german_df1 <- data.frame(effect_of = character(0), effect_on = character(0), log2FC = numeric(0), variable = character(0))

for (strain in levels(norm_pairs$Strain)) {
  focal_strain <- norm_pairs %>% filter(Strain == strain)
  co_variables <- focal_strain %>% filter(col_type == "co") %>% pull('variable') %>% unique()
  co_variables1 <- focal_strain %>% filter(col_type =="co")
  mono_variables <- focal_strain %>% filter(col_type == "mono")
  
  for (i in co_variables) {
    pair_name <- pairs_info %>% filter(variable == i, Strain != strain) %>% pull('Strain')
    data_frame <- bind_rows(mono_variables, focal_strain %>% filter(variable == i))
    
    average_mono <- mean(mono_variables$Normalised_counts)
    
    average_mono <- rep(average_mono, length(data_frame$Normalised_counts[data_frame$col_type=="mono"]))
    
    #log2FC <- log2((data_frame$Normalised_counts[data_frame$col_type=="co"])/(average_mono))
    
    log2FC <- log2((data_frame$Normalised_counts[data_frame$col_type=="mono"])/(average_mono))
    
    subet <- subset(co_variables1, co_variables1$variable==i)
    
    german_df1 <- bind_rows(
      german_df1,
      data.frame(
        effect_of = mono_variables$Strain,
        effect_on = mono_variables$Strain,
        log2FC = log2FC,
        variable = mono_variables$variable))
    deduped.data1 <- unique(german_df1[ , 1:4 ])
    
  }
}

german_df_new <- rbind(deduped.data, deduped.data1)

write.csv(german_df_new, "strains_interactions1.csv")

# Coordinates ----
# I  transformed this dataframe in excel to have two separate columns
# with the effect of A on B and B on A - symmetric.

cartesian_coord <- read.csv("cartesian_coordinates.csv")

# Filter one symmetry axis

vector = cartesian_coord %>% filter(if_else(A_on_B < 0, B_on_A/A_on_B > 1, B_on_A/A_on_B < 1))

# Plot cartesian coordinates 

pdf("cartesian_coord.pdf")
ggplot(cartesian_coord, aes(A_on_B, B_on_A)) +
  coord_cartesian(xlim = c(-8.5, 3), ylim = c(-8.5, 3)) +
  geom_point() +
  #scale_color_manual(values = colors_code) +
  #scale_fill_manual(values = colors_code) +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_reverse() +
  geom_abline(slope=-1, intercept = 0, lty = 2)
dev.off()

# Convert cartesian coordinates into polar coordinates

cart2polar <- function(x, y) {
  #data.frame(r = sqrt(x^2 + y^2), theta = deg(atan2(y, x)))
  data.frame(r = sqrt(x^2 + y^2), theta = circular(deg(atan2(y, x)),units="degrees") ) # <- I'm transforming "theta" directly to a cricular class
}

cart2pol_res <- cart2polar(vector$A_on_B, vector$B_on_A)

# Transform theta based on line y=-x

cart2pol_res$theta <- (45+(cart2pol_res$theta))

# Plot polar coordinates (all replicates)

pdf("polar_plot1.pdf")
ggplot(cart2pol_res, aes(x=theta, y=r))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=theta, yend=r), size = 0.3) +
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45))
dev.off()

polar_dataframe <- data.frame(A= vector$A, B = vector$B, pair_infor = vector$variable, A_on_B= vector$A_on_B, B_on_A = vector$B_on_A, radius = cart2pol_res$r, angle = cart2pol_res$theta)


# Add same/different species info

for (i in 1:nrow(polar_dataframe)){    
  for (j in 1:nrow(pairs_info)){ 
    if(polar_dataframe[i,3]==as.character(pairs_info[j,1])){  # if variable match
      polar_dataframe[i,8]<-pairs_info[j,4]
      print(i)
      break
    }
    
  }
}

colnames(polar_dataframe)[8] <- "Status"

write.csv(polar_dataframe, "polar_dataframe.csv")

# Create polar plots divided by conspecific/allospecific pairs

subset_Same <- subset(polar_dataframe, polar_dataframe$Status=="SameSDP")

pdf("samespecies_polar.pdf")
ggplot(subset_Same, aes(x=angle, y=radius))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=angle, yend=radius), size = 0.3) +
  ggtitle("Same Species")+
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45))
dev.off()


subset_Different <- subset(polar_dataframe, polar_dataframe$Status=="DifferentSDP")

pdf("different_species_polar.pdf")
ggplot(subset_Different, aes(x=angle, y=radius))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=angle, yend=radius), size = 0.3) +
  ggtitle("Different Species")+
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45))
dev.off()

# Polar statistics analysis ----

# Calculate means for all replicates between categories

polar_dataframe$angle.vector <- as.vector(polar_dataframe$angle)
polar_dataframe$radians <- conversion.circular(polar_dataframe$angle, units = "radians")
polar_dataframe$rads.vector <- as.vector(polar_dataframe$radius)
subsetMI_allrep <- subset(polar_dataframe, angle.vector < -45)
subsetMI_allrep_same <- subset(subsetMI_allrep, subsetMI_allrep$Status=="SameSDP")
subsetMI_allrep_different <- subset(subsetMI_allrep, subsetMI_allrep$Status=="DifferentSDP")

mean(subsetMI_allrep_same$angle.vector)
sd(subsetMI_allrep_same$angle.vector)
mean(subsetMI_allrep_different$angle.vector)
sd(subsetMI_allrep_different$angle.vector)

mean(subsetMI_allrep_same$radius)
mean(subsetMI_allrep_different$radius)

# Calculate the following stats 
# theta  = circular mean (average direction calculated from all replicates)
# rho    = mean resultant length (spread metric in the range 0-1, 0 = spread, 1 = concentrated in one angle; calculated using unit vector lenght =1 )
# kappa  = concentration (if 0=uniform; if large, concentrated around mu; variance (sigma squared) = 1/kappa) 

library(CircStats)
library(dplyr)
library(circular)
polar_dataframe$pair_infor <- as.factor(polar_dataframe$pair_infor)

polar_stats <- polar_dataframe %>% 
  dplyr::group_by(pair_infor) %>%
  dplyr::summarize(theta = mean(angle,na.rm=TRUE), 
                   rho = rho.circular(angle,na.rm=TRUE), 
                   kappa = est.kappa(angle,bias=TRUE), 
                   A_B = mean(A_on_B), 
                   B_A = mean(B_on_A),
  )
polar_stats$pair_infor <- as.factor(polar_stats$pair_infor)

cart2polar <- function(x, y) {
  data.frame(r = sqrt(x^2 + y^2), theta = circular(deg(atan2(y, x)),units="degrees") ) # <- I'm transforming "theta" directly to a circular class
}

tempeh <- cart2polar(polar_stats$A_B, polar_stats$B_A)
tempeh$pair_infor <- polar_stats$pair_infor
tmp1 <- left_join(polar_stats, tempeh, by = "pair_infor")  

library(plyr)
inforpairs <- read.csv(file="pairs_info2.txt", sep = "\t")
pairs_info <- read.csv(file="pairs_info.txt", sep = "\t")

pairs_info$pair_infor <- as.factor(inforpairs$pair_infor)
polar_final <- plyr::join(tmp1,inforpairs,by = "pair_infor") # Add info of strains

# Same plot by category as before but focusing on MI and using mean values...

pdf("polar_plot_means.pdf")
ggplot(polar_final, aes(x=theta.x, y=r))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=theta.x, yend=r), size = 0.3) +
  ggtitle("All pairs")+
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45)) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6, 8))
dev.off()

subsetMI1 <- subset(polar_final, theta.x < -45)
subsetMI1_same <- subset(subsetMI1, subsetMI1$Status=="SameSDP")
subsetMI1_different <- subset(subsetMI1, subsetMI1$Status=="DifferentSDP")

pdf("subsetMI1_same.pdf")
ggplot(subsetMI1_same, aes(x=theta.x, y=r))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=theta.x, yend=r), size = 0.3) +
  ggtitle("Same species")+
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45))+
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6, 8))
dev.off()

pdf("subsetMI1_different.pdf")
ggplot(subsetMI1_different, aes(x=theta.x, y=r))+ 
  coord_polar(start = 135/360*2*pi, theta = "x", direction = 1, clip = "on") +
  geom_segment(aes(y=0, xend=theta.x, yend=r), size = 0.3) +
  ggtitle("Different species")+
  theme_linedraw() +
  scale_x_continuous(limits = c(-180,180), breaks = c(0,45,90,-90,-45))+
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6, 8))
dev.off()


# Define functions ----
# This is the R code in Landler et al 2019 BMC Ecol https://doi.org/10.1186/s12898-019-0246-8
# 2. R code for the self-written functions and an example of applying the tests to a circular data set – (input data in degrees!!!)

#setting function for HR infinity (the originally proposed HR test referred to as HR∞ in the paper)
HermansRassonT <- function(sample){
  n <- length(sample)
  total <- 0
  for (i in 1:n){
    for (j in 1:n){total <- total + abs(sin(sample[i]-sample[j]))}}
  T <- abs((n/pi)-(total/(2*n)))
  return(T)}

HermansRassonPDeg <- function(sample){
  univals <- 9999
  n <- length(sample)
  sample<-circular(sample)
  sample <- rad(sample)
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- matrix(rcircularuniform(n,
                                     control.circular=list(units="radians")))
    testset[f] <- HermansRassonT(data1)}
  Tsample <- HermansRassonT(sample)
  counter <- 1
  for(j in 1:univals){if(testset[j]>Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)
}

#setting function for HR (a new version of the HR test referred to as HR in the paper)
HermansRasson2T <- function(sample){
  n <- length(sample) 
  total <- 0
  for (i in 1:n){
    for (j in 1:n){ total <- total + abs(abs(sample[i]-sample[j])-pi)-
      (pi/2)
    total <- total - (2.895*(abs(sin(sample[i]-sample[j]))-(2/pi)))}}
  T <- total/n
  return(T)}

HermansRasson2PDeg <- function(sample){
  univals <- 9999
  n <- length(sample)
  sample<-circular(sample)
  sample <- rad(sample)
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- matrix(rcircularuniform(n,
                                     control.circular=list(units="radians")))
    testset[f] <- HermansRasson2T(data1)}
  Tsample <- HermansRasson2T(sample)
  counter <- 1
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}

#Setting function for Pycke test
PyckeT <- function(sample){
  n <- length(sample)
  total <- 0
  for (i in 2:n){
    for (j in 1:(i-1)){
      numerator <- cos(sample[i]-sample[j])-sqrt(0.5)
      denominator <- 1.5-(2*sqrt(0.5)*cos(sample[i]-
                                            sample[j]))
      total <- total + (numerator/denominator)}}
  T <- (2*total)/(n-1)
  return(T)}
PyckePDeg <- function(sample){ 
  univals <- 9999
  n <- length(sample)
  sample<-circular(sample)
  sample <- rad(sample)
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- matrix(rcircularuniform(n,
                                     control.circular=list(units="radians")))
    testset[f] <- PyckeT(data1)}
  Tsample <- PyckeT(sample)
  counter <- 1
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}

# Pycke test applied ----

polar_tests <- matrix(ncol=5,nrow=length(polar_final$pair_infor))
colnames(polar_tests) <- c("pair_infor", "rayleigh", "HermansRassonInf", "HermansRasson2PDeg", "Pyke")
polar_final$pair_infor <- as.factor(polar_final$pair_infor)
polar_final$theta.y=NULL

for (var in as.numeric(levels(polar_final$pair_infor) )) {
  #for (var in c(1,2,3,4,5) ) {
  print(var) 
  datas <- polar_dataframe %>% 
    filter(pair_infor == var) 
  ray <- rayleigh.test(datas$angle.vector)
  hr1 <- HermansRassonPDeg(datas$angle.vector)
  hr2 <- HermansRasson2PDeg(datas$angle.vector)
  pyk <- PyckePDeg(datas$angle.vector)
  polar_tests[var,] <- c(var, ray$p.value, hr1, hr2, pyk)
}
polar_tests <- as.tibble(polar_tests[c(1:66),])
polar_tests$pair_infor <- as.factor(polar_tests$pair_infor) 
# <- This table contains the p-values of each test to reject the null hypothesis of uniformity

# Merge all aggregated data into a single datafrmae
tmp <- merge(as.data.frame(polar_final), polar_tests, by.x = "pair_infor")
infotrab <- as_tibble(read.table(file="pairs_info2.txt",header=TRUE))
colnames(infotrab) <- c("pair_infor","Strain","SDP","type")
tmp2 <- infotrab %>% 
  dplyr::select(pair_infor,type) %>%
  distinct()
polar_aggregated <- merge(tmp, tmp2, by.x = "pair_infor")

# Add a factor to indicate whether the direction is singificant
# tests are : "rayleigh" "HermansRassonInf"   "HermansRasson2PDeg" "Pyke"

polar_aggregated <- polar_aggregated %>% mutate(signif = if_else(HermansRasson2PDeg < 0.05, "sig", "not", missing = NULL ))
polar_aggregated %>% mutate(signif = if_else(Pyke < 0.05, "sig", "not", missing = NULL ))

write.csv(polar_aggregated, "polar_aggregated.csv")

# rho vs theta plot

pdf("significance_plot.pdf")
ggplot(polar_aggregated, aes(x=theta.x, y=rho, shape = signif))+
  scale_x_continuous(limits = c(-90,90), breaks = c(0,45,90,135,180,-90,-45,-135))+
  geom_point(alpha = 1/4, size=3)+
  geom_vline(xintercept = c(-90,-45,45,90), linetype="dotted", colour="gray", size=0.3)+
  geom_hline(yintercept = c(0.5,0.75), linetype="dotted", colour="blue", size=0.3)+
  theme(legend.position="bottom", legend.box = "vertical ") +
  labs(y= "rho", x = "theta")
dev.off()

# Correlation with phylogenetic distance ----

phylo_dist = read.csv("phyl_dist_R.csv") #table containing phylogenetic distances per pair of strain, values obtained with script = phylodist.R and also jaccard distances for accessory genome obtained with script = jaccard_dist_accessory.R
polar_dataframe <- read.csv("polar_dataframe.csv")
polar_dataframe$X <- NULL
#polar_dataframe <- data.frame(A= vector$A, B = vector$B, pair_infor = vector$variable, A_on_B= vector$A_on_B, B_on_A = vector$B_on_A, radius = cart2pol_res$r, angle = cart2pol_res$theta)
# Add phylo_dist information to polar dataframe


for (i in 1:nrow(polar_dataframe)){    
  for (j in 1:nrow(phylo_dist)){ 
    if(polar_dataframe[i,1]==as.character(phylo_dist[j,1])&polar_dataframe[i,2]==as.character(phylo_dist[j,2])){  # if variable match
      polar_dataframe[i,8]<-phylo_dist[j,3]
      polar_dataframe[i,9]<-phylo_dist[j,4]
      print(i)
      break
    }
    
  }
}

colnames(polar_dataframe)[8] <- "phylo_dist"
colnames(polar_dataframe)[9] <- "AG_overlap"

# Plot correlations radius/angle - phylo_dist

polar_dataframe$phylo_dist <- as.numeric(polar_dataframe$phylo_dist)
ggplot(polar_dataframe, aes(x = phylo_dist, y = radius)) +
  geom_point() + 
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  ylab("radius") +
  xlab("Core−genes−based phylogenetic distance") +
  theme_bw() +
  theme(legend.position = "none")

polar_dataframe$phylo_dist <- as.numeric(polar_dataframe$phylo_dist)
ggplot(polar_dataframe, aes(x = phylo_dist, y = angle)) +
  geom_point() + 
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  ylab("angle") +
  xlab("Core−genes−based phylogenetic distance") +
  theme_bw() +
  theme(legend.position = "none")

# Plot correlations radius/angle - AG_overlap

polar_dataframe$phylo_dist <- as.numeric(polar_dataframe$phylo_dist)
ggplot(polar_dataframe, aes(x = AG_overlap, y = radius)) +
  geom_point() + 
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  ylab("radius") +
  xlab("AG overlap") +
  theme_bw() +
  theme(legend.position = "none")

polar_dataframe$phylo_dist <- as.numeric(polar_dataframe$phylo_dist)
ggplot(polar_dataframe, aes(x = AG_overlap, y = angle)) +
  geom_point() + 
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  ylab("angle") +
  xlab("AG overlap") +
  theme_bw() +
  theme(legend.position = "none")


# For means

polar_aggregated <- read.csv("polar_aggregated.csv")
polar_aggregated$X <- NULL
for (i in 1:nrow(phylo_dist)){    
  for (j in 1:nrow(polar_aggregated)){ 
    if(phylo_dist[i,5]==as.character(polar_aggregated[j,1])){  # if variable match
      phylo_dist[i,6]<-polar_aggregated[j,2]
      phylo_dist[i,7]<-polar_aggregated[j,7]
      print(i)
      break
    }
    
  }
}
colnames(phylo_dist)[6] <- "theta.mean"
colnames(phylo_dist)[7] <- "radius.mean"

pairs_info <- read.table("pairs_info.txt", header=TRUE, check.names = FALSE)

for (i in 1:nrow(phylo_dist)){    
  for (j in 1:nrow(pairs_info)){ 
    if(phylo_dist[i,5]==(pairs_info[j,1])){  # if variable match
      phylo_dist[i,8]<-pairs_info[j,4]
      print(i)
      break
    }
  }
}
colnames(phylo_dist)[8] <- "Status"

# angle/radius vs phylogenetic distance

polar_dataframe$phylo_dist <- as.factor(polar_dataframe$phylo_dist)
phylodist_subset <- subset(phylo_dist, theta.mean < -45)
ggplot(phylodist_subset, aes(x = distance, y = theta.mean)) +
  geom_point(aes(colour=Status)) + 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  #scale_colour_manual(values = mycolors) +
  ylab("Angle-mean") +
  xlab("Core−genes−based phylogenetic distance") +
  theme_bw() +
  theme(legend.position = "none") #+
#scale_y_continuous(limits = c(-90,-45)) 
dev.off()

polar_dataframe$pair_infor <- as.factor(polar_dataframe$phylo_dist)
phylodist_subset <- subset(phylo_dist, theta.mean < -45)
ggplot(phylodist_subset, aes(x = distance, y = radius.mean)) +
  geom_point(aes(colour=Status)) + 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  #scale_colour_manual(values = mycolors) +
  ylab("Radius-mean") +
  xlab("Core−genes−based phylogenetic distance") +
  theme_bw() +
  theme(legend.position = "none") #+
  #scale_y_continuous(limits = c(-90,-45)) 

# angle/radius vs AG_overlap

polar_dataframe$phylo_dist <- as.factor(polar_dataframe$phylo_dist)
phylodist_subset <- subset(phylo_dist, theta.mean < -45)
ggplot(phylodist_subset, aes(x = AG_jaccard, y = theta.mean)) +
  geom_point(aes(colour=Status)) + 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  #scale_colour_manual(values = mycolors) +
  ylab("Angle-mean") +
  xlab("AG_overlap") +
  theme_bw() +
  theme(legend.position = "none") #+
#scale_y_continuous(limits = c(-90,-45)) 
dev.off()

polar_dataframe$pair_infor <- as.factor(polar_dataframe$phylo_dist)
phylodist_subset <- subset(phylo_dist, theta.mean < -45)
ggplot(phylodist_subset, aes(x = AG_jaccard, y = radius.mean)) +
  geom_point(aes(colour=Status)) + 
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  stat_smooth(method = lm) +
  stat_poly_eq(formula = x~y,
               aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
               parse = TRUE) +
  #scale_colour_manual(values = mycolors) +
  ylab("Radius-mean") +
  xlab("AG_overlap") +
  theme_bw() +
  theme(legend.position = "none") #+
#scale_y_continuous(limits = c(-90,-45)) 


# Calculate distance from 50-50 ----

calc50 <- (0.5-norm_pairs$value)

calc50_dataframe <- data.frame(pair=norm_pairs$variable, value = norm_pairs$value, calc50 = calc50)

toDelete <- seq(1, nrow(calc50_dataframe), 2)
calc50_dataframe <- calc50_dataframe[ toDelete ,]

# Add status
for (i in 1:nrow(calc50_dataframe)){    
  for (j in 1:nrow(pairs_info)){ 
    if(calc50_dataframe[i,1]==(pairs_info[j,1])){  # if variable match
      calc50_dataframe[i,4]<-pairs_info[j,4]
      print(i)
      break
    }
  }
}

colnames(calc50_dataframe)[4] <- "Status"

calc50_dataframe <- calc50_dataframe[calc50_dataframe$Status != "Mono",]

pdf("hist_50-50.pdf")

ggplot(calc50_dataframe, aes(x=calc50, fill = Status)) + 
  geom_histogram(aes(y=..density..), binwidth=0.1, position="identity", alpha = 0.5) +
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_density(alpha=.1) +
  theme_bw() +
  xlab("Distance to 50%") +
  ylab("Density") +
  theme(legend.position = "none")

dev.off()

subset_same50 <- subset(calc50_dataframe, calc50_dataframe$Status=="SameSDP")
subset_different50 <- subset(calc50_dataframe, calc50_dataframe$Status=="DifferentSDP")

ggplot(subset_same50, aes(x=calc50)) + 
  geom_histogram(aes(y=..density..), binwidth=0.1, position="identity", alpha = 0.5) +
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_density(alpha=.1) +
  theme_bw() +
  xlab("Distance to 50%") +
  ylab("Density") +
  theme(legend.position = "none")

ggplot(subset_different50, aes(x=calc50)) + 
  geom_histogram(aes(y=..density..), binwidth=0.1, position="identity", alpha = 0.5) +
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  geom_density(alpha=.1) +
  theme_bw() +
  xlab("Distance to 50%") +
  ylab("Density") +
  theme(legend.position = "none")

shapiro.test(subset_same50$calc50)
shapiro.test(subset_different50$calc50)



  
# Calculate limit of detection ----

data <- read.table("invivo_pairwise.txt", h=T, check.names = FALSE) # updated dataset without contaminated samples

lod <- melt(data, id.vars = "Strain")
start <- 0
end <- (unlist(gregexpr('_', lod$variable))[2] - 1)
lod$variable <- substr(lod$variable, start, end)
lod$variable <- sub("_$","",lod$variable)

lod_sum <- data.frame(aggregate(lod$value, by=list(variable=lod$variable), FUN=sum))
colnames(lod_sum)[2] <- "total_reads"

# add CFUs

CFUs <- read.csv("CFUs.csv")

for (i in 1:nrow(lod_sum)){    
  for (j in 1:nrow(CFUs)){ 
    if(lod_sum[i,1]==as.character(CFUs[j,1])){  # if variable match
      lod_sum[i,3]<-CFUs[j,2]
      print(i)
      break
    }
    
  }
}

colnames(lod_sum)[3] <- "CFUs"

lod_sum <- na.omit(lod_sum)
lod_sum$total_reads <- as.character(lod_sum$total_reads)

lod_sum$total_reads <- as.numeric(lod_sum$total_reads)

lod_sum$LOD <- (lod_sum$CFUs / lod_sum$total_reads)

write.csv(lod_sum, "LOD.csv")