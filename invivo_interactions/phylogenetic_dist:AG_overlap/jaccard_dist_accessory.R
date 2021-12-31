setwd("~/Documents/University/Manuscripts/Brochet_strain-interactions/Datasets/invivo_interactions")

library(ade4)

jaccard_test <- read.csv("jaccard_test.csv")
jaccard_test1 <- read.csv("jaccard_test_no_head.csv", header=FALSE)
jaccard_test1[,1] <- NULL
d <- dist.binary(jaccard_test1, method = 1, diag = FALSE, upper = FALSE) #method 1 is Jaccard index (1901) S3 coefficient of Gower & Legendre
hc <- hclust(d)               # apply hierarchical clustering 
plot(hc, labels=jaccard_test$Strain)    # plot the dendrogram

d <- as.matrix(d)
write.csv(d, "jaccard_results.csv")
