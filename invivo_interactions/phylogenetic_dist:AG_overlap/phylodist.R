setwd("~/Documents/University/Manuscripts/Brochet_strain-interactions/Datasets/invivo_interactions")

library(seqinr)
require(ape)
require(ggtree)
require(ggplot2)

alignment20052010 <- read.alignment(file="CoreGenesAlignment.fasta", format="fasta")

alignment20052010$seq #showing the sequence of the file

alignment20052010$nb #number of sequence

alignment20052010$nam #name of the sequence

printMultipleAlignment <- function(alignment, chunksize=60)
{
  # this function requires the Biostrings package
  require("Biostrings")
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  starts <- seq(1, alignmentlen, by=chunksize)
  n <- length(starts)
  # get the alignment for each of the sequences:
  aln <- vector()
  lettersprinted <- vector()
  for (j in 1:numseqs)
  {
    alignmentj <- alignment$seq[[j]]
    aln[j] <- alignmentj
    lettersprinted[j] <- 0
  }
  # print out the alignment in blocks of 'chunksize' columns:
  for (i in 1:n) { # for each of n chunks
    for (j in 1:numseqs)
    {
      alnj <- aln[j]
      chunkseqjaln <- substring(alnj, starts[i], starts[i]+chunksize-1)
      chunkseqjaln <- toupper(chunkseqjaln)
      # Find out how many gaps there are in chunkseqjaln:
      gapsj <- countPattern("-",chunkseqjaln) # countPattern() is from Biostrings package
      # Calculate how many residues of the first sequence we have printed so far in the alignment:
      lettersprinted[j] <- lettersprinted[j] + chunksize - gapsj
      print(paste(chunkseqjaln,lettersprinted[j]))
    }
    print(paste(' '))
  }
}

printMultipleAlignment(alignment20052010, 60)

virusdist <- dist.alignment(alignment20052010)
mydf<-as.data.frame(as.matrix(virusdist))
write.csv(mydf,file = "filename.csv")

head(virusdist)

rootedNJtree <- function(alignment, theoutgroup, type)
{
  # load the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    myrootedtree <- root(mytree, outgroup, r=TRUE)
    return(myrootedtree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
  # bootstrap the tree
  myboot <- boot.phylo(myrootedtree, mymat, makemytree)
  # plot the tree:
  plot.phylo(myrootedtree, type="p", show.tip.label = TRUE, show.node.label = TRUE, tip.color = "black", font=3)  # plot the rooted phylogenetic tree
  nodelabels(myboot,cex=0.10)          # plot the bootstrap values
  myrootedtree$node.label <- myboot   # make the bootstrap values be the node labels
  return(myrootedtree)
}

firm5DNAtree <- rootedNJtree(alignment20052010, "ESL0186", type="DNA")

#Saving the tree in Newrick type. 
write.tree(firm5DNAtree, "firm5DNAtree.tre")

require("gridExtra")
tree <- read.tree("firm5DNAtree.tre")
#View the tree in normal format
ggtree(tree, ladderize = FALSE, branch.length = "none", layout="circular")+geom_point(aes(shape=isTip, color=isTip), size=3)+geom_text(aes(x=branch, label=label, angle=angle), size=3.5, color="red", vjust=-0.3)+ggtitle("firm5DNAtree_circ")

