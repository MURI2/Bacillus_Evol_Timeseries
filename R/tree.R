rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")


package.list <- c('ape', 'seqinr', 'phylobase', 'adephylo', 'geiger', 'picante', 'stats', 'RColorBrewer', 'caper', 'phylolm', 'pmc', 'ggplot2', 'tidyr', 'dplyr', 'phangorn', 'pander') 
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}


# Read Alignment File {seqinr}
read.aln <- read.alignment(file = "./data/tree/Task2_16S_rRNA.afa", format = "fasta")  
# Convert Alignment File to DNAbin Object {ape}
p.DNAbin <- as.DNAbin(read.aln) 
# Identify Base Pair Region of 16S rRNA Gene to Visuzlize
window <- p.DNAbin[, 100:1500] 
# Command to Visusalize Sequence Alignment {ape}
image.DNAbin(window, cex.lab = 0.50) 
# Optional Code Adds Grid to Help Visualize Rows of Sequences 
#grid(ncol(window), nrow(window), col = "lightgrey") 





ml.bootstrap <- read.tree("./data/tree/RAxML_bipartitions.T20")
par(mar = c(1,1,2,1) + 0.1)
plot.phylo(ml.bootstrap, type = "phylogram", direction = "right", show.tip.label=TRUE,
           use.edge.length = FALSE, cex = 0.6, label.offset = 1, main = "Maximum Likelihood with Support Values")
add.scale.bar(cex = 0.7)
nodelabels(ml.bootstrap$node.label, font = 2, bg = "white", frame = "r", cex = 0.5)
