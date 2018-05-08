rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")

#strains <- c('B', 'C', 'D', 'F', 'J', 'P', 'S')
strains <- c('P')
strains.vec <- vector(mode="list", length=6)
names(strains.vec) <- c('B', 'C', 'D', 'F', 'J', 'P')
strains.vec[[1]] <- 'Bacillus'; strains.vec[[2]] <- 'Caulobacter'
strains.vec[[3]] <- 'Deinococcus'; strains.vec[[4]] <- 'Pedobacter'; 
strains.vec[[5]] <- 'Janthinobacterium';strains.vec[[6]] <- 'Pseudomonas'
package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR', 'indicspecies')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}

for (strain in strains) {
  sdata <- c("data/gene_by_sample_long/gene_by_sample_long_", strain, ".txt")
  classes <- as.character(c("numeric", "character", "numeric", "numeric", "character", 
               "character", "numeric", "numeric", "character", "numeric"))
  df <- read.table(paste(sdata, collapse = ''), sep = "\t", header = TRUE, colClasses = classes)
  df_sub <- select(df, LocusTag, Line, Muts)
  gene_by_pop <- spread(df_sub, Line, Muts)
  pop_by_gene <- t(subset(gene_by_pop, select = -c(LocusTag) ))
  pop_by_gene.no0 <- pop_by_gene[rowSums(pop_by_gene[,-1]) != 0,]
  pop_by_gene.db <- vegdist(pop_by_gene.no0, method = "bray", upper = TRUE, diag = TRUE)
  fish.pcoa <- cmdscale(pop_by_gene.db, eig = TRUE, k = 3) 
  
  explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
  explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
  explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/D100_CDS_", strain, ".png"), collapse = ''),
      width = 1200, height = 1200, res = 96*2)
  # Define Plot Parameters
  par(mar = c(5, 5, 1, 2) + 0.1)
  # Initiate Plot
  plot(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
       xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
       ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
       pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)
  
  # Add Axes
  axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  abline(h = 0, v = 0, lty = 3)
  box(lwd = 2)
  
  # Add Points & Labels
  cols <- c()
  treats <- c()
  for (x in rownames(fish.pcoa$points)){
    if (grepl("frequency_L0", x)){
      treats <- c(treats, "1")
      cols <- c(cols, "#87CEEB")
    } else if ( grepl("frequency_L1", x)) {
      treats <- c(treats, "10")
      cols <- c(cols, "#FFA500")
    } else if (grepl("frequency_L2", x)) {
      treats <- c(treats, "100")
      cols <- c(cols, "#FF6347")
    }
  }
  #plot
  points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
         pch = 19, cex = 3, bg = "gray", col = cols)
  ordiellipse(fish.pcoa, treats, conf = 0.95)
  dev.off()
  
  # run PCoA for gene_by_pop that includes intergenic regions
  sdata.all <- c("data/gene_by_sample/", strain, "/sample_by_gene.txt")
  df.all <- read.table(paste(sdata.all, collapse = ''), sep = "\t", header = TRUE, row.names = 1)
  df.all.no0 <- df.all[rowSums(df.all[,-1]) != 0,]
  df.all.db <- vegdist(df.all.no0, method = "bray", upper = TRUE, diag = TRUE)
  df.all.pcoa <- cmdscale(df.all.db, eig = TRUE, k = 3) 
  explainvar1 <- round(df.all.pcoa$eig[1] / sum(df.all.pcoa$eig), 3) * 100
  explainvar2 <- round(df.all.pcoa$eig[2] / sum(df.all.pcoa$eig), 3) * 100
  explainvar3 <- round(df.all.pcoa$eig[3] / sum(df.all.pcoa$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/D100_all_", strain, ".png"), collapse = ''),
      width = 1200, height = 1200, res = 96*2)
  # Define Plot Parameters
  par(mar = c(5, 5, 1, 2) + 0.1)
  # Initiate Plot
  plot(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
       xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
       ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
       pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)
  
  # Add Axes
  axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  abline(h = 0, v = 0, lty = 3)
  box(lwd = 2)
  
  # Add Points & Labels
  cols <- c()
  treats <- c()
  for (x in rownames(df.all.pcoa$points)){
    if (grepl("frequency_L0", x)){
      treats <- c(treats, "1")
      cols <- c(cols, "#87CEEB")
    } else if ( grepl("frequency_L1", x)) {
      treats <- c(treats, "10")
      cols <- c(cols, "#FFA500")
    } else if (grepl("frequency_L2", x)) {
      treats <- c(treats, "100")
      cols <- c(cols, "#FF6347")
    }
  }
  
  times <- c()
  for (x in rownames(fish.pcoa$points)){
    if (grepl("D100", x)){
      times <- c(times, "100")
    } else if ( grepl("D200", x)) {
      times <- c(times, "200")
    } else if (grepl("D300", x)) {
      times <- c(times, "300")
    }
  }
  #plot
  points(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2],
         pch = 19, cex = 3, bg = "gray", col = cols)
  ordiellipse(df.all.pcoa, treats, conf = 0.95)
  dev.off()
 
  # run PCoA for gene_by_pop that includes intergenic regions and filters by genes with multiple hits

  sdata.all <- c("data/gene_by_sample/", strain, "/sample_by_gene_multiple_hits.txt")
  df.all <- read.table(paste(sdata.all, collapse = ''), sep = "\t", header = TRUE, row.names = 1)
  df.all.no0 <- df.all[rowSums(df.all[,-1]) != 0,]
  df.all.db <- vegdist(df.all.no0, method = "bray", upper = TRUE, diag = TRUE)
  df.all.pcoa <- cmdscale(df.all.db, eig = TRUE, k = 3) 
  print(vegdist(df.all.pcoa$points, method = "euclidean", upper = TRUE, diag = TRUE))
  explainvar1 <- round(df.all.pcoa$eig[1] / sum(df.all.pcoa$eig), 3) * 100
  explainvar2 <- round(df.all.pcoa$eig[2] / sum(df.all.pcoa$eig), 3) * 100
  explainvar3 <- round(df.all.pcoa$eig[3] / sum(df.all.pcoa$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/D100_all_mult_", strain, ".png"), collapse = ''),
      width = 1200, height = 1200, res = 96*2)
  # Add Points & Labels
  cols <- c()
  treats <- c()
  for (x in rownames(df.all.pcoa$points)){
    if (grepl("frequency_L0", x)){
      treats <- c(treats, "1")
      cols <- c(cols, "#87CEEB")
    } else if ( grepl("frequency_L1", x)) {
      treats <- c(treats, "10")
      cols <- c(cols, "#FFA500")
    } else if (grepl("frequency_L2", x)) {
      treats <- c(treats, "100")
      cols <- c(cols, "#FF6347")
    }
  }
  
  times <- c()
  for (x in rownames(df.all.pcoa$points)){
    if (grepl("D100", x)){
      times <- c(times, "100")
    } else if ( grepl("D200", x)) {
      times <- c(times, "200")
    } else if (grepl("D300", x)) {
      times <- c(times, "300")
    }
  }
  
  # Define Plot Parameters
  par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
  #par(mar = c(5, 5, 1, 2) + 0.1)
  # Initiate Plot
  strain.string <- toString(strains.vec[strain])
  plot(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], 
       main=bquote(~bolditalic(.(strain.string))),
       xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
       xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
       ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
       pch = 16, cex = 2.5, type = "n", cex.lab = 1.8, 
       cex.axis = 1.2, axes = FALSE, cex.main=1.75)
  
 
  # Add Axes
  axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  abline(h = 0, v = 0, lty = 3)
  box(lwd = 2)
  
  #plot
  points(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2],
         pch = 19, cex = 3, bg = "gray", col = cols)
  ordiellipse(df.all.pcoa, treats, conf = 0.95)
  text(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], labels=times, cex= 0.7)
  #op <- par(cex = 1.2)
  legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
         col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
         pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
  dev.off()
  print(strain)
  # First we calculate the relative abundances of each species at each site
  df.all.no0.REL <- df.all.no0
  for(i in 1:nrow(df.all.no0)){
    df.all.no0.REL[i, ] = df.all.no0[i, ] / sum(df.all.no0[i, ])
  } 
  #gene.corr <- add.spec.scores(df.all.pcoa, df.all.no0.REL, method = "cor.scores")$cproj
  #corrcut  <- 0.7       # user defined cutoff
  #imp.spp  <- gene.corr[abs(gene.corr[, 1]) >= corrcut | abs(gene.corr[, 2]) >= corrcut, ]
  #print(imp.spp)
  # Permutation Test for Species Abundances Across Axes
  #fit <- envfit(df.all.pcoa, df.all.no0.REL, perm = 999)
  #fit.slice <- c(which(fit$vectors$pvals < 0.05))
  #print(fit$vectors)
  #print(fit$vectors$r[fit.slice])
  
  # Run PERMANOVA with adonis function
  #print(adonis(df.all.no0 ~ treats, method = "bray", permutations = 999))
  # Indicator Value
  #print(rownames(df.all.no0))
  #print(treats)
  #indval <- multipatt(df.all.no0, cluster = treats, func = "IndVal.g", control = how(nperm=1000)) 
  #summary(indval)
  
  # environment preference of mutations in each gene
  #df.all.no0.rel <- decostand(df.all.no0, method = "total")
  #phi <- multipatt(df.all.no0.rel, cluster = treats, func = "r.g", control = how(nperm=1000)) 
  #summary(phi)
}



