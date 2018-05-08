rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")

#strains <- c('B', 'C', 'D', 'F', 'J', 'P', 'S')
strains <- c('B')
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
  # run PCoA for gene_by_pop containing G scores for all time points
  gene_by_pop <- c("data/gene_by_sample/", strain, "/sample_by_gene_Gscore.txt")
  df.all <- read.table(paste(gene_by_pop, collapse = ''), sep = "\t", header = TRUE, row.names = 1)
  df.all.no0 <- df.all[rowSums(df.all[,-1]) != 0,]
  df.all.db <- vegdist(df.all.no0, method = "bray", upper = TRUE, diag = TRUE)
  df.all.pcoa <- cmdscale(df.all.db, eig = TRUE, k = 14) 
  explainvar1 <- round(df.all.pcoa$eig[1] / sum(df.all.pcoa$eig), 3) * 100
  explainvar2 <- round(df.all.pcoa$eig[2] / sum(df.all.pcoa$eig), 3) * 100
  explainvar3 <- round(df.all.pcoa$eig[3] / sum(df.all.pcoa$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/pcoa_", strain, ".png"), collapse = ''),
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
  for (x in rownames(df.all.pcoa$points)){
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
  text(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], labels=times, cex= 0.7)
  dev.off()
  
  # run PCoA for gene_by_pop containing G scores for just day 100
  df.all.no0.100 <- df.all.no0[grep('D100', rownames(df.all.no0)), ]
  df.all.db.100 <- vegdist(df.all.no0.100, method = "bray", upper = TRUE, diag = TRUE)
  df.all.pcoa.100 <- cmdscale(df.all.db.100, eig = TRUE, k = 10) 
  explainvar1 <- round(df.all.pcoa.100$eig[1] / sum(df.all.pcoa.100$eig), 3) * 100
  explainvar2 <- round(df.all.pcoa.100$eig[2] / sum(df.all.pcoa.100$eig), 3) * 100
  explainvar3 <- round(df.all.pcoa.100$eig[3] / sum(df.all.pcoa.100$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/pcoa_D100_", strain, ".png"), collapse = ''),
      width = 1200, height = 1200, res = 96*2)
  # Define Plot Parameters
  par(mar = c(5, 5, 1, 2) + 0.1)
  # Initiate Plot
  plot(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
       xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
       ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
       pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)
  
  # Add Axes
  axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  abline(h = 0, v = 0, lty = 3)
  box(lwd = 2)
  
  # Add Points & Labels
  cols.D100 <- c()
  treats.D100 <- c()
  for (x in rownames(df.all.pcoa.100$points)){
    if (grepl("frequency_L0", x)){
      treats.D100 <- c(treats.D100, "1")
      cols.D100 <- c(cols.D100, "#87CEEB")
    } else if ( grepl("frequency_L1", x)) {
      treats.D100 <- c(treats.D100, "10")
      cols.D100 <- c(cols.D100, "#FFA500")
    } else if (grepl("frequency_L2", x)) {
      treats.D100 <- c(treats.D100, "100")
      cols.D100 <- c(cols.D100, "#FF6347")
    }
  }
  

  #plot
  points(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2],
         pch = 19, cex = 3, bg = "gray", col = cols.D100)
  ordiellipse(df.all.pcoa.100, treats.D100, conf = 0.95)
  dev.off()
  
  if(strain == 'S'){
    days <- c('D100')
  } else {
    days <- c('D100', 'D200', 'D300')
  }
  for (day in days){
    to.keep <- c()
    df.all.no0.day <- df.all.no0[grep(day, rownames(df.all.no0)), ]
    treats.day <- c()
    for (x in rownames(df.all.no0.day)){
      if (grepl("frequency_L0", x)){
        treats.day <- c(treats.day, "1")
      } else if ( grepl("frequency_L1", x)) {
        treats.day <- c(treats.day, "10")
      } else if (grepl("frequency_L2", x)) {
        treats.day <- c(treats.day, "100")
      }
    }
    #print(day)
    # Run PERMANOVA with adonis function
    #print(adonis(df.all.no0.day ~ treats.day, method = "bray", permutations = 999))
    # Indicator Value
    #indval <- multipatt(df.all.no0.day, cluster = treats.day, func = "IndVal.g", control = how(nperm=1000)) 
    #summary(indval)
  }
  # calculate distance between vector time-points 
  transfer.times <- c('0', '1', '2')
  reps <- c('1', '2', '3', '4', '5')
  days <- c(100, 200)
  fileName <- paste("data/euclidean_distance/", strain, ".txt", sep = "")
  sink(fileName)
  header <- paste("Line","Time1", "Time2", "Distance", sep = "\t")
  cat(header)
  cat("\n")
  for (transfer.time in transfer.times){
    for (rep in reps){
      for (day in days){
        v1 <- paste("frequency_L", transfer.time, strain, rep, '_D', toString(day), sep = "")
        v2 <- paste("frequency_L", transfer.time, strain, rep, '_D', toString(day+100), sep = "")
        v1.exist <- any(row.names(df.all.pcoa$points) == v1)
        v2.exist <- any(row.names(df.all.pcoa$points) == v2)
        if (v1.exist && v2.exist){
          dist.v1.v2 <- dist(rbind(df.all.pcoa$points[v1,], df.all.pcoa$points[v2,]))[1]
          line.out <- paste(paste("L", transfer.time, strain, rep, sep="") , toString(day), toString(day+100), toString(dist.v1.v2), sep = "\t")
          cat(line.out)
          cat("\n")
        }
      }
    }
  }
  sink()
}
