rm(list = ls())
getwd()
setwd("~/GitHub/Bacillus_Evol_Timeseries")


package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR', 'indicspecies', 
                  'viridis', 'png', 'grid', 'scales', 'data.table')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}


######
# run PCoA for gene_by_pop containing G scores for all time points
######
gene_by_pop <- c("data/pool_pop_seq/gene_by_pop_delta.txt")
df.all <- read.table(paste(gene_by_pop, collapse = ''), sep = "\t", header = TRUE, row.names = 1)
df.all.no0 <- df.all[rowSums(df.all[,-1]) != 0,]
df.all.db <- vegdist(df.all.no0, method = "bray", upper = TRUE, diag = TRUE)
df.all.pcoa <- cmdscale(df.all.db, eig = TRUE, k = 14) 
explainvar1 <- round(df.all.pcoa$eig[1] / sum(df.all.pcoa$eig), 3) * 100
explainvar2 <- round(df.all.pcoa$eig[2] / sum(df.all.pcoa$eig), 3) * 100
explainvar3 <- round(df.all.pcoa$eig[3] / sum(df.all.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)
# Plot Eigenvalues
plot(df.all.pcoa$eig, xlab = "PCoA Axis", ylab = "Eigenvalue", 
     las = 1, cex.lab = 1.5, pch = 16)

# Add Expectation based on Kaiser-Guttman criterion and Broken Stick Model
abline(h = mean(df.all.pcoa$eig), lty = 2, lwd = 2, col = "blue")
b.stick <- bstick(29, sum(df.all.pcoa$eig))
lines(1:29, b.stick, type = "l", lty = 4, lwd = 2, col = "red")

# Add Legend
legend("topright", legend = c("Avg Eigenvalue", "Broken-Stick"), 
       lty = c(2, 4), bty = "n", col = c("blue", "red"))



shape <- c()
for (x in rownames(df.all.pcoa$points)){
  if (grepl("B", x)){
    shape <- c(shape, 16)
  } else if ( grepl("S", x)) {
    shape <- c(shape, 1)
  } 
}

# Add Points & Labels
cols <- c()
treats <- c()
for (x in rownames(df.all.pcoa$points)){
  if (grepl("L0", x)){
    treats <- c(treats, "1")
    cols <- c(cols, "#87CEEB")
  } else if ( grepl("L1", x)) {
    treats <- c(treats, "10")
    cols <- c(cols, "#FFA500")
  } else if (grepl("L2", x)) {
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

png(filename = paste(c("figs/pcoa_B_S.png"), collapse = ''),
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
#plot
points(0, 0, pch = 16, cex = 3, bg = 'gray', col = 'gray', lwd = 4)
points(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2],
       pch = shape, cex = 3, bg = "gray", col = alpha(cols, 0.8), lwd =4)
ellipse.cols <- c("#87CEEB", "#FFA500", "#FF6347", "#87CEEB", "#FFA500", "#FF6347")
ordiellipse(df.all.pcoa, paste0(shape, treats), col = ellipse.cols, lwd = 2, conf = 0.95)
text(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], labels=times, cex= 0.7)
dev.off()



######
# run PCA for gene_by_pop containing G scores for all time points
######
df.all.no0.H <- decostand(df.all.no0, method = 'hellinger')
###### worry about this later.....

#########
#PERMANOVA
#########
B_S.category <- c()
for (x in rownames(df.all.no0)){
  if (grepl("B", x)){
    B_S.category <- c(B_S.category, 'B')
  } else if ( grepl("S", x)) {
    B_S.category <- c(B_S.category, 'S')
  } 
}

treats.perm <- c()
for (x in rownames(df.all.no0)){
  if (grepl("L0", x)){
    treats.perm <- c(treats.perm, 'L0')
  } else if ( grepl("L1", x)) {
    treats.perm <- c(treats.perm, 'L1')
  } else if (grepl("L2", x)) {
    treats.perm <- c(treats.perm, 'L2')
  }
}

perm <- adonis(df.all.no0 ~ B_S.category * treats.perm, method = "bray", permutations = 999)

interact.beta <- perm$coef.sites[4,]
dim(t(perm$coef.sites) *  perm$model.matrix)
t(perm$coef.sites) *  perm$model.matrix
perm$coef.sites * t(perm$model.matrix)
write.table(t(perm$coef.sites) *  perm$model.matrix, "data/betas.txt", sep="\t")
perm
# plot coefficients 


###########################
## Indicator Value Analysis
###########################
test.indval <- cbind.data.frame(B_S.category, treats.perm)
test.indval$merge <- paste( test.indval$treats.perm, test.indval$B_S.category, sep="")

indval <- multipatt(df.all.no0, cluster = test.indval$merge, func = "IndVal.g", control = how(nperm=999)) 
summary(indval)

#########################
#Mean centroid dispersion
########################
beta.disp <- betadisper(d = df.all.db, group = test.indval$merge)

get.euc.dist.2D <- function(beta.disp, groups, axes_number){
  centroids <- beta.disp$centroids[,1:axes_number]
  posistions <- beta.disp$vectors[,1:axes_number]
  # get euclidean distancs from first three axes
  eucs <- c()
  pop.name <- c()
  treatment <- c()
  for(i in groups) {
    centroid.i <- centroids[i, ]
    positions.i <- posistions[rownames(posistions) %like% paste("_", i, sep = ""), ]
    for(j in 1:nrow(positions.i)) {
      position.j <- positions.i[j,]
      sample <- rownames(positions.i)[j]
      euc.dist.j <- dist(rbind(position.j, centroid.i))
      eucs <- c(eucs, euc.dist.j)
      pop.name <- c(pop.name, sample)
      treatment <- c(treatment, i)
    }
  }
  pop.euc <- cbind(pop.name, as.numeric(eucs), treatment) 
  return(data.frame(pop.euc))
}


euc.mat.m <- get.euc.dist.2D(beta.disp, unique(test.indval$merge), 3)
euc.mat.m$V2 <- as.numeric(as.character(euc.mat.m$V2))



png(filename = paste(c("figs/boxplot_B_S.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
boxplot(euc.mat.m$V2 ~ euc.mat.m$treatment)
dev.off()

