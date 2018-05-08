################################################################################
#                                                                              #
# ParEvol Source Code: Functions for Parallel Evolution Analyses               #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Will Shoemaker                                                   #
#                                                                              #
# Created: 2018/02/07                                                          #
#                                                                              #
# Last update:  2018/02/07                                                     #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code provides numerous functions to be used in the analysis of   #
#        mutation data post breseq anlaysis                                    #
#        The ideas in the code are the products of a project between           #
#        Will Shoemaker and Jay Lennon as well as communication                #
#        with Mario Muscarella                                                 #
################################################################################

make.null.matrix <- function(df, num.muts, pmf){
  # make an empty matrix using the dimensions of the pop_by_gene matrix
  df.null <- matrix(data = 0, nrow = nrow(df), ncol = ncol(df), byrow = FALSE,
                    dimnames = list(rownames(df), colnames(df)))
  for (i in seq_along(num.muts)){
    muts <- num.muts[i]
    sample.muts <- sample(colnames(df), size = muts, replace = TRUE, prob = pmf)
    for(j in sample.muts){
      df.null[i,j] <- df.null[i,j] + 1
    }
  }
  return(df.null)
}


make.G.matrix <- function(df, pmf){
  df.G <- matrix(data = 0, nrow = nrow(df), ncol = ncol(df), byrow = FALSE,
                 dimnames = list(rownames(df), colnames(df)))
  # iterate row 
  for(i in 1:nrow(df)) {
    row <- df[i,]
    mut.total <- sum(row)
    E <- mut.total * pmf
    G_i <- 2*row*log(row/E)
    G_i <- replace(G_i, is.na(G_i), 0)
    # make negative values zero
    G_i <- replace(G_i, which(G_i < 0), 0)
    df.G[i,] <- G_i
  }
  return(df.G)
}


# multiplicity score
make.m.matrix <- function(df, sizes){
  df.m <- matrix(data = 0, nrow = nrow(df), ncol = ncol(df), byrow = FALSE,
                 dimnames = list(rownames(df), colnames(df)))
  # iterate row 
  mean.size <- mean(sizes)
  mean.over.size <- mean.size / sizes
  for(i in 1:nrow(df)) {
    row <- df[i,]
    m_i <- mean.over.size * row
    mean.m <- sum(row) / length(row)
    #print(m_i/mean.m)
    #print(length(mean.m))
    #print(row * log(m_i/mean.m))
    M_i <- row * log(m_i/mean.m)
    M_i <- replace(M_i, which(M_i < 0), 0)
    M_i <- replace(M_i, is.na(M_i), 0)
    #df.m[i,] <- m_i
    df.m[i,] <- M_i
  }
  return(df.m)
}


get.euc.dist.2D <- function(beta.disp){
  centroids.2 <- beta.disp$centroids[,1:3]
  posistions.2 <- beta.disp$vectors[,1:3]
  # get euclidean distancs from first two axes
  eucs <- c()
  pop.name <- c()
  treat.strain.names <- c()
  for(i in 1:nrow(centroids.2)) {
    centroid.i <- centroids.2[i,]
    treat.strain <- rownames(centroids.2)[i]
    pos.treat.strain <- posistions.2[grep(treat.strain, rownames(posistions.2)), ]
    for(j in 1:nrow(pos.treat.strain)) {
      dist.j <- pos.treat.strain[j,]
      sample <- rownames(pos.treat.strain)[j]
      euc.dist <- dist(rbind(dist.j, centroid.i))
      eucs <- c(eucs, euc.dist)
      pop.name <- c(pop.name, sample)
      treat.strain.names <- c(treat.strain.names, treat.strain)
    }
  }
  pop.euc <- cbind(pop.name, eucs, treat.strain.names) 
  return(data.frame(pop.euc))
}


sim.euc.dist.G <- function(df, pop.muts, pmf, df.G.groups, iter){
  mat <- data.frame(matrix(NA, nrow=0, ncol=length(unique(df.G.groups))))
  for(i in 1:iter) {
    df.null <- make.null.matrix(df.merge, pop.muts, pmf)
    df.null.G <- make.G.matrix(df.null, pmf)
    df.null.G.no0 <- df.null.G[rowSums(df.null.G[,-1]) != 0,]
    df.null.G.no0.db <- vegdist(df.null.G.no0, method = "bray", upper = TRUE, diag = TRUE)
    # get groups, we can use the same group labels
    #df.G.null.groups <- substr(c(sapply(strsplit(rownames(df.null.G.no0),"_"), `[`, 2)), 1, 3)
    beta.disp.null <- betadisper(d = df.null.G.no0.db, group = df.G.groups)
    euc.mat.null <- get.euc.dist.2D(beta.disp.null)
    euc.mean.null <- aggregate(eucs ~ treat.strain.names, euc.mat.null, FUN = function(x) mean(as.numeric(as.character(x))))
    colnames(mat) <- euc.mean.null$treat.strain.names
    if(i==1) {
      names(mat) <- euc.mean.null$treat.strain.names
    }
    mat <- rbind(mat, euc.mean.null$eucs)
  }
  return(mat)
}


sim.euc.dist.m <- function(df, pop.muts, sizes, df.m.groups, iter){
  mat <- data.frame(matrix(NA, nrow=0, ncol=length(unique(df.m.groups))))
  for(i in 1:iter) {
    df.null <- make.null.matrix(df.merge, pop.muts, pmf)
    df.null.m <- make.m.matrix(df.null, sizes)
    df.null.m.no0 <- df.null.m[rowSums(df.null.m[,-1]) != 0,]
    df.null.m.no0.db <- vegdist(df.null.m.no0, method = "bray", upper = TRUE, diag = TRUE)
    # get groups, we can use the same group labels
    #df.G.null.groups <- substr(c(sapply(strsplit(rownames(df.null.G.no0),"_"), `[`, 2)), 1, 3)
    beta.disp.null <- betadisper(d = df.null.m.no0.db, group = df.m.groups)
    euc.mat.null <- get.euc.dist.2D(beta.disp.null)
    euc.mean.null <- aggregate(eucs ~ treat.strain.names, euc.mat.null, FUN = function(x) mean(as.numeric(as.character(x))))
    colnames(mat) <- euc.mean.null$treat.strain.names
    if(i==1) {
      names(mat) <- euc.mean.null$treat.strain.names
    }
    mat <- rbind(mat, euc.mean.null$eucs)
  }
  return(mat)
}


#make.M.matrix
