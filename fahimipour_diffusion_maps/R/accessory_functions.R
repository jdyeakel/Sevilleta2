## load some libraries
if(!require(tidyverse)) install.packages('tidyverse'); library(tidyverse)
if(!require(magrittr)) install.packages('magrittr'); library(magrittr)
if(!require(reshape2)) install.packages('reshape2'); library(reshape2)
if(!require(spatstat)) install.packages('spatstat'); library(spatstat)
if(!require(RColorBrewer)) install.packages('RColorBrewer'); library(RColorBrewer)
if(!require(parallelDist)) install.packages('parallelDist'); library(parallelDist)

##---
## functions we'll need for later
##---
## turn continuous vec to colors
map2color <- function(x, pal, limits = NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


## standardize columns func
norm.mat <- function(mat){
  ms <- apply(mat, 2, function(col) (col - mean(col)) / sd(col))
  ms
}

## compute euclid dist and standrdize to similarities
get.euc <- function(mat, n.threads = 1, alt = 'euclidean'){
  if(n.threads > 1){
    eu <- parDist(mat, method = alt, threads = n.threads) %>% as.matrix()
    ieu <- 1 / eu
    diag(ieu) <- 0
  }
  
  if(n.threads == 1){
    eu <- dist(mat, method = alt) %>% as.matrix()
    ieu <- 1 / eu
    diag(ieu) <- 0
  }

  ieu
}

## function to threshold the normalized distance mat
threshold <- function(mat, top_k = 10){
  
  thr <- mat
  tnr <- nrow(thr)
  
  ## set similarities outside of the top k to 0
  for(i in 1:tnr){
    
    ## rank entries in each row in reverse order (so largest value == rank 1), 
    ## and set entries that are outside of 1:k to 0
    thr[i, !rank(-thr[i, ], ties.method = 'random') %in% 1:top_k] <- 0
    
  }
  
  for(i in 1:tnr){
    for(j in 1:tnr){
      if(thr[i, j] == 0 & thr[j, i] != 0){
        thr[i, j] <- thr[j, i]
      }
    }
  }
  
  thr
  
}

## function to calculate norm laplacian
get.laplac <- function(mat){
  
  L <- -mat
  S <- rowSums(mat)
  nL <- L / S
  diag(nL) <- 1
  nL
  
}
