##---
## This script demonstrates the diffusion map procedure for
## a subset of the metabolic networks analyzed by
## Fahimipour and Gross (in review).
##---
## load libraries
library(dplyr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(tibble)

## set working directory
setwd('/path/to/this/directory')

## point to the subdir containing metabolic network edge lists
net_loc <- './data'

## get a list of files
file_list <- list.files(path = net_loc, full.names = T)

## load all the metabolic networks into a list
nets <- list()

## loop through files
for(k in 1:length(file_list)){
  
  ## read the network
  temp <- read.table(file_list[k], sep = '\t', header = T)
  
  ## store in our list
  nets[[k]] <- temp
  
}

## bind our edge lists together
all_nets <- do.call('rbind', nets)

## retain the reaction and species columns
all_nets <- all_nets[, names(all_nets) %in% c('rxn', 'taxon')]

## add a dummy unity vector for casting into wide format
all_nets$dummy <- rep(1)

## cast into wide form
m <- dcast(all_nets, formula = taxon ~ rxn, value.var = 'dummy', fill = 0) %>%
  column_to_rownames(., var = 'taxon') %>%
  as.matrix()

##
## end data preprocessing
##


##---
## diffusion mapping starts here
##---
## load our defined functions for diffusion map steps
source('./R/accessory_functions.R')

## set seed
set.seed(777)

## normalize
nm <- m %>% norm.mat()

## make affinity matrix
aff <- nm %>%
  get.euc(., n.threads = 1) %>% 
  threshold(., top_k = 10)

## compute laplacian
Lij <- aff %>% 
  get.laplac()

## smallest keig vectors
eig <- eigen(Lij)

## get eigenvalues
evl <- eig$values %>%
  Re() %>%
  round(., 10)

## get eigenvectors
evc <- eig$vectors %>%
  Re() %>%
  round(., 10)

## create objects containing diffusion variables
for(d in 1:length(evl)){
  assign(paste('dim', d, sep = '_'),
         Re(evc[, rank(evl) == (d + 1)])
  )
}

## specify number of variables you want to retain
k_eig <- 10

## merge in array
dat <- do.call(mapply, c(FUN = cbind, mget(paste0("dim_", 1:(k_eig))))) %>%
  t() %>%
  as.data.frame()

## add labels
colnames(dat) <- paste(paste('dim', 1:(k_eig), sep = '_'))
rownames(dat) <- rownames(Lij)

##
## you now have a data frame called 'dat' that contains
## the first 10 diffusion variables describing major
## variation in this subset of 100 genomes.
## We can explore these quickly by visualizing
## the first two dimensions as new coordinates.
##
ggplot(dat, aes(x = dim_1, y = dim_2, fill = dim_1)) +
  theme_classic() +
  xlab('Variable 1') +
  ylab('Variable 2') +
  geom_hline(yintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_vline(xintercept = 0, size = 0.25, colour = '#bdbdbd') +
  geom_point(shape = 21, size = 1) +
  scale_fill_gradientn(colours = brewer.pal(11, 'RdBu'), guide = F) +
  theme(aspect.ratio = 1)
  






