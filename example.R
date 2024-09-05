####################################
# Get B0/B1 as func of editing rate
# Also show how to run simulations 
# note: simulations can be slow,
# only run locally for v. small trees
####################################
#--- load req packages and funcs
library(ape)
library(TreeSimGM)
library(tidyverse)
library(Rcpp)
sourceCpp("src/funcs.cpp")
source("src/model.R")
library(phangorn)
options(warn=-1)

library(ggplot2); 
theme_set(theme_minimal(base_size = 14)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))) 
#--- Set-up
n <- 5 # number of generations (small for demonstration)
ell <- 1/(n+1) # under synchronous division
chars <- c("a","b","c","d","e","f","g","h","i","j",
	"k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v",
	"w", "x", "y", "z", "$", "%","@", "&","(", ")","#",
	"!","^","?","<",">","/","{","}","[","]", "A", "B","C", "D", "E",
    "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P",
    "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
j <- length(chars) # num unique edits
q <- 1/j # under uniform insertion probs
k <- c(5,7,9) # number target sites per tape
m <- c(10,30) # number of tape copies per cell
lambda <- seq(0, 20, by=0.2)
nsim <- 100 # num simulation repeats



#--- Get B0 & B1 curves
pars <- crossing(lambda=lambda, k=k,ell=ell,m=m,n=n,q=q, d=1-ell)

res <- pars
res$B0 <- pars %>% pmap_dbl(pfull_0)
res$B1 <- pars %>% pmap_dbl(pfull_1)
res_long <- res %>% pivot_longer(cols=c(B0, B1))

#--- plot
g1 <- ggplot(res_long, aes(x=lambda,y=value, col=name))+
    geom_line(linewidth=1.1)+
    facet_grid(k~m, labeller=label_both)+
    labs(x="editing rate", y="probability", col="")

######################################################
# Verify with simulations (small tree)
######################################################
# Get tree
tree <- generate_tree(alpha=1, beta=200, n=n, ell=ell/2)
# check actual vs assumed min branch length 
ell_obs <- get_min_branch(tree)
print(paste("actual min branch length:", ell_obs))
print(paste("assumed min branch length:", ell))

# Get true distance matrix
true_dists <- cophenetic.phylo(tree)

# RUN SIMS
lambda <- seq(1,20, by=2) #editing rates
pars <- crossing(i=1:nsim, lambda=lambda, k=k, m=m)
sim_res <- pars %>% pmap_dfr(., ~get_RF_score(..1,tree, true_dists, ..3,..2,..4,chars))
sim_res$lambda <- pars$lambda

# get proportion of exact sims
sim_res <- sim_res %>%
    group_by(k,m,lambda) %>%
    summarize(value=sum(sim_dist==0)/n())
sim_res$name <- "simulated"


#--plot together
g2 <- g1 + geom_point(data=sim_res,shape=4)
