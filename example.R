load_res <- FALSE
save_figs <- FALSE
source("setup.R")


###################
# Uniform tree
##################
# 2^16 cells, sample 1000, exp division times
n <- 10 # approx. number of generations
alpha <- 1
beta <- 100
trees <- sim.taxa(numbsim=1, n=2^n, waitsp=paste0("rweibull(",beta,",",alpha,")"), waitext="rexp(0)")

tree <- makeNodeLabel(trees[[1]])

# scale branch lengths
d <- max(node.depth.edgelength(tree))+tree$root.edge
tree$edge.length <- tree$edge.length/d
tree$root.edge <- tree$root.edge/d

# get min branch length
ell <- min(min(tree$edge.length), tree$root.edge)


#####################
# Parameters
#####################
k <- 12 # num target sites
j <- length(all_chars)
rates <- seq(5,20, by=2) #mutation rates
nsim <- 3 # num times repeat experiment
m <- 20 # copy num

params <- crossing(nsim=nsim, k=k, lambda=rates, m=m, ell=ell, j=j)

####################
# Simulate & Analyse
####################
# will take a few mins
res_unif <- params %>% pmap_dfr(simulate_and_score, tree=tree, true_dists=cophenetic.phylo(tree))

############################
############################
############################

###################
# Generate one tree
##################
# 2^16 cells, sample 1000, exp division times
n <- 16 # approx. number of generations
samp <- 1000
frac <- samp/2^n
nsim <- 1
# note: takes a few mins
trees <- sim.bd.taxa(samp, nsim, 1, 0, frac=frac, complete=FALSE)
tree <- makeNodeLabel(trees[[1]])

# scale branch lengths
d <- max(node.depth.edgelength(tree))+tree$root.edge
tree$edge.length <- tree$edge.length/d
tree$root.edge <- tree$root.edge/d

# get min branch length
ell <- min(min(tree$edge.length), tree$root.edge)
params <- crossing(nsim=nsim, k=k, lambda=rates, m=m, ell=ell, j=j)
####################
# Simulate & Analyse
####################
# will take a few mins
res <- params %>% pmap_dfr(simulate_and_score, tree=tree, true_dists=cophenetic.phylo(tree))

