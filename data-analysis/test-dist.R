source("../analysis/setup.R")

# generate tree & simulate
tree <- generate_tree(1, 100, 9, 0.001)
ell <- get_min_branch(tree)

true_dist <- cophenetic.phylo(tree)

# params
k <- 5
m <- 20
lambda_vec <- rep(6, m)
j <- 4
chars <- all_chars[1:j]
# run sim
res <- simulate_barcodes(tree, k, lambda_vec, m, chars)
# Build distance matrix from all barcodes
dists <- get_distance_sequential(res)
dists <- max(dists)-dists	
diag(dists) <- Inf
alpha <- min(dists)

#---- compare trees (true accuracy)
rownames(dists) <- colnames(dists) <- tree$tip.label
upgma_tree <- upgma(dists)	
# get RF dist
RF_dist <- dist.topo(tree, upgma_tree, method="PH85")


