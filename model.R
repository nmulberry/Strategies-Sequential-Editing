min_expected_trips <- function(lambda, k,n,q){
	expected_trips(1-2/n,1-1/n, lambda,k,q)
}

min_shared_edits <- function(lambda, k,n,q){
	min_edits_on_branch(1-2/n,1-1/n, lambda,k,q)
}

# Note: doesn't work do not use
min_variance_trips <- function(lambda, k,n,q,E0){
	variance_trips(1-2/n,1-1/n, lambda,k,q,E0)
}

get_resolvable_trips <- function(k,lambda,n,m,q,eta){
	#d <- seq(0, 1-2/n, by=1/n)
	d <- 1-2/n
	p0 <- (p_edits_on_branch(0, lambda,k, 1/n, 1-2/n))^m

	#p00 <- mapply(p_trips_0, x=0, du=d, dv=d+1/n, lambda=lambda, k=k, q=q)
	#p0 <- dbinom(0, m, 1-p00) #prob 0 barcodes out of m have more than 0 edits 
	#dmax <- d[max(which(p0 >= eta))]
	#if (all((p0 >= eta) == FALSE)) {return(dmax=0)}
	#ntrip <- choose(2^(n*(d+2/n)),3)	
	#trips
	#tt <- ntrip*p0
	return(log(max(p0)))
}

min_edits_on_branch <- function(du, dv, lambda, k, q){
	ell <- dv-du
	x <- seq(0, k, by=1)
	res <- sum(x*sapply(x, p_edits_on_branch, lambda, k, ell, du))
	return(res) 
}

p_edits_on_branch <- function(x,lambda, k, ell, d){
	k_vec <- seq(0,k, by = 1)
	p_new <- mapply(trunc_poi_R, x, k-k_vec, lambda*ell)
	p_hist <- mapply(trunc_poi_R, k_vec, k, lambda*d)
	return(sum(p_new*p_hist)) 
}


trunc_poi_R <- function(n, max_n, rate){
	if (n > max_n) {
		p_x <- 0
	} else if (n==max_n) {
			p_x <- 1-sum(sapply(0:(max_n-1), dpois, rate)) 
	} else {
			p_x <- dpois(n, rate)
	}
	return(p_x)	
}

##-- SIMULATE AND SCORE TREES---


## Simulating barcode evolution on trees
# x: branch, l: branch length
# lambda: editing rate
# sample_p: dist for character sub
# chars: list of chars
simulate_barcode <- function(rate, chars, sampling,tree,root){
	# first, mutate along root edge	
	new_root <- colouring(root, tree$root.edge/tree$age, rate, sampling, chars)	
	mut_tree <- rTraitMult(tree, 
        	colouring, 
        	root.value=new_root, 
        	ancestor=TRUE, 
        	lambda=rate, 
        	sample_p=sampling, 
        	chars=chars)
   	ntip <- tree$Nnode+1
  	Y <- mut_tree$x1[1:ntip]
return(Y)
}

colouring <- function(x,l, lambda,sample_p,chars){
	p1 <- stringr::str_split(x,"")[[1]]
    t <- 0
    M <- length(p1)
    while (t <= l){
        time_to_mut <- rexp(1,lambda)
        t <- t + time_to_mut	
        if (p1[M] == "0" && t <= l) {
            pos <- min(which(p1 == "0"))
            new <- sample(chars, size=1, prob=sample_p) 
            p1[pos] <- new
        } 
    }
    return(paste(p1, collapse=""))
}


#colouring <- function(x,l,lambda,sample_p, chars){
#	p1 <- stringr::str_split(x, "")[[1]]
#	M <- length(p1)
#	n_edits <- rpois(1, lambda*l) # num edits
#	i <- 1
#	while (i <= n_edits && p1[M] == "0"){
#		pos <- min(which(p1 == "0"))
#		p1[pos] <- sample(chars, size=1, prob=sample_p)
#		i <- i+1
#	}
#	return(paste(p1, collapse="")) 
#}


simulate_barcodes <- function(tree, k, lambda_vec, m, chars){
	# sim bars over tree
	root <- paste0(rep("0",k), collapse='')
    sample_p <- rep(1/length(chars), length(chars))
	res <- lapply(lambda_vec, simulate_barcode, chars, sample_p,tree,root)
	return(res)
}

# returns dataframe with simulated trip score and upgma dist
get_RF_score <- function(i,tree,true_dists,k,lambda,m,chars){
	# simulate
	if (length(lambda)==m){
		lambda_vec <- lambda
	} else {
		lambda_vec <- rep(lambda, m)
	}
	res <- simulate_barcodes(tree, k, lambda_vec, m, chars)
	# Build distance matrix from all barcodes
	dists <- get_distance_sequential(res)
	dists <- max(dists)-dists	
	rownames(dists) <- colnames(dists) <- tree$tip.label
	upgma_tree <- upgma(dists)	
	# get RF dist
	#dist <- TreeDistance(tree, upgma_tree)
	dist <- dist.topo(tree, upgma_tree, method="PH85")
	# observed triplets
    trips <- proportion_triplets(true_dists, dists)
	#print(paste("rate:", lambda, "i:", i, "done"))
	return(data.frame(k=k,i=i, m=m, triplets=trips, sim_dist=dist))
}

# returns ultrametric tree
# alpha, beta: shape/scale of Weibull branching times
# n: number generations (2^n tips)
generate_tree <- function(alpha, beta, n){
    dev_trees <- sim.taxa(numbsim=1, n=2^n, waitsp=paste0("rweibull(",beta,",",alpha,")"), waitext="rexp(0)")
    tree <- dev_trees[[1]]
    # make tree ultrametric   
	tree <- makeNodeLabel(dev_trees[[1]])
	dd <-max(node.depth.edgelength(tree))
	tree$edge.length <- tree$edge.length/dd 
	# and scale root edge
	tree$root.edge <- tree$root.edge/dd
	return(tree)
}

generate_twostage_tree <- function(alpha, beta, n1, n, ratio){
	
    dev_trees <- sim.taxa(numbsim=1, n=2^n, waitsp=paste0("rweibull(",beta,",",alpha,")"), waitext="rexp(0)")
    tree <- dev_trees[[1]]
	tree <- makeNodeLabel(dev_trees[[1]])
	# RE SCALE

	depth <- node.depth(tree)
	tree$edge.length <- sapply(1:nrow(tree$edge), function(i) {
		end_node <- tree$edge[i,2]
		if (depth[end_node] >= 2^n1){
			l <- tree$edge.length[i]*ratio
		} else {
			l <- tree$edge.length[i]
		} 
	})		
	# also fix stem
	tree$root.edge <- tree$root.edge*ratio

	# fix tip lengths--make ultrametric
	tree <- force.ultrametric(tree) # tree should be close to ultra
	# rescale to age 1
	dd <- max(node.depth.edgelength(tree))
	tree$edge.length <- tree$edge.length/dd 
	tree$root.edge <- tree$root.edge/dd
	return(tree)

}


# Simulate process and return df with different scores
simulate_and_score <- function(nsim, k, lambda, m,n,j, tree, true_dists){

	chars <- all_chars[1:j]
	df <- 1:nsim %>% map_dfr(get_RF_score, tree, true_dists, k, lambda, m, chars)
	df$j <- j
	df$n <- n
	return(df)
}



read_results <- function(file,dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	model <- name[[1]][1]
	n <- name[[1]][2]
	m <- name[[1]][3]
	j <- name[[1]][4]
	k <- name[[1]][5]
	k <- strsplit(k, split="\\.")[[1]][1]
	df$model <- model
	df$n <- as.numeric(n)
	df$m <- as.numeric(m)
	df$j <- as.numeric(j)
	df$k <- as.numeric(k)
	return(df)
}

get_character_matrix <- function(){
	# get tree
	tree <- generate_tree(1, 200, n)
	# run simulation
	res <- simulate_barcodes(tree, k, lambda, m, chars)
	# process results:
	# first get priors on insertion freqs
	freqs <- lapply(res, get_insert_freqs)	
	freqs <- do.call(rbind, freqs)
	prior <- data.frame(chars=1:length(chars))
	prior$freqs <- rowMeans(freqs)	
	# next, remove nonsequentially identical edits
	clean_res <- lapply(res, remove_homo_edits)
	# now create character matrix & combine bars
	# n_ij = mutation in ith cell at position j
	# each char becomes an integer, 0 becomes -1 ("misssing")
	tot_res <- do.call(paste0, clean_res)
	tot_res_list <- strsplit(tot_res, "")
	M <- unlist(tot_res_list)
	# map each character to an integer	
	M2 <- match(M, chars)
	M2[M=="*"] <- -1
	M2[M=="0"] <- 0
	M2 <- matrix(M2, ncol=m*k, byrow=TRUE)		
	rownames(M2) <- paste0("c", seq(1:2^n))
	colnames(M2) <- paste0("p", seq(1:(m*k)))

}

get_insert_freqs <- function(bar){
	char <- unlist(strsplit(bar, ''))
	char <- char[char != "0"]
	# count
	counts <- table(char)
	counts <- counts/sum(counts)
	counts <- data.frame(counts)
	return(counts$Freq)	
}
