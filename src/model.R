pfull_0 <- function(lambda,k,ell,m,n,q,d){
    p0 <- ptrip_0(d, lambda, k, ell, m,n,0)
    p00 <- 1-(n-2)*p0
    return(max(0,p00))
}


pfull_1 <- function(lambda,k,ell,m,n,q, dmax){
    # first max p0 over T
    if (dmax > 1-ell) {dmax = 1-ell}
    res <- optim(dmax,ptrip_1, lambda=lambda,
        k=k,ell=ell,m=m,n=n,q=q, method="L-BFGS", lower=ell, upper=dmax,
        control=list(fnscale=-1))
    pmax <- res$value
    p0 <- 1-(choose(n,3)*pmax)
    return(max(0, p0))
}



ptrip_0 <- function(d, lambda,k,ell,m,n,q){
	p00 <- prob_seq_edits_in_grp(0,d, d+ell,lambda,k,0)^m
    return(p00)
}

ptrip_1 <- function(d,lambda,k,ell,m,n,q){
    p0 <- prob_trip_approx(lambda,k,ell,d,m,q)
    return(p0)
}

## fix m, find largest possible d (if it exists)
## fixed rate, m
p1_max_d <- function(k,ell,lambda, m, n, q,eps){
    d_max <- 1-ell
    
	f1 <- function(dneg, k,ell,lambda, n, eps,m){
	    d <- -dneg
	    res <- pfull_1(lambda, k, ell, m,n, q, d)
	return(res-eps)
	} 
	## check end points
	p0_end <- f1(-d_max, k, ell,lambda,n, eps, m)
	p0_start <- f1(0, k, ell, lambda, n, eps, m)
	if (p0_end >= 0){
		return(d=d_max)
	} else if (p0_start < 0) {
		return(d=0)
	} else {
	## find root closest to tips (d_max)
	root <- -1*uniroot(f1, c(-d_max, 0), k=k, ell=ell, lambda=lambda, n=n,eps=eps, m=m)$root
	return(root)
	 
	}
}

prob_trip_approx <- function(lambda, k, ell,d,m,q){
    mu_in <- m*expected_in_grp(d,d+ell,lambda,k,q)
    mu_out<- m*expected_out_grp(d,lambda,k,q)
    var_in <- m*var_in_grp(d,d+ell,lambda,k,q)
    var_out <- m*var_out_grp(d,lambda,k,q)
    ## normal approx
    mu_tot <- mu_in-mu_out
    var_tot <- var_in+var_out
    p_approx <- pnorm(1, mu_tot, sqrt(var_tot))
    return(p_approx)
}

prob_tripR_full <- function(lambda,k,ell,d,m,n,q){
    # pre compute pin and pout
    f_in <- sapply(0:k, prob_seq_edits_in_grp, d, d+ell, lambda, k, q)
    f_out <- sapply(0:k, prob_on_branch_homo, d,lambda,k,q)
    p0 <- 0
    mm <- m*k
    #in >= 1
    psum_in <- sapply(1:(mm), prob_allcomb, f_in,  k, m)
    psum_out <- sapply(0:(mm), prob_allcomb, f_out, k, m)
    psum_out_strict_less <- sapply(1:mm, function(i) sum(psum_out[1:i]))

    # put together
    
    p <- sum(psum_in*(psum_out_strict_less)^(n/2))
    
    return(p^(n-2))
}



prob_allcomb <- function(x,f,k,m){
   
    # get all combinations of (m-1) bars  
    table <- expand.grid(rep(list(0:k), m-1)) # all combinations
    row_sums <- rowSums(table)
    table <- table[row_sums <= x & row_sums >= x-k,] 
    #row_sums <- rowSums(table)
    ## get prod f...f(m-1)
    
    prod_f <- apply(table, 1, function(r) prod(f[r+1])) 
    f_a <- apply(table,1, function(r) f[x-sum(r)+1]) 

    return(sum(prod_f*f_a))
}

prob_zero <- function(lambda, k, ell, d, m, q){
    p0 <- prob_seq_edits_in_grp(0,d, d+ell, lambda, k, q)
    return(p0^m)
}

opt_p0_fix_d <- function(k,ell,eps,n, d2){
    res <- optimize(prob_zero, c(1,1/ell), tol=0.00001, k=k, ell=ell, d=d2,m=1,q=0)
    # find min m s.t. ptrip < eps
    lambda <- res$minimum
    p01 <- res$objective
    m <- log((1-eps)/(n-2))/log(p01)
    ###
    m <- ceiling(m)
    p01 <- pfull_0(lambda,k,ell,m,n,0,d2)
    return(data.frame(m=m, lambda=lambda, p0=p01))
}


## fix m, find largest possible d (if it exists)
## NB: rates <= 1000
opt_p0_fix_m <- function(k,ell,eps,m){
    d_max <- 1-ell
	f <- function(dneg, k,ell,eps,m){
		res <- optimize(ptrip_d, c(1, 100), tol=0.00001, k=k,ell=ell, d=-1*dneg, m=m)
		return(res$objective-eps)
	}    
	## check end points
	p0_end <- f(-d_max, k, ell, eps, m)
	p0_start <- f(0, k, ell, eps, m)
	if (p0_end <= 0){
		return(d=d_max)
	} else if (p0_start > 0) {
		return(d=0)
	} else {
	## find root closest to tips (d_max)
	root <- -1*uniroot(f, c(-d_max, 0), k=k, ell=ell, eps=eps, m=m)$root
	return(root) 
	}
}
## fix m, find largest possible d (if it exists)
## fixed rate, m
est_max_d <- function(k,ell,lambda, m, eps){
    d_max <- 1-ell
	f <- function(dneg, k,ell,eps,m){
	d <- -dneg
	res <- ptrip_d(lambda, k, ell, d, m)	
	return(res-eps)
	} 
	## check end points
	p0_end <- f(-d_max, k, ell, eps, m)
	p0_start <- f(0, k, ell, eps, m)
	if (p0_end <= 0){
		return(d=d_max)
	} else if (p0_start > 0) {
		return(d=0)
	} else {
	## find root closest to tips (d_max)
	root <- -1*uniroot(f, c(-d_max, 0), k=k, ell=ell, eps=eps, m=m)$root
	return(root)
	 
	}
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
	new_root <- colouring(root, tree$root.edge, rate, sampling, chars)	
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
	dist <- dist.topo(tree, upgma_tree, method="PH85")
	return(data.frame(k=k,i=i, m=m, sim_dist=dist))
}


generate_twostage_tree <- function(alpha, beta, n1, n, ratio){
	
    dev_trees <- sim.taxa(numbsim=1, n=2^n, waitsp=paste0("rweibull(",beta,",",alpha,")"), waitext="rexp(0)")
    tree <- dev_trees[[1]]
	tree <- makeNodeLabel(dev_trees[[1]])
	# RE SCALE

	depth <- node.depth(tree)
	tree$edge.length <- sapply(1:nrow(tree$edge), function(i) {
		end_node <- tree$edge[i,2]
		if (depth[end_node] > 2^(n-n1)){
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
simulate_and_score <- function(nsim, k, lambda, m,ell,j, tree, true_dists){
	chars <- all_chars[1:j]
	df <- 1:nsim %>% map_dfr(get_RF_score, tree, true_dists, k, lambda, m, chars)
	df$j <- j
	df$ell <- ell
	df$lambda1 <- unique(lambda)[1]
    df$lambda2 <- unique(lambda)[2]
    df$n <- Ntip(tree)
    return(df)
}

get_min_branch <- function(tree){
	#tree2 <- drop.tip(tree, tree$tip.label, trim.internal=FALSE)
	ee <- min(tree$edge.length[tree$edge[,2] > Ntip(tree)])
	return(min(ee, tree$root.edge))
}

generate_test_tree <- function(alpha, beta, n){
    dev_trees <- sim.taxa(numbsim=1, n=2^n, waitsp=paste0("rweibull(",beta,",",alpha,")"), waitext="rexp(0)")
    tree <- dev_trees[[1]]
    # make tree ultrametric   
    tree <- makeNodeLabel(dev_trees[[1]])
    dd <-max(node.depth.edgelength(tree))+tree$root.edge
    tree$edge.length <- tree$edge.length/dd 
    # and scale root edge
    tree$root.edge <- tree$root.edge/dd     
    return(tree)
}

generate_tree <- function(alpha, beta, n, ell){
	tree <- generate_test_tree(alpha,beta,n)
    ell2 <- get_min_branch(tree)
    it <- 1
    while (ell2 < ell && it < 10){
        tree <- generate_test_tree(alpha, beta, n)
        ell2 <- get_min_branch(tree)
        it <- it +1
    }
    return(tree)
}


prune_tree <- function(tree, rho){
	ntips <- length(tree$tip.label)
	n_delete <- sum(rbinom(ntips, 1, 1-rho))
	tips_to_delete <- sample(1:ntips, n_delete)
	tree <- drop.tip(tree, tips_to_delete)
	return(tree)
}

