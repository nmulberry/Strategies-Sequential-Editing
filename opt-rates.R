
source("setup.R")

## ---OPTIMISE

ell_vec <-  c(1e-1, 1e-2, 1e-3, 1e-4, 5e-5, 1e-5, 1e-7)
k_vec <- seq(1,20, by=1)
q_vec <- c(1, 1/2, 1/3, 1/100)
pars <- crossing(k=k_vec, ell=ell_vec, q=q_vec)

opt_E_over_lambda <- function(k,ell,q){
    # run optim
    res <- optim(k, min_expected_trips, k=k,n=1/ell,q=q,
        method="L-BFGS", lower=1, upper=k*n, control=list(fnscale=-1))
    return(data.frame(lambda=res$par, 
        E=res$value,
        k=k, ell=ell, q=q))
}

res <- pars %>% pmap_dfr(opt_E_over_lambda)
res$test <- res$k/(1-res$ell)
## -----PLOT
res$ell <- factor(res$ell)
res$q <- factor(res$q)
ggplot(res, aes(x=k,y = lambda, col=ell, group=ell))+
    geom_line()+
    facet_wrap(~q)




### test
j <- 10
m <- c(5,10)
k <- c(12)
nsim <- 50
n <- 10
params <- crossing(nsim=nsim, k=k, lambda=c(1,2), m=m, ell=1/n, j=j)
params <- params %>% mutate(lambda = case_when(lambda==1~k, TRUE~10.39))

tree <- generate_tree(alpha=1, beta=200, n=n)
true_dists <- cophenetic.phylo(tree)
res <- params %>% pmap_dfr(simulate_and_score, tree=tree, true_dists=true_dists)

ggplot(res, aes(x=k, y=log(triplets+1), col=m, group=m))+
    geom_boxplot()+
    facet_wrap(~lambda)
