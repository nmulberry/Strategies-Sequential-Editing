source("../analysis/setup.R")
#----------------------------#
test_q_estimator <- function(it, k, j, tree){
    chars <- all_chars[1:j]    
    true_q <- 1/j
    eps <- 0.9 #desired accuracy
    n <- length(tree$tip.label)
    opt_res <- opt_p0_fix_d(k,ell,eps,n,1-ell)
    m <- opt_res$m
    # goal rate:
    lambda_goal <- opt_res$lambda
    # actual rate: (dispersed about lambda_goal)
    lambda <- rgamma(m, shape=10, scale=lambda_goal/10)
    ## simulate (one exp)
    dat <- simulate_barcodes(tree, k, lambda, m, chars)

    ## get estimates per tape 
    res <- lapply(dat, mean_homo)
    res <- as.data.frame(do.call(rbind,res))
    colnames(res) <- c("h", "N")

    # per tape estimates
    res$q_est <- res$h/res$N

    # combine
    #tot_est <- sum(res$h)/sum(res$N)
    res$it <- it
    res$true_q <- 1/j
    res$k <- k
    return(res)
}

# generate tree
tree <- generate_tree(1, 100, 9, 0.001)
ell <- get_min_branch(tree)

# Run over k, q
num_it <- 100
pars <- crossing(it=1:num_it, k=c(5,13), j=c(4,8,16,32,64))
res <- pars %>% pmap_dfr(test_q_estimator, tree=tree)
#---------------------------------#
#-- plot one iteration, hi/lo q
res <- res %>% mutate(k_lab = paste(k, "sites"))%>%
group_by(true_q, k_lab, it) %>% mutate(pooled_q=sum(h)/sum(N)) %>% ungroup()
g0 <- ggplot(filter(res, it==1), aes(x=q_est))+
    geom_histogram(col="black", fill="lightblue")+
    facet_grid(k_lab~true_q, scales="free")+
    geom_vline(aes(xintercept=true_q), col="red")+
    geom_vline(aes(xintercept=pooled_q), col="red", linetype="dashed")+
    labs(y="", x="Estimated q per tape")

ggsave("estimate_q_by_tape.pdf", height=4, width=9)

#-- plot pooled estimates
res_pooled <- res %>% 
    group_by(k_lab,true_q, it) %>%
    summarize(pooled_q = sum(h)/sum(N))

g1 <- ggplot(res_pooled, aes(x=true_q, y=pooled_q,
    group=true_q))+
    geom_boxplot()+
    geom_abline(slope=1, intercept=0, col="red", linetype="dashed")+
    facet_wrap(~k_lab)+
    labs(x="True q", y="Pooled estimate q")
ggsave("pooled_estimates_q_100it.pdf", height=4, width=8)
