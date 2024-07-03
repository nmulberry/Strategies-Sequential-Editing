source("setup.R")


if (FALSE) {


nsim <- 500
chars <- all_chars
k <- c(5,7)
lambda <- seq(1,30, by=1)
m <- c(10,30)
i <- 1:nsim
j <- c(2,4,16)
d <- 0.8
ell <- c(0.1, 0.01)
pars0 <- crossing(i=i, lambda=lambda,j=j, m=m,k=k,d=d, ell=ell)
##===========================#
## test diff
run_colour_model_diff <- function(i,lambda,j,m,k,d,ell){
    chars <- all_chars[1:j]
    root <- paste0(rep("0", k), collapse='')
    sample_p <- rep(1/length(chars), length(chars))
    s_in <- 0
    s_outa <- 0
    s_outb <- 0
    for (i in 1:m){
        res0 <- colouring(root, d, lambda, sample_p, chars)
        resc <- colouring(res0, 1.-d, lambda, sample_p, chars)
        resv <- colouring(res0, ell, lambda, sample_p, chars)
        resa <- colouring(resv, 1-d-ell, lambda, sample_p, chars)
        resb <- colouring(resv, 1-d-ell, lambda, sample_p, chars)
        # count
        s_in <- s_in + count_sequential_matches(resa,resb)
        s_outa <- s_outa + count_sequential_matches(resa,resc)
        s_outb <- s_outb + count_sequential_matches(resb, resc)
     }

    return(s_in-s_outa)
}

res <- pars0
res$num_edits <- res %>% pmap_dbl(run_colour_model_diff)

res22 <- res %>%
    group_by(lambda,j,m,k,d,ell)%>%
    summarize(ptrip=sum(num_edits > 0)/n())%>%
    ungroup()



res22 <- res22 %>% mutate(q=1/j) %>%
mutate(p01 = pmap_dbl(., ~ ptrip_d2(..1,..4,..6,..5,..3,2,..8)), 
    p02 = pmap_dbl(., ~ ptrip_d(..1,..4,..6,..5,..3,2,..8)))

res22 <- res22 %>%
    pivot_longer(cols=starts_with("p"))%>%
    mutate(name=case_when(name=="ptrip"~"Simulated", name=="p01"~"Approx", name=="p02"~"P0"))
    

ggplot(filter(res22), aes(lambda, y=value, group=name, col=name, linetype=name))+
    geom_line()+
    scale_colour_manual(values=c(Simulated="black", Approx="red", Full="blue"))+
    scale_linetype_manual(values=c(Simulated="solid", Approx="dashed", Full="dotted"))+
    facet_grid(m+j~k+ell, labeller=label_both)+theme(legend.position="bottom")+
    labs(y="Probability resolve (a,b|c)",col="", linetype="", x=expression(lambda))

ggsave("onetrip-sim-all.pdf")


ggplot(filter(res22, ell==0.1, k==5), aes(lambda, y=value, group=name, col=name, linetype=name))+
    geom_line()+
    scale_colour_manual(values=c(Simulated="black", Approx="red", Full="blue"))+
    scale_linetype_manual(values=c(Simulated="solid", Approx="dashed", Full="dotted"))+
    facet_grid(m~j, labeller=label_both)+theme(legend.position="bottom")+
    labs(y="Probability resolve (a,b|c)",col="", linetype="", x=expression(lambda))

ggsave("onetrip-sim-small.pdf")
}


#### trees
test_p_over_trees <- function(n, beta, ell){
    tree <- generate_tree(alpha=1, beta=beta, n=n, ell=ell)
    ell<- get_min_branch(tree)
    true_dists <- cophenetic.phylo(tree)
    ## run sims
    pars <- crossing(nsim=100, k=c(1,5), lambda = seq(1,30, by=1), m=c(10,50),ell=ell, j=c(4,64))
    res <- pars %>% pmap_dfr(simulate_and_score, tree=tree, true_dists=true_dists)
    res$lambda <- res$lambda1
    res_summ <- res %>%
        group_by(k,m,j,ell,lambda) %>%
        summarize(trip_score = sum(triplets==0)/n(),
            upgma_score = sum(sim_dist==0)/n())
    res_summ$n <- 2^n
    res_summ$q <- 1/res_summ$j
    res_summ$p0 <- res_summ %>% pmap_dbl(., ~pfull_0(..5,..1,..4,..2,..8,..9))
    res_summ$p1 <- res_summ %>% pmap_dbl(., ~pfull_1(..5,..1,..4,..2,..8,..9))
    res_summ$p2 <- res_summ %>% pmap_dbl(., ~pfull_2(..5,..1,..4,..2,..8,..9))
    res_summ$beta <- beta
    return(res_summ)
}



pars <- crossing(n=c(4), beta=c(10, 100), ell=c(0.01)) 
res <- pars %>% pmap_dfr(test_p_over_trees)


res$k <- factor(res$k)
gg <- ggplot(res, aes(x=lambda))+
    geom_line(aes(y=upgma_score), linewidth=1.05)+
    geom_line(aes(y=p0),linewidth=1.05, linetype="dashed", col="red")+
    geom_line(aes(y=p2), linewidth=1.05, alpha=1, linetype="dashed", col="lightblue")+
    geom_line(aes(y=p1), linewidth=1.05, alpha=1, linetype="dashed",col="pink")+
    facet_grid(m+beta~j+k, labeller=label_both)+
    labs(y="Reconstruction Probability", x=expression(lambda))

ggsave("all_test_trees.pdf")
