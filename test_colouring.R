source("setup.R")
nsim <- 1000
chars <- all_chars
k <- 9
root <- paste0(rep("0", k), collapse='')
j <- c(1,2,10,20)
#j <-1
q <- 1/j
lambda <- seq(1,16,by=5)
#lambda <- 11
i <- 1:nsim
pars0 <- crossing(i=i, lambda=lambda,j=j)
save_figs <- TRUE
fig_dir <- "./"
#====================================#
# -------MORE TESTS------------------#
#====================================#
## test homo
run_colour_model_homo <- function(i,lambda,j){
    chars <- all_chars[1:j]
    sample_p <- rep(1/length(chars), length(chars))
    resa <- colouring(root, 1.0, lambda, sample_p, chars)
    resb <- colouring(root, 1.0, lambda, sample_p, chars)
    # count
    s <- count_sequential_matches(resa,resb)
    return(s)
}

res <- pars0
res$num_edits <- res %>% pmap_dbl(run_colour_model_homo)

gg1_diffs <- ggplot(res, aes(x=num_edits))+geom_histogram(col="blue", alpha=0.3, fill="blue")
## compare to prob distn
pars <- crossing(x=0:k, lambda=lambda,k=k,d=0.0, q=q)
pars$p <- pars %>% pmap_dbl(corr_prob_homo)
pars$j <- 1/pars$q

# put together
res2 <- res %>% group_by(lambda, num_edits,j) %>%
    summarize(sim_freq=n()/nsim)%>%
    rename(x=num_edits)

res2 <- merge(res2, pars)

gg3_homo <- ggplot(res2)+
   geom_col(aes(x=x, y=p),col="red", alpha=0.3, fill="red")+
   geom_col(aes(x=x, y=sim_freq),col="blue", alpha=0.3, fill="blue")+
   facet_grid(j~lambda, labeller=label_bquote(rows=j:.(j), cols=lambda:.(lambda)))+
  labs(x="# Shared Edits", y="")

if (save_figs){
    ggsave(paste0(fig_dir, "/prob_homoplasy.pdf"), height=6, width=8)
}
#=======================================#
## test outgrp
run_colour_model_outgrp <- function(i,lambda,j){
    chars <- all_chars[1:j]
    sample_p <- rep(1/length(chars), length(chars))
    resu <- colouring(root, 0.2, lambda, sample_p, chars)
    resa <- colouring(resu, 0.8, lambda, sample_p, chars)
    resb <- colouring(resu, 0.8, lambda, sample_p, chars)
    # count
    s <- count_sequential_matches(resa,resb)-count_sequential_matches(resa, resu)
    return(s)
}

res <- pars0
res$num_edits <- res %>% pmap_dbl(run_colour_model_outgrp)

gg1_diffs <- ggplot(res, aes(x=num_edits))+geom_histogram(col="blue", alpha=0.3, fill="blue")
## compare to prob distn
pars <- crossing(x=0:k,d=0.2,lambda=lambda, k=k, q=q)
pars$p <- pars %>% pmap_dbl(prob_on_branch_homo)
pars$j <- 1/pars$q

# put together
res2 <- res %>% group_by(lambda, num_edits,j) %>%
    summarize(sim_freq=n()/nsim)%>%
    rename(x=num_edits)

res2 <- merge(res2, pars)

gg3_outgrp <- ggplot(res2)+
   geom_col(aes(x=x, y=p),col="red", alpha=0.3, fill="red")+
   geom_col(aes(x=x, y=sim_freq),col="blue", alpha=0.3, fill="blue")+
   facet_grid(j~lambda, labeller=label_bquote(rows=j:.(j), cols=lambda:.(lambda)))+
  labs(x="# Shared Edits", y="")

if (save_figs){
    ggsave(paste0(fig_dir, "/prob_outgrp.pdf"), height=6, width=8)
}

## test ingrp
run_colour_model_ingrp <- function(i,lambda,j, du, dv){
    chars <- all_chars[1:j]
    sample_p <- rep(1/length(chars), length(chars))
    resu <- colouring(root, du, lambda, sample_p, chars)
    resv <- colouring(resu, dv-du, lambda, sample_p, chars)
    resa <- colouring(resv, 1-dv, lambda, sample_p, chars)
    resb <- colouring(resv, 1-dv, lambda, sample_p, chars)
    # count
    s <- count_sequential_matches(resa,resb)-count_sequential_matches(resa,resu)
    return(s)
}

res <- pars0
du <- 0.6
dv <- du + 0.1
res$du <- du
res$dv <- dv
res$num_edits <- res %>% pmap_dbl(run_colour_model_ingrp)




gg1_diffs <- ggplot(res, aes(x=num_edits))+geom_histogram(col="blue", alpha=0.3, fill="blue")
## compare to prob distn

pars <- crossing(x=0:k,du=du, dv=dv,lambda=lambda, k=k, q=q)
pars$p <- pars %>% pmap_dbl(prob_seq_edits_in_grp)
pars$j <- 1/pars$q

# put together
res2 <- res %>% group_by(lambda, num_edits,j) %>%
    summarize(sim_freq=n()/nsim)%>%
    rename(x=num_edits)

res2 <- merge(res2, pars)

gg3_ingrp <- ggplot(res2)+
   geom_col(aes(x=x, y=p),col="red", alpha=0.3, fill="red")+
   geom_col(aes(x=x, y=sim_freq),col="blue", alpha=0.3, fill="blue")+
   facet_grid(j~lambda, labeller=label_bquote(rows=j:.(j), cols=lambda:.(lambda)))+
  labs(x="# Shared Edits", y="")

if (save_figs){
    ggsave(paste0(fig_dir, "/prob_ingrp.pdf"), height=6, width=8)
}
#=======================================#
## test num edits
run_colour_model <- function(i,rate){
    chars <- all_chars
    sample_p <- rep(1/length(chars), length(chars))
    resa <- colouring(root, 1.0, rate, sample_p, chars)
    # count number edits
    s <- count_sequential_matches(resa,resa)
    return(s)
}

res <- crossing(i=i, rate=lambda)
res$num_edits <- res %>% pmap_dbl(run_colour_model)

gg1_diffs <- ggplot(res, aes(x=num_edits))+geom_histogram(col="blue", alpha=0.3, fill="blue")
## compare to prob distn
pars <- crossing(n=0:k,max_n=k, rate=lambda)
pars$p <- pars %>% pmap_dbl(trunc_poi)
pars$x <- pars$n

# put together
res2 <- res %>% group_by(rate, num_edits) %>%
    summarize(sim_freq=n()/nsim)%>%
    rename(x=num_edits)

res2 <- merge(res2, pars)

gg3_base <- ggplot(res2)+
   geom_col(aes(x=x, y=p),col="red", alpha=0.3, fill="red")+
   geom_col(aes(x=x, y=sim_freq),col="blue", alpha=0.3, fill="blue")+
   facet_grid(.~rate, labeller=label_bquote(cols=lambda:.(rate)))+
  labs(x="Number edits", y="")

if (save_figs){
    ggsave(paste0(fig_dir, "/prob_trunc_poi.pdf"), height=5, width=10)
}

#############################
