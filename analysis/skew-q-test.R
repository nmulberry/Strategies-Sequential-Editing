##--- Optimal Rates---##
# REDO WITH SKEWED Q DISTRIBUTION
# on small trees (can run locally)
set.seed(124)
source("setup.R")


## generate a skewed distribution
skewed_edits <- rlnorm(20)
skewed_edits <- skewed_edits/sum(skewed_edits)
q <- sum(skewed_edits^2)

# now ensure that 1/q is approx. integer
q_desired <- 1/floor(1/q)

# fudge the dist
skewed_edits[21] <- sqrt(q_desired-q)
skewed_edits <- sort(skewed_edits)

n_skew <- length(skewed_edits)
skew_chars <- all_chars[1:n_skew]
# check
q <- sum(skewed_edits^2)
n_unif <- floor(1/q)
unif_edits <- rep(q, n_unif) 
skew_chars <- all_chars[1:n_skew]
unif_chars <- all_chars[1:n_unif]


#==========#
#simple plot
dists_df <- data.frame(x=skew_chars, Skewed=skewed_edits)
unif_df <- data.frame(x=unif_chars, Uniform=unif_edits) 
df_tot <- left_join(dists_df, unif_df)



####################################
#Set up
n <- 6 # number of generations (small for demonstration)
ell <- 1/(n+1) # under synchronous division

k <- c(5,9) # number target sites per tape
m <- c(10,30,50) # number of tape copies per cell
lambda <- seq(1, 20, by=0.2)
nsim <- 50 # num simulation repeats

#--- Get Binf & Bq curves
pars <- crossing(lambda=lambda, k=k,ell=ell,m=m,n=2^n,q=q, d=1-ell)

res <- pars
res$Binf <- pars %>% pmap_dbl(pfull_0)
res$Bq <- pars %>% pmap_dbl(pfull_1)
res_long <- res %>% pivot_longer(cols=c(Binf, Bq))

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
pars <- crossing(i=1:nsim, k=k, lambda=lambda, m=m)
sim_res_unif <- pars %>% pmap_dfr(., get_RF_score, 
    chars=unif_chars, edit_dist=unif_edits, method="UPGMA",
    tree=tree, true_dists=true_dists)

sim_res_unif$lambda <- pars$lambda
sim_res_unif$type <- "Uniform"
print("re-running with skewed distribution")
# redo for skewed dist
sim_res_skew <- pars %>% pmap_dfr(., get_RF_score, 
    chars=skew_chars, edit_dist=skewed_edits, method="UPGMA",
    tree=tree, true_dists=true_dists)

sim_res_skew$lambda <- pars$lambda
sim_res_skew$type <- "Skewed"

sim_res <- rbind(sim_res_skew, sim_res_unif)


# get proportion of exact sims
sim_res <- sim_res %>%
    group_by(k,m,lambda, type) %>%
    summarize(value=sum(sim_dist==0)/n())
sim_res$name <- "Simulated triplet score"
#-----------------------#
# plot vals

gg0 <- ggplot(df_tot %>% pivot_longer(cols=c(Skewed, Uniform)),
    aes(x=x, y=value))+
    geom_hline(yintercept=1/12, col="gray")+
    geom_col(col="black", fill="gray")+facet_wrap(~name)+
#    scale_fill_manual(values=c("gray", "gold"))+
    theme(legend.position="none", axis.text.x=element_blank())+
    scale_y_continuous(expand=c(0,0), breaks=c(1/6, 1/12, 1/24), labels=c("", "1/12", ""))+
    labs(y="", x="Edit")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "#f0f0f0"),
          strip.text = element_text(face = "bold"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    geom_text(aes(x=6, y=0.12, label="q=1/12"), col="gray")

gg <- ggplot()+
    geom_line(data=res_long, linewidth=1, aes(x=lambda, y=value, col=name), alpha=0.75)+
    geom_point(data=sim_res, size=2,aes(x=lambda, y=value, shape=type), col="black")+
#    geom_point(data=sim_res, size=1, aes(x=lambda, y=value, col=type, shape=type, fill=type))+
    facet_grid(k~m, labeller=label_both)+
    scale_y_continuous(breaks=c(0.5, 1.0))+
    scale_color_manual(values=c("lightblue", "red", "gray", "gold"))+
    scale_shape_manual(values=c(2,8))+
    labs(x=expression(lambda), y="Accuracy", col="", type="")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "#f0f0f0"),
          strip.text = element_text(face = "bold"))



cowplot::plot_grid(
    cowplot::plot_grid(gg0, ggplot()+theme_minimal()),
    gg, nrow=2, rel_heights=c(1,2), labels=c("a", "b")
)

