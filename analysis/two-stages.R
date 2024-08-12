source("setup.R")
if (FALSE){
k <- 5
j <- 64
ratio <- 0.1
eps <- 0.9
n <- 10
n1 <- 6
m <- 40
nsim <- 10
############################
run_over_alpha <- function(ratio){
tree <- generate_twostage_tree(alpha=1, beta=200, n1=6, n=n, ratio=ratio)
true_dists <- cophenetic.phylo(tree)
ell1 <- get_min_branch(tree)
ell2 <- ell1/ratio #approx
# get lambda1,m1
print(n1*ell1)
res1 <- opt_p0_fix_d(k,ell1, eps, 2^n1, n1*ell1)
res2 <- opt_p0_fix_d(k, ell2, eps, 2^n, 1-ell2)

print(res1)
m1 <- min(res1$m,m)
m2 <- m-m1

return(m1)
}

#lambda_vec_two <- c(rep(res1$lambda,m1),rep(res2$lambda, m2))
#lambda_vec_single <- rep(res1$lambda, m)

#res_two <- simulate_and_score(nsim, k, lambda_vec_two, m, ell1,j,tree, true_dists)
#print("done 2 rate")
#res_single <- simulate_and_score(nsim, k, lambda_vec_single, m, ell1, j, tree, true_dists)
}
###############################################
###############################################
#read output

##--------------------##
## LOAD DATA IF AVAIL
##--------------------##
load_res <- TRUE
read_unif_results <- function(file, dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	ratio <- name[[1]][1]
    df$ratio <- as.numeric(ratio)
    df$n1 <- as.numeric(name[[1]][2])
    df$sim_dist <- as.numeric(df$sim_dist)
    return(df)
}




if (load_res){
    files <- list.files("output/two-stage")
    if (length(files) > 0){
        res0 <- map_dfr(files, read_unif_results, dir="output/two-stage")
    } else {
        print("no results")
    }
}



res <- res0
res <- filter(res, ratio != 0.08, ratio!= 0.09)
res <- res %>%
    group_by(k,m,q,n1,ratio,type) %>%
    mutate(score=sum(sim_dist==0)/n())%>%
    ungroup()

res <- res %>% mutate(type=case_when(type=="single" ~ "One rate", TRUE~"Two rates"),
    n1 = case_when(n1==4~ "Phase shift after 4 generations", TRUE ~ "Phase shift after 6 generations"))
pal <- brewer.pal(n = 4, "Set1")[1:4]
g1 <- ggplot(res, aes(x=ratio, y=sim_dist, group=interaction(type, ratio, n1), fill=type))+
    geom_boxplot()+
    facet_wrap(~n1)+
    labs(y="RF Distance",x=expression("Ratio of early:late divisions (" *alpha* ")"), fill="")+
    scale_fill_manual(values=pal[1:2])


res2 <- res %>% pivot_longer(cols=c("lambda1", "lambda2"))%>%
    mutate(name=case_when(name=="lambda1" ~ "First rate", TRUE~ "Second rate"))
g2 <- ggplot(filter(res2,type=="Two rates"), aes(x=ratio, y=value, fill=name))+
    geom_col(col="black", position="dodge")+facet_wrap(~n1)+
    labs(y=expression(lambda),x=expression("Ratio of early:late divisions (" *alpha* ")"), fill="")+
    scale_fill_manual(values=pal[3:4])+
theme(axis.title.y = element_text(margin = margin(t = 0, r = 14, b = 0, l = 0)))
cowplot::plot_grid(g1,g2, nrow=2, labels=c("a","b"))
ggsave("multistage-growth.png", height=6, width=10)
