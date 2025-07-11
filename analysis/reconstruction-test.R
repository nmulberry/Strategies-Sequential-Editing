##-----SET-UP-----##
j <- 8
chars <- all_chars[1:j]
q <- 1/j
#---parameter regime---#
k <- c(5,9)
m <- c(10,20,30)
Ng <- 7
n <- 2^Ng # num gen
ell <- 1/(Ng+1.5) # slightly smaller than 1/Ng

##---- GENERATE DATA--------##
# don't run over all rates
lambda <- seq(1,20, by=1)
tree <- generate_test_tree(1,200,Ng)

##--- TEST ON UPGMA, NJ-----## 
nsim <- 20
pars <- crossing(nsim=nsim, k=k, lambda=lambda, m=m, ell=ell, j=j)
true_dists <- cophenetic.phylo(tree)

print("starting simulations")
sim_res <- pars %>%
    pmap_dfr(., simulate_and_score, tree=tree, true_dists=true_dists, chars=chars, method="ALL")
print("done simulations!")
    

sim_res <- sim_res %>%
    pivot_longer(cols=c("upgma", "nj", "triplets"),
    names_to="method", values_to="score")

## SCORE (note: distances not directly comparable across tree reconstruction vs trips)
sim_res <- sim_res %>%
    group_by(lambda=lambda1,m,k,j,method) %>%
    summarize(score=sum(score==0)/n())


gg <- ggplot(sim_res, aes(x=lambda, y=score, col=method, shape=method))+
    geom_point(col="black")+geom_line()+
    facet_grid(m~k, labeller=label_both)+
    labs(x=expression(paste("Editing rate (", lambda, ")")), 
        y="Simulated Accuracy", col="Method")+
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "#f0f0f0"),
                  strip.text = element_text(face = "bold"))+
    scale_colour_manual(values=c("gold", "purple", "orange"))

