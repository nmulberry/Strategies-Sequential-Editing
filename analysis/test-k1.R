source("setup.R")
##--------------------##
## LOAD DATA IF AVAIL
##--------------------##
load_res <- TRUE
read_unif_results <- function(file, dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	model <- name[[1]][1]
	n <- name[[1]][2]
	m <- name[[1]][3]
	j <- name[[1]][4]
	k <- name[[1]][5]
	df$n <- as.numeric(n)
	df$m <- as.numeric(m)
	df$j <- as.numeric(j)
	df$k <- as.numeric(k)
    df$sim_dist <- as.numeric(df$sim_dist)
    df$model <- model
    return(df)
}


if (load_res){
    files <- list.files("output/test_k1")
    if (length(files) > 0){
        res0 <- map_dfr(files, read_unif_results, dir="output/test_k1")
    } else {
        print("no results")
    }
}



#res <- filter(res0, m < 70)

res <- res0
res$lambda <- res$lambda1


# Get p0

res <- res %>% group_by(n) %>% mutate(ell=min(ell)) %>% ungroup()%>% 
    mutate(d = 1-ell, n= 2^n, q=(1/j))  



#res <- readRDS("uniform_simulations.RDS")

# count prop trips exact
res_trunc <- res %>%
	group_by(lambda,k,ell,m,n,q,model,j) %>%
	summarize(score = sum(sim_dist == 0)/n())
res_trunc$dmax <- 1-res_trunc$ell
res_trunc$p0 <- pmap_dbl(res_trunc, ~pfull_0(..1,..2,..3,..4,..5,..6,..10))
res_trunc$p1 <- pmap_dbl(res_trunc, ~pfull_1(..1,..2,..3,..4,..5,..6,..10))

res <- merge(res, res_trunc)


### PLOT HEATMAPS & TRAJECTORIES

res4 <- res %>% group_by(k,m,n,j,ell,model,lambda) %>% summarize(f50=mean(sim_dist),
	fmin=min(sim_dist),fmax=max(sim_dist))
res4$m <- factor(res4$m)	

gg_traj <- ggplot(filter(res4),
    aes(x=lambda, y=f50, col=m, group=m))+
    geom_line(linewidth=1.25)+ 
    geom_ribbon(aes(x=lambda,ymin=fmin,ymax=fmax, fill=m), col=NA, alpha=0.2)+
    facet_grid(ell~j, labeller=label_both)+
    labs(x=expression(lambda), y="Simulated RF Distance")+
    scale_fill_brewer(palette="Dark2")+
    scale_colour_brewer(palette="Dark2")+
	theme(legend.margin=margin(0,0,0,20))



res2 <- res %>% 
	pivot_longer(cols=c("p0", "p1", "score")) %>%
	mutate(name=case_when(name=="p0" ~ "B0", name=="p1"~"B1", name=="score"~"Simulated"))

gg_test <- ggplot(res2, aes(x=lambda, y=value, col=name, linetype=name))+
	geom_line(linewidth=1.2)+
	facet_grid(m~j, labeller=label_both)+
	scale_colour_manual(values=c("lightblue", "red", "black"))+
    scale_linetype_manual(values=c("dashed", "dashed", "solid"))+
	labs(x=expression(lambda), y="Reconstruction probability", col="", linetype="")+
	theme(legend.position="bottom")

ggsave("k1-results.pdf",height=8, width=10)



## BARS
res <- filter(res, j %in% c(4,64))
res_upgma <- res %>%
	group_by(k,n,j,m,lambda,model) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.9) %>%
	group_by(k,n,j,m,model) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_p0 <- res %>%
	group_by(k,n,j,m,model) %>%
	filter(p0 >= 0.9) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))






