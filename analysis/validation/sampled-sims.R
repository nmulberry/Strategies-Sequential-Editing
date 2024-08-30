##--------------------##
## LOAD DATA IF AVAIL
##--------------------##
load_res <- TRUE
read_samp_results <- function(file, dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	model <- name[[1]][1]
	n_gen <- name[[1]][2]
	m <- name[[1]][3]
	j <- name[[1]][4]
	k <- name[[1]][5]
    df$n_gen <- as.numeric(n_gen)
    df$m <- as.numeric(m)
	df$j <- as.numeric(j)
	df$k <- as.numeric(k)
    df$model <- model
    return(df)
}




if (load_res){
    files <- list.files("output/test_samp")

    if (length(files) > 0){
        res0 <- map_dfr(files, read_samp_results, dir="output/test_samp")
    } else {
        print("no results")
    }
}



#res <- filter(res0, m < 70)

res <- res0
res$lambda <- res$lambda1


# Get p0

res <- res %>% group_by(n_gen) %>% mutate(ell=1/(n_gen+1), n= max(n)) %>% ungroup() %>%  
    mutate(d = 1-ell, q=(1/j))  



#res <- readRDS("uniform_simulations.RDS")

# count prop trips exact
res_trunc <- res %>%
	group_by(lambda,k,ell,m,n,q,n_gen,j) %>%
	summarize(score = sum(sim_dist == 0)/n())
res_trunc$dmax <- 1-res_trunc$ell
res_trunc$p0 <- pmap_dbl(res_trunc, ~pfull_0(..1,..2,..3,..4,..5,..6, ..10))
res_trunc$p1 <- pmap_dbl(res_trunc, ~pfull_1(..1,..2,..3,..4,..5,..6,..10))

res <- merge(res, res_trunc)


### PLOT HEATMAPS & TRAJECTORIES

res4 <- res %>% group_by(k,m,n,j,ell,n_gen,lambda) %>% summarize(f50=mean(sim_dist),
	fmin=min(sim_dist),fmax=max(sim_dist))
res4$k <- factor(res4$k)	
#res4 <- filter(res4, lambda <= 20)



## BARS
res <- filter(res, j %in% c(4,64))

res <- res %>% mutate(k= case_when(k==5~ "5 sites", k==9~"9 sites", k==13~"13 sites"),
	j = case_when(j==4~"Low Diversity", j==64~"High Diversity"),
	n_gen=case_when(n_gen == 10 ~ "10 generations", n_gen==11~ "11 generations", n_gen==12~"12 generations"))
res$n_gen <- factor(res$n_gen, levels=c("10 generations", "11 generations", "12 generations"))
res$k <- factor(res$k, levels=c("5 sites", "9 sites", "13 sites"))

res_upgma <- res %>%
	group_by(k,n,j,m,lambda,n_gen,ell) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.9) %>%
	group_by(k,n,j,m,n_gen) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_p0 <- res %>%
	group_by(k,n,j,m,n_gen) %>%
	filter(p0 >= 0.9) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))
res_p1 <- res %>%
	group_by(k,n,j,m,n_gen) %>%
	filter(p1 >= 0.9) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_upgma$type <- "Simulated (UPGMA)"
res_p0$type <- "B0"
res_p1$type <- "B1"




gg3 <- ggplot() +
  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=3) +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=1, width=0, alpha=0.75) +
  geom_errorbar(data=filter(res_p1), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
  facet_grid(n_gen~j+k) +
  scale_x_continuous(limits=c(1,30), breaks=c(5,10,15,20,25,30), labels=c("5","","15","","","30"))+
  labs(x=expression("Editing rate (" * lambda * ")"), col="", shape="",
	y = "Tape copies per cell (m)") +
  scale_color_manual(values=c("lightblue","red", "black"))+
	theme(legend.position="bottom")

ggsave("sampled_bars.pdf", height=6, width=9)












