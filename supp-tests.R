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
    df$model <- model
    return(df)
}




if (load_res){
    files <- list.files("output/test_unif")

    if (length(files) > 0){
        res0 <- map_dfr(files, read_unif_results, dir="output/test_unif")
    } else {
        print("no results")
    }
}



res <- res0
res$lambda <- res$lambda1



# Get p0

res <- res %>% group_by(n) %>% mutate(ell=min(ell)) %>% ungroup()%>% 
    mutate(d = 1-2*ell, n= 2^n, q=(1/j)) %>%  
    mutate(p0 = pmap_dbl(., ~ptrip_d2(..13,..1,..8,..14,..3,..11,..15)))




# count prop trips exact
res <- res %>%
	group_by(k, lambda, n, j, m, model) %>%
	mutate(score = sum(sim_dist == 0)/n()) 


res3 <- res %>% 
    pivot_longer(cols=c("score", "sim_dist", "p0"), 
                names_to="Objective", 
                values_to="value") %>%
	group_by(k,lambda, m,n,j, model, Objective) %>%
    summarize(f50=median(value), fmin=min(value), fmax=max(value))

## rename
res3 <- res3 %>% 
    mutate(Objective=case_when(
        Objective=="score" ~ "Simulations", 
        Objective=="sim_dist"~"RF", 
        TRUE ~ "Theory"))
res3$k <- factor(res3$k)
res3$m <- factor(res3$m)



## BARS
j_star <- 16
res2 <- filter(res, j==j_star)

res_upgma <- res2 %>%
	group_by(k,n,j,m,lambda,p0,model) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.95) %>%
	group_by(k,n,j,m,model) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_p0 <- res2 %>%
	group_by(k,n,j,m,model) %>%
	filter(p0 >= 0.95) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_upgma$type <- "Simulations"
res_p0$type <- "Theory"

## also get lambda which minimises p0
res2 <- res2 %>% group_by(k, model, m,n) %>% filter(p0 == max(p0))




gg3 <- ggplot() +
  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=3) +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=3, width=0, alpha=0.5) +
  geom_point(data=res2, aes(x=lambda, y=m), shape=8, size=2)+  
  facet_grid(k~n, labeller=label_both) +
  labs(x=expression(lambda), col="", shape="") +
  scale_color_manual(values=c("black", "red"))+
	theme(legend.position="bottom")

ggsave(paste0("test_unif_j_", j_star, ".pdf"))













