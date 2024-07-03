source("setup.R")
if (FALSE){
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
    mutate(d = 1-2*ell, n= 2^n, q=(1/j))  


}

res <- readRDS("uniform_simulations.RDS")

# count prop trips exact
res_trunc <- res %>%
	group_by(lambda,k,ell,m,n,q,model,j) %>%
	summarize(score = sum(sim_dist == 0)/n())

res_trunc$p0 <- pmap_dbl(res_trunc, ~pfull_0(..1,..2,..3,..4,..5,..6))
res_trunc$p1 <- pmap_dbl(res_trunc, ~pfull_1(..1,..2,..3,..4,..5,..6))
#res_trunc$p2 <- pmap_dbl(res_trunc, ~pfull_2(..1,..2,..3,..4,..5,..6))

res <- merge(res, res_trunc)


### PLOT HEATMAPS & TRAJECTORIES

res4 <- res %>% group_by(k,m,n,j,ell,model,lambda) %>% summarize(f50=mean(sim_dist),
	fmin=min(sim_dist),fmax=max(sim_dist))
res4$k <- factor(res4$k)	
res4 <- filter(res4, lambda <= 20)

gg_traj <- ggplot(filter(res4,m==10,n==256),
    aes(x=lambda, y=f50/1000, col=k, group=k))+
    geom_line(linewidth=1.25)+# scale_y_continuous(labels=c("0.0", "","","","1.0"))+
    geom_ribbon(aes(x=lambda,ymin=fmin/1000,ymax=fmax/1000, fill=k), col=NA, alpha=0.2)+
    facet_grid(m~j, labeller=label_both)+
    labs(x=expression(lambda), y="Simulated RF Distance")+
    scale_fill_brewer(palette="Dark2")+
    scale_colour_brewer(palette="Dark2")+
	theme(legend.margin=margin(0,0,0,20))

gg_trips_grid <- ggplot(filter(res,n==256, lambda <= 20), aes(x=lambda, y=m, fill=score, col=score))+
	geom_raster()+
	scale_fill_viridis(option="plasma", breaks=c(0,0.5, 0.95))+scale_colour_viridis(option="plasma", breaks=c(0,0.5,0.95))+
	stat_contour(aes(z=p0),col="lightblue", linewidth=1.25, breaks=c(0.9))+
	stat_contour(aes(z=p1),col="red", linewidth=1.25, breaks=c(0.9))+
#	stat_contour(aes(z=p2),col="lightblue", linewidth=1.25, breaks=c(0.9))+
	labs(x=expression(lambda), y="# barcodes (m)",col="Sim", fill="Sim")+
	scale_x_continuous(expand=c(0,0))+
	scale_y_continuous(expand=c(0,0))+
	facet_grid(k~j, labeller=label_both)+
 	theme(
        legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(0,0,0,10),
		legend.key.height=unit(0.8,"cm"))


gg_trips_grid <- gg_trips_grid + guides(fill = guide_colourbar(ticks=FALSE, 
                                       frame.colour = "black",
                                       frame.linewidth = 1))



res2 <- res %>% 
	pivot_longer(cols=c("p0", "p1")) %>%
	mutate(name=case_when(name=="p0" ~ "B0", name=="p1"~"B1"))

gg_test <- ggplot(filter(res2, n==256, m==50), aes(x=lambda, y=value, col=name))+
	geom_line(linewidth=1.2)+
	facet_grid(k~j, labeller=label_both)+
	scale_colour_manual(values=c("lightblue", "red"))+
	labs(x=expression(lambda), y="reconstruction prob", col="Predicted Accuracy")+
	theme(legend.position="bottom")


library(cowplot)
legend2 <- cowplot::get_legend(gg_test)

plot_grid(gg_traj, 
	gg_trips_grid,
nrow=2,labels=c("a","b"), rel_heights=c(1,1.8))

ggsave("uniform-results.pdf",height=8, width=10)






## BARS
res <- filter(res, k %in% c(7,9))
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
res_p1 <- res %>%
	group_by(k,n,j,m,model) %>%
	filter(p1 >= 0.9) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_upgma$type <- "Simulated"
res_p0$type <- "B0"
res_p1$type <- "B1"




gg3 <- ggplot() +
  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=3) +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=1, width=0, alpha=0.75) +
  geom_errorbar(data=filter(res_p1), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
  facet_grid(n~k+j, labeller=label_both) +
  scale_x_continuous(limits=c(1,20))+
  labs(x=expression(lambda), col="", shape="") +
  scale_color_manual(values=c("lightblue","red", "black"))+
	theme(legend.position="bottom")

ggsave("uniform_bars.pdf", height=6, width=9)












