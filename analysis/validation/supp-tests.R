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



#res <- filter(res0, m < 70)

res <- res0
res$lambda <- res$lambda1


# Get p0

res <- res %>% group_by(n) %>% mutate(ell=min(ell)) %>% ungroup()%>% 
    mutate(d = 1-2*ell, n= 2^n, q=(1/j))  



#res <- readRDS("uniform_simulations.RDS")

# count prop trips exact
res_trunc <- res %>%
	group_by(lambda,k,ell,m,n,q,model,j) %>%
	summarize(score = sum(sim_dist == 0)/n())
res_trunc$dmax <- 1-res_trunc$ell
res_trunc$p0 <- pmap_dbl(res_trunc, ~pfull_0(..1,..2,..3,..4,..5,..6, ..10))
res_trunc$p1 <- pmap_dbl(res_trunc, ~pfull_1(..1,..2,..3,..4,..5,..6,..10))

res <- merge(res, res_trunc)


### PLOT HEATMAPS & TRAJECTORIES

res4 <- res %>% group_by(k,m,n,j,ell,model,lambda) %>% summarize(f50=mean(sim_dist),
	fmin=min(sim_dist),fmax=max(sim_dist))
res4$k <- factor(res4$k)	
#res4 <- filter(res4, lambda <= 20)
res4 <- res4 %>% mutate(j_lab = case_when(j==4 ~ "Low diversity", j==64 ~ "High diversity"))%>% mutate(m_lab="")
gg_traj <- ggplot(filter(res4,m==10,n==128, k !=13, j != 16),
    aes(x=lambda, y=f50/1000, col=k, group=k))+
    geom_line(linewidth=1.25)+ scale_y_continuous(breaks=c(0, 0.125,0.25),labels=c("", "",""))+
    geom_ribbon(aes(x=lambda,ymin=fmin/1000,ymax=fmax/1000, fill=k), col=NA, alpha=0.2)+
    facet_grid(j_lab~m_lab)+
    labs(x=expression("Editing rate (" * lambda * ")"), y="", col="# sites (k)", fill="# sites (k)")+
    scale_fill_brewer(palette="Dark2")+
    scale_colour_brewer(palette="Dark2")+
	theme(legend.position="inside", legend.position.inside=c(0.17,0.9))+
	ggtitle("Simulated RF Distance")


res_cont <- filter(res, m==10)
res_cont$score <- 0
res_cont$p0 <- 0
res_cont$p1 <- 0
res_cont$m <- -10
res_cont <- rbind(res, res_cont)

res_cont <- res_cont %>% mutate(j_lab=case_when(j==4 ~ "Low diversity", j==64 ~ "High diversity"))%>%
	mutate(k_lab = case_when(k==5~ "5 sites", k==9 ~ "9 sites"))


gg_trips_grid <- ggplot(filter(res_cont,n==128, k!=13, j!= 16), aes(x=lambda, y=m, fill=score, col=score))+
	geom_raster()+
	scale_fill_viridis(option="plasma", breaks=c(0,0.5, 1.0))+scale_colour_viridis(option="plasma", breaks=c(0,0.5,1.0))+
	stat_contour(aes(z=p0),col="lightblue", linewidth=1.25, breaks=c(0.9))+
	stat_contour(aes(z=p1),col="red", linewidth=1.25, breaks=c(0.9))+
	labs(x=expression("Editing rate (" * lambda * ")"), y="Tape copies per cell (m)",col="", fill="")+
	scale_x_continuous(expand=c(0,0))+
	scale_y_continuous(expand=c(0,0))+
	facet_grid(j_lab~k_lab)+
 	theme(
		legend.position = "right",
        legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(0,0,0,10),
		legend.key.height=unit(1,"cm"))+
	ggtitle("Simulated vs Theoretical Accuracy")

gg_trips_grid <- gg_trips_grid + guides(fill = guide_colourbar(ticks.colour="black", 
                                       frame.colour = "black",
                                       frame.linewidth = 1))



res2 <- res %>% 
	pivot_longer(cols=c("p0", "p1")) %>%
	mutate(name=case_when(name=="p0" ~ "B0", name=="p1"~"B1"))

gg_test <- ggplot(filter(res2, n==128, m==50), aes(x=lambda, y=value, col=name))+
	geom_line(linewidth=1.2)+
	facet_grid(k~j, labeller=label_both)+
	scale_colour_manual(values=c("lightblue", "red"))+
	labs(x=expression(lambda), y="reconstruction prob", col="Predicted 0.9 Accuracy")+
	theme(legend.position="top")

ggsave("test-legend.pdf")
library(cowplot)
legend2 <- cowplot::get_legend(gg_test)

plot_grid(gg_trips_grid, 
	gg_traj,
nrow=1,labels=c("a","b"), rel_widths=c(2,1))

ggsave("uniform-results.png",height=7, width=10)














## BARS
res <- filter(res, j %in% c(4,64))
res <- res %>% mutate(k= case_when(k==5~ "5 sites", k==9~"9 sites", k==13~"13 sites"),
	j = case_when(j==4~"Low diversity", j==64~"High diversity"),
	n=case_when(n==128~"128 tips", n==256 ~ "256 tips", n==512~"512 tips", n==1024~"1024 tips"))
res$n <- factor(res$n, levels=c("128 tips", "256 tips", "512 tips", "1024 tips"))
res$k <- factor(res$k, levels=c("5 sites", "9 sites", "13 sites"))


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
  facet_grid(n~j+k) +
  scale_x_continuous(limits=c(1,30), breaks=c(5,10,15,20,25,30), labels=c("5","","15","","","30"))+
  labs(x=expression("Editing rate (" * lambda * ")"), col="", shape="",
	y = "Tape copies per cell (m)") +
  scale_color_manual(values=c("lightblue","red", "black"))+
	theme(legend.position="bottom")

ggsave("uniform_bars.png", height=6, width=9)












