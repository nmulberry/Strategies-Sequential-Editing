source("setup.R")

##--------------------##
## LOAD DATA IF AVAIL
##--------------------##
load_res <- TRUE
read_test_results <- function(file, dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	model <- name[[1]][1]
	m <- name[[1]][2]
	j <- name[[1]][3]
	k <- name[[1]][4]
	k <- strsplit(k, split="\\.")[[1]][1]
	df$model <- model
	df$m <- as.numeric(m)
	df$j <- as.numeric(j)
	df$k <- as.numeric(k)
	return(df)
}




if (load_res){
    files <- list.files("output/test_runs")

    if (length(files) > 0){
        res0 <- map_dfr(files, read_test_results, dir="output/test_runs")
    } else {
        print("no results")
    }
}


res <- res0
res$lambda <- res$lambda1


#res <- res %>% group_by(model) %>% mutate(ell=min(ell))%>%ungroup()

# Get p0
res$d <- 1-2*res$ell
res$q <- 1/res$j
res <- res %>% mutate(p0 = pmap_dbl(., ~ptrip_d2(..13,..1,..8,..14,..3,..11,..15)))




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
        Objective=="score" ~ "UPGMA", 
        Objective=="sim_dist"~"RF", 
        TRUE ~ "P0"))
res3$k <- factor(res3$k)
res3$m <- factor(res3$m)



## FULL PLOTS
gg_all_A <- ggplot(filter(res3, model=="A",j==64), aes(x=lambda, y=f50, col=k, group=k))+
	geom_line(linewidth=1.25)+
	geom_ribbon(aes(x=lambda, ymin=fmin, ymax=fmax, fill=k), col=NA, alpha=0.2)+
	facet_grid(Objective~m, scales="free_y", labeller=labeller(.cols=label_both))+
	labs(x=expression(lambda),y="")+
    theme(legend.position="bottom")+
    scale_colour_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Dark2")
ggsave("tests_modelA_j64.pdf")



gg_all_B <- ggplot(filter(res3,model=="B", j==64), aes(x=lambda, y=f50, col=k, group=k))+
	geom_line(linewidth=1.25)+
	geom_ribbon(aes(x=lambda, ymin=fmin, ymax=fmax, fill=k), col=NA, alpha=0.2)+
	facet_grid(Objective~m, scales="free_y", labeller=labeller(.cols=label_both))+
	labs(x=expression(lambda),y="")+
    theme(legend.position="bottom")+
    scale_colour_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Dark2")
ggsave("tests_modelB_j64.pdf")

gg_all_C <- ggplot(filter(res3,model=="C",j==64), aes(x=lambda, y=f50, col=k, group=k))+
	geom_line(linewidth=1.25)+
	geom_ribbon(aes(x=lambda, ymin=fmin, ymax=fmax, fill=k), col=NA, alpha=0.2)+
	facet_grid(Objective~m, scales="free_y", labeller=labeller(.cols=label_both))+
	labs(x=expression(lambda),y="")+
    theme(legend.position="bottom")+
    scale_colour_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Dark2")
ggsave("tests_modelC_j64.pdf")



####heatmap?
## higher res 

gg_trips_grid <- ggplot(filter(res3, j==64, model=="C", Objective=="P0", k==5), aes(x=lambda, y=m, fill=f50, col=f50))+
	geom_raster()+
	scale_fill_viridis(option="plasma", breaks=c(0,0.5, 0.9))+scale_colour_viridis(option="plasma", breaks=c(0,0.5,0.9))+
	labs(x=expression(lambda), y="# Barcodes (m)",col="P0", fill="P0")+
	scale_x_continuous(expand=c(0,0))+
	scale_y_discrete(expand=c(0,0))+
 	theme(
        legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(0,0,0,18),
		legend.key.height=unit(0.8,"cm"))


gg_trips_grid <- gg_trips_grid + guides(fill = guide_colourbar(ticks=FALSE, 
                                       frame.colour = "black",
                                       frame.linewidth = 1))





gg_upgma <- ggplot(filter(res3, j==64,model=="C", m %in% c(10,20,50), Objective=="RF", k==5),
	aes(x=lambda, y=f50/100, col=m, group=m))+
	geom_line(linewidth=1.25)+
	geom_ribbon(aes(x=lambda, ymin=fmin/100, ymax=fmax/100, fill=m), col=NA, alpha=0.2)+
#	facet_grid(.~k, labeller=label_both)+
	scale_colour_brewer(palette="Dark2")+
	scale_fill_brewer(palette="Dark2")+
	labs(y="RF distance", x=expression(lambda), fill="m", col="m")+theme(legend.position="right")


## BARS

res2 <- filter(res, j==64)



res_upgma <- res2 %>%
	group_by(k,n,j,m,lambda,p0,model,ell) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.975) %>%
	group_by(k,n,j,m,model,ell) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))%>%
	mutate(ell=round(ell,3))

res_p0 <- res2 %>%
	group_by(k,n,j,m,model,ell) %>%
	filter(p0 >= 0.975) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))%>%
	mutate(ell=round(ell,3))

res_upgma$type <- "UPGMA"
res_p0$type <- "P0"

## also get lambda which minimises p0
res2 <- res2 %>% group_by(k, model, m) %>% filter(p0 == min(p0))%>% ungroup()




gg3 <- ggplot() +
  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=3) +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=3, width=0, alpha=0.5) +
  
  facet_grid(k~ell+n, labeller=label_both) +
  labs(x=expression(lambda), col="", shape="") +
  scale_color_manual(values=c("red", "black"))+
	theme(legend.position="bottom")+
	ggtitle("95% Reconstruction Accuracy")

## put together

library(cowplot)
plot_grid(plot_grid(gg_trips_grid, gg_upgma, nrow=2,labels=c("A","B"), rel_heights=c(1,1.2)),
	gg3, nrow=1, labels=c("","C"),rel_widths=c(1,1.5))

ggsave("test-results-full.pdf")








