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



# Get p0
res$q <- 1/res$j

## for now..
res <- filter(res, model != "A")


# count prop trips exact
res <- res %>%
	group_by(k, lambda, ell, n, j, m,q,model) %>%
	summarize(score = sum(sim_dist == 0)/n()) 

res$p0 <- res %>% pmap_dbl(., ~ptrip_min(..2,..1,..3,..6,..4,..7))

res2 <- res


res_upgma <- res2 %>%
	dplyr::filter(score >= 0.95) %>%
	group_by(k,n,j,m,model,ell) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))%>%
	mutate(ell=round(ell,3))

res_p0 <- res2 %>%
	group_by(k,n,j,m,model,ell) %>%
	filter(p0 >= 0.95) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))%>%
	mutate(ell=round(ell,3))

res_upgma$type <- "UPGMA"
res_p0$type <- "P0"

## also get lambda which minimises p0
res2 <- res2 %>% group_by(k, model, m) %>% filter(p0 == min(p0))%>% ungroup()



gg3 <- ggplot() +
  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=1.5) +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=3, width=0, alpha=0.5) +
  facet_grid(ell+n~k+j, labeller=label_both) +
  labs(x=expression(lambda), col="", shape="") +
  scale_color_manual(values=c("red", "black"))+
	theme(legend.position="bottom")+
	scale_x_continuous(limits=c(1,20))

ggsave("synch-asynch-comparison.pdf")






