#source("setup.R")
res <- res0 %>% filter(model=="synchr")
## get pars from runs
lambda <- unique(res$lambda)
k <- as.numeric(unique(res$k))
n <- as.numeric(unique(res$n))
j <- as.numeric(unique(res$j))


res$q <- 1/as.numeric(res$j)
## Get expected trips and rho_0
pars2 <- crossing(lambda=lambda, k=k, n=n, q=1/j)
pars2$E0 <- pars2 %>% pmap_dbl(min_expected_trips) 

res <- res %>% mutate(n = as.numeric(n), m=as.numeric(m), j=as.numeric(j))
res2 <- merge(res,pars2)



res2 <- res2 %>%
    mutate(En = m*E0,
    p0 = pmap_dbl(., ~ get_resolvable_trips(..1,..2,..3,..6,..4, eta=0995)))
#	p0 = exp(-m*E0^2/(2*k^2)))


##=======================##
## PLOT RATES ----------##
##======================##


res_trips <- res2 %>%
	group_by(k,n,j,m,lambda, En) %>%
	summarize(p_trips = sum(triplets==0)/n()) %>%
	dplyr::filter(p_trips >= 0.95) %>%
	group_by(k,n,j,m) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))	
res_upgma <- res2 %>%
	group_by(k,n,j,m,lambda, En) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.95) %>%
	group_by(k,n,j,m) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))
res_E <- res2 %>%
	group_by(k,n,j,m) %>%
	filter(En == max(En))
res_p <- res2 %>%
	group_by(k,n,j,m) %>%
	filter(p0 == min(p0))


res_trips$type <- "Triplets"
res_upgma$type <- "UPGMA"
res_E$type <- "E[Delta]"
res_p$type <- "p[0]"


gg_bars <- ggplot() +
  geom_errorbar(data=filter(res_upgma,n==10), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=3) +
  geom_errorbar(data=filter(res_trips,n==10), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=3, width=0, alpha=0.5) +
  facet_grid(k ~ j, labeller=label_both) +
  labs(x=expression(lambda), col="", shape="") +
  scale_color_manual(values=c("red", "black"))+
	theme(legend.position="bottom")+
	geom_point(data=filter(res_E, n==10),
		aes(x=lambda, y=m), shape=4)


gg2_bars <- ggplot() +
  geom_errorbar(data=filter(res_upgma,n==10, j%in% c(5,20,40)), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type), width=3) +
  geom_errorbar(data=filter(res_trips,n==10, j%in%c(5,20,40)), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),linewidth=3, width=0, alpha=0.5) +
  facet_grid(k ~ j, labeller=label_both) +
  labs(x=expression(lambda), col="", shape="", y = "# Barcodes (m)") +
  scale_color_manual(values=c("red", "black"))+
	theme(legend.position="top")+
	geom_point(data=filter(res_E,j%in%c(5,20,40),n==10),
		aes(x=lambda, y=m), shape=4)

### save
#ggsave(paste0(fig_dir, "/opt-rates-bars-grid.pdf"), height=9, width=11)
