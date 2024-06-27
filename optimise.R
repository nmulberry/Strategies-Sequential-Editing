source("setup.R")

## --- PART1 : OPTIMISEi
m <- 30
k <- c(5,7,9)
ell <- c(0.1, 0.03, 0.01)
q <- c(0.01, 0.05, 0.1, 0.3)
lambda <- seq(1, 20, by=1)
pars <- crossing(lambda=lambda, k=k,ell=ell, m=m,q=q)
res0 <- pars
## max wrt d
run_over_d <- function(lambda,k,ell,m,q){
	d <- seq(0, 1-2*ell, by=ell)
	pars <- data.frame(lambda=lambda,k=k,ell=ell,d=d, m=m,q=q)
	pars$p <- pars %>% pmap_dbl(prob_trip_approx)
	return(pars)
}


opt_over_d <- function(lambda,k,ell,m,q){
	res <- optim(0.5, prob_trip_approx_d, lambda=lambda, k=k, ell=ell,
	 m=m, q=q, method="L-BFGS", lower=0, upper=1-2*ell, control=list(fnscale=-1))	
	return(data.frame(d=res$par, p=res$value, k=k,m=m,ell=ell, q=q, lambda=lambda))
}



res <- pars %>% pmap_dfr(run_over_d)

## look at p(d)
res1 <- filter(res, ell==0.1)

gg1 <- ggplot(res1, aes(x=lambda, y=d, fill=p,col=p))+
	geom_tile()+
	facet_grid(k~q, labeller=label_both)+
	scale_fill_viridis()+scale_colour_viridis()

ggsave("optim_over_d_ell_0.1.pdf")
res2 <- filter(res, ell==0.01)

gg1 <- ggplot(res2, aes(x=lambda, y=d, fill=p,col=p))+
	geom_tile()+
	facet_grid(k~q, labeller=label_both)+
	scale_fill_viridis()+scale_colour_viridis()+
	labs(x=expression(lambda))+scale_x_continuous(expand=c(0,0))+
	scale_y_continuous(expand=c(0,0))+ 	
	theme(
       legend.key.width = unit(0.4, "cm"),
       legend.margin = margin(0,0,0,18),
       legend.key.height=unit(0.8,"cm"))

gg1 <- gg1 + guides(fill = guide_colourbar(ticks=FALSE, 
                                       frame.colour = "black",
                                       frame.linewidth = 1))


ggsave("optim_over_d_ell_0.01.pdf")

## look at max d wrt lambda,q,k (m? shouldn't matter?)
res3 <- res %>% group_by(ell,lambda,k,m,q) %>% 
	filter(p==max(p)) %>%
	filter(d==max(d)) %>% ungroup() #just one obs per grp
res3$k <- factor(res3$k)
gg3 <- ggplot(res3, aes(x=lambda, y=d, col=k, group=k))+
	geom_line(linewidth=1.2)+
	facet_grid(ell~q, labeller=label_both)+
	scale_colour_brewer(palette="Set2")+
	labs(x=expression(lambda))

ggsave("optim_over_d_summ.pdf")

## check results from optim
res_test <- pars %>% pmap_dfr(opt_over_d)
res_test$k <- factor(res_test$k)
gg3 <- ggplot(res_test, aes(x=lambda, y=d, col=k, group=k))+
	geom_line(linewidth=1.2)+
	facet_grid(ell~q, labeller=label_both)+
	scale_colour_brewer(palette="Set2")+
	labs(x=expression(lambda))



if(FALSE){
#### look at full prob
m <- 5
k <- c(5)
ell <- c(0.01)
q <- c(0.8)
lambda <- seq(1, 10, by=1)
pars <- crossing(lambda=lambda, k=k,ell=ell, m=m,q=q)
res0 <- pars
## max wrt d
opt_over_d2 <- function(lambda,k,ell,m,q){
	d <- seq(0, 1-2*ell, by=4*ell)
	pars0 <- data.frame(lambda=lambda,k=k,ell=ell,d=d, m=m,q=q)
	pars <- pars0
	pars$pfull <- pars0 %>% pmap_dbl(prob_tripR_full)
	pars$papprox <- pars0 %>% pmap_dbl(prob_trip_approx)
	return(pars)
}


res <- pars %>% pmap_dfr(opt_over_d2)
res <- res %>% pivot_longer(cols=c("pfull", "papprox"))


## look at p(d)
res1 <- filter(res, ell==0.01)

gg4 <- ggplot(res1, aes(x=lambda, y=d, fill=value,col=value))+
	geom_tile()+
	facet_grid(.~name, labeller=label_both)+
	scale_fill_viridis()+scale_colour_viridis()

ggsave("optim_over_d_full.pdf")}
