source("setup.R")

lambda <- seq(1,30, by=1)
q <- c(0)
n <- c(2^10)
m <- c(10,30,50)
k <- c(5,7,9)
ell <- 0.1
d <- 0.8
pars <- crossing(lambda=lambda, k=k, ell=ell, d=d,m=m,n=n,q=q)

res <- pars
res$p0 <- pars %>% pmap_dbl(pfull_0)
res$p <- pars %>% pmap_dbl(ptrip_d2)
res <- res %>% pivot_longer(cols=c("p0", "p"))


gg <- ggplot(res,
    aes(x=lambda, y=value, linetype=name, col=name))+
    geom_line()+
    facet_grid(k~m, labeller=label_both)

ggsave("test-p.pdf")
