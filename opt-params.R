source("setup.R")

k <- c(5,7,9,12)
ell <- c(0.01,0.1)
eps <- 0.9
n <- c(10,100,1000,2000)
q <- c(0,0.03, 0.1)

df <- crossing(k=k,ell=ell,eps=eps,n=n,q=q)

get_m_over_pars <- function(k,ell,eps,n,q){
#######################
## first find opt lambda, m for p0
# should be close to opt lambda as long as q relatively small
res <- opt_p0_fix_d(k,ell,eps,n)
m0 <- res$m
## now look trip bound
# first test
m_max <- 1000
p <- pfull_1(res$lambda, k, ell, res$m, n, q)
pend <- pfull_1(res$lambda, k, ell, m_max, n, q)
if (p > eps){m1 <- res$m}
else if (pend < eps) {m1 <- NA}
else {
    # use root finder
    f <- function(m,lambda,k,ell,n,q){
        return(pfull_1(lambda,k,ell,m,n,q)-eps)
    }
    r <- uniroot(f, c(res$m, m_max), 
        lambda=res$lambda, k=k, ell=ell, n=n, q=q)

    m1 <- ceiling(r$root)
}

return(data.frame(k=k,ell=ell, n=n, q=q, eps=eps, m=c(m0,m1), type=c("B0", "B1")))

}


res <- df %>% pmap_dfr(get_m_over_pars)

res$q <- factor(res$q)
gg1 <- ggplot(res, aes(x=n, y=m, group=interaction(type,q),col=type, linetype=q))+
    geom_line()+geom_point()+scale_x_continuous(breaks=c(0,1000,2000))+
    labs(x="Number tips (n)", y="Predicted # barcodes for 90% accuracy", col="")+
    facet_grid(ell~k, labeller=label_both, scales="free")+
    scale_colour_manual(values=c("B1"="red", "B0"="lightblue"))+
    theme(legend.position="bottom")

ggsave("min_m_eps_0.9.pdf", width=10)
