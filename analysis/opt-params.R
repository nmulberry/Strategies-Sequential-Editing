source("setup.R")

k <- c(5,9,13)
ell <- seq(0.05, 0.1, by=0.01)
#ell <- 0.1079
eps <- 0.9
n <- c(100,1000,5000)
#n <- 256
j <- 1/64
q <- c(0,0.25)

df <- crossing(k=k,ell=ell,eps=eps,n=n,q=q)

get_m_over_pars <- function(k,ell,eps,n,q){
#######################
## first find opt lambda, m for p0
# should be close to opt lambda as long as q relatively small
res <- opt_p0_fix_d(k,ell,eps,n, 1-ell)
m0 <- res$m
## now look trip bound
# first test
m_max <- 1000
p <- pfull_1(res$lambda, k, ell, res$m, n, q,1-ell)

pend <- pfull_1(res$lambda, k, ell, m_max, n, q,1-ell)
if (p > eps){m1 <- res$m}
else if (pend < eps) {m1 <- NA}
else {
    # use root finder
    f <- function(m,lambda,k,ell,n,q){
        return(pfull_1(lambda,k,ell,m,n,q, 1-ell)-eps)
    }
    r <- uniroot(f, c(res$m, m_max), 
        lambda=res$lambda, k=k, ell=ell, n=n, q=q)

    m1 <- ceiling(r$root)
}

return(data.frame(k=k,ell=ell, n=n, q=q, eps=eps, m=c(m0,m1), type=c("B0", "B1")))

}


res <- df %>% pmap_dfr(get_m_over_pars)
res$q <- factor(res$q)

res <- res %>% mutate(q_lab = case_when(q==0 ~ "No chance homoplasy", q==0.25~"Low diversity"))%>%
    mutate(k_lab = case_when(k==5~"5 sites", k==9~"9 sites", k==13 ~ "13 sites"))%>%
    mutate(n_lab = case_when(n==100 ~ "100 tips", n==1000~ "1000 tips", n==5000~ "5000 tips"))

res$k_lab <- factor(res$k_lab, levels=c("5 sites", "9 sites", "13 sites"))
res$n_lab <- factor(res$n_lab, levels=c("100 tips", "1000 tips", "5000 tips"))


gg1 <- ggplot(filter(res, type=="B0", q==0), aes(x=ell, y=m, group=interaction(type,q),col=type, linetype=q))+
    geom_line()+geom_point()+
    labs(x="Resolution", y="Predicted # tapes  for 90% accuracy", col="")+
    facet_grid(n_lab~k_lab)+
    scale_colour_manual(values=c("B1"="red", "B0"="lightblue"))+
    theme(legend.position="bottom")

ggsave("min_m_eps_0.9_largek.pdf", width=10)


gg1 <- ggplot(res, aes(x=ell, y=log10(m), group=interaction(type,q),col=type, linetype=q_lab))+
    geom_line()+geom_point()+
    labs(x="Resolution", y="Predicted # tapes for 90% accuracy (log10)", col="", linetype="")+
    facet_grid(n_lab~k_lab)+
    scale_colour_manual(values=c("B1"="red", "B0"="lightblue"))+
    scale_linetype_manual(values=c("No chance homoplasy" = "solid", "Low diversity" = "dashed"))+
    theme(legend.position="bottom")
ggsave("min_m_eps_0.9.pdf", width=10)
