source("setup.R")

## compare approx & full prob
# wrt q, m, k,...(fix lambda?max over lambda?)

pars <- crossing(lambda=c(4,5,6,7), k=c(5,7), ell = 0.1, d=0.8, 
    m = c(3,5,7), n = 2, q=seq(0,0.1,by=0.005)) 
res <- pars
res$eps <- pars%>% pmap_dbl(eps)
res <- res %>% group_by(k,m,n,q) %>% summarize(eps = max(eps))

res$m<-factor(res$m)
ggplot(res, aes(x=q, y=eps, linetype=m))+
    geom_line()+
    geom_line(aes(x=q, y=q))+
    facet_grid(.~k, labeller=label_both) 
