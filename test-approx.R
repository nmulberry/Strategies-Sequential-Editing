source("setup.R")
lambda <- seq(1, 30, by=1)
k <- c(5)

q <- seq(0, 1, by=0.2)
m <- c(10,50,100)

d <- 0.8
ell <- 0.1

df <- crossing(lambda=lambda, k=k, ell=ell, d=d, m=m,q=q)

res <- df
res$p_approx <- df %>% pmap_dbl(prob_trip_approx)
res$p0 <- df %>% pmap_dbl(prob_zero)
res <- res %>% pivot_longer(cols=c("p_approx", "p0"))

test <- ggplot(res, 
    aes(x=lambda, y=value, col=name))+
    geom_line()+
    facet_grid(k+m~q, labeller=label_both)+
    labs(x=expression(lambda), y="", col="")+
    scale_colour_brewer(palette="Set1")
ggsave("test-approx1.pdf", height=7,width=10)
## test over different values of k?

## also test against full distribution (small m)
if (FALSE){
q <- c(0.25,0.5,0.75)
m <- c(5,7)
k <- 5

df <- crossing(lambda=lambda, k=k, ell=ell,d=d, m=m, q=q)
res <- df
res$p_approx <- df %>% pmap_dbl(prob_trip_approx)
res$p_full <- df %>% pmap_dbl(prob_tripR_full)
res <- res %>% pivot_longer(cols=c("p_approx", "p_full"))
test2 <- ggplot(res, 
    aes(x=lambda, y=value, col=name))+
    geom_line()+
    facet_grid(k+m~q, labeller=label_both)+
    labs(x=expression(lambda), y="", col="")+
    scale_colour_brewer(palette="Set1")
ggsave("test-approx2.pdf", height=7,width=10)
}

## test when q=0:
k <- c(5,7,9,12)
m <- c(1,10,100,1000)
q <- 0
d <- 0.8
ell <- 0.1

df <- crossing(lambda=lambda, k=k, ell=ell, d=d, m=m,q=q)

res <- df
res$p_approx <- df %>% pmap_dbl(prob_trip_approx)
res$p0 <- df %>% pmap_dbl(prob_zero)
res <- res %>% 
    group_by(k,m,q,d,ell) %>%
    mutate(err = norm(p_approx-p0, type="2"))



res2 <- res %>% pivot_longer(cols=c("p_approx", "p0"))

res2$k <- factor(res2$k)
test3 <- ggplot(res2, 
    aes(x=lambda, y=value, col=k, linetype=name))+
    geom_line(linewidth=1.1)+
    facet_grid(~m, labeller=label_both)+
    labs(x=expression(lambda), y="", col="k", linetype="")+
    scale_colour_brewer(palette="Set1")
ggsave("test-approx-q0.pdf")


test2 <- ggplot(res, aes(x=m, y=err))+
    geom_point()+geom_line()+
    facet_wrap(~k)
