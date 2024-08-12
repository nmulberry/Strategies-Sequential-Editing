source("setup.R")
lambda <- seq(1, 30, by=1)


## test when q=0:
k <- c(1,5,13)
m <- c(10,100)
q <- 0.0
ell <- c(0.1,0.005)
n <- c(3,100,1000) 

df <- crossing(lambda=lambda, k=k, ell=ell, m=m,n=n,q=q)
df$dmax <- 1-df$ell

res <- df
res$B1 <- df %>% pmap_dbl(pfull_1)
res$B0 <- df %>% pmap_dbl(., ~pfull_0(..1,..2,..3,..4,..5,..6,..7))



res2 <- res %>% pivot_longer(cols=c("B1", "B0"))

res2$k <- factor(res2$k)

test3 <- ggplot(filter(res2), 
    aes(x=lambda, y=value, col=name, linetype=name))+
    geom_line(linewidth=1.1)+ facet_grid(m+ell~n+k, labeller=label_both)+
    labs(x=expression(lambda), y="", col="", linetype="")+
    scale_linetype_manual(values=c("B0"="solid", "B1"="twodash"))+
    scale_colour_manual(values=c("B1"="red", "B0"="lightblue"))+
	theme(legend.position="bottom")
ggsave("test-approx-q0.png", width=10, height=9)

########################
# TEST 2: TRIPLET PROB
########################
nsim <- 100
chars <- all_chars
k <- c(5,13)
lambda <- seq(1,30, by=1)
m <- c(10,30)
i <- 1:nsim
j <- c(4,16,64)
d <- 0.8
ell <- c(0.1, 0.01)
pars0 <- crossing(i=i, lambda=lambda,j=j, m=m,k=k,d=d, ell=ell)
##===========================#
## test diff
run_colour_model_diff <- function(i,lambda,j,m,k,d,ell){
    chars <- all_chars[1:j]
    root <- paste0(rep("0", k), collapse='')
    sample_p <- rep(1/length(chars), length(chars))
    s_in <- 0
    s_outa <- 0
    s_outb <- 0
    for (i in 1:m){
        res0 <- colouring(root, d, lambda, sample_p, chars)
        resc <- colouring(res0, 1.-d, lambda, sample_p, chars)
        resv <- colouring(res0, ell, lambda, sample_p, chars)
        resa <- colouring(resv, 1-d-ell, lambda, sample_p, chars)
        resb <- colouring(resv, 1-d-ell, lambda, sample_p, chars)
        # count
        s_in <- s_in + count_sequential_matches(resa,resb)
        s_outa <- s_outa + count_sequential_matches(resa,resc)
        s_outb <- s_outb + count_sequential_matches(resb, resc)
     }

    return(s_in-s_outa)
}

res <- pars0
res$num_edits <- res %>% pmap_dbl(run_colour_model_diff)

res22 <- res %>%
    group_by(lambda,j,m,k,d,ell)%>%
    summarize(ptrip=sum(num_edits > 0)/n())%>%
    ungroup()%>%
	mutate(ptrip=1-ptrip) #prob not resolve



res22 <- res22 %>% mutate(q=1/j) %>%
mutate(p0 = pmap_dbl(., ~ ptrip_0(..5,..1,..4,..6,..3,3,..8)), 
    p1 = pmap_dbl(., ~ ptrip_1(..5,..1,..4,..6,..3,3,..8)))

res22 <- res22 %>%
    pivot_longer(cols=starts_with("p"))%>%
    mutate(name=case_when(name=="ptrip"~"Simulated", name=="p0"~"B0", name=="p1"~"B1",
	name=="p2"~"p_approx_split"))
    
res22$name <- factor(res22$name, levels=c("Simulated", "B1", "B0"))
ggplot(filter(res22), aes(lambda, y=value, group=name, col=name, linetype=name))+
    geom_line()+
    scale_linetype_manual(values=c(Simulated="solid", "B0"="dashed","B1" = "dashed"))+
	scale_colour_manual(values=c(Simulated="black", "B1"="red", "B0" = "lightblue"))+
    facet_grid(m+j~k+ell, labeller=label_both)+theme(legend.position="bottom")+
    labs(y="Prob not resolve triplet",col="", linetype="", x=expression(lambda))

ggsave("onetrip-sim-all.png", width=10, height=9)


