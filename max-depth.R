

read_results <- function(file,dir){
	df <- readRDS(paste0(dir, "/", file))
    return(df)
}
files <- list.files("output/max-ell-d")

res0 <- map_dfr(files, read_results, dir="output/max-ell-d")
res0 <- res0 %>% dplyr::filter(m != 5)
## correlation between trip score and ell (not n)
ggplot(res0, aes(x=ell, y=max_dist, col=n))+
    geom_point()+
    facet_grid(n~m)




## predicted max dist
pred_p0 <- function(ell,m,k,lambda,eps){
    d <- 1-2*ell
    p0 <- ptrip_d(k, k,ell,d,m)
    return(p0)

}
eps <- -13

res <- res0
res$p0 <- res %>% pmap_dbl(., ~ pred_p0(..9,..3,..1,..1,eps))
