##--- Optimal Rates---##
# read in data

read_results <- function(file,dir){
	df <- readRDS(paste0(dir, "/", file))
	name <- strsplit(file, split="_")
	return(df)
}


if (T){
    files <- list.files("output/test_unif/")

    if (length(files) > 0){
        res0 <- map_dfr(files, read_results, dir="output/test_unif")
    } else {
        print("no results")
    }
}


# compare to bounds

res2 <- res0 %>%
    rename(lambda=lambda1) %>%
    dplyr::select(-c(lambda2)) %>%
    dplyr::filter(j != 16)%>%
    mutate(j_lab = case_when(j==4 ~ "Low diversity", j==64 ~ "High diversity"))

res2$n <- as.numeric(res2$n)

res_upgma <- res2 %>%
	group_by(k,ell,j,m,lambda,n, j_lab) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_p0 <- res2 %>%
	group_by(k,ell,j,m,lambda,n,j_lab) %>%
	nest()%>%
	summarize(B0=pfull_0(lambda,k,ell,m,n,q=1/j,d=1-ell))%>%
	filter(B0 >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab) %>%
	summarize(lambda_min=min(lambda), lambda_max=max(lambda))

res_pq <- res2 %>%
	group_by(k,ell,j,m, lambda,n,j_lab) %>%
	nest() %>%
	summarize(Bq= pfull_1(lambda,k,ell,m,n,q=1/j,d=1-ell))%>%
	filter(Bq >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_trip <- res2 %>%
	group_by(k,ell,j,m,lambda,n, j_lab) %>%
	summarize(p_trip = sum(triplets==0)/n()) %>%
	dplyr::filter(p_trip >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))


res_upgma$type <- "UPGMA"
res_trip$type <- "Simulated triplet score"
res_pq$type <- "Bq"
res_p0$type <- "Binf"


library(ggh4x)
n_labeller <- function(n) {
  paste0("n: ", n, " tips")
}
j_labeller <- function(j_lab) {
  paste0(j_lab, " diversity")
}
k_labeller <- function(k) {
  paste0("k: ", k, " sites")
}

facet_labeller <- labeller(
  n = n_labeller,
  j = j_labeller,
  k = k_labeller)


gg3 <- ggplot() +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
  geom_errorbar(data=filter(res_pq), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
#  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
#                 width=2) +
  geom_errorbar(data=filter(res_trip), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=6) +
  facet_nested(n~j_lab+k, labeller=facet_labeller) +
  theme(strip.placement = "outside",
          strip.background = element_rect(fill = "#f0f0f0"),
                  strip.text = element_text(face = "bold"))+
labs(x=expression(paste("Editing rate (", lambda, ")")), col="90% Accuracy:",
  y="Tape copies per cell (m)") +
  scale_color_manual(values=c("lightblue", "red", "black"))+
	theme(legend.position="bottom")

######################
# SUBSAMPLED TREES



if (T){
    files <- list.files("output/test_samp/")
    if (length(files) > 0){
        res0 <- map_dfr(files, read_results, dir="output/test_samp")
    } else {
        print("no results")
    }
}


# compare to bounds

res2 <- res0 %>%
    rename(lambda=lambda1) %>%
    dplyr::select(-c(lambda2)) %>%
    dplyr::filter(j != 16)%>%
    mutate(j_lab = case_when(j==4 ~ "Low diversity", j==64 ~ "High diversity"))%>%
    mutate(gen=case_when(ell < 0.08 ~ 12, ell < 0.085 ~ 11, TRUE ~ 10))

res2$n <- as.numeric(res2$n)

res_upgma <- res2 %>%
	group_by(k,ell,j,m,lambda,n, j_lab,gen) %>%
	summarize(p_upgma = sum(sim_dist==0)/n()) %>%
	dplyr::filter(p_upgma >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab, gen) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_p0 <- res2 %>%
	group_by(k,ell,j,m,lambda,n,j_lab, gen) %>%
	nest()%>%
	summarize(B0=pfull_0(lambda,k,ell,m,n,q=1/j,d=1-ell))%>%
	filter(B0 >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab, gen) %>%
	summarize(lambda_min=min(lambda), lambda_max=max(lambda))

res_pq <- res2 %>%
	group_by(k,ell,j,m, lambda,n,j_lab, gen) %>%
	nest() %>%
	summarize(Bq= pfull_1(lambda,k,ell,m,n,q=1/j,d=1-ell))%>%
	filter(Bq >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab,gen) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))

res_trip <- res2 %>%
	group_by(k,ell,j,m,lambda,n, j_lab, gen) %>%
	summarize(p_trip = sum(triplets==0)/n()) %>%
	dplyr::filter(p_trip >= 0.95) %>%
	group_by(k,ell,j,m,n, j_lab, gen) %>%
	summarize(lambda_min = min(lambda), lambda_max = max(lambda))


res_upgma$type <- "UPGMA"
res_trip$type <- "Simulated triplet score"
res_pq$type <- "Bq"
res_p0$type <- "Binf"


ng_labeller <- function(n) {
  paste0(n, " generations")
}
j_labeller <- function(j_lab) {
  paste0(j_lab, " diversity")
}
k_labeller <- function(k) {
  paste0("k: ", k, " sites")
}

facet_labeller <- labeller(
  gen = ng_labeller,
  j = j_labeller,
  k = k_labeller)


gg4 <- ggplot() +
  geom_errorbar(data=filter(res_p0), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
  geom_errorbar(data=filter(res_pq), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                linewidth=2, width=0, alpha=0.75) +
#  geom_errorbar(data=filter(res_upgma), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
#                 width=2) +
  geom_errorbar(data=filter(res_trip), aes(xmin=lambda_min, xmax=lambda_max, y=m, color=type),
                width=6) +
  facet_nested(gen~j_lab+k, labeller=facet_labeller) +
  theme(strip.placement = "outside",
          strip.background = element_rect(fill = "#f0f0f0"),
                  strip.text = element_text(face = "bold"))+
labs(x=expression(paste("Editing rate (", lambda, ")")), col="90% Accuracy:",
  y="Tape copies per cell (m)") +
  scale_color_manual(values=c("lightblue", "red", "black"))+
	theme(legend.position="bottom")




