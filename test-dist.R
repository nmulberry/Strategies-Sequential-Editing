source("setup.R")


lambda <- 5
k <- 5
ell <- 0.1
d <- 0.8
x <- seq(0,k)
p_in <- sapply(x, prob_seq_edits_in_grp,
    d, d+ell, lambda, k, 0)

# plot
plot(x,p_in)
