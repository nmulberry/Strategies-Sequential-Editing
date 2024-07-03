source("setup.R")

lambda <- 5
k <- 5
du <- 0.2
dv <- 0.6
x <- 0:k
p_x <- sapply(x, prob_seq_edits_in_grp, du, dv, lambda, k, 0)

test_x <- sapply(x, trunc_poi, k, )

