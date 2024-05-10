library(ape)
library(TreeSimGM)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(forcats)
library(ggplot2); 
theme_set(theme_bw(base_size = 14))
library(Rcpp)
sourceCpp("funcs.cpp")
source("model.R")
library(phytools)
library(viridis)
library(phangorn)
fig_dir <- "../report/Figures"
save_figs <- TRUE
asynchronous <- TRUE
load_res <- TRUE

all_chars <- c("a","b","c","d","e","f","g","h","i","j",
	"k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v",
	"w", "x", "y", "z", "$", "%","@", "&","(", ")","#",
	"!","^","?","<",">","/","{","}","[","]")

##--------------------##
## LOAD DATA IF AVAIL
##--------------------##
if (load_res){
    files <- list.files("output")

    if (length(files) > 0){
        res0 <- map_dfr(files, read_results, dir="output")
    } else {
        print("no results")
    }
}
