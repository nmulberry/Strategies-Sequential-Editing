library(ape)
library(TreeSimGM)
library(TreeSim)
library(stringr)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(forcats)
library(ggplot2); 
theme_set(theme_minimal(base_size = 12)+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))) 
library(Rcpp)
sourceCpp("funcs.cpp")
source("model.R")
library(phytools)
library(viridis)
library(phangorn)
library(RColorBrewer)
#fig_dir <- "../report/Figures"
#save_figs <- TRUE
#asynchronous <- TRUE
#load_res <- TRUE

all_chars <- c("a","b","c","d","e","f","g","h","i","j",
	"k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v",
	"w", "x", "y", "z", "$", "%","@", "&","(", ")","#",
	"!","^","?","<",">","/","{","}","[","]", "A", "B","C", "D", "E",
    "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P",
    "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")

