res <- readRDS("output/two-stage/twostage.RDS")

ggplot(res, aes(x=ratio, y=triplets, col=type, group=interaction(type,ratio)))+
    geom_boxplot()+
    facet_wrap(~m)

