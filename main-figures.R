source("setup.R")
library(cowplot)

## FIGURE 1

source("test.R")
source("lambda-bars.R")

plot_grid(
    plot_grid(ggplot+theme_minimal(),
    plot_grid(gg_trips_grid, gg_upgma,align="v", ncol=1),
    rel_heights=c(1,9,9), ncol=1),
    gg2_bars,
    ncol=2,labels = "AUTO",  rel_widths = c(1, 1.5))
