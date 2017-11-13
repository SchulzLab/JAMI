benchmark <- read.csv("/local/home/mlist/Projects/JAMI/data/benchmark.csv")
View(benchmark)
benchmark$label <- paste(benchmark$Genes,"genes with",
                         prettyNum(benchmark$Interactions, big.mark = ","), "triplets")

library(ggplot2)
ggplot(benchmark,
       aes(x = Threads,
           y = Time / 60,
           color = as.factor(Software),
           shape = as.factor(Software))) +
    geom_line() +
    theme_bw(base_size = 18) +
    facet_wrap(~label) +
    geom_point(size = 5) +
    labs(shape = "Software", color = "Software",
         x = "Number of threads",
         y = "Time in minutes") +
    scale_color_manual(values = c("#FF8C00", "#8A2BE2")) +
    scale_x_continuous(trans = "log2", breaks = c(c(1,2,4,8,16)))

