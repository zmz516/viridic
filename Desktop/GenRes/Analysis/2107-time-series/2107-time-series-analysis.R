# load packages -----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)


# import data and summarise ------------------------------------------------

pfu <- read_csv("data-raw/pfu.csv")

summary <- pfu |> 
  group_by(culture, time) |> 
  summarise(mean_pfu = mean(PFU, na.rm = T), .groups = "drop")


# plotting ----------------------------------------------------------------

plot <- ggplot(summary, aes(x = time, y = mean_pfu, colour = culture)) +
  geom_line()+
  geom_point()+
  scale_y_log10(labels = label_scientific())+
  labs(x = "Time (hours)",
       y = "PFU/mL",
       title = "Lysogen PFU Time Series",
       colour = "Replicate")+
  theme_bw()

ggsave("plots/plot.png", plot = plot, width = 7, height = 5, dpi = 300)
