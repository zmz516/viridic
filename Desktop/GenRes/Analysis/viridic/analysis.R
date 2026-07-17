# load packages -----------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)

# importing ---------------------------------------------------------------
cluster_table <- read_tsv("data-raw/VIRIDIC_cluster_table.tsv")

sim_dist_table <- read_tsv("data-raw/VIRIDIC_sim-dist_table.tsv")


# heatmap -----------------------------------------------------------------
install.packages("pheatmap")
library(pheatmap)

heatmap_matrix <- as.matrix(sim_dist_table[,-1]) # exclude the first column which contains the genome names
rownames(heatmap_matrix) <- sim_dist_table$genome # paste the values from the genome column to the row names of the matrix
