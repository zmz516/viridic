
# Lysogen MOI 0.005 -------------------------------------------------------

# time series data for  infection of 009 host
# with lysogenic wild type phage plasmid
# OD readings for bacterial population 
# PFU for free virions



# load packages -----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)

# # importing CFU ---------------------------------------------------------------
# 
# import_cfu <- read_csv("data-raw/0508_OD.csv")
# 
# import_cfu <- import_cfu |> 
#   mutate(across(3:4, as.numeric))
# 
# # extract the LB blanks
# # normalize to the blanks = od_norm
# 
# blank_od <- import_cfu |> 
#   filter(sample == "LB") |> 
#   summarise(mean_od = mean(od, na.rm = TRUE)) |> 
#   pull(mean_od)
# 
# data <- import_cfu |> 
#   filter(sample !="LB") |> 
#   mutate(od_norm = pmax(0, od - blank_od))
# 
# # OD readings in the plate reader need to be corrected to spectrophotometer readings
# 
#   # function to convert plate reader OD to spectrophotometer OD = od_corr
#   calculate_od_corr <- function(od_norm) {0.178 * (od_norm^2) + 2.40 * od_norm}
#   
#   # apply the function
#   data <- data |> 
#     mutate(od_corr = calculate_od_corr(od_norm))
# 
# # lysate samples were measured directly in the spectrophotometer so do not need to be corrected
# # add these to the data frame
# 
# lysate_od <- read_csv("data-raw/lysate_OD.csv")
# 
# lysate_od <- lysate_od |> 
#   mutate(across(3:4, as.numeric))
# 
# data <- data |> 
#   bind_rows(data, lysate_od)
# 
# # calcualte CFU from the corrected od readings
# 
#   # function to calculate CFU
#   calculate_log_cfu <- function(od_corr) {1.092 * od_corr + 7.732}
# 
#   # apply the functions
#   data <- data |> 
#     mutate(logCFU = calculate_log_cfu(od_corr),
#            cfu = 10^logCFU)
# 
# # summarise ---------------------------------------------------------------
# 
# cfu_sum <- data |> 
#   group_by(sample, time) |> 
#   summarise(mean_cfu = mean(cfu, na.rm = TRUE),
#             cfu_sd = sd(cfu, na.rm = T),
#             cfu_n = n(),
#             cfu_se= cfu_sd / sqrt(cfu_n),
#             .groups = "drop")
# 
# 
# # plotting CFU vs time ----------------------------------------------------
# 
# cfu_time <- ggplot(cfu_sum, aes(x = time, y = mean_cfu, colour = sample))+
#   geom_line()+
#   geom_point()+
#   geom_errorbar(aes(ymin = mean_cfu - cfu_se,
#                     ymax = mean_cfu + cfu_se),
#                 width = 0.2,
#                 colour = "black", 
#                 alpha = 0.5)+
#   scale_y_log10(labels = label_scientific())+
#   scale_colour_manual(name = "Sample",
#                       limits = c("host","lysogen","lysate"),
#                       labels = c("host" = "Host 009",
#                                  "lysogen" = "Lysogen",
#                                  "lysate" = "Lysate"),
#                       values = c("host" = "orange2",
#                                  "lysogen" = "lightblue3",
#                                  "lysate" = "grey"))+
#     scale_x_continuous(breaks = seq(0, max(data$time, na.rm = TRUE), by = 6))+
#   labs(x = "Time (hours)",
#        y = "CFU/mL")+
#   theme_bw()
# cfu_time  
# 
# ggsave("plots/cfu_time.png", plot = cfu_time, width = 7, height = 5, dpi = 300)
# 
# 
# import pfu -----------------------------------------------------------

import_pfu <- read_csv("data-raw/0508_pfu.csv")

import_pfu <- import_pfu |>
  mutate(across(3:5, as.numeric))

# calcualte PFU
calculate_pfu <- function(plaque, dilution) {plaque / (0.01 * dilution)}

pfu <- import_pfu |>
  mutate(pfu = calculate_pfu(plaque, dilution)) |>
  drop_na(plaque, dilution)

# summarise ---------------------------------------------------------------

pfu_sum <- pfu |>
  group_by(sample, time) |>
  summarise(mean_pfu = mean(pfu, na.rm = TRUE),
            pfu_sd = sd(pfu, na.rm = TRUE),
            pfu_n = n(),
            pfu_se = pfu_sd / sqrt(pfu_n),
            .groups = "drop")


# plotting CFU vs TIME ----------------------------------------------------

pfu_time <- pfu_sum |>
  ggplot(aes(x = time, y = mean_pfu, colour = sample))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin = mean_pfu - pfu_se,
                    ymax = mean_pfu + pfu_se),
                width = 0.2,
                colour = "black",
                alpha = 0.5)+
  scale_y_log10(labels = label_scientific())+
  scale_colour_manual(name = "Sample",
                      limits = c("host","lysogen","lysate"),
                      labels = c("host" = "Host 009",
                                 "lysogen" = "Lysogen",
                                 "lysate" = "Lysate"),
                      values = c("host" = "orange2",
                                 "lysogen" = "lightblue3",
                                 "lysate" = "grey"))+
  labs(x = "Time (hours)",
       y = "PFU/mL")+
  theme_bw()
pfu_time

ggsave("plots/pfu_time.png", plot = pfu_time, width = 7, height = 5, dpi = 300)


# pfu vs cfu vs time -------------------------------------------------------------------

# combined <- cfu_sum |> 
#   left_join(pfu_sum, by = c("sample", "time"))
# 
# 
# # scaling the axes and plotting together
# primary_min <- 1e7
# primary_max <- 1e12
# sec_min <- 1e1
# sec_max <- 1e6
# 
# scale_factor <- primary_min / sec_min  
# 
# pfucfuplot <- ggplot(data = combined, aes(x = time))+
#   geom_line(aes(y = mean_cfu, colour = sample, linetype = "CFU"))+
#   geom_point(aes(y = mean_cfu, colour = sample), shape = 16, size = 2)+
#   geom_errorbar(aes(ymin = mean_cfu - cfu_se,
#                     ymax = mean_cfu + cfu_se),
#                 width = 0.2,
#                 colour = "black", 
#                 alpha = 0.5)+
#   
#   
#   geom_line(aes(y = mean_pfu * scale_factor, colour = sample, linetype = "PFU"))+
#   geom_point(aes(y = mean_pfu * scale_factor, colour = sample), shape = 17, size = 2)+
#   geom_errorbar(aes(ymin = (mean_pfu*scale_factor) - (pfu_se*scale_factor),
#                     ymax = (mean_pfu*scale_factor) + (pfu_se*scale_factor)),
#                 width = 0.2,
#                 colour = "black", 
#                 alpha = 0.5)+
#   
#   scale_y_log10(name = "Colony Density (CFU/mL)",
#                 limits = c(primary_min, primary_max),
#                 sec.axis = sec_axis(transform = ~./ scale_factor ,
#                                     name = "Phage Concentration (PFU/mL)"))+
#   labs(x = "Time (Hours)",
#        colour = "Sample")+
#   
#   scale_colour_manual(name = "Sample",
#                       limits = c("host","lysogen","lysate"),
#                       labels = c("host" = "Host 009",
#                                  "lysogen" = "Lysogen",
#                                  "lysate" = "Lysate"),
#                       values = c("host" = "orange3",
#                                  "lysogen" = "lightblue4",
#                                  "lysate" = "grey"))+
#   
#   scale_linetype_manual(name = "Measurement",
#                         values = c("CFU" = "dotdash", "PFU" = "solid"))+
#   theme_bw()
# 
# 
# pfucfuplot 
# ggsave("plots/pfu_cfu.png", plot = pfucfuplot, width = 7, height = 5, dpi = 300)
# 

# it is better to show OD rather than CFU ---------------------------------

import_od<- read_csv("data-raw/0508_OD.csv")

import_od <- import_od |> 
  mutate(across(3:4, as.numeric))

# extract the LB blanks
# normalize to the blanks = od_norm

blank_od <- import_od |> 
  filter(sample == "LB") |> 
  summarise(mean_od = mean(od, na.rm = TRUE)) |> 
  pull(mean_od)

od_data <- import_od |> 
  filter(sample !="LB") |> 
  mutate(od_norm = pmax(0, od - blank_od))

# OD readings in the plate reader need to be corrected to spectrophotometer readings

# function to convert plate reader OD to spectrophotometer OD = od_corr
calculate_od_corr <- function(od_norm) {0.178 * (od_norm^2) + 2.40 * od_norm}

# apply the function
od_data <- od_data |> 
  mutate(od_corr = calculate_od_corr(od_norm))

# lysate od were already normalised in the spectrophotometer, don;t need to be normalised or correlated

lysate_od <- read_csv("data-raw/lysate_OD.csv")

lysate_od <- lysate_od |> 
  rename(od_norm = od_corr) |> 
  mutate(across(3:4, as.numeric))

od_data <- od_data |> 
  bind_rows(od_data, lysate_od) |> 
  mutate(od_corr = if_else(sample == "lysate", od_norm, od_corr))


# summarise  OD ---------------------------------------------------------------

od_sum <- od_data |> 
  group_by(sample, time) |> 
  summarise(mean_od = mean(od_norm, na.rm = TRUE),
            od_sd = sd(od_norm, na.rm = T),
            od_n = n(),
            od_se= od_sd / sqrt(od_n),
            .groups = "drop")


# plotting OD vs time ----------------------------------------------------

# better to plot the plate reader OD as these are more accurate and don't go over 1



OD_time <- ggplot(od_sum, aes(x = time, y = mean_od, colour = sample))+
  geom_line(linewidth = 1.2)+
  geom_point()+
  geom_errorbar(aes(ymin = mean_od - od_se,
                    ymax = mean_od + od_se),
                width = 0.2,
                colour = "black", 
                alpha = 0.5)+
  scale_colour_manual(name = "Sample",
                      limits = c("host","lysogen","lysate"),
                      labels = c("host" = "Host 009",
                                 "lysogen" = "Lysogen",
                                 "lysate" = "Lysate"),
                      values = c("host" = "orange2",
                                 "lysogen" = "lightblue3",
                                 "lysate" = "grey"))+
  scale_x_continuous(breaks = seq(0, max(od_sum$time, na.rm = TRUE), by = 6))+
  labs(x = "Time (hours)",
       y = "OD600")+
  theme_bw()
OD_time  

ggsave("plots/OD_time.png", plot = OD_time, width = 7, height = 5, dpi = 300)



# plotting OD vs PFU vs time ----------------------------------------------

# od_pfu <- od_sum |>
#   left_join(pfu_sum, by = c("sample", "time"))


# # scaling the axes and plotting together
# primary_min <- 1e-3
# primary_max <- 1e2
# sec_min <- 1e1
# sec_max <- 1e6
# 
# scale_factor <- primary_min / sec_min


# od_pfu_plot <- ggplot(data = od_pfu, aes(x = time))+
#   geom_line(aes(y = mean_od, group = sample, linetype = "OD"), linewidth = 1)+
#   geom_point(aes(y = mean_od, colour = sample), size = 2)+
#   geom_errorbar(aes(ymin = mean_od - od_se,
#                     ymax = mean_od + od_se),
#                 width = 0.2,
#                 colour = "black",
#                 alpha = 0.5)+
# 
# 
#   geom_line(aes(y = mean_pfu * scale_factor, colour = sample, linetype = "PFU"), linewidth = 1.2)+
#   geom_point(aes(y = mean_pfu * scale_factor, colour = sample), size = 2)+
#   geom_errorbar(aes(ymin = (mean_pfu*scale_factor) - (pfu_se*scale_factor),
#                     ymax = (mean_pfu*scale_factor) + (pfu_se*scale_factor)),
#                 width = 0.2,
#                 colour = "black",
#                 alpha = 0.5)+
# 
#   scale_y_log10(name = "OD600",
#                 sec.axis = sec_axis(transform = ~./ scale_factor ,
#                                     name = "Phage Concentration (PFU/mL)"))+
#   labs(x = "Time (Hours)",
#        colour = "Sample")+
# 
#   scale_colour_manual(name = "Sample",
#                       limits = c("host","lysogen","lysate"),
#                       labels = c("host" = "Host 009",
#                                  "lysogen" = "Lysogen",
#                                  "lysate" = "Lysate"),
#                       values = c("host" = "red",
#                                  "lysogen" = "blue",
#                                  "lysate" = "grey"))+
#     scale_linetype_manual(name = "Measurement",
#                           values = c("OD" = "dotted", "PFU" = "solid"))+
# 
#   theme_bw()


# ^ the scales are all wank
# try again


od_pfu <- od_sum |>
  left_join(pfu_sum, by = c("sample", "time"))


od_pfu_log <- od_pfu |>  
  mutate(
    # Transforms 0 PFU -> log10(1) = 0, avoiding -Inf
    log_pfu_clean = log10(mean_pfu + 1),
    
    # Delta method standard error & standard deviation (using offset x + 1)
    log_pfu_se = ifelse(mean_pfu > 0, pfu_se / ((mean_pfu + 1) * log(10)), 0),
    log_pfu_sd = ifelse(mean_pfu > 0, pfu_sd / ((mean_pfu + 1) * log(10)), 0),

    # Asymmetric Bounds incorporating the +1 offset
    # Set lower bound to 0 (log10(1)) if mean_pfu < pfu_se, or NA if preferred
    log_pfu_lower = ifelse(mean_pfu >= pfu_se, log10(mean_pfu - pfu_se + 1), log10(1)),
    log_pfu_upper = log10(mean_pfu + pfu_se + 1))


scale_factor <- max(od_pfu_log$log_pfu_clean, na.rm = TRUE) / max(od_pfu_log$mean_od, na.rm = TRUE)


# ----

od_pfu_plot <-ggplot(data = od_pfu_log, aes(x = time))+
  geom_line(aes(y = mean_od, colour = sample, linetype = "OD"), linewidth = 1)+
  geom_point(aes(y = mean_od, colour = sample), size = 2)+
  geom_errorbar(aes(ymin = mean_od - od_se,
                    ymax = mean_od + od_se),
                width = 0.2,
                colour = "black",
                alpha = 0.5)+
  
  geom_line(aes(y = log_pfu_clean / scale_factor, colour = sample, linetype = "PFU"), linewidth = 1.2) +
  geom_point(aes(y = log_pfu_clean / scale_factor, colour = sample), size = 2)+
  geom_errorbar(aes(ymin = log_pfu_lower / scale_factor,
                ymax = log_pfu_upper / scale_factor),
                width = 0.2,
                colour = "black",
                alpha = 0.5)+
  
  
  scale_y_continuous(name = expression("Cell Growth (" * OD[600] * ")"),
                     sec.axis = sec_axis(~ . * scale_factor, name = expression("Free Virions (" * log[10] * " PFU/mL)"))) +
  scale_x_continuous(name = "Time Post Infection (hours)",
                     limits = c(0,max(od_pfu_log$time, na.rm = TRUE)),
                     breaks = seq(0, max(od_pfu_log$time, na.rm = TRUE), by = 2 ))+
  
  scale_colour_manual(name = "Sample",
                      limits = c("host","lysogen","lysate"),
                      labels = c("host" = "Host 009",
                                 "lysogen" = "Lysogen",
                                 "lysate" = "Lysate"),
                      values = c("host" = "red",
                                 "lysogen" = "blue",
                                 "lysate" = "grey"))+
  scale_linetype_manual(name = "Measurement",
                        values = c("OD" = "dotted", "PFU" = "solid"))+
  labs(title = "MOI = 0.005")+
  theme_bw()+
  theme(plot.title = element_text(face = "bold"))

od_pfu_plot
ggsave("plots/od_pfu.png", plot = od_pfu_plot, width = 7, height = 5, dpi = 300)


# stunning


# plot for presentation ---------------------------------------------------

pres_plot <- od_pfu_log |> filter(sample == "host") |> 
  ggplot( aes(x = time))+
  geom_line(aes(y = log_pfu_clean, linetype = "PFU"), linewidth = 1.2, colour = "red") +
  geom_point(aes(y = log_pfu_clean), size = 2, colour = "red")+
  geom_errorbar(aes(ymin = log_pfu_lower,
                    ymax = log_pfu_upper),
                width = 0.2,
                colour = "black",
                alpha = 0.5)+
  
  geom_line(aes(y = mean_od * scale_factor, linetype = "OD"), linewidth = 1, colour = "black")+
  geom_point(aes(y = mean_od* scale_factor), size = 2, colour = "black")+
  geom_errorbar(aes(ymin = (mean_od - od_se)* scale_factor,
                    ymax = (mean_od + od_se)* scale_factor),
                width = 0.2,
                colour = "black",
                alpha = 0.5)+
  

  scale_y_continuous(name = expression("Free Virions (" * log[10] * " PFU/mL)"),
                     sec.axis = sec_axis(~ . * scale_factor, name = expression("Cell Growth (" * OD[600] * ")"))) +
  
  scale_x_continuous(name = "Time Post Infection (hours)",
                     limits = c(0,max(od_pfu_log$time, na.rm = TRUE)),
                     breaks = seq(0, max(od_pfu_log$time, na.rm = TRUE), by = 2 ))+
  
  scale_linetype_manual(name = "Measurement",
                        values = c("OD" = "dotted", "PFU" = "solid"))+
  labs(title = "Virus : Host = 0.005")+
  theme_bw()+
  theme(legend.position = "none",
    plot.title = element_text(face = "bold"),
        axis.title.y = element_text(colour = "red"),
        axis.title.y.right = element_text(colour = "black"))

pres_plot      
ggsave("plots/pres_plot.png", plot = pres_plot, width = 5, height = 4, dpi = 300)

