#!/usr/bin/Rscript

pacman::p_load(tidyverse, ggplot2, patchwork, parallel)

## Run on local machine
model_name <- "M3"
base_model_name <- "Omicron20"
dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")
stopifnot(file.exists(dir_model_data))

## Load data
values_change_IDR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_change_DSR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_traveler_weight <- c(1/1000, 1/100, seq(0.1, 0.9, 0.2), 1)
# values_strategies <- c(expand.grid(c("T3", "T3a", "T7", "T7a", "T11", "T11a"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))
values_strategies <- c(expand.grid(c("T3", "T7", "T11"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))

df_all_values <- list.files(dir_model_data, "df_all_values_M3_\\d+.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}
lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
df_vline <- tibble(lag=lags_intro)
df_vline$date_decimal <- floor(df_vline$lag)/365.25 + data_fitting$dates$date_start_d
df_vline$date <- date(date_decimal(df_vline$date_decimal))
df_vline$variant <- c("Omicron BA.1", "Omicron BA.2")

# Below is for visualizing the EDT results
source("scripts/model_simulation/helper/plot_fig_3.R")
plot_fig_3(dir_rst, dir_model_data, df_all_values, df_vline, ncores=64)
