#!/usr/bin/Rscript

pacman::p_load(tidyverse, ggplot2, patchwork, parallel, doParallel)

## Run on local machine
model_name <- "M4"
base_model_name <- "Omicron20"
dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")

n_cores <- 64

## Load data
values_traveler_weight <- c(1/1000, 1/100, seq(0.1, 0.9, 0.2), 1)
# values_strategies <- c(expand.grid(c("T2", "T2a", "T7", "T7a", "T11", "T11a"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))
values_strategies <- c(expand.grid(c("T2", "T7", "T11"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))
values_id_unit_intro <- seq(0, 28, 1)
df_all_values <- list.files(dir_model_data, "df_all_values_M4_\\d{1,2}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

df_all_values_fig_4abcd <- df_all_values %>% filter(!grepl("+", scenario, fixed=TRUE)) %>% filter(!grepl("T3-", scenario, fixed=TRUE)) %>% filter(!grepl("T3a-", scenario, fixed=TRUE))
df_all_values_fig_S5 <- df_all_values %>% filter(grepl("+", scenario, fixed=TRUE) | grepl("STp_TW_ID_T3-H", scenario, fixed=TRUE))

data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

## Below is for visualizing the results
source("scripts/model_simulation/helper/plot_fig_4abcd.R")

list_p_fig_4abcd <- plot_fig_4abcd(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values_todo = df_all_values_fig_4abcd,
  data_fitting = data_fitting,
  bootstrap_each_param_pair=1024*10,
  n_cores = n_cores
)

source("scripts/model_simulation/helper/plot_fig_S5.R")

list_p_fig_S5 <- plot_fig_S5(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values_todo = df_all_values_fig_S5,
  data_fitting = data_fitting,
  bootstrap_each_param_pair=1024*10,
  n_cores = n_cores
)

## simulate the data for Fig 5abc
list_p_fig_4abcd <- readRDS(paste0(dir_rst, "list_p_fig_4abcd.rds"))
ggsave(list_p_fig_4abcd$p_fig_4abcd, filename = paste0(dir_rst, "fig_4.jpg"), width = 20, height = 25, units = "in", dpi = 300)
ggsave(list_p_fig_4abcd$p_fig_4abcd, filename = paste0(dir_rst, "fig_4.pdf"), width = 20, height = 25, units = "in")

source("scripts/model_simulation/helper/prepare_simulation_fig_5abc.R")

## Manually select the best parameter set, based on the results of Fig 4abcd
best_param_set <- list_p_fig_4abcd$df_EDT_fig_4bcd_percentile %>%
  filter((parval_2=="P2-H" & parval_1=="0.1%") | (parval_2=="T2-R" & parval_1=="70%")) %>% 
  mutate(parval_1 = as.numeric(gsub("%","",as.character(parval_1)))/100) 

df_all_values_new <- prepare_simulation_fig_5abc(
  best_param_set = best_param_set,
  base_model_name = base_model_name,
  values_change_IDR = c(1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 1),
  values_change_DSR = c(1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 1),
  values_id_unit_intro = seq(0, 28, 1),
  data_fitting = data_fitting,
  check_redo = TRUE
)
