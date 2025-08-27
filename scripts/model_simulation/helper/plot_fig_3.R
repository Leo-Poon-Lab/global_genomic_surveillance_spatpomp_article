source("scripts/data_processing/helper/color_scheme.R")
library(scales) # Add this line to import the scales package

# Fig. 3a: Differences in earliest detection dates in four scenarios comparing to the base scenario. Four scenarios: (1) Increase detection capacity; (2) Increase overall sequencing capacity; (3) Adjust Traveler weight (SP) for imported cases; (4) Same resources but increase sequencing capacity in travel hubs. ([ggstatsplot](https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html)). For BA.1 and BA.2 separately.## upper panel: BA.1; lower panel: BA.2

plot_fig_3a <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores=64
){
  df_all_values_fig_3a <- df_all_values %>% filter(fig=="fig_3a")

  # load simulation data and get EDT
  df_EDT <- mclapply(seq_along(df_all_values_fig_3a$scenario), function(i){
    model_name_i <- paste0("M3_", df_all_values_fig_3a$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    df_sims$scenario <- df_all_values_fig_3a$scenario[i]
    return(df_sims)
  }, mc.cores = ncores) %>% bind_rows()

  df_EDT_base <- df_EDT %>% filter(scenario=="sequencing_1")
  df_EDT_base_hline <- df_EDT_base %>% group_by(variant, case_type) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(round(mean(lag), 4), 4)) %>% mutate(scenario="base")
  df_EDT_base_hline$lineage_type <- paste0(df_EDT_base_hline$variant, " (", tolower(df_EDT_base_hline$case_type), ")")
  print(df_EDT_base_hline)

  df_EDT$lineage_type <- paste0(df_EDT$variant, " (", tolower(df_EDT$case_type), ")")

  df_EDT$params <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[1]
  })
  df_EDT$params <- factor(df_EDT$params, levels = c("detection", "sequencing", "TW"), labels = c("Changing diagnostics", "Changing sequencing", "Traveler weight (all regions)"))
  df_EDT$parval <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[2]
  })
  df_EDT$parval[df_EDT$parval=="NA" & df_EDT$params=="Changing diagnostics"] <- "1"
  df_EDT$parval[df_EDT$parval=="NA" & df_EDT$params=="Changing sequencing"] <- "1"
  df_EDT$parval <- as.numeric(df_EDT$parval)
  df_EDT$parval <- factor(df_EDT$parval, levels = sort(unique(df_EDT$parval)), labels = gsub(".0%", "%", scales::percent(sort(unique(df_EDT$parval)), accuracy = 0.1)))

  df_medians <- df_EDT %>% group_by(parval, lineage_type, variant, params) %>% summarise(median_lag = median(lag), max_lag = boxplot.stats(lag)[[1]][5])
  df_EDT %>% filter(params=="Changing sequencing") %>% filter(parval=="0.1%") %>% filter(variant=="Omicron BA.2") %>% .$lag %>% boxplot.stats()
  df_medians %>% filter(params=="Changing sequencing") %>% filter(parval=="0.1%") %>% filter(variant=="Omicron BA.2")

  df_EDT_base_hline_note <- df_EDT_base_hline %>% summarise(note = paste0("Baseline:\n", paste0(round(lag_median, 1), " (", substr(case_type, 1, 3), ")", collapse = ", ")))

  # plot
  p_fig_3a <- ggplot(df_EDT) + 
    geom_hline(data=df_EDT_base_hline, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
    geom_boxplot(aes(x=parval, y=lag, color=lineage_type), alpha=0.5, outlier.shape=NA) +
    geom_text(data=df_medians, aes(x=parval, y=max_lag, label=round(median_lag, 1), color=lineage_type), vjust=-0.2, position=position_dodge(width=0.75)) +
    geom_text(data=df_EDT_base_hline_note, aes(x="50%", label=note), y=110, hjust=0, vjust=0, color="darkred") +
    facet_grid(rows = vars(variant), cols = vars(params), scales="free_x", space="free_x") +
    scale_color_manual(name="", values=colors_lineage_type) +
    scale_linetype_manual(name="", values=c("dashed", "solid"))+
    # scale_y_continuous(limits=c(0, ylimit)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "bottom"
    )+
    ylab("Identification lag (days)")+
    guides(color = guide_legend(nrow = 1)) +
    NULL
  ggsave(paste0(dir_rst, "fig_3a.jpg"), p_fig_3a, width=12, height=6)
  ggsave(paste0(dir_rst, "fig_3a.pdf"), p_fig_3a, width=12, height=6)
  writexl::write_xlsx(df_EDT, paste0(dir_rst, "fig_3a_source.xlsx"))
  writexl::write_xlsx(df_EDT_base_hline, paste0(dir_rst, "fig_3a_base_hline.xlsx"))
  return(p_fig_3a)
}

# Fig. 3b - 3g: 2-D density plots of the frequencies of achieving lowest detection lags, for 2 parameters. 3b: Traveler weight with detection capacity (BA.1/BA.2); 3c: Traveler weight with sequencing capacity (BA.1/BA.2); 3d: Traveler weight with travel hubs (BA.1/BA.2)

plot_fig_3bcdef <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores
){
  name_figs_all <- c("fig_3b", "fig_3c", "fig_3d", "fig_3e", "fig_3f")
  df_all_values_todo <- df_all_values %>% filter((fig %in% name_figs_all) | (scenario=="sequencing_1"))

  # load simulation data and get EDT
  df_EDT <- mclapply(seq_along(df_all_values_todo$scenario), function(i){
    model_name_i <- paste0("M3_", df_all_values_todo$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    df_sims$scenario <- df_all_values_todo$scenario[i]
    df_sims$fig <- df_all_values_todo$fig[i]
    return(df_sims)
  }, mc.cores = ncores) %>% bind_rows()

  df_EDT_base <- df_EDT %>% filter(scenario=="sequencing_1")
  df_EDT_base_median <- df_EDT_base %>% group_by(variant, case_type) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(mean(lag), 4)) %>% mutate(scenario="base")
  df_EDT_base_median$lineage_type <- paste0(df_EDT_base_median$variant, " (", tolower(df_EDT_base_median$case_type), ")")

  range_color_1 <- c(1,85)
  range_color_2 <- c(8,45)

  list_p_fig_3bcdef <- lapply(seq_along(name_figs_all), function(i){
    # i=1
    title_sub_fig <- c("b", "c", "d", "e", "f")[i]

    range_color_i <- if(i %in% 1:3){
      range_color_1
    } else{
      range_color_2
    }

    name_fig_i <- name_figs_all[i]
    df_EDT_fig_i <- df_EDT %>% filter(fig==name_fig_i)
    levels_to_use <- c("detection", "sequencing", "TW", "STp", "STt")
    if(title_sub_fig %in% c("e", "f")){
      labels_to_use <- c("Changing diagnostics", "Changing sequencing", "Traveler weight (travel hubs)", "Strategies (per capita)", "Strategies (total)")
    } else {
      labels_to_use <- c("Changing diagnostics", "Changing sequencing", "Traveler weight (all regions)", "Strategies (per capita)", "Strategies (total)")
    }
    df_EDT_fig_i <- df_EDT_fig_i %>% separate(scenario, into=c("params_1", "params_2", "parval_1", "parval_2"), sep="_") %>% mutate(
      params_1 = factor(params_1, levels = levels_to_use, labels = labels_to_use),
      parval_1 = as.character(parval_1),
      params_2 = factor(params_2, levels = levels_to_use, labels = labels_to_use),
      parval_2 = as.character(parval_2)
    )

    if(grepl("Traveler weight", df_EDT_fig_i$params_2[1])){
      # switch the order of params_1 and params_2 if params_2 is Traveler weight
      df_EDT_fig_i <- df_EDT_fig_i %>% mutate(
        tmp_params = params_1,
        tmp_parval = parval_1,
        params_1 = params_2,
        parval_1 = parval_2,
        params_2 = tmp_params,
        parval_2 = tmp_parval
      ) %>% select(-tmp_params, -tmp_parval)
    }

    if(any(grepl("T3-R", df_EDT_fig_i$parval_2))){
      df_EDT_fig_i <- df_EDT_fig_i %>% filter(parval_2 %in% values_strategies)
      if(df_EDT_fig_i$params_2[1]=="Strategies (per capita)"){
        values_strategies_for_use <- gsub("^T", "P", values_strategies)
      } else{
        values_strategies_for_use <- values_strategies
      }
      values_strategies_new <- gsub("-R", "-R^0", values_strategies_for_use, fixed = TRUE)
      values_strategies_new <- gsub("-H", "-R^50", values_strategies_new, fixed = TRUE)
      values_strategies_new <- gsub("-M", "-R^90", values_strategies_new, fixed = TRUE)

      df_EDT_fig_i$parval_2_new <- factor(df_EDT_fig_i$parval_2, levels = values_strategies, labels = values_strategies_new)
      df_EDT_fig_i$parval_2 <- factor(df_EDT_fig_i$parval_2, levels = values_strategies, labels = values_strategies_for_use)
    } else {
      df_EDT_fig_i$parval_2_new <- df_EDT_fig_i$parval_2
    }

    list_p_fig_i <- lapply(c("Omicron BA.1", "Omicron BA.2"), function(this_variant){
      # this_variant = "Omicron BA.1"
      df_EDT_fig_i_j <- df_EDT_fig_i %>% filter(variant==this_variant) # filter by variant

      # Check and modify parval_1
      if(!any(grepl("[TP]3-R", df_EDT_fig_i_j$parval_1))){
        df_EDT_fig_i_j$parval_1 <- as.numeric(as.character(df_EDT_fig_i_j$parval_1))
        df_EDT_fig_i_j$parval_1 <- factor(df_EDT_fig_i_j$parval_1, levels = sort(unique(df_EDT_fig_i_j$parval_1)), labels = gsub(".0%", "%", scales::percent(sort(unique(df_EDT_fig_i_j$parval_1)), accuracy = 0.1)))
      }
      
      # Check and modify parval_2
      if(!any(grepl("[TP]3-R", df_EDT_fig_i_j$parval_2))){
        df_EDT_fig_i_j$parval_2 <- as.numeric(as.character(df_EDT_fig_i_j$parval_2))
        df_EDT_fig_i_j$parval_2 <- factor(df_EDT_fig_i_j$parval_2, levels = sort(unique(df_EDT_fig_i_j$parval_2)), labels = gsub(".0%", "%", scales::percent(sort(unique(df_EDT_fig_i_j$parval_2)), accuracy = 0.1)))
      }

      set.seed(2024)
      df_EDT_fig_i_j <- df_EDT_fig_i_j %>% group_by(id_sim, variant, params_1, params_2, parval_1, parval_2, parval_2_new) %>% arrange(lag, desc(value)) %>% slice_min(lag) %>% slice_max(value) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if two case types have with the same lowest lag, we will then choose the one with highest number of sequenced cases, if still ties, randomly sample one.
      df_EDT_fig_i_j_median <- df_EDT_fig_i_j %>% group_by(variant, params_1, params_2, parval_1, parval_2, parval_2_new) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(mean(lag), 4)) %>% ungroup() # showing the median lag for each tile (heatmap background)
      write_csv(df_EDT_fig_i_j_median, paste0(dir_rst, name_fig_i, "_", gsub(" ", "_", this_variant), "_heatmap_source.csv"))
      df_EDT_fig_i_j_earliest <- df_EDT_fig_i_j %>% group_by(variant, params_1, params_2, parval_1, parval_2, parval_2_new) %>% mutate(lag_median=round(median(lag), 4), value_median=round(median(value), 4)) %>% ungroup() %>% group_by(variant, id_sim) %>% slice_min(lag_median) %>% ungroup() # choose lowest lag_median, if there are ties, keep all.
      write_csv(df_EDT_fig_i_j_earliest, paste0(dir_rst, name_fig_i, "_", gsub(" ", "_", this_variant), "_label_source.csv"))
      df_EDT_fig_i_j_min <- df_EDT_fig_i_j_earliest %>% ungroup() %>% group_by(parval_1, parval_2, parval_2_new) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_sim=n()) %>% ungroup() %>% mutate(percent_Community=round(N_Community/n_sim*100), percent_imported=round(N_imported/n_sim*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag, stratified by community and imported cases, for each tile (label/text)

      midpoint_value <- df_EDT_base_median %>% filter(variant==this_variant) %>% pull(lag_median) %>% min()

      p_fig_i_j <- ggplot()+
        geom_tile(data=df_EDT_fig_i_j_median, aes(x=parval_1, y=parval_2_new, fill=lag_median))+
        geom_text(data=df_EDT_fig_i_j_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2_new, label=percent_Community), color="black", nudge_x = -0.2)+
        geom_label(data=df_EDT_fig_i_j_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2_new, label=percent_imported), color="black", nudge_x = 0.2, alpha=0.9, label.padding=unit(0.1, "lines"))+
        geom_text(data=df_EDT_fig_i_j_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2_new, label=percent_Community), color="black")+
        geom_label(data=df_EDT_fig_i_j_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2_new, label=percent_imported), color="black", alpha=0.9, label.padding=unit(0.1, "lines"))+
        xlab(df_EDT_fig_i$params_1[1])+
        ylab(df_EDT_fig_i$params_2[1])+
        scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color_i)+
        scale_size_continuous(range=c(2, 3))+
        facet_grid(rows = vars(variant), scales="free_x", space="free_x")+
        theme_minimal()+
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.key.height = unit(0.3, 'cm'),
          # axis.text.x = element_text(angle = 15),
          legend.position = "bottom"
          )+
        guides(size="none")+
        NULL

      if(any(grepl("[TP]3-R", df_EDT_fig_i_j$parval_2))){
        p_fig_i_j <- p_fig_i_j+
          scale_y_discrete(labels = str2expression(as.character(df_EDT_fig_i_j_median$parval_2_new)))
      } else {
        p_fig_i_j <- p_fig_i_j+
          scale_y_discrete(labels = as.character(df_EDT_fig_i_j_median$parval_2))
      }

      return(p_fig_i_j)
      })

    p_fig_i <- (list_p_fig_i[[1]] + ggtitle(title_sub_fig)) + list_p_fig_i[[2]] + plot_layout(ncol=1, axes = "keep", guides = "keep") & theme(legend.position = "bottom")
    ggsave(paste0(dir_rst, name_fig_i, ".jpg"), p_fig_i, width=6, height=6)
    ggsave(paste0(dir_rst, name_fig_i, ".pdf"), p_fig_i, width=6, height=6)
    writexl::write_xlsx(df_EDT_fig_i, paste0(dir_rst, name_fig_i, "_source.xlsx"))
    return(p_fig_i)
  })

  p_fig_3bcd <- wrap_plots(list_p_fig_3bcdef[[1]], list_p_fig_3bcdef[[2]], list_p_fig_3bcdef[[3]], nrow=1)
  ggsave(paste0(dir_rst, "fig_3bcd.jpg"), p_fig_3bcd, width=29.7-3, height=21-5, scale=0.5)
  ggsave(paste0(dir_rst, "fig_3bcd.pdf"), p_fig_3bcd, width=29.7-3, height=21-5, scale=0.5)
  
  p_fig_3ef <- wrap_plots(list_p_fig_3bcdef[[4]], list_p_fig_3bcdef[[5]], nrow=1)
  ggsave(paste0(dir_rst, "fig_3ef.jpg"), p_fig_3ef, width=21, height=29.7*0.45, scale=0.55)
  ggsave(paste0(dir_rst, "fig_3ef.pdf"), p_fig_3ef, width=21, height=29.7*0.45, scale=0.55)

  p_fig_3bcdef <- wrap_plots(p_fig_3bcd, p_fig_3ef, nrow=2, byrow=TRUE, heights = c(7,7))
  ggsave(paste0(dir_rst, "fig_3bcdef.jpg"), p_fig_3bcdef, width=12, height=(12*297/210))
  ggsave(paste0(dir_rst, "fig_3bcdef.pdf"), p_fig_3bcdef, width=12, height=(12*297/210))
  return(p_fig_3bcdef)
}

plot_fig_3 <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores
){
  p_fig_3a <- plot_fig_3a(dir_rst, dir_model_data, df_all_values, df_vline, ncores)+ggtitle("a")
  p_fig_3bcdef <- plot_fig_3bcdef(dir_rst, dir_model_data, df_all_values, df_vline, ncores)

  p_fig_3 <- p_fig_3a + p_fig_3bcdef + plot_layout(ncol=1, heights = c(0.8, 2))
  ggsave(paste0(dir_rst, "fig_3.jpg"), p_fig_3, width=15, height=20)
  ggsave(paste0(dir_rst, "fig_3.pdf"), p_fig_3, width=15, height=20)

}
