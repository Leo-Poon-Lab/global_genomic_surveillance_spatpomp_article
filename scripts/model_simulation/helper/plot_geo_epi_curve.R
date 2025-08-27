library(geofacet)
library(ggflags)
library(dplyr)
library(tidyr)
library(shadowtext)
source("scripts/model_fitting/helper/sim_cases.R")

plot_geo_epi_curve <- function(
  data_quants,
  df_meas,
  events,
  lineages_for_events,
  colors_events,
  colors_labels,
  geo_grid=geo_grid,
  add_flags=FALSE,
  add_text=TRUE
){
  lineages_for_events <- gsub(" ", "_", lineages_for_events)
  names(colors_events) <- gsub(" ", "_", names(colors_events))

  df_meas <- aggregating_measurements(df_meas)
  df_meas$date <- as.Date(lubridate::date_decimal(df_meas$date_decimal))
  if(any(names(df_meas)=="daily_cases")){
    names(df_meas)[names(df_meas)=="daily_cases"] <- "Cases"
  }

  data_quants$code <- data_quants$unitname

  df_meas_plot <- df_meas %>% select(date, code, loc_name, all_of(lineages_for_events)) %>% pivot_longer(cols = c(-date, -code, -loc_name), names_to = "lineage", values_to = "Cases")

  # single line plot
  p <- ggplot(data = df_meas_plot %>% filter(Cases>1)) +
    geom_line(aes(x = date, y = Cases, group = lineage, color = lineage), linewidth = 0.7, linetype="solid", alpha=0.9) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10(  # Plot on log-scale
      limits = c(1, 10^8.75),
      breaks = c(1, 10, 10^3, 10^5, 10^7),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = c(0, 0)
    ) +
    scale_color_manual(values = colors_events, labels = colors_labels) +
    NULL

  # ribbon bands
  for (i in seq_along(events)){
    this_event <- events[i]
    data_quants_tmp <- data_quants %>% 
      dplyr::select(code, date, quantile_level, all_of(this_event)) %>% 
      filter(date!=min(data_quants$date)) %>% 
      mutate_at(vars(all_of(this_event)), ~ifelse(.<1, 1, .)) %>%
      pivot_wider(names_from = quantile_level, values_from = this_event, values_fn = mean)
    p <- p +
      geom_ribbon(data = data_quants_tmp, aes(x = date, ymin = get("0.05"), ymax = get("0.95")), fill = colors_events[i], alpha = 0.3) +
      NULL
  }

  # To quantify the model fit, we need to calculate the Pearson's r, MAE and RMSE for each region and each event
  df_summary_fit <- tibble(
    code = character(),
    event = character(),
    Lineage = character(),
    pearson_r = numeric(),
    mae = numeric(),
    rmse = numeric())

  for (i in seq_along(events)){
    this_event <- events[i]
    this_lineage <- lineages_for_events[i]
    for (j in seq_along(geo_grid$code)){
      this_code <- geo_grid$code[j]

      df_values_actual <- df_meas_plot %>%
        filter(code == this_code, lineage == this_lineage) %>%
        select(date, Cases)

      df_values_predicted <- data_quants %>%
        filter(quantile_level == 0.5, code == this_code) %>%
        select(date, all_of(this_event))

      df_values_combined <- left_join(df_values_actual, df_values_predicted, by = "date")
      stat_pearson <- cor(df_values_combined$Cases, df_values_combined[[this_event]], use = "complete.obs")
      stat_mae <- mean(abs(df_values_combined$Cases - df_values_combined[[this_event]]), na.rm = TRUE)
      stat_rmse <- sqrt(mean((df_values_combined$Cases - df_values_combined[[this_event]])^2, na.rm = TRUE))

      # Store the statistics in the summary dataframe
      df_summary_fit <- bind_rows(df_summary_fit, tibble(
        code = this_code,
        event = this_event,
        Lineage = this_lineage,
        pearson_r = stat_pearson,
        mae = stat_mae,
        rmse = stat_rmse
      ))
    }
  }

  write_csv(df_summary_fit, paste0(dir_rst, "model_fit_summary.csv"))

  # text of region
  df_meas_plot_text <- df_meas_plot %>% filter(date==min(df_meas_plot$date), lineage=="Cases")
  df_meas_plot_text <- left_join(df_meas_plot_text, geo_grid %>% select(code, name), by="code")
  p <- p +
    geom_label(data = df_meas_plot_text, aes(x = date, y = 10^8.25, label = name), hjust = 0, vjust = 0.5, size = 2, nudge_x = 5, color=colors_spatial_units[order(names(colors_spatial_units))]) +
    NULL

  df_meas_plot_text$code_flags <- tolower(df_meas_plot_text$code)
  if(add_flags){
    p <- p +
      geom_flag(data = df_meas_plot_text %>% filter(!grepl("others", code)), aes(x = date+165, y = 10^8, country = code_flags), size = 5) +
      NULL
  }

  df_summary_fit$event <- factor(df_summary_fit$event, levels = events, labels = lineages_for_events)
  if(add_text){
    # add pearson_r, mae and rmse in df_summary_fit. For each grid, show three lines for these three statistics respectively. In each line, the numbers for six different events are shown using different color per colors_events. use scientific notation for the numbers.
    
    # Create text annotations for model fit statistics
    df_summary_fit_text <- df_summary_fit %>%
      select(code, event, pearson_r, mae, rmse) %>%
      pivot_longer(cols = c(pearson_r, mae, rmse), names_to = "metric", values_to = "value") %>%
      mutate(
      value_text = case_when(
        metric == "pearson_r" ~ sprintf("%.2f", value),
        TRUE ~ sprintf("%.1e", value)
      )
      ) %>%
      group_by(code, metric) %>%
      arrange(match(event, events)) %>%
      summarise(
      x_pos = case_when(
        metric[1] == "pearson_r" ~ min(df_meas_plot$date),
        metric[1] == "mae" ~ min(df_meas_plot$date) + 30,
        metric[1] == "rmse" ~ min(df_meas_plot$date) + 80
      ),
      y_pos = 10^7.8,
      .groups = "drop"
      )

    # Add metric labels
    df_metric_labels <- df_summary_fit_text %>%
      select(code, metric, x_pos, y_pos) %>%
      distinct() %>%
      mutate(
      label = case_when(
        metric == "pearson_r" ~ "r",
        metric == "mae" ~ "MAE",
        metric == "rmse" ~ "RMSE"
      )
      )

    # Add individual colored values
    df_metric_values <- df_summary_fit %>%
      select(code, event, pearson_r, mae, rmse) %>%
      pivot_longer(cols = c(pearson_r, mae, rmse), names_to = "metric", values_to = "value") %>%
      mutate(
      value_text = case_when(
        metric == "pearson_r" ~ sprintf("%.2f", value),
        TRUE ~ {
        exponent <- floor(log10(value))
        mantissa <- value / 10^exponent
        sprintf("%.1f%s10^%g", mantissa, "%*%", exponent)
        }
      ),
      x_pos = case_when(
        metric == "pearson_r" ~ min(df_meas_plot$date),
        metric == "mae" ~ min(df_meas_plot$date) + 30,
        metric == "rmse" ~ min(df_meas_plot$date) + 80
      )
      ) %>%
      group_by(code, metric) %>%
      mutate(
      y_offset = (seq_along(event) - 1) * 0.35,
      y_pos = 10^(7 - y_offset)
      ) %>%
      ungroup()
    
    df_metric_values$value_text <- gsub("%*%10^0", "", df_metric_values$value_text, fixed = TRUE) # Remove the "10^0" part from the text
    df_metric_values$value_text <- gsub("%*%10^1", "%*%10", df_metric_values$value_text, fixed = TRUE) 

    p <- p +
      geom_text(data = df_metric_labels,
          aes(x = x_pos, y = y_pos, label = label),
          hjust = 0, vjust = 1, size = 2, color = "grey30", fontface = "bold") +
      geom_shadowtext(data = df_metric_values %>% filter(!grepl("NaN", value_text)),
          aes(x = x_pos, y = y_pos, label = value_text, color = event),
          hjust = 0, vjust = 0, size = 1.8, parse = TRUE, show.legend = FALSE,
          bg.colour = "white", bg.r = 0.1) +
      scale_color_manual(values = colors_events)
  }

  p <- p + 
    facet_geo(~code, grid = geo_grid, label="name") +
    theme_minimal() +
    theme(
      # Top-right position
      legend.position = "bottom",
      # Elements within a guide are placed one next to the other in the same row
      legend.direction = "horizontal",
      # Different guides are stacked vertically
      legend.box = "vertical",
      # No legend title
      legend.title = element_blank(),
      # Light background color
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 30, 20, 30),
      # # Customize the title. Note the new font family and its larger size.
      # plot.title = element_text(
      #   margin = margin(0, 0, -100, 0), 
      #   size = 26, 
      #   family = "KyivType Sans", 
      #   face = "bold", 
      #   vjust = 0, 
      #   color = "grey25"
      # ),
      plot.caption = element_text(size = 11),
      # Remove titles for x and y axes.
      axis.title = element_blank(),
      # Specify color for the tick labels along both axes 
      axis.text = element_text(color = "grey40"),
      # Specify face and color for the text on top of each panel/facet
      # strip.text = element_text(color = "grey20")
      strip.background = element_blank(),
      strip.text = element_blank()
    ) + 
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 2)), 
           fill = guide_legend(nrow = 1))

}

get_geo_grid <- function(
  num_country=86,
  regions_code_3,
  other_codes
  ){
  # get geo grid
  if(num_country == 86){
    data_geo_grid <- read_csv("https://raw.githubusercontent.com/hafen/grid-designer/master/grids/world_86countries_grid.csv")
  } else if(num_country == 192){
    data_geo_grid <- read_csv("https://raw.githubusercontent.com/hafen/grid-designer/master/grids/world_countries_grid1.csv")
    names(data_geo_grid)[2] <- "code"
  } else{
    stop("num_country must be either 86 or 192")
  }
  data_geo_grid_sub <- data_geo_grid %>% filter(code %in% regions_code_3)
  if(any(!regions_code_3 %in% data_geo_grid$code)){
    code_to_add <- regions_code_3[!regions_code_3 %in% data_geo_grid$code]
    data_geo_grid_sub <- bind_rows(data_geo_grid_sub, tibble(code=code_to_add, name=code_to_add, row=max(data_geo_grid$row)+1, col=max(data_geo_grid$col)+seq_along(code_to_add)))
  }
  if(all(!is.na(other_codes))){
    data_geo_grid_sub <- bind_rows(data_geo_grid_sub, tibble(code=other_codes, name=other_codes, row=max(data_geo_grid$row)+2, col=max(data_geo_grid$col)+seq_along(other_codes)))
  }
  ## change row to be sequential
  old_row <- sort(unique(data_geo_grid_sub$row))
  new_row <- 1:length(old_row)
  data_geo_grid_sub <- data_geo_grid_sub %>% mutate(row = new_row[match(row, old_row)])
  ## change col to be sequential
  old_col <- sort(unique(data_geo_grid_sub$col))
  new_col <- 1:length(old_col)
  data_geo_grid_sub <- data_geo_grid_sub %>% mutate(col = new_col[match(col, old_col)])

  return(data_geo_grid_sub)
}