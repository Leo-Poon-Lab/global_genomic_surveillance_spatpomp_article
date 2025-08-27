source("scripts/data_processing/install_prerequisite.R")
pacman::p_load(tidyverse, countrycode, janitor, padr, hablar, ggflags, ggbump, lubridate, wesanderson, shadowtext, sf, rnaturalearth, BBmisc, patchwork)

source("scripts/data_processing/helper/color_scheme.R")

# Fig 1a: plot the bump graph for different ranked factors
# Fig 1b: global transmission map linking to the bump graph
# Fig 1c: plot the epi curve with the geo_grid

plot_fig_1 <- function(
  model_name,
  model_sims,
  model_quants,
  dir_rst
){
  # read data
  data_measurements <- readRDS(paste0("results/model_data/data_measurements_", model_name, ".rds")) # for reported cases, imported/community sequenced cases
  data_fitting <- readRDS(paste0("results/model_data/data_fitting_", model_name, ".rds")) # for travel volume, population, infected cases
  
  # Fig 1a. data
  ## Fig 1a. left panel: global map
  spdf_world <- ne_countries()
  # spdf_world <- spdf_world@data
  all_codes_2 <- sort(data_fitting$country_under_investigation)
  all_position <- lapply(seq_along(all_codes_2), function(i){
    # Get geospatial data
    if(all_codes_2[i]=="AF_others"){
      return(c(28, 11))
    } else if(all_codes_2[i]=="AS_others"){
      return(c(99, 10))
    } else if(all_codes_2[i]=="EU_others"){
      return(c(15, 62))
    } else if(all_codes_2[i]=="NA_others"){
      return(c(-87, 20))
    } else if(all_codes_2[i]=="OC_others"){
      return(c(161, -22))
    } else if(all_codes_2[i]=="SA_others"){
      return(c(-66, -20))
    } else if(all_codes_2[i]=="HK"){
      return(c(114.17, 22.29))
    }
    data_country <- spdf_world[which(spdf_world$iso_a2_eh == all_codes_2[i]), ]
    point <- c(data_country$label_x[1], data_country$label_y[1])
  })
  df_pos <- tibble(code=all_codes_2, lon=sapply(all_position, function(x) {x[1]}), lat=sapply(all_position, function(x) {x[2]}))
  df_pos$code_flags <- tolower(df_pos$code)

  ## map data
  world <- rnaturalearthdata::countries50 %>% st_as_sf() 
  world$code <- world$iso_a2
  world <- left_join(world %>% select(-continent), data_fitting$df_new_code, "code")

  ## Fig 1a. right panel: bump graph
  ## ranked columns: Population, travel volume, infected cases, reported cases, imported/community sequenced cases;
  df_plot_sub1.1 <- data_fitting$travel$mov_mat %>% 
    filter(date>=data_fitting$dates$date_start, date<=data_fitting$dates$date_end) %>%
    group_by(code_o) %>% 
    summarise(travel_volume_outflow = sum(flow_all)) %>%
    dplyr::rename(code=code_o)
  df_plot_sub1.2 <- data_fitting$travel$mov_mat %>% 
    filter(date>=data_fitting$dates$date_start, date<=data_fitting$dates$date_end) %>%
    group_by(code_d) %>% 
    summarise(travel_volume_inflow = sum(flow_all)) %>%
    dplyr::rename(code=code_d)
  df_plot_sub1 <- left_join(df_plot_sub1.1, df_plot_sub1.2, by="code") %>% 
    mutate(travel_volume = travel_volume_outflow + travel_volume_inflow) %>% 
    select(code, travel_volume)

  df_plot_sub2 <- data_fitting$df_covar %>% 
    select(date_decimal, code, population, inf_cuml_mean) %>% 
    filter(date_decimal>=data_fitting$dates$date_start_d, date_decimal<=data_fitting$dates$date_end_d) %>% 
    group_by(code) %>% 
    summarise(
      population = max(population),
      infected_new = max(inf_cuml_mean) - min(inf_cuml_mean)
      )

  df_plot_sub3 <- data_measurements %>% 
    dplyr::rename(report_cases = daily_cases) %>%
    filter(date_decimal>=data_fitting$dates$date_start_d, date_decimal<=data_fitting$dates$date_end_d) %>% 
    group_by(code, date_decimal) %>%
    transmute(
      imported_sequences = sum(c_across(matches("C_._i_.+"))),
      community_sequences = sum(c_across(matches("C_._c"))),
      report_cases
      ) %>% 
    ungroup() %>% 
    group_by(code) %>%
    summarise(
      imported_sequences = sum(imported_sequences),
      community_sequences = sum(community_sequences),
      report_cases = sum(report_cases)
    )

  df_plot_1a_right <- left_join(df_plot_sub1, df_plot_sub2, by="code") %>% left_join(df_plot_sub3, by="code")
  df_plot_1a_right <- df_plot_1a_right %>% transmute(
    code,
    population,
    travel_volume = travel_volume,
    infected_new = infected_new,
    report_cases = report_cases,
    community_sequences = community_sequences,
    imported_sequences = imported_sequences
    )
  df_plot_1a_right <- df_plot_1a_right %>% 
    transmute(
      code,
      rank_population = rank(-population, ties.method = "random"),
      rank_travel_volume = rank(-travel_volume, ties.method = "random"),
      rank_infected_new = rank(-infected_new, ties.method = "random"),
      rank_report_cases = rank(-report_cases, ties.method = "random"),
      rank_community_sequences = rank(-community_sequences, ties.method = "random"),
      rank_imported_sequences = rank(-imported_sequences, ties.method = "random")
    ) %>% 
    pivot_longer(cols = c(rank_travel_volume, rank_population, rank_infected_new, rank_imported_sequences, rank_community_sequences, rank_report_cases), names_to = "rank_var", values_to = "rank")

  df_plot_1a_right$rank_var <- factor(
    df_plot_1a_right$rank_var, 
    levels = c(
      "rank_travel_volume", "rank_population", "rank_infected_new", "rank_report_cases", "rank_community_sequences", "rank_imported_sequences"
    ),
    labels = c(
      "Travel volume", "Population", "Infected cases", "Reported cases", "Sequenced cases\n(Community)", "Sequenced cases\n(Imported)"
    )
  )

  df_plot_1a_right$code_flags <- tolower(df_plot_1a_right$code)
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "eu_others"] <- "Europe\n(others)"
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "af_others"] <- "Africa\n(others)"
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "as_others"] <- "Asia\n(others)"
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "sa_others"] <- "South America\n(others)"
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "na_others"] <- "North America\n(others)"
  df_plot_1a_right$code_flags[df_plot_1a_right$code_flags == "oc_others"] <- "Oceania\n(others)"

  df_plot_1a_right_link <- df_plot_1a_right %>% group_by(rank_var) %>% 
    mutate(
      # bump_y_start = -90,
      bump_y_start = normalize(-rank, range = c(-303, 90), method = "range"),
      # bump_x_start = normalize(-rank, range = c(-180, 180), method = "range")
      bump_x_start = 200
    )
  df_plot_1a_right_link <- left_join(df_plot_1a_right_link, df_pos %>% select(code, lon, lat), "code")

  p_fig_1a_left <- ggplot(data = world) + 
    geom_sf(data = world %>% filter(grepl("others", new_code)), aes(fill = new_code), lwd = 0, colour = "white", alpha=0.3) + 
    geom_sf(data = world %>% filter(!grepl("others", new_code)), aes(fill = new_code), lwd = 0, colour = "white", alpha=0.9) + 
    scale_fill_manual(values = c(colors_spatial_units), na.value = "grey90") +
    geom_sigmoid(data = df_plot_1a_right_link %>% filter(rank_var == "Travel volume", !grepl("others", code)), 
      aes(x = lon, y = lat, xend = bump_x_start, yend = bump_y_start, color = code), 
      alpha = 0.9, smooth = 12, size = 1, direction = "x") +
    geom_sigmoid(data = df_plot_1a_right_link %>% filter(rank_var == "Travel volume", grepl("others", code)), 
      aes(x = lon, y = lat, xend = bump_x_start, yend = bump_y_start, color = code), 
      alpha = 0.3, smooth = 12, size = 1, direction = "x") +
    geom_flag(
      data = df_pos %>% filter(!grepl("others", code_flags)),
      aes(x = lon, y = lat, country = code_flags),
      size = 6
    )+
    geom_label(data = df_pos %>% filter(grepl("others", code_flags)), aes(x = lon, y = lat, label = code, color=code), size = 2.3, alpha=0.9) + 
    scale_color_manual(values = colors_spatial_units)+
    scale_size_identity()+
    theme_void() +
    coord_sf(xlim = c(-150, 180), ylim = c(-60, 80), expand = TRUE)+
    theme(
      legend.position = "none"
      ) + 
    NULL
  ggsave(paste0(dir_rst, "/fig_1a_left.pdf"), plot=p_fig_1a_left, width=8*1.2, height=8, dpi=400)
  ggsave(paste0(dir_rst, "/fig_1a_left.jpg"), plot=p_fig_1a_left, width=8*1.2, height=8, dpi=400)

  ## try to match Fig 1a left and 1a right
  p_fig_1a_right <- df_plot_1a_right %>% 
    ggplot(aes(rank_var, rank, group = code, color = code, fill = code)) +
    geom_bump(size = 1.5, lineend = "round") + 
    # big flag at the two ends
    geom_flag(
      data = df_plot_1a_right %>% filter(rank_var %in% c("Travel volume", "Sequenced cases\n(Imported)"), !grepl("others", code_flags)),
      aes(country = code_flags),
      size = 8
    )+
    geom_label(
      data = df_plot_1a_right %>% filter(rank_var %in% c("Travel volume", "Sequenced cases\n(Imported)"), grepl("others", code_flags)),
      aes(label=code_flags),
      fill="white",
      size=2
    )+
    # small flag in the middle
    geom_flag(
      data = df_plot_1a_right %>% filter(!rank_var %in% c("Travel volume", "Sequenced cases\n(Imported)"), !grepl("others", code_flags)),
      aes(country = code_flags),
      size = 6
    )+
    geom_label(
      data = df_plot_1a_right %>% filter(!rank_var %in% c("Travel volume", "Sequenced cases\n(Imported)"), grepl("others", code_flags)),
      aes(label=code_flags),
      fill="white",
      size=1.75
    )+
    scale_y_reverse(breaks = 1:29, position = "right")+
    scale_color_manual(values = colors_spatial_units) +
    scale_fill_manual(values = colors_spatial_units) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "transparent"),
      text = element_text(color = "black"),
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      # axis.text.x=element_blank(),
      axis.title = element_blank(),
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      axis.ticks.length = unit(0, "pt")
      )+
    # coord_flip()+
    NULL
  ggsave(paste0(dir_rst, "/fig_1a_right.jpg"), plot=p_fig_1a_right, width=8*1.4, height=8, dpi=400)
  ggsave(paste0(dir_rst, "/fig_1a_right.pdf"), plot=p_fig_1a_right, width=8*1.4, height=8, dpi=400)

  p_fig_1a_leftright <- p_fig_1a_left + p_fig_1a_right + plot_layout(ncol=2, widths = c(1, 1.5))
  ggsave(paste0(dir_rst, "/fig_1a_leftright.pdf"), plot=p_fig_1a_leftright, width=8*2, height=8, dpi=400)

  write_csv(df_plot_1a_right, paste0(dir_rst, "/fig_1a_right.csv"))

  # Fig 1b. data
  ## Transmission map showing the spread routes of Omicron BA.1 and BA.2 variants
  ## separate for 3 sub maps, for the first, the second and the third 2-weeks after emergence.
  lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
  dates_intro <- floor(lags_intro)/365.25+min(model_sims$date_decimal)

  df_sims <- model_sims %>%
    mutate(
      C_1_i_infected_new = rowSums(select(., starts_with("C_1_i_infected_origin_"))),
      C_1_i_detected_new = rowSums(select(., starts_with("C_1_i_detected_origin_"))),
      C_1_i_sequenced_new = rowSums(select(., starts_with("C_1_i_o_"))),
      C_2_i_infected_new = rowSums(select(., starts_with("C_2_i_infected_origin_"))),
      C_2_i_detected_new = rowSums(select(., starts_with("C_2_i_detected_origin_"))),
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_o_"))),
      C_3_i_infected_new = rowSums(select(., starts_with("C_3_i_infected_origin_"))),
      C_3_i_detected_new = rowSums(select(., starts_with("C_3_i_detected_origin_"))),
      C_3_i_sequenced_new = rowSums(select(., starts_with("C_3_i_o_"))
      ),
      C_4_i_infected_new = rowSums(select(., starts_with("C_4_i_infected_origin_"))),
      C_4_i_detected_new = rowSums(select(., starts_with("C_4_i_detected_origin_"))),
      C_4_i_sequenced_new = rowSums(select(., starts_with("C_4_i_o_"))
      )
    )

  transform_to_df_spread <- function(
    df_sims,
    date_lower,
    date_upper,
    col_i_infected_new,
    col_i_infected_origin_pattern,
    variant_name,
    df_pos,
    point_diameter,
    square_side_length,
    x_adjust,
    spread_direction
  ){
    ## average the imported cases in simulation results
    df_spread <- df_sims %>% select(date_decimal, `.id`, unitname, all_of(col_i_infected_new), matches(col_i_infected_origin_pattern)) %>% filter(date_decimal>=date_lower, date_decimal<date_upper) %>% filter(get(col_i_infected_new)>=1) %>% select(-all_of(col_i_infected_new)) %>% pivot_longer(cols = matches(col_i_infected_origin_pattern), names_to = "origin", values_to = "value") %>% separate(origin, into=c("variant", "origin"), sep="infected_origin_") %>% mutate(variant=variant_name, origin=gsub("\\D", "", origin)) %>% mutate(origin = factor(origin, levels=0:28, labels = data_fitting$country_under_investigation)) %>% filter(value>0) %>% group_by(date_decimal, origin, unitname, variant) %>% summarise(value=mean(value)) %>% ungroup()

    if(nrow(df_spread)==0){
      return(tibble(
        origin=character(0),
        unitname=character(0),
        variant=character(0),
        value=numeric(0),
        ori_lon=numeric(0),
        ori_lat=numeric(0),
        des_lon=numeric(0),
        des_lat=numeric(0),
        color=NA,
        id=NA,
        des_lat_new=numeric(0),
        des_lon_new=numeric(0),
        ori_lat_new=numeric(0),
        ori_lon_new=numeric(0)
      ))
    }

    ## aggregate the cases for the whole period
    df_spread <- df_spread %>% group_by(origin, unitname, variant) %>% summarise(value=sum(value)) %>% ungroup()

    df_spread <- left_join(df_spread, df_pos %>% transmute(origin=code, ori_lon=lon, ori_lat=lat))
    df_spread <- left_join(df_spread, df_pos %>% transmute(unitname=code, des_lon=lon, des_lat=lat))
    # df_spread$day <- round((df_spread$date_decimal-date_lower)*365.25)
    
    df_spread$color <- colors_lineage[variant_name]
    df_spread <- spread_point_locations(df_spread, df_pos, point_diameter=point_diameter, square_side_length=square_side_length, x_adjust=x_adjust, direction=spread_direction)
    return(df_spread)
  }

  spread_point_locations <- function(
    df_spread,
    df_pos,
    point_diameter = 0.1,
    square_side_length= 0.3,
    x_adjust = 0.1,
    direction = "right"
  ){
    direction <- tolower(direction)
    stopifnot(direction %in% c("right", "left"))
    nrow_max <- ceiling(square_side_length/point_diameter)
    
    positions_all <- df_spread %>% select(origin, unitname) %>% unlist() %>% unique()
    df_spread_new <- df_spread %>% mutate(id=1:nrow(df_spread))
    df_spread_new$des_lat_new <- NA
    df_spread_new$des_lon_new <- NA
    df_spread_new$ori_lat_new <- NA
    df_spread_new$ori_lon_new <- NA

    # for each unit, we generate new positions for points
    lapply(seq_along(positions_all), function(i){
      code_i <- positions_all[i]
      df_origin_i <- df_spread_new %>% filter(origin==code_i)
      df_dest_i <- df_spread_new %>% filter(unitname==code_i)
      freq_origin_i <- nrow(df_origin_i)
      freq_dest_i <- nrow(df_dest_i)
      freq_all_i <- freq_origin_i + freq_dest_i

      ## get the central position
      lon_i <- df_pos$lon[df_pos$code==code_i]
      lat_i <- df_pos$lat[df_pos$code==code_i]

      ## **Dest points** should be placed on the inner right side of the central position
      if(freq_dest_i>0){
        ncol_dest_points <- ceiling(freq_dest_i/nrow_max)
        y_new_dest_points <- seq(from=0, by=point_diameter, length.out=min(c(nrow_max, freq_dest_i)))
        ### center and adjust the latitude position
        y_new_dest_points <- y_new_dest_points - mean(y_new_dest_points) + lat_i
        ### adjust the longitude position
        x_new_dest_points <- seq(from=0, by=point_diameter, length.out=ncol_dest_points)
        if(direction == "right"){
          x_new_dest_points <- lon_i + x_adjust + x_new_dest_points
        } else {
          x_new_dest_points <- lon_i - x_adjust - x_new_dest_points
        }
        df_xy_new_dest <- expand.grid(des_lat_new=y_new_dest_points, des_lon_new=x_new_dest_points)[1:freq_dest_i,]
        ### sort by the origin latitude
        df_dest_i <- df_dest_i %>% arrange(desc(ori_lat), ori_lon)
        ### assign the new positions
        df_dest_i$des_lat_new <- df_xy_new_dest$des_lat_new
        df_dest_i$des_lon_new <- df_xy_new_dest$des_lon_new
        df_dest_i <- df_dest_i %>% arrange(id)
        df_spread_new$des_lat_new[df_spread_new$id %in% df_dest_i$id] <<- df_dest_i$des_lat_new
        df_spread_new$des_lon_new[df_spread_new$id %in% df_dest_i$id] <<- df_dest_i$des_lon_new
      }

      ## **Origin points** should be placed on the outer right side
      if(freq_origin_i>0){
        ncol_origin_points_i <- ceiling(freq_origin_i/nrow_max)
        y_new_origin_points <- seq(from=0, by=point_diameter, length.out=min(c(nrow_max, freq_origin_i)))
        ### center and adjust the latitude position
        y_new_origin_points <- y_new_origin_points - mean(y_new_origin_points) + lat_i
        ### adjust the longitude position
        x_new_origin_points <- seq(from=0, by=point_diameter, length.out=ncol_origin_points_i)
        if(direction == "right"){
          if(freq_dest_i>0){
            x_new_origin_points <- max(c(df_xy_new_dest$des_lon_new, lon_i), na.rm = T) + x_adjust + x_new_origin_points
          } else{
            x_new_origin_points <- lon_i + x_adjust + x_new_origin_points
          }
        } else {
          if(freq_dest_i>0){
            x_new_origin_points <- min(c(df_xy_new_dest$des_lon_new, lon_i), na.rm = T) - x_adjust - x_new_origin_points
          } else{
            x_new_origin_points <- lon_i - x_adjust - x_new_origin_points
          }
        }
        df_xy_new_origin <- expand.grid(ori_lat_new=y_new_origin_points, ori_lon_new=x_new_origin_points)[1:freq_origin_i,]
        ### sort by the destination latitude
        df_origin_i <- df_origin_i %>% arrange(des_lat, des_lon)
        ### assign the new positions
        df_origin_i$ori_lat_new <- df_xy_new_origin$ori_lat_new
        df_origin_i$ori_lon_new <- df_xy_new_origin$ori_lon_new
        df_origin_i <- df_origin_i %>% arrange(id)
        df_spread_new$ori_lat_new[df_spread_new$id %in% df_origin_i$id] <<- df_origin_i$ori_lat_new
        df_spread_new$ori_lon_new[df_spread_new$id %in% df_origin_i$id] <<- df_origin_i$ori_lon_new
      }
      return(NULL)
    }) 
    return(df_spread_new)
  }

  plot_spread_map <- function(
    df_spread_1,
    df_spread_2,
    curvature,
    size,
    title
  ){
    line_size=size
    point_size=size*0.6

    ggplot(data = world) + 
    geom_sf(data = world %>% filter(is.na(new_code)), lwd = 0, colour = "white", alpha=0.2) + 
    geom_sf(data = world %>% filter(!is.na(new_code)), lwd = 0.1, colour = "grey30", alpha=0.8) + 
    geom_curve(data = df_spread_1, aes(x = ori_lon_new, y = ori_lat_new, xend = des_lon_new, yend = des_lat_new, color = color, alpha = value), size = line_size, curvature = curvature) +
    geom_curve(data = df_spread_2, aes(x = ori_lon_new, y = ori_lat_new, xend = des_lon_new, yend = des_lat_new, color = color, alpha = value), size = line_size, curvature = -curvature) +
    geom_flag(
      data = tibble(code = unique(c(df_spread_1$origin, df_spread_1$unitname, df_spread_2$origin, df_spread_2$unitname))) %>% left_join(df_pos, "code") %>% filter(!grepl("others", code_flags)),
      aes(x = lon, y = lat+3, country = code_flags),
      size = 6
    )+
    geom_label(data = df_pos %>% filter(grepl("others", code_flags)), aes(x = lon, y = lat+3, label = code), size = 2.3, alpha=0.9) + 
    geom_point(data = df_spread_1, aes(x = ori_lon_new, y = ori_lat_new, color = color), size = point_size, shape = 1) + # circle for origin
    geom_point(data = df_spread_1, aes(x = des_lon_new, y = des_lat_new, color = color), size = point_size, shape = 2) + # triangle for destination
    geom_point(data = df_spread_2, aes(x = des_lon_new, y = des_lat_new, color = color), size = point_size, shape = 2) + 
    geom_point(data = df_spread_2, aes(x = ori_lon_new, y = ori_lat_new, color = color), size = point_size, shape = 1) +
    # geom_curve(data = df_spread_1 %>% filter(origin!="ZA"), aes(x = ori_lon, y = ori_lat, xend = des_lon, yend = des_lat), alpha = alpha*1.2, color="black", curvature = curvature, size=0.3) +
    # geom_curve(data = df_spread_2 %>% filter(origin!="ZA"), aes(x = ori_lon, y = ori_lat, xend = des_lon, yend = des_lat), alpha = alpha*1.2, color="black", curvature = -curvature, size=0.3) +
    coord_sf(
      xlim = c(
        min(c(df_spread_1$ori_lon_new, df_spread_1$des_lon_new, df_spread_2$ori_lon_new, df_spread_2$des_lon_new)),
        max(c(df_spread_1$ori_lon_new, df_spread_1$des_lon_new, df_spread_2$ori_lon_new, df_spread_2$des_lon_new))
      ),
      ylim = c(
        min(c(df_spread_1$ori_lat_new, df_spread_1$des_lat_new, df_spread_2$ori_lat_new, df_spread_2$des_lat_new)),
        max(c(df_spread_1$ori_lat_new, df_spread_1$des_lat_new, df_spread_2$ori_lat_new, df_spread_2$des_lat_new))
        ),
      expand = TRUE)+
    scale_alpha(name="Number of cases", range = c(0.1, 1))+
    scale_color_identity(name = "Variant",
      breaks = c(df_spread_1$color[1], df_spread_2$color[1]),
      labels = c(df_spread_1$variant[1], df_spread_2$variant[1]),
      guide = "legend") +
    theme_void() +
    theme(
      legend.position = "bottom"
      ) + 
    guides(
      color = guide_legend(override.aes = list(alpha = 1, lwd = 2)),
    )+
    ggtitle(title)+
    NULL
  }

  df_spread_BA1_first_1w <- transform_to_df_spread(
    df_sims,
    date_lower=dates_intro[1],
    date_upper=dates_intro[1]+7*1/365.25,
    col_i_infected_new="C_2_i_infected_new",
    col_i_infected_origin_pattern="C_2_i_infected_origin.+_new$",
    variant_name="Omicron BA.1",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "right"
  )
  df_spread_BA2_first_1w <- transform_to_df_spread(
    df_sims,
    dates_intro[2],
    dates_intro[2]+7*1/365.25,
    "C_3_i_infected_new",
    "C_3_i_infected_origin.+_new$",
    "Omicron BA.2",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "left"
  )
  if(nrow(df_spread_BA1_first_1w)>0 | nrow(df_spread_BA2_first_1w)>0){
    p_fig_1b_top <- plot_spread_map(df_spread_1=df_spread_BA1_first_1w, df_spread_2=df_spread_BA2_first_1w, 0.4, size=1, "The first week after emergence")
    ggsave(paste0(dir_rst, "/fig_1b_top.jpg"), plot=p_fig_1b_top, width=10, height=10, dpi=400)
    ggsave(paste0(dir_rst, "/fig_1b_top.pdf"), plot=p_fig_1b_top, width=10, height=10, dpi=400)
  } else {
    file.remove(paste0(dir_rst, "/fig_1b_top.jpg"))
    file.remove(paste0(dir_rst, "/fig_1b_top.pdf"))
  }

  df_spread_BA1_second_1w <- transform_to_df_spread(
    df_sims,
    date_lower=dates_intro[1]+7*1/365.25,
    date_upper=dates_intro[1]+7*2/365.25,
    col_i_infected_new="C_2_i_infected_new",
    col_i_infected_origin_pattern="C_2_i_infected_origin.+_new$",
    variant_name="Omicron BA.1",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "right"
  )
  df_spread_BA2_second_1w <- transform_to_df_spread(
    df_sims,
    dates_intro[2]+7*1/365.25,
    dates_intro[2]+7*2/365.25,
    "C_3_i_infected_new",
    "C_3_i_infected_origin.+_new$",
    "Omicron BA.2",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "left"
  )
  if(nrow(df_spread_BA1_second_1w)>0 | nrow(df_spread_BA2_second_1w)>0){
    p_fig_1b_middle <- plot_spread_map(df_spread_1=df_spread_BA1_second_1w, df_spread_2=df_spread_BA2_second_1w, 0.4, size=0.5, "The second week after emergence")
    ggsave(paste0(dir_rst, "/fig_1b_middle.jpg"), plot=p_fig_1b_middle, width=10, height=10, dpi=400)
    ggsave(paste0(dir_rst, "/fig_1b_middle.pdf"), plot=p_fig_1b_middle, width=10, height=10, dpi=400)
  } else {
    file.remove(paste0(dir_rst, "/fig_1b_middle.jpg"))
    file.remove(paste0(dir_rst, "/fig_1b_middle.pdf"))
  }

df_spread_BA1_first_2w <- transform_to_df_spread(
    df_sims,
    date_lower=dates_intro[1],
    date_upper=dates_intro[1]+7*2/365.25,
    col_i_infected_new="C_2_i_infected_new",
    col_i_infected_origin_pattern="C_2_i_infected_origin.+_new$",
    variant_name="Omicron BA.1",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "right"
  )
  df_spread_BA2_first_2w <- transform_to_df_spread(
    df_sims,
    dates_intro[2],
    dates_intro[2]+7*2/365.25,
    "C_3_i_infected_new",
    "C_3_i_infected_origin.+_new$",
    "Omicron BA.2",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "left"
  )
  if(nrow(df_spread_BA1_first_2w)>0 | nrow(df_spread_BA2_first_2w)>0){
    p_fig_1b_middle <- plot_spread_map(df_spread_1=df_spread_BA1_first_2w, df_spread_2=df_spread_BA2_first_2w, 0.4, size=0.5, "The 1st to 2nd week after emergence")
    ggsave(paste0(dir_rst, "/fig_1b_top_middle.jpg"), plot=p_fig_1b_middle, width=10, height=10, dpi=400)
    ggsave(paste0(dir_rst, "/fig_1b_top_middle.pdf"), plot=p_fig_1b_middle, width=10, height=10, dpi=400)
  } else {
    file.remove(paste0(dir_rst, "/fig_1b_top_middle.jpg"))
    file.remove(paste0(dir_rst, "/fig_1b_top_middle.pdf"))
  }

  df_spread_BA1_third_2w <- transform_to_df_spread(
    df_sims,
    dates_intro[1]+7*2/365.25,
    dates_intro[1]+7*4/365.25,
    "C_2_i_infected_new",
    "C_2_i_infected_origin.+_new$",
    "Omicron BA.1",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "right"
  )
  df_spread_BA2_third_2w <- transform_to_df_spread(
    df_sims,
    dates_intro[2]+7*2/365.25,
    dates_intro[2]+7*4/365.25,
    "C_3_i_infected_new",
    "C_3_i_infected_origin.+_new$",
    "Omicron BA.2",
    df_pos,
    point_diameter=3/4,
    square_side_length=3,
    x_adjust=1,
    spread_direction = "left"
  )
  if(nrow(df_spread_BA1_third_2w)>0 | nrow(df_spread_BA2_third_2w)>0){
    p_fig_1b_bottom <- plot_spread_map(df_spread_1=df_spread_BA1_third_2w, df_spread_2=df_spread_BA2_third_2w, 0.4, size=0.5, "The 3rd to 4th weeks after emergence")
    ggsave(paste0(dir_rst, "/fig_1b_bottom.jpg"), plot=p_fig_1b_bottom, width=10, height=10, dpi=400)
    ggsave(paste0(dir_rst, "/fig_1b_bottom.pdf"), plot=p_fig_1b_bottom, width=10, height=10, dpi=400)
  } else {
    file.remove(paste0(dir_rst, "/fig_1b_bottom.jpg"))
    file.remove(paste0(dir_rst, "/fig_1b_bottom.pdf"))
  }

  # Fig 1c. 
  plot_fig_1c(model_name, model_sims, model_quants, dir_rst)
  return(NULL)

}

plot_fig_1c <- function(
  model_name,
  model_sims,
  model_quants,
  dir_rst
){
  data_measurements <- readRDS(paste0("results/model_data/data_measurements_", model_name, ".rds")) # for reported cases, imported/community sequenced cases
  cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

  # Fig 1c. data
  df_meas <- data_measurements
  df_meas$date <- as.Date(lubridate::date_decimal(df_meas$date_decimal))
  df_meas <- left_join(df_meas, cross_check_table %>% select(code, loc_name))
  df_meas$loc_name[is.na(df_meas$loc_name)] <- df_meas$code[is.na(df_meas$loc_name)]

  model_quants <- model_quants %>% filter(date_decimal>model_quants$date_decimal[1])
  model_quants$code <- model_quants$unitname
  model_quants <- left_join(model_quants, cross_check_table %>% select(code, loc_name), by="code")
  model_quants$loc_name[is.na(model_quants$loc_name)] <- model_quants$code[is.na(model_quants$loc_name)]

  source("scripts/model_fitting/helper/geo_grid.R")
  layout_geo_grid_29 <- geo_grid_29()
  source("scripts/model_simulation/helper/plot_geo_epi_curve.R")
  p_fig_1c <- plot_geo_epi_curve(
    data_quants=model_quants,
    df_meas=df_meas,
    events=c("daily_cases", "daily_deaths", "C_1_all", "C_2_all", "C_3_all", "C_4_all"),
    lineages_for_events=c("Cases", "Deaths", "Delta", "Omicron BA_one", "Omicron BA_two", "Others"),
    colors_events=c("Cases"="grey", "Deaths"="chocolate4", colors_lineage_old),
    colors_labels=c("All cases (diagnostic)", "All deaths (diagnostic)", "Delta (sequenced)", "Omicron BA.1 (sequenced)", "Omicron BA.2 (sequenced)", "Others (sequenced)"),
    geo_grid=layout_geo_grid_29,
    add_flags=TRUE
    )
  ggsave(paste0(dir_rst, "/fig_1c.jpg"), plot=p_fig_1c, width=10.5, height=10.5, dpi=400)
  ggsave(paste0(dir_rst, "/fig_1c.pdf"), plot=p_fig_1c, width=10.5, height=10.5, dpi=400)
}
