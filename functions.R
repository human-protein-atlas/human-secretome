


make_dendrogram <- function(var1, var2, value, method = "spearman"){
  if(method == "spearman") {
    data.frame(var1, var2, value) %>%
      spread(key = var1, value = value) %>%
      column_to_rownames("var2") %>% 
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram() 
  } else if(method == "euclidean") {
    data.frame(var1, var2, value) %>%
      spread(key = var1, value = value) %>%
      column_to_rownames("var2") %>% 
      as.matrix() %>% 
      t() %>% 
      na.omit() %>%
      dist() %>%
      hclust("average")  %>%
      as.dendrogram() 
  } else if(method == "ward.D2") {
    data.frame(var1, var2, value) %>%
      spread(key = var1, value = value) %>%
      column_to_rownames("var2") %>% 
      as.matrix() %>%
      t() %>%
      na.omit() %>%
      dist() %>%
      hclust("ward.D2") %>%
      as.dendrogram() 
  }
  
}

get_dendrogram_segments <- function(dendrogram) {
  dendrogram %>%
    ggdendro::dendro_data() %$%
    left_join(segments, labels, by = c("xend" = "x", "yend" = "y"))
}

range_scale_manual <- function(x, xmax, xmin, span, dodge = 0) { 
  ((x - xmin)/(xmax - xmin))*
    (span[2] - span[1]) + span[1] + 
    dodge
}
range_scale <- function(x, span) { 
  xmin = min(x)
  xmax = max(x)
  ((x - xmin)/(xmax - xmin))*
    (span[2] - span[1]) + span[1]
}
ggdendroheat <- function(x, y, value, fill_factor = NA, show.legend = T, xdendrogram = T, ydendrogram = T,
                         x_margin = 0.2, y_margin = 0.2, xaxis_order = NULL, yaxis_order = NULL, 
                         xdendrogram_color = NULL, ydendrogram_color = NULL,
                         range_scale_x = F, range_scale_y = F, x_clustering_method = "spearman", y_clustering_method = "spearman",
                         under_lim_tile_color = NA, under_lim = 1){
  
  x_dendrogram <- 
    make_dendrogram(x, y, value, method = x_clustering_method)
  x_dendrogram_segments <- get_dendrogram_segments(x_dendrogram)
  
  y_dendrogram <- 
    make_dendrogram(y, x, value, method = y_clustering_method)
  y_dendrogram_segments <- get_dendrogram_segments(y_dendrogram)
  
  g <- 
    tibble(x, y, value, 
           under_limit = value < under_lim,
           fill_factor = fill_factor) 
  
  if(range_scale_x){
    g <- g %>%
      group_by(x) %>%
      mutate(value = range_scale(value, c(0,1))) %>%
      ungroup()
  }
  
  if(range_scale_y){
    g <- g %>%
      group_by(y) %>%
      mutate(value = range_scale(value, c(0,1))) %>%
      ungroup()
  }
  
  if(is.null(xaxis_order)){ 
    g <- g %>%
      mutate(x = factor(x, levels = labels(x_dendrogram))) 
  } else {
    if(xdendrogram) warning("displaying dendrogram and giving explicit axis order is not recommended! Axis and dendrogram will show conflicting order.")
    g <- g %>%
      mutate(x = factor(x, levels = xaxis_order)) 
  }
  
  if(is.null(yaxis_order)){ 
    g <- g %>%
      mutate(y = factor(y, levels = labels(y_dendrogram))) 
  } else {
    if(ydendrogram) warning("displaying dendrogram and giving explicit axis order is not recommended! Axis and dendrogram will show conflicting order.")
    g <- g %>%
      mutate(y = factor(y, levels = yaxis_order)) 
  }
  
  
  if(is.na(fill_factor[1])) {
    g <- g %>%
      ggplot() +
      geom_tile(aes(x, y, fill = value), show.legend = show.legend) + 
      geom_tile(data = filter(g, under_limit), aes(x, y), fill = under_lim_tile_color, show.legend = F)
  } else {
    g <- g %>%
      ggplot() +
      geom_tile(aes(x, y, fill = fill_factor), show.legend = show.legend) + 
      geom_tile(data = filter(g, under_limit), aes(x, y), fill = under_lim_tile_color, show.legend = F)
  }
  
  
  xp <- 0
  yp <- 0
  
  
  x_factors <- length(labels(x_dendrogram))
  y_factors <- length(labels(y_dendrogram))
  
  if(xdendrogram) {
    
    if(!is.null(xdendrogram_color)) {
      g <- g + 
        geom_segment(data = x_dendrogram_segments %>%
                       mutate(label = xdendrogram_color[match(x_dendrogram_segments$label, names(xdendrogram_color))]), 
                     aes(x = x, 
                         xend = xend,
                         y = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(y_factors, y_factors/(1 - y_margin)), 
                                                dodge = 0.5), 
                         yend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(y_factors, y_factors/(1 - y_margin)), 
                                                   dodge = 0.5),
                         color = label))
    } else{
      g <- g + 
        geom_segment(data = x_dendrogram_segments, 
                     aes(x = x, 
                         xend = xend,
                         y = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(y_factors, y_factors/(1 - y_margin)), 
                                                dodge = 0.5), 
                         yend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(y_factors, y_factors/(1 - y_margin)), 
                                                   dodge = 0.5)))
    }
    
    
    xp <- c(x_factors + x_margin)
  } 
  
  if(ydendrogram) {
    if(!is.null(ydendrogram_color)) {
      g <- g + 
        geom_segment(data = y_dendrogram_segments %>%
                       mutate(label = ydendrogram_color[match(y_dendrogram_segments$label, names(ydendrogram_color))]), 
                     aes(x = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(x_factors, x_factors/(1 - x_margin)), 
                                                dodge = 0.5), 
                         xend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(x_factors, x_factors/(1 - x_margin)), 
                                                   dodge = 0.5),
                         y = x, 
                         yend = xend,
                         color = label))
    } else{
      g <- g + 
        geom_segment(data = y_dendrogram_segments, 
                     aes(x = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(x_factors, x_factors/(1 - x_margin)), 
                                                dodge = 0.5), 
                         xend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(x_factors, x_factors/(1 - x_margin)), 
                                                   dodge = 0.5),
                         y = x, 
                         yend = xend))
      
    }
    
    
    
    yp <- c(y_factors + y_margin)
  }
  
  g + 
    annotate("point", x = xp, y = yp, color = NA) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8), 
          panel.grid = element_blank())
  
  
}


Aasa_heatmap <- function(df, xcol, ycol, value_col, xband_color_cols, 
                         x_clustering_method = "euclidean",
                         y_clustering_method = "euclidean",
                         x_margin = 0.2,
                         y_margin = 0.2,
                         x_dendrogram_nudge = 0,
                         y_dendrogram_nudge = 0,
                         value_title = "value") {
  
  
  g <- 
    df %>%
    select(setNames(c(xcol, ycol, value_col, xband_color_cols),
                    c("x", "y", "value", paste0("color_", 1:length(xband_color_cols))))) 
  
  x_dendrogram <- 
    g %$%
    make_dendrogram(x, y, value, method = x_clustering_method)
  x_dendrogram_segments <- get_dendrogram_segments(x_dendrogram)
  x_factors <- length(labels(x_dendrogram))
  
  y_dendrogram <- 
    g %$%
    make_dendrogram(y, x, value, method = y_clustering_method)
  y_dendrogram_segments <- get_dendrogram_segments(y_dendrogram)
  y_factors <- length(labels(y_dendrogram))
  
  g <- 
    g %>% 
    mutate(x = factor(x, levels = labels(x_dendrogram)), 
           y = factor(y, levels = labels(y_dendrogram)))
  
  p <- g %>%
    ggplot() +
    geom_tile(aes(x, y, fill = value)) 
  
  p <- 
    p + {
      data <- 
        left_join(x_dendrogram_segments,
                  g %>% 
                    select(x, 
                           setNames(paste0("color_", 1:length(xband_color_cols)), 1:length(xband_color_cols))) %>%
                    unique() %>%
                    gather(key = "i", value = "color", -x),
                  by = c("label" = "x")) %>% 
        filter(!is.na(label)) %>%
        mutate(i = as.integer(i))
      
      
      geom_rect(data = data ,
                aes(xmin = x - 0.5, 
                    xmax = x + 0.5,
                    ymin = y_factors + 0.5 + (x_dendrogram_nudge / length(xband_color_cols)) * (i - 1),
                    ymax = y_factors + 0.5 + (x_dendrogram_nudge / length(xband_color_cols)) * (i)), 
                fill = data$color)
    }
  
  p +
    geom_segment(data = x_dendrogram_segments , 
                 aes(x = x, 
                     xend = xend,
                     y = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                            span = c(y_factors + x_dendrogram_nudge, y_factors/(1 - y_margin)), 
                                            dodge = 0.5), 
                     yend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                               span = c(y_factors + x_dendrogram_nudge, y_factors/(1 - y_margin)), 
                                               dodge = 0.5))) + 
    geom_segment(data = y_dendrogram_segments, 
                 aes(x = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                            span = c(x_factors, x_factors/(1 - x_margin)), 
                                            dodge = 0.5), 
                     xend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                               span = c(x_factors, x_factors/(1 - x_margin)), 
                                               dodge = 0.5),
                     y = x, 
                     yend = xend)) + 
    scale_fill_viridis_c(name = value_title, option = "inferno", direction = -1) +
    scale_color_manual(values = secreted_to_colors)+
    theme(axis.text.x = element_blank(), panel.background = element_rect(fill = "white"))
}

chord_classification <- function(from, to, sizes, grid.col, groups = rep(1, length(from)), 
                                 plot.order = c(unique(from), unique(to)), line_expansion = 10000, size_labels = F, link_alpha = 0.5){
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, 
                 order = plot.order, transparency = 1 - link_alpha)
  
  if(size_labels) {
    for(i in 1:nrow(chord)) {
      value <- chord$value[i]
      if(is.null(value)) value <- chord$value1[1]
      x1 <- chord$x1[i] - value / 2
      x2 <- chord$x2[i] - value / 2
      to_ <- chord$cn[i]
      from_ <- chord$rn[i]
      circos.text(x = x1, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = from_, niceFacing = T)
      circos.text(x = x2, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = to_, niceFacing = T)
    }
  }
  
  
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)
    
    width <- strwidth(sector.name)*line_expansion
    
    circos.segments(x0 = mean(xlim), x1 = mean(xlim), 
                    y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment, 
                    sector.name)
    
    circos.segments(x0 = mean(xlim) - width/2, x1 = mean(xlim) + width/2, 
                    y0 = mean(ylim) - 0.2 + adjustment, y1 = mean(ylim) - 0.2 + adjustment, 
                    sector.name) 
    
    circos.text(mean(xlim), mean(ylim) + adjustment, sector.name, niceFacing = TRUE, facing = "bending")
  }, bg.border = NA)
  
  circos.clear()
}



Aasa_facet_heatmap <- function(df, xcol, ycol, value_col, xband_color_cols, facet_col, 
                               x_clustering_method = "euclidean",
                               y_clustering_method = "euclidean",
                               x_margin = 0.2,
                               y_margin = 0.2,
                               x_dendrogram_nudge = 0,
                               y_dendrogram_nudge = 0,
                               value_title = "value") {
  
  
  g <- 
    df %>%
    select(setNames(c(xcol, ycol, value_col, facet_col, xband_color_cols),
                    c("x", "y", "value", "facet", paste0("color_", 1:length(xband_color_cols))))) %>%
    group_by(facet) %>%
    mutate(facet_genenum = paste0(facet, " (n=", length(unique(x)), ")")) %>%
    ungroup() %>%
    mutate(facet = factor(facet_genenum, levels = loc_levels_n))
  
  facets <- unique(g$facet)
  x_dendrograms <- 
    lapply(facets, FUN = function(f) {
      g %>% 
        filter(facet == f) %$%
        make_dendrogram(x, y, value, method = x_clustering_method)
    })
  
  x_dendrogram_segments <- 
    lapply(x_dendrograms, FUN = function(d) {
      get_dendrogram_segments(d)
    })
  
  x_labels <- 
    lapply(x_dendrograms, FUN = labels) 
  
  x_factors <- 
    sapply(x_labels, FUN = length) %>% 
    setNames(facets)
    
  x_dendrogram_segments_df <- 
    x_dendrogram_segments %>% 
    set_names(facets) %>%
    plyr::ldply(.id = "facet") %>% 
    as.tibble()
  
  x_labels_df <- 
    x_labels %>% 
    {tibble(x = unlist(.),
            facet = rep(facets, x_factors))}
  
  y_factors <- length(unique(g$y))
  
  g <- 
    g %>% 
    mutate(x = factor(x, levels = x_labels_df$x)) 
  
  p <- g %>%
    ggplot() +
    geom_tile(aes(x, y, fill = value)) +
    facet_wrap(~ facet, scales = "free_x")
  
  p <- 
    p + {
      data <- 
        left_join(x_dendrogram_segments_df,
                  g %>% 
                    select(x, 
                           setNames(paste0("color_", 1:length(xband_color_cols)), 1:length(xband_color_cols))) %>%
                    unique() %>%
                    gather(key = "i", value = "color", -x),
                  by = c("label" = "x")) %>% 
        filter(!is.na(label)) %>%
        mutate(i = as.integer(i))
      
      
      geom_rect(data = data ,
                aes(xmin = x - 0.5, 
                    xmax = x + 0.5,
                    ymin = y_factors + 0.5 + (x_dendrogram_nudge / length(xband_color_cols)) * (i - 1),
                    ymax = y_factors + 0.5 + (x_dendrogram_nudge / length(xband_color_cols)) * (i)), 
                fill = data$color)
    }
  
  p +
    geom_segment(data = x_dendrogram_segments_df %>%
                   group_by(facet) %>% 
                   mutate(y_ = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                 span = c(y_factors + x_dendrogram_nudge, y_factors/(1 - y_margin)), 
                                                 dodge = 0.5),
                          yend_ = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                    span = c(y_factors + x_dendrogram_nudge, y_factors/(1 - y_margin)), 
                                                    dodge = 0.5)), 
                 aes(x = x, 
                     xend = xend,
                     y = y_, 
                     yend = yend_)) + 
   
    scale_fill_viridis_c(name = value_title, option = "inferno", direction = -1) +
    scale_color_manual(values = secreted_to_colors)+
    theme(axis.text.x = element_blank(), panel.background = element_rect(fill = "white"))
}

