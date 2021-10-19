# Project and Code Overview ====
#
# Project:      Bioclassification
# Contact:      e-mail: John.O'brien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z
#               (2021) Modelling demersal fish and benthic invertebrate assemblages 
#               in support of marine conservation planning. Ecol Appl
#
# Overview:
# Functions to assist in analyses of distribution of fish and benthic invertebrate
# assemblages
#

# Grid Filter ====

# Function that creates grid within study area of given resolution 
# and excludes grid cells that intersect with land

GridFilter <- function(regionshape, resol = 1, landshape){
  
  require(raster)
  require(fasterize)
  require(sf)
  
  grid <- raster(regionshape)
  res(grid) <- resol
  grid <- fasterize(regionshape, grid)
  gridpolygon <- rasterToPolygons(grid) %>% 
    st_as_sf() %>% 
    transmute(GridID = c(1:length(layer)))
  # Remove grid cells that intersect with land
  Grid.sel <- st_intersects(gridpolygon, landshape, sparse = TRUE) %>%
    lengths(.) == 0
  gridpolygonNoLand <- gridpolygon[Grid.sel,] %>% 
    mutate(GridID = 1:length(GridID))
}

# ggplot colourtree ====

# displays colour-coded dendrogram with ggplot

plot_ggdendro <- function(dendro_data, palette, region_tag, tag_position, ytick_args,
                          ylab_args, legend_position, legend_justification){
  
  require(ggplot2)
  
  dplot <- ggplot(dendro_data) + 
    theme_classic() +
    geom_segment(aes(x = x, 
                     y = y, 
                     xend = xend, 
                     yend = yend, 
                     colour = col, 
                     size = lwd),
                 lineend = 'square',
                 linejoin = 'bevel') +
    scale_size_identity() +
    scale_color_identity('', 
                         guide = 'legend', 
                         breaks = palette$assigned, 
                         labels = palette$name) +
    guides(col = guide_legend(override.aes = list(size = 2))) +
    scale_y_continuous(name = 'Simpson Dissimilarity', limits = c(0,1), expand = c(0,0)) +
    labs(tag = region_tag) +
    theme(legend.text = element_text(size = 12),
          legend.position = legend_position,
          legend.justification = legend_justification,
          legend.direction = 'vertical',
          legend.key = element_rect(fill = 'white', colour = 'white'),
          axis.text.y = ytick_args,
          axis.title.y = ylab_args,
          axis.line.x = element_line(colour = 'black'),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.tag = element_text(face = 'bold', size = 12),
          plot.tag.position = tag_position,
          panel.border = element_rect(fill = NA, size = 1), 
          panel.background = element_rect(fill = 'grey85'))
}

# ggplot colour grid ====

# displays colour-coded map of assemblages assigned to grid cells with ggplot

plot_colourgrid <- function(col_grid, boundaries_poly, land_poly, 
                          label_position, expand_scale, xlim, ylim,
                          x_breaks, y_breaks, expansion_const,
                          region_tag, tag_position, ytext_angle,
                          scalebar_position, scalebar_width, scalebar_ypad 
                          ){
  
  require(ggplot2)
  
  gridPlot <- ggplot() +
    geom_sf(data = col_grid, 
            aes(fill = assigned),
            col = NA,
            show.legend = F,
            na.rm = T) + 
    scale_fill_identity() +
    geom_sf(data = boundaries_poly,
            fill = NA, 
            col = 'black') +
    geom_sf(data = land_poly, 
            fill = 'grey30',
            col = 'grey20', 
            lwd = 0.25) +
    coord_sf(label_graticule = label_position,
             expand = expand_scale,
             xlim = xlim,
             ylim = ylim) +
    scale_x_continuous(name = '',
                       breaks = x_breaks,
                       expand = expansion_const) +
    scale_y_continuous(name = '',
                       breaks = y_breaks,
                       expand = expansion_const) +
    annotate('text', 
             label = paste0("bold(", region_tag, ")"),
             parse = T,
             size = 6,
             x = tag_position[1],
             y = tag_position[2]) +
    annotation_scale(bar_cols = c('grey40','grey10'),
                     text_cex = 1,
                     line_col = 'grey10',
                     location = scalebar_position,
                     width_hint = scalebar_width,
                     pad_y = unit(scalebar_ypad, 'cm')) +
    theme(panel.background = element_rect(fill = '#e0e0e0'),
          panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
          plot.margin = unit(c(0,0,0,0), 'lines'),
          panel.border = element_rect(fill = NA, colour = 'black'),
          axis.text = element_text(size = 12.5),
          axis.text.y.right = element_text(hjust = 0.5, angle = ytext_angle))
  
}

# ggplot colour map ====

# displays colour-coded map of predicted assemblage distributions with ggplot

plot_colourmap <- function(col_map, boundaries_poly, land_poly, uncertainty,
                           label_position, expand_scale, xlim, ylim,
                           x_breaks, y_breaks, expansion_const,
                           region_tag, tag_position, ytext_angle,
                           scalebar_position, scalebar_width, scalebar_ypad 
){
  
  require(ggplot2)
  
  mapPlot <- ggplot() +
    geom_sf(data = col_map, 
              aes(fill = assigned),
              col = NA,
              show.legend = F,
              na.rm = T) + 
    scale_fill_identity() +
    geom_sf(data = boundaries_poly,
            fill = NA, 
            col = 'black') +
    geom_sf(data = land_poly, 
            fill = 'grey30',
            col = 'grey20', 
            lwd = 0.25) +
    geom_sf(data = uncertainty, 
            col = 'black', 
            fill = NA,
            show.legend = F) +
    coord_sf(label_graticule = label_position,
             expand = expand_scale,
             xlim = xlim,
             ylim = ylim) +
    scale_x_continuous(name = '',
                       breaks = x_breaks,
                       expand = expansion_const) +
    scale_y_continuous(name = '',
                       breaks = y_breaks,
                       expand = expansion_const) +
    annotate('text', 
             label = paste0("bold(", region_tag, ")"),
             parse = T,
             size = 6,
             x = tag_position[1],
             y = tag_position[2]) +
    annotation_scale(bar_cols = c('grey40','grey10'),
                     text_cex = 1,
                     line_col = 'grey10',
                     location = scalebar_position,
                     width_hint = scalebar_width,
                     pad_y = unit(scalebar_ypad, 'cm')) +
    theme(panel.background = element_rect(fill = '#e0e0e0'),
          panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
          plot.margin = unit(c(0,0,0,0), 'lines'),
          panel.border = element_rect(fill = NA, colour = 'black'),
          axis.text = element_text(size = 12.5),
          axis.text.y.right = element_text(hjust = 0.5, angle = ytext_angle))
  
}

# ggplot colour boxplots

# illustrates variation in environmental predictors within and across
# assemblage types as colour-coded boxplots

boxplot_colour <- function(.data, x, y, colour_var, ylim = NULL, ytitle, xlabs) {
  
  require(ggplot2)
  
  x = enquo(x)
  y = enquo(y)
  colour_var = enquo(colour_var)
  
  ggplot(.data, aes(x = !!x, y = !!y, fill = !!colour_var)) +
    geom_boxplot(show.legend = F) +
    theme_classic() +
    scale_fill_identity() +
    labs(x = NULL, y = eval(expression(ytitle))) +
    scale_x_discrete(labels = xlabs) +
    scale_y_continuous(limits = ylim) +
    theme(text = element_text(size = 14), 
          axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12, "pt"), 
                                      hjust = 0.5),
          axis.ticks.length = unit(2, "mm"))
 
}

# arranges 4 colour-coded boxplots into 2 X 2 array

boxplot_4panel <- function(tl, tr, bl, br, bioregion) {
  
  require(ggplot2)
  require(gridExtra)
  require(gtable)
  
  # Coerce ggplot objects to graphical objects
  
  tl <- tl + theme(axis.text.x = element_blank())
  tl <- ggplotGrob(tl)
  tr <- tr + theme(axis.text.x = element_blank())
  tr <- ggplotGrob(tr)
  bl <- bl + theme(axis.text.x = element_text(angle = 50, hjust = 1))
  bl <- ggplotGrob(bl)
  br <- br + theme(axis.text.x = element_text(angle = 50, hjust = 1))
  br <- ggplotGrob(br)
  
  #Arrange grobs and plot in 2 X 2 array
  
  c1 <- rbind(tl, bl, size = "first") #bind/align plot elements of column 1
  c2 <- rbind(tr, br, size = "first") #bind/align plot elements of column 2
  boxplot <- cbind(c1, c2, size = "first") #bind/align plot elements of both rows
  
  # write to file
  
  tiff(filename = paste0("Output/", bioregion, "_EnvVariation.tiff"), 
       width = 8, 
       height = 8, 
       units = "in", 
       res = 300, 
       compression = "lzw")
  
  plot(boxplot)
  
  dev.off()
  
}
