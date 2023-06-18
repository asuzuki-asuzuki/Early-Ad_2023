


aes_string2 <- function(...){
  args <- lapply(list(...), function(x) sprintf("`%s`", x))
  do.call(ggplot2::aes_string, args)
}


plot_spat_point_layer_ggplot = function(ggobject,
                                        sdimx = NULL,
                                        sdimy = NULL,
                                        cell_locations_metadata_selected,
                                        cell_locations_metadata_other,
                                        cell_color = NULL,
                                        color_as_factor = T,
                                        cell_color_code = NULL,
                                        cell_color_gradient = c('blue', 'white', 'red'),
                                        gradient_midpoint = NULL,
                                        gradient_limits = NULL,
                                        select_cell_groups = NULL,
                                        select_cells = NULL,
                                        point_size = 2,
                                        point_alpha = 1,
                                        point_border_col = 'lightgrey',
                                        point_border_stroke = 0.1,
                                        show_cluster_center = F,
                                        show_center_label = T,
                                        center_point_size = 4,
                                        center_point_border_col = 'black',
                                        center_point_border_stroke = 0.1,
                                        label_size = 4,
                                        label_fontface = 'bold',
                                        show_other_cells = T,
                                        other_cell_color = 'lightgrey',
                                        other_point_size = 1,
                                        show_legend = TRUE

) {

  ## specify spatial dimensions first
  if(is.null(sdimx) | is.null(sdimy)) {

    warning("plot_method = ggplot, but spatial dimensions for sdimx and/or sdimy are not specified. \n
            It will default to the 'sdimx' and 'sdimy' ")
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  ## ggplot object
  pl = ggobject

  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, sdimy),
                                   color = other_cell_color, show.legend = F, size = other_point_size, alpha = point_alpha)
  }


  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor


  # cell color default
  if(is.null(cell_color)) {

    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                   aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21,
                                   fill = cell_color, size = point_size,
                                   stroke = point_border_stroke, color = point_border_col,
                                   alpha = point_alpha)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      cell_locations_metadata_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy, fill = 'temp_color'),
                                     show.legend = show_legend, shape = 21,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)

    }

  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)

    } else {

      class_cell_color = class(cell_locations_metadata_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {
        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = cell_locations_metadata_selected[[cell_color]]
          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          cell_locations_metadata_selected[[cell_color]] = limit_numeric_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21,
                                       size = point_size,
                                       color = point_border_col,
                                       stroke = point_border_stroke,
                                       alpha = point_alpha)



      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = stats::median(get('sdimx')),
                                                                      center_2 = stats::median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke,
                                       alpha = point_alpha)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', fill = cell_color),
                                         color = center_point_border_col, stroke = center_point_border_stroke,
                                         size = center_point_size, shape = 21,
                                         alpha = point_alpha)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface)
        }

      }

      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(cell_locations_metadata_selected[[cell_color]])
        }

		print(gradient_limits)
        pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
#                                                 mid = cell_color_gradient[[2]],
                                                 high = cell_color_gradient[[2]],
                                                 limit = gradient_limits,   midpoint = gradient_midpoint)

      }
    }
  }
  return(pl)
}




## USE
MyspatPlot2D_single = function(gobject,
                             show_image = F,
                             gimage = NULL,
                             image_name = 'image',
                             sdimx = 'sdimx',
                             sdimy = 'sdimy',
                             spat_enr_names = NULL,
                             cell_color = NULL,
                             color_as_factor = T,
                             cell_color_code = NULL,
                             cell_color_gradient = c('blue', 'white', 'red'),
                             gradient_midpoint = NULL,
                             gradient_limits = NULL,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             point_shape = c('border', 'no_border', 'voronoi'),
                             point_size = 3,
                             point_alpha = 1,
                             point_border_col = 'black',
                             point_border_stroke = 0.1,
                             show_cluster_center = F,
                             show_center_label = F,
                             center_point_size = 4,
                             center_point_border_col = 'black',
                             center_point_border_stroke = 0.1,
                             label_size = 4,
                             label_fontface = 'bold',
                             show_network = F,
                             spatial_network_name = 'Delaunay_network',
                             network_color = NULL,
                             network_alpha = 1,
                             show_grid = F,
                             spatial_grid_name = 'spatial_grid',
                             grid_color = NULL,
                             show_other_cells = T,
                             other_cell_color = 'lightgrey',
                             other_point_size = 1,
                             other_cells_alpha = 0.1,
                             coord_fix_ratio = NULL,
                             title = NULL,
                             show_legend = T,
                             legend_text = 8,
                             legend_symbol_size = 1,
                             background_color = 'white',
                             vor_border_color = 'white',
                             vor_max_radius = 200,
                             vor_alpha = 1,
                             axis_text = 8,
                             axis_title = 8,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'spatPlot2D_single'
) {


  ## giotto image ##
  if(show_image == TRUE) {
    if(!is.null(gimage)) gimage = gimage
    else if(!is.null(image_name)) {
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exists')
    }
  }


  ## point shape ##
  point_shape = match.arg(point_shape, choices = c('border', 'no_border', 'voronoi'))

  ## get spatial cell locations
  cell_locations  = gobject@spatial_locs

  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = select_spatialGrid(gobject, spatial_grid_name)
  } else {
    spatial_grid = NULL
  }

  ## get cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cell_metadata
  }

  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]

    # if specific cells are selected
    # cell_locations_metadata = cell_locations_metadata_selected

  } else if(is.null(select_cells)) {

    cell_locations_metadata_selected = cell_locations_metadata
    cell_locations_metadata_other = NULL

  }


  # data.table and ggplot variables
  sdimx_begin = sdimy_begin = sdimx_end = sdimy_end = x_start = x_end = y_start = y_end = NULL


  ### create 2D plot with ggplot ###
  #cat('create 2D plot with ggplot \n')


  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_bw()


  ## plot image ##
  if(show_image == TRUE & !is.null(gimage)) {
    pl = plot_spat_image_layer_ggplot(ggplot = pl,
                                      gobject = gobject,
                                      gimage = gimage,
                                      sdimx = sdimx,
                                      sdimy = sdimy)
  }

  ## plot spatial network
  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) network_color = 'red'
    pl <- pl + ggplot2::geom_segment(data = spatial_network, aes(x = sdimx_begin, y = sdimy_begin,
                                                                 xend = sdimx_end, yend = sdimy_end),
                                     color = network_color, size = 0.5, alpha = network_alpha)
  }


  ## plot spatial grid
  if(!is.null(spatial_grid) & show_grid == TRUE) {
    if(is.null(grid_color)) grid_color = 'black'
    pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end,
                                                           ymin = y_start, ymax = y_end),
                                  color = grid_color, fill = NA)
  }

  ## plot point layer
  if(point_shape == 'border') {
    pl = plot_spat_point_layer_ggplot(ggobject = pl,
                                      sdimx = sdimx,
                                      sdimy = sdimy,
                                      cell_locations_metadata_selected = cell_locations_metadata_selected,
                                      cell_locations_metadata_other = cell_locations_metadata_other,
                                      cell_color = cell_color,
                                      color_as_factor = color_as_factor,
                                      cell_color_code = cell_color_code,
                                      cell_color_gradient = cell_color_gradient,
                                      gradient_midpoint = gradient_midpoint,
                                      gradient_limits = gradient_limits,
                                      select_cell_groups = select_cell_groups,
                                      select_cells = select_cells,
                                      point_size = point_size,
                                      point_alpha = point_alpha,
                                      point_border_stroke = point_border_stroke,
                                      point_border_col = point_border_col,
                                      show_cluster_center = show_cluster_center,
                                      show_center_label = show_center_label,
                                      center_point_size = center_point_size,
                                      center_point_border_col = center_point_border_col,
                                      center_point_border_stroke = center_point_border_stroke,
                                      label_size = label_size,
                                      label_fontface = label_fontface,
                                      show_other_cells = show_other_cells,
                                      other_cell_color = other_cell_color,
                                      other_point_size = other_point_size,
                                      show_legend = show_legend)
  } else if(point_shape == 'no_border') {
    pl = plot_spat_point_layer_ggplot_noFILL(ggobject = pl,
                                             sdimx = sdimx,
                                             sdimy = sdimy,
                                             cell_locations_metadata_selected = cell_locations_metadata_selected,
                                             cell_locations_metadata_other = cell_locations_metadata_other,
                                             cell_color = cell_color,
                                             color_as_factor = color_as_factor,
                                             cell_color_code = cell_color_code,
                                             cell_color_gradient = cell_color_gradient,
                                             gradient_midpoint = gradient_midpoint,
                                             gradient_limits = gradient_limits,
                                             select_cell_groups = select_cell_groups,
                                             select_cells = select_cells,
                                             point_size = point_size,
                                             point_alpha = point_alpha,
                                             show_cluster_center = show_cluster_center,
                                             show_center_label = show_center_label,
                                             center_point_size = center_point_size,
                                             label_size = label_size,
                                             label_fontface = label_fontface,
                                             show_other_cells = show_other_cells,
                                             other_cell_color = other_cell_color,
                                             other_point_size = other_point_size,
                                             show_legend = show_legend)

  } else if(point_shape == 'voronoi') {

    pl = plot_spat_voronoi_layer_ggplot(ggobject = pl,
                                        sdimx = sdimx,
                                        sdimy = sdimy,
                                        cell_locations_metadata_selected = cell_locations_metadata_selected,
                                        cell_locations_metadata_other = cell_locations_metadata_other,
                                        cell_color = cell_color,
                                        color_as_factor = color_as_factor,
                                        cell_color_code = cell_color_code,
                                        cell_color_gradient = cell_color_gradient,
                                        gradient_midpoint = gradient_midpoint,
                                        gradient_limits = gradient_limits,
                                        select_cell_groups = select_cell_groups,
                                        select_cells = select_cells,
                                        point_size = point_size,
                                        point_alpha = point_alpha,
                                        show_cluster_center = show_cluster_center,
                                        show_center_label = show_center_label,
                                        center_point_size = center_point_size,
                                        label_size = label_size,
                                        label_fontface = label_fontface,
                                        show_other_cells = show_other_cells,
                                        other_cell_color = other_cell_color,
                                        other_point_size = other_point_size,
                                        background_color = background_color,
                                        vor_border_color = vor_border_color,
                                        vor_max_radius = vor_max_radius,
                                        vor_alpha = vor_alpha,
                                        show_legend = show_legend)

  }



  ## adjust theme settings
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                            legend.title = element_blank(),
                            legend.text = element_text(size = legend_text),
                            axis.title = element_text(size = axis_title),
                            axis.text = element_text(size = axis_text),
                            panel.grid = element_blank(),
                            panel.background = element_rect(fill = background_color))

  ## change symbol size of legend
  if(color_as_factor == TRUE) {
    if(point_shape %in% c('border', 'voronoi')) {
      pl = pl + guides(fill = guide_legend(override.aes = list(size = legend_symbol_size)))
    } else if(point_shape == 'no_border') {
      pl = pl + guides(color = guide_legend(override.aes = list(size = legend_symbol_size)))
    }
  }


  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  # provide x, y and plot titles
  if(is.null(title)) title = cell_color
  pl <- pl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates', title = title)


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }


  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }

}



