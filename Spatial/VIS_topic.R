library(scatterpie)
vizAllTopics <- function(theta, Pos,
                         topicOrder=seq(ncol(theta)),
                         topicCols=rainbow(ncol(theta)),
                         groups = NA,
                         group_cols = NA,
                         r = max(0.4, max(Pos)/nrow(Pos)*4),
                         lwd = 0.01,
                         showLegend = TRUE,
                         plotTitle = NA,
                         overlay = NA) {
  
  ## check that theta and Pos are either data.frames or matrices
 # if( !is.matrix(theta) & !is.data.frame(theta) ){
 #   stop("`theta` must be a matrix or data.frame.")
 # }
  if( !is.matrix(Pos) & !is.data.frame(Pos) ){
    stop("`Pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  
  if( (any(!colnames(Pos) %in% c("x", "y")) == TRUE) | (dim(Pos)[2] != 2) ){
    stop("`Pos` must have exactly 2 columns named `x` and `y`.")
  }
  
  # pixel cell-type distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  colnames(theta_ordered) <- paste0("X", colnames(theta_ordered))
  
  # ensure that `theta` and `Pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta_ordered), rownames(Pos))
  pixels <- rownames(theta_ordered)[which(rownames(theta_ordered) %in% pixels)]
  
  # add columns "x", "y" with document Positions from `Pos`
  theta_ordered_Pos <- merge(data.frame(theta_ordered),
                             data.frame(Pos), by=0)
  rownames(theta_ordered_Pos) <- theta_ordered_Pos[,"Row.names"]
  ## make sure pixels in the original order before the merge
  theta_ordered_Pos <- theta_ordered_Pos[pixels,]
  
  # first column after merge is "Row.names", last two are "x" and "y"
  # problem is that data frame will replace "-" and " " with "."
  topicColumns <- colnames(theta_ordered_Pos)[2:(dim(theta_ordered_Pos)[2]-2)]
  
  # color of piechart groups (lines of piechart):
  if (is.na(groups[1])) {
    groups <- rep("0", dim(theta_ordered_Pos)[1])
    theta_ordered_Pos$Pixel.Groups <- groups
  } else {
    theta_ordered_Pos$Pixel.Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c("0" = "gray")
  }
  
  message("Plotting scatterpies for ", dim(theta_ordered_Pos)[1], " pixels with ", length(topicColumns),
      " cell-types...this could take a while if the dataset is large.", "\n")
  
  if (!is.na(overlay[1])){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
      ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_Pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
        ) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_Pos,
                                  cols = topicColumns,
                                  legend_name = "CellTypes") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.Position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
}


