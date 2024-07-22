gg_ordiplot4 <- function(ord, groups, groups2, scaling = 1, conf=NULL, show.groups="all", spiders = FALSE, pt.size = 3, plot=TRUE) {
  library(ggordiplots)
  library(vegan)
  x <- y <- cntr.x <- cntr.y <- Group <- NULL
  groups <- as.factor(groups)
  groups2 = as.factor(groups2)
  if (show.groups[1]=="all") {
    show.groups <- as.vector(levels(groups))
    
  }
  
  # Get site coordinates to plot.
  #df_ord <- ord$points #scaling=scaling, choices=choices)
  names<- rownames(ord)
  #df_ord<- df_ord %>% as.data.frame() #%>% tibble::column_to_rownames(var = "SampleID")
  #axis.labels <- ord$ProportionExplained
  df_ord <- data.frame(x=ord[ ,1], y=ord[ ,2], Group=groups, Group2= groups2)
  rownames(df_ord)<-names

  # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
  
  
  # Make a data frame for the spiders.
  df_spiders <- df_ord[-4]
  df_spiders$cntr.x <- NA
  df_spiders$cntr.y <- NA
  for (g in show.groups) {
    df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
  }
  df_spiders <- df_spiders[ , c(3,4,5,1,2)]
  df_spiders <- df_spiders[order(df_spiders$Group), ]
  df_spiders <- df_spiders[df_spiders$Group %in% show.groups, ]

  # Make basic ggplot with ellipses.
  xlab <- paste("NMDS1")
  ylab <- paste("NMDS2")
  plt <- ggplot2::ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group, shape=Group2), size = pt.size) +
    xlab(xlab) + ylab(ylab)
  

  # Add spiders.
  if (spiders == TRUE) {
    plt <- plt + geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE)
  }
  
  plt <- plt + coord_fixed(ratio=1)
  
  # Plot?
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, df_mean.ord=df_mean.ord, df_spiders=df_spiders, plot=plt))
}
