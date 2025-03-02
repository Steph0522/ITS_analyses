gg_ordiplot5 <- function(ord, groups, groups2, scaling = 1, choices = c(1,2), kind = c("sd", "se", "ehull"), conf=NULL, show.groups="all", ellipse = TRUE, label = FALSE, hull = FALSE, spiders = FALSE, pt.size = 3, plot=TRUE) {
  library(ggordiplots)
  library(vegan)
  x <- y <- cntr.x <- cntr.y <- Group <- NULL
  groups <- as.factor(groups)
  groups2 = as.factor(groups2)
  if (show.groups[1]=="all") {
    show.groups <- as.vector(levels(groups))
    
  }
  
  # Get site coordinates to plot.
  df_ord <- ord$Vectors %>% column_to_rownames(var = "SampleID") #scaling=scaling, choices=choices)
  axis.labels <- ord$ProportionExplained[choices]
  df_ord <- data.frame(x=df_ord[ , 1], y=df_ord[ , 2], Group=groups, Group2= groups2)
  
  # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
  
  # Get parameters from the ordiellipse function.
  if (is.null(conf)) {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", label = label)
  } else {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", conf = conf, label = label)
  }
  library(vegan)
  # Get points to plot for the ellipses.
  df_ellipse <- data.frame()
  for(g in show.groups) {
   df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                            vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
  }
  colnames(df_ellipse) <- c("x", "y", "Group")
  df_ellipse <- df_ellipse[ , c(3,1,2)]
  
  # Make data frame for hulls.
  rslt.hull <- vegan::ordihull(ord, groups = groups, scaling = scaling, choices = choices, show.groups = show.groups, draw = "none")
  df_hull <- data.frame()
  df_temp <- data.frame()
  for (g in show.groups) {
    x <- rslt.hull[[g]][ , 1]
    y <- rslt.hull[[g]][ , 2]
    Group <- rep(g, length(x))
    df_temp <- data.frame(Group = Group, x=x, y=y)
    df_hull <- rbind(df_hull, df_temp)
  }
  
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
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  plt <- ggplot2::ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group, shape=Group2), size = pt.size) +
    xlab(xlab) + ylab(ylab)
  
  # Add ellipses.
  if (ellipse == TRUE) {
    plt <- plt + geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add labels.
  if (label == TRUE) {
    plt <- plt + geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE)
  }
  
  # Add hulls.
  if (hull == TRUE) {
    plt <- plt + geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add spiders.
  if (spiders == TRUE) {
    plt <- plt + geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE)
  }
  
  plt <- plt + coord_fixed(ratio=1)
  
  # Plot?
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, df_mean.ord=df_mean.ord, df_ellipse=df_ellipse, df_hull=df_hull, df_spiders=df_spiders, plot=plt))
}
