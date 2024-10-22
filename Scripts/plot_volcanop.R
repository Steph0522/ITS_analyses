library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

plot_volcanop <- function(x, taxa, cutoff.pval = 0.05) {
  # Asegúrate de que 'wi.ep' esté definido
  if (!"wi.ep" %in% colnames(x)) {
    stop("El dataframe debe contener la columna 'wi.ep'.")
  }
  
  # Calcula p.add y all.p
  p.add <- min(x$wi.ep[x$wi.ep > 0]) / 10
  called <- x$wi.ep <= cutoff.pval
  x$all.p <- x$wi.ep + p.add
  
  # Crea una nueva columna para el log10 de p-values
  x$log_pvalue <- -1 * log10(x$all.p)
  
  # Procesa los datos
  x_processed <- x %>%
    rownames_to_column(var = "Feature.ID") %>%
    inner_join(taxa, by = "Feature.ID") %>%
    mutate(taxa = case_when(
      str_detect(Taxon, "g__") ~ str_extract(Taxon, "(?<=g__)[^_;]+"),
      str_detect(Taxon, "f__") ~ str_extract(Taxon, "(?<=f__)[^_;]+"),
      str_detect(Taxon, "c__") ~ str_extract(Taxon, "(?<=c__)[^_;]+"),
      str_detect(Taxon, "o__") ~ str_extract(Taxon, "(?<=o__)[^_;]+"),
      TRUE ~ NA_character_
    ))
  
  # Filtra solo los significativos
  significant_x <- x_processed %>%
    filter(all.p < cutoff.pval)
  
  # Filtra los 3 taxones con los mayores effect sizes (rojos y azules)
  top_red_taxa <- significant_x %>%
    filter(diff.btw > 0) %>%
    top_n(3, diff.btw) %>%
    dplyr::select(diff.btw, log_pvalue, taxa)
  
  top_blue_taxa <- significant_x %>%
    filter(diff.btw < 0) %>%
    top_n(3, abs(diff.btw)) %>%
    dplyr::select(diff.btw, log_pvalue, taxa)
  
  # Combina los resultados
  top_taxa <- bind_rows(top_red_taxa, top_blue_taxa)
  
  # Crea el gráfico
  ggplot(x, aes(x = diff.btw, y = log_pvalue)) +
    geom_point(data = x[!called, ], color = "gray", size = 3, shape = 19) +  # Puntos no significativos en gris
    geom_point(data = x[called, ], aes(color = ifelse(diff.btw < 0, "blue", "red")), size = 3, shape = 19) +  # Significativos en azul o rojo
    geom_vline(xintercept = c(-1.5, 1.5), color = 'black', linetype = 'dashed') +
    geom_hline(yintercept = -1 * log10(cutoff.pval), color = 'black', linetype = 'dashed') +
    labs(x = expression("Median Log"[2]~" Difference"), 
         y = expression("-1 * Median Log"[10]~" p value")) +
    theme_minimal() +
    scale_color_manual(values = c("blue" = "blue", "red" = "red")) +  # Define los colores manualmente
    # Agregar anotaciones para los taxones seleccionados
    geom_text(data = top_taxa, aes(x = diff.btw, y = log_pvalue, label = taxa), 
              vjust = -0.5, color = "black", fontface = "italic", size=2) +  # Texto en cursiva
    annotate("text", x = min(x$diff.btw)+1, y = 0.1, 
             label = colnames(x)[grep("rab.win", colnames(x))[1]], 
             color = "black", size = 3, vjust = 0) +
    annotate("text", x = max(x$diff.btw)-1, y = 0.1, 
             label = colnames(x)[grep("rab.win", colnames(x))[2]], 
             color = "black", size = 3, vjust = 0) +
    theme(legend.position = "none")
}

# Uso de la función
 plot_volcano(x = x, taxa = taxa)

