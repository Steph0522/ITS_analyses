#if (!requireNamespace("devtools", quietly=TRUE))
 # install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggmsa")

library(ggmsa)
#protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
sequences <- "~/Documents/ITS_corredor/seqs_check.fasta"

png("align.png",width=12,height=4, units = "in", res=1200)
ggmsa(sequences, start =0, end = 680, char_width = 0.5, seq_name = T, font = NULL, color = "Chemistry_NT") + geom_seqlogo(color = "Chemistry_NT") + geom_msaBar()
dev.off()
