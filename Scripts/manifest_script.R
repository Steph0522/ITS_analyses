setwd("/home/yendi/Documents/ITS_corredor/all_seqs/joined_trimmed_raw/joined_trimmed/")

library(qiime2R)
#?qiime2R::write_q2manifest()

write_q2manifest("/home/yendi/Documents/ITS_corredor/all_seqs/joined_trimmed_raw/joined_trimmed//manifest.txt",
                 "/home/yendi/Documents/ITS_corredor/all_seqs/joined_trimmed_raw/joined_trimmed/",
                 extension=".fastq.gz", paired=FALSE)

#update.packages("ggplot2")
