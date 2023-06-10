# Analisis ITS 
# Cliente: Curso de microbioma para R
#if (!requireNamespace("BiocManager", quietly=TRUE))
 # install.packages("BiocManager")
#BiocManager::install("dada2")
#----Llamado de librerias----
# BiocManager::install("dada2", force = T)
library(dada2)
library(ShortRead)
library(Biostrings)

#----Deteccion de primers----

# Ajuste de rutas y directorio de trabajo
# ruta que contiene los archivos fastq
ruta_fastq <- "/home/yendi/Documents/ITS_corredor/ITS_analysis/all_seqs/joined_trimmed/joined_trimmed/"

# Crear carpeta para guardar resultados de ITS y ajusta directorio de trabajo
ruta_resultados_ITS <- "/home/yendi/Documents/ITS_corredor/ITS_analysis/all_seqs/joined_trimmed/joined_trimmed/resultados/"
setwd(ruta_resultados_ITS)

fnFs <- sort(list.files(ruta_fastq, pattern = "_R1_001.fastq.gz", full.names = TRUE))
#fnRs <- sort(list.files(ruta_fastq, pattern = "_R2_001.fastq.gz", full.names = TRUE))
fnFs
# Se obtuvieron los primers del articulo: Zhang2022
FWD <- "AACTTTYRRCAAYGGATCWCT"  
REV <- "AGCCTCCGCTTATTGATATGCTTAART" 

allOrients <- function(primer) {
  # Crear orientaciones de los primers
  # Biostrings funciona con objetos DNAString en lugar de vectores de caracteres
  require(Biostrings)
  dna <- DNAString(primer) 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convertir de nuevo a vector de caracteres
}

# Determinar todas las orientaciones posibles de los primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Crear un directorio nuevo con filtrado 
fnFs.filtN <- file.path(ruta_resultados_ITS, "filtN", basename(fnFs)) 
#fnRs.filtN <- file.path(ruta_resultados_ITS, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN,  maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
  # Cuenta el numero de lecturas en donde se encontraron primers
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Crear tabla para la deteccion de los primers
primeres_detected <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                           #FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                           REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))
                           #REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

primeres_detected

#write.csv(primeres_detected, "primeres_detected.csv")

#remove primer
#cutadapt <- "/usr/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
#system2(cutadapt, args = "--version") # Run shell commands from R


#path.cut <- file.path(ruta_fastq, "cutadapt")
#if(!dir.exists(path.cut)) dir.create(path.cut)
#fnFs.cut <- file.path(path.cut, basename(fnFs))
#fnRs.cut <- file.path(path.cut, basename(fnRs))

#FWD.RC <- dada2:::rc(FWD)
#REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
#R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
#R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
#for(i in seq_along(fnFs)) {
 # system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
  #                           "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
   #                          fnFs.filtN[i], fnRs.filtN[i], "--minimum-length=1")) # input files
#}
#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                   #     primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                    #                                                               fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
#
#
#----Pipeline desde secuencias cortadas o limpias----

# Dado que en la secuencias no se detectaron primers, podemos continuar. 
# Trabajar con secuencias cortadas
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(ruta_fastq, pattern = "_R1_001.fastq.gz", full.names = TRUE))
#cutRs <- sort(list.files(ruta_fastq, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][3]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#----Analisis forward y reverse----
plotQualityProfile(cutFs[1:2])
#plotQualityProfile(cutRs[1:2])


#-----Filtrado y corte----------
filtFs <- file.path(ruta_fastq, "filtered", basename(cutFs))
#filtRs <- file.path(path.cut, "filtered", basename(cutRs))
filtFs

names(filtFs) <- sample.names
#names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs,  minLen = 50,
                     maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                     compress = T, multithread = T)
head(out)

#----Modelos de error----
errF <- learnErrors(filtFs, multithread = TRUE)
#errR <- learnErrors(filtRs, multithread = TRUE)

#png("error_model.png", 
 #   units = "in",
  #  height = 7,
   # width = 10,
    #res = 300)
plotErrors(errF, nominalQ = TRUE)
#dev.off()


#----Dereplicar----
derepFs <- derepFastq(filtFs, verbose = TRUE)
#derepRs <- derepFastq(filtRs, verbose = TRUE)


#----Inferencia de ASVs----
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = "pseudo")
#dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#----Pareamiento de secuencias----
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#----Tabla de secuencias----
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

#----Quitar quimeras----
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = T, 
                                    verbose = T)
dim(seqtab_nochim)

#---------Denoising stats-----------
getN <-function(x) sum(getUniques(x))
# ?getUniques
stats <- cbind(out, sapply(dadaFs, getN), 
               rowSums(seqtab_nochim))

colnames(stats) <- c("input", "filtered", "denoisedF",  "nonchim")

stats
#write.csv(stats, "denoising-stats-its.csv")

#----------Asignacion taxonomica-----------
ruta_clasificador <- "/home/yendi/Downloads/1E25CA4CC30A31C2E2B8CB2C89824C83D080A7F5A62E6263A0E95B37C6628067/sh_general_release_dynamic_29.11.2022_dev.fasta" 
taxa <- assignTaxonomy(seqtab_nochim, ruta_clasificador, multithread = TRUE, tryRC = TRUE)

# Visualizar lo que se genero despues de la asignacion
taxa_print <- taxa
rownames(taxa_print) <- NULL

head(taxa_print)
dim(taxa_print)

# Exportar objetos generados durante el preprocesamiento 
save( errF, dadaFs, seqtab_nochim, taxa,
     file = "ITS_dada2_results_joined_trimmed.RData")
write.csv(taxa, "~/Documents/ITS_corredor/ITS_analysis/all_seqs/joined_trimmed//taxonomy.csv")
write.csv(seqtab_nochim, "~/Documents/ITS_corredor/ITS_analysis/all_seqs/joined_trimmed/table.csv")
write.csv(stats, "~/Documents/ITS_corredor/ITS_analysis/all_seqs/joined_trimmed/stats.csv")


