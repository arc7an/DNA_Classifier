library(ape)
library(reldna)
library(parallel)
source("/home/artem/work/2018/classifier_on_other_genomes/calculate_sequence_phys_profile.R")
# e.coli_U00096.2_char <- read.GenBank(access.nb = 'U00096.2', as.character = T)$U00096.2
# e.coli_U00096.2_string <- paste0(e.coli_U00096.2_char, collapse = '')[[1]]
# 
# 
# 
# shifted_by <- 250
# extended_e.coli_U00096.2_string <- paste0(
#   substr(e.coli_U00096.2_string, nchar(e.coli_U00096.2_string) - shifted_by, nchar(e.coli_U00096.2_string)),
#   e.coli_U00096.2_string,
#   substr(e.coli_U00096.2_string, 0, shifted_by)
# )
# 
# extended_ecoli_string <- extended_e.coli_U00096.2_string
parallel_calculate_ep <- function(genome_file_location, part=1, strand="forward"){
  genome_string <- unlist(read.fasta(genome_file_location, as.string = F, seqonly = T))
  extended_genome_string <- paste0(substr(genome_string, (nchar(genome_string) - 400), nchar(genome_string)), genome_string, substr(genome_string, 1, 400), sep="")
  #pseudo_tss <- 1:250000
  
  no_cores = detectCores() - 1
  cl <- makeCluster(no_cores)
  
  
  # calculate_EP_on_interval <- function(tss_position, extended_string, shifted_by, zout) {
  #   p <- lseqspline1D(substr(extended_string, tss_position-250, tss_position+150), bound=c(50, 350), ref=251 )
  #   return (p$mpot[p$x %in% zout])
  # }
  
  # clusterEvalQ(cl, 'extended_ecoli_string')
  clusterEvalQ(cl, 'extended_genome_string')
  clusterEvalQ(cl, {library(reldna)})
  clusterExport(cl, "calculate_EP_on_interval")
  clusterExport(cl, "get_forward_substring")
  clusterExport(cl, "get_reverse_substring")
  
  #res <- parSapply(cl, X = pseudo_tss, FUN = function(x) calculate_EP_on_interval(x, extended_ecoli_string, 250, -480:239))
  #stopCluster(cl)
  
  #save(res, file = '/home/jane/Документы/Misha/mol_a_2018/parsapply_ep_ecoli1strand_1_.rda')
  
  bounds <- c(seq(1, nchar(extended_genome_string), 250000), nchar(extended_genome_string))
  for (i in part:(length(bounds) - 1)){
    pseudo_tss <- bounds[i]:bounds[i+1] + 400
    print(length(pseudo_tss))
    res <- parSapply(cl, X = pseudo_tss, FUN = function(x) calculate_EP_on_interval(x, unlist(strsplit(extended_genome_string, '')), c(250, 150), -480:239, strand))
    
   # assign( paste0('sliced_ep_', min(pseudo_tss), '_', max(pseudo_tss)), res)
    save( res, file = paste0('/home/artem/work/2018/classifier_on_other_genomes/calculated_EP_ecoli/',  paste0('sliced_ep_', min(pseudo_tss), '_', max(pseudo_tss)), '.rda'))
    rm(list = 'res')
    gc(verbose = FALSE)
  }
  stopCluster(cl)
}