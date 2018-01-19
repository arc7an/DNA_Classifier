if(!require(ape)){stop('Code requires library "ape".\n')}
if(!require(dplyr)){stop('Code requires library "dplyr".\n')}

# used functions
printConfusionMatrix <- function(dirName) {
  setwd(dirName)
  files <- dir(dirName)
  output <- c()
  outputFactor <- c()
  for (f in files) {
    if (lapply(strsplit(f, '_'), tail, 1) == "matrix.Rdata"){
      load(f)
      spl <- strsplit(f, '_')
      res <- c(
        confusionMatrix_promoter_vs_other$overall[1],
        confusionMatrix_promoter_vs_other$byClass
      )
      factorRes <- c(
        paste(spl[[1]][1], "vs Promoter"),
        spl[[1]][4],
        spl[[1]][2],
        spl[[1]][3],
        spl[[1]][5]
      )
      names(factorRes)[1]  <- 'TypeOfModel'
      names(factorRes)[2]  <- 'LearningMethod'
      names(factorRes)[3]  <- 'NumberOfComponents'
      names(factorRes)[4]  <- 'TrainingProportion'
      names(factorRes)[5]  <- 'Seed'

      output <- rbind(output, res)
      outputFactor <- rbind(outputFactor, factorRes)
    }

  }
  rownames(output) <- NULL
  rownames(outputFactor) <- NULL
  return(data.frame(output, outputFactor))
}

filter_best_model <- function(models_folder, classifier_type, ...){
  conf_matrices <- printConfusionMatrix(models_folder)
  filtered <- conf_matrices[conf_matrices$TypeOfModel == classifier_type,]
  return(arrange(filtered, ... = ...)[1:10,])
}

sliding_window <- function(input, len, by=2){
  starts <- seq(1,length(input),by)
  tt <- lapply(starts, function(y) input[y:(y+(len-1))])
  sapply(tt, function(x) x[!is.na(x)])
}

# find best models by Accuracy
# best_150_50 <- lapply("Nonpromoter vs Promoter", function(x) filter_best_model('/home/jane/r_ML/interval_variation/interval_150_50_roc/', x, -Accuracy))
# bounds of genomes slices of e.coli
#bounds_e_coli <- c(1,  250001,  500001, 750001, 1000001, 1250001, 1500001, 1750001, 2000001, 2250001, 2500001, 2750001, 3000001, 3250001, 3500001, 3750001, 4000001, 4250001, 4500001, 4640176)
# load pca data set that was used for model training. Variable name: prcomp_res
load('/home/artem/work/2016/iteb/pca_data_set_200_150_50.Rdata')

# go to directory of calculated EP slices of the genome
setwd("/home/artem/work/2018/classifier_on_other_genomes/calculated_EP_ecoli/")
epSlices <- list.files()
# load profiles of dynamic characteristics and GC and cut the end for proper cancatenation
load("/home/artem/work/2018/classifier_on_other_genomes/drive-download-20171201T095355Z-001/df_whole_ecoli_dynamic_characteristics.rda")
df_whole_ecoli_dynamic_characteristics <- df_whole_ecoli_dynamic_characteristics[1:(length(df_whole_ecoli_dynamic_characteristics[,1]) - 193),] #75 193

make_prediction_for_slice <- function(ep_slice, model, strand=1) {
# load chosen model Nonpropmoter-Promoter and adjust apriori probabilities
#  load("/home/artem/work/2017/grant_R_code_and_data/interval_150_50/Nonpromoter_50_0.9_nb_1116.Rdata")
  gc(verbose = F)
  load(model)
  fit_model$finalModel$apriori[1] <- 0.9999
  fit_model$finalModel$apriori[2] <- 0.0001
  bounds <- as.integer(strsplit(strsplit(ep_slice, ".r")[[1]][1], "_")[[1]][3:4])
  print("Bounds:")
  print(bounds)
  # load slice of EP profile, variable name - res
  load(ep_slice)
  res <- scale(res)
  # only for end EP slice
  # res <- res[,1:(ncol(res) - 294)]
# split data for one strand
  if (strand == 1){
    df_whole_ecoli_dynamic_characteristics <- df_whole_ecoli_dynamic_characteristics[bounds[1]:bounds[2], c(1,3,5)]
    # for right end only
    # df_whole_ecoli_dynamic_characteristics <- df_whole_ecoli_dynamic_characteristics[bounds[1]:(bounds[2] - 294), c(1,3,5)]
  }
  else {
    df_whole_ecoli_dynamic_characteristics <- df_whole_ecoli_dynamic_characteristics[bounds[1]:bounds[2], c(2,4,6)]
    # for right end only
    # df_whole_ecoli_dynamic_characteristics <- df_whole_ecoli_dynamic_characteristics[bounds[1]:(bounds[2] - 294), c(2,4,6)]
    res <- res[nrow(res):1,]
  }

  sliced_e01 <- sliding_window(df_whole_ecoli_dynamic_characteristics[,1], 201, 1)
  sliced_e01 <- t(do.call(rbind, sliced_e01))
  sliced_d1 <- sliding_window(df_whole_ecoli_dynamic_characteristics[,2], 201, 1)
  sliced_d1 <- t(do.call(rbind, sliced_d1))
  sliced_gc1 <- sliding_window(df_whole_ecoli_dynamic_characteristics[,3], 201, 1)
  sliced_gc1 <- t(do.call(rbind, sliced_gc1))

  rm(df_whole_ecoli_dynamic_characteristics)
  gc(verbose = F)

  combined <- rbind(sliced_e01, sliced_d1, res, sliced_gc1)

  rm(sliced_e01, sliced_d1, sliced_gc1, res)
  gc(verbose = F)

  combined_pca <- t(combined) %*% prcomp_res$rotation - prcomp_res$center
  rm(combined)
  gc(verbose = F)

  prediction <- predict(fit_model, combined_pca, type="prob")

  # print("Amount of Promoters with probability > 50%:")
  # print(length(prediction[prediction$Promoter > 0.5,2]))
  # print("Amount of Promoters with probability > 75%:")
  # print(length(prediction[prediction$Promoter > 0.75,2]))
  # print("Amount of Promoters with probability > 90%:")
  # print(length(prediction[prediction$Promoter > 0.9,2]))

  rm(combined_pca)
  gc(verbose = F)
  # filename <- cat("/home/artem/work/2018/MCE/predictions/prediction_", bounds[1], "_", bounds[2], ".Rdata", sep="")
  save(prediction, file = paste0("/home/artem/work/2018/MCE/predictions/reverse_strand/",  paste0("prediction_", bounds[1], "_", bounds[2]), '.Rdata'))
  rm(prediction)
  rm(combined_pca)
  gc(verbose = F)
  # save(prediction, file=filename)
}

# predictions_list <- c()

# #for (i in 1:length(epSlices)) {
# for (i in c(1)) {
#   prediction <- make_prediction_for_slice(epSlices[i], "/home/artem/work/2017/grant_R_code_and_data/interval_150_50/Nonpromoter_50_0.9_nb_1116.Rdata", 1)
#   predictions_list <- c(predictions_list, prediction)
# }

# predictions_result <- lapply(epSlices, function(slice) make_prediction_for_slice(slice, "/home/artem/work/2018/MCE/Nonpromoter_50_0.9_nb_1116.Rdata", 1))
lapply(epSlices[17], function(slice) make_prediction_for_slice(slice, "/home/artem/work/2018/MCE/Nonpromoter_50_0.9_nb_1116.Rdata", 2))
