#setting working directory where the script, 2 octave scripts and .Rdata files are stored
setwd('/home/artem/work/2016/iteb/')

# removing saved working space
rm(list=ls())

#loading R libraries
library(R.matlab)
library(seqinr)
library(pca3d)
library(rgl)
library(factoextra)
## calculation of dynamical properties (activation energy and size of open states) from within R using 2 octave script
#writing in R E.coli K12 MG1655 genome from fasta file (must be copied to the working directory)


window_size <- 200
interval <- c(150, 50)

ecoli<-unlist(read.fasta('e.coli_U00096.2.fasta', as.string = F, seqonly = T))

#tranformation genome into character since for function below usage
ecoli_char<-unlist(strsplit(ecoli, ''))
#creating function for dynamical properties calculation

#loading data on sequences of different types (promoters, non-promoters, genes, islands, and lowscore) from .Rdata files (must be copied separetely)
load('spline_dataset_pro.Rdata')
load('spline_dataset_notpro.Rdata')
load('spline_dataset_gen.Rdata')
load('spline_dataset_isl.Rdata')
load('dataset_lowscore.Rdata')
dataset_lowscore <- dataset_lowscore[1:2000]

# extracting data on all promoters and on experimentaly found ones - including previosely calculated electrostatic potential profiles

calculate_dynamic_characteristics <- function(seq, interval_size) {
  if (missing(seq))
    stop("Need to specify sequence (as a vector of chars)")
  
  if (missing(interval_size))
    stop("Need to specify interval size")
  
  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  
  seq<-toupper(seq)
  # seq<-c(seq, seq[2:(interval_size)])
  
  a <- 3.4*10^(-10)
  I <- c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K <- c(227, 155, 220, 149)*10^(-20)
  V <- c(2.09, 1.43, 3.12, 2.12)*10^(-20)
  tau <- c(127, 99, 140, 84)
  
  csA<-cumsum(seq=='A') 
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')
  
  countA <- csA[interval_size:length(csA)]-c(0, csA[1:(length(csA)-interval_size)])
  countT <- csT[interval_size:length(csT)]-c(0, csT[1:(length(csT)-interval_size)])
  countG <- csG[interval_size:length(csG)]-c(0, csG[1:(length(csG)-interval_size)])
  countC <- csC[interval_size:length(csC)]-c(0, csC[1:(length(csC)-interval_size)])
  
  M<-cbind(countA, countT, countG, countC)/interval_size
  M_comp<-cbind(countT, countA, countC, countG)/interval_size
  M_comp<-apply(t(M_comp), 1, rev) 
  Is<-as.numeric(M%*%I)#! numeric conversion
  Ks<-as.numeric(M%*%K)
  Vs<-as.numeric(M%*%V)
  
  E01<-(8*(Ks*Vs)^0.5)* 6E23 / 4184
  d1<-((Ks*a^2)/Vs)^(0.5)/a;
  c1<-(Ks*a^2/Is)^0.5
  m1<-E01/c1/6.011E-26
  taus1<-as.numeric(M%*%tau) #!as.numeric conversion
  gc = M[,3] + M[,4]
  
  Is<-as.numeric(M_comp%*%I)#! numeric conversion
  Ks<-as.numeric(M_comp%*%K)
  Vs<-as.numeric(M_comp%*%V)
  
  E02<- 8*(Ks*Vs)^0.5  * 6E23 / 4184;
  d2<-((Ks*a^2)/Vs)^(0.5)/a;
  c2<-(Ks*a^2/Is)^0.5;
  m2<-E02/c2/6.011E-26;
  taus2<-as.numeric(M_comp%*%tau)
  
  dynchars_return<-list(E01=E01, d1=d1, c1=c1, m1=m1, taus1=taus1, gc=gc, E02=E02, d2=d2, c2=c2, m2=m2, taus2=taus2)
  return(dynchars_return)
}

prepare_sequences <- function(dataset, experimental_flag = F, promoter_names = c()) {
  tss_vector <- c()
  strands_vector <- c()
  ep_vector <- c()
  names_vector <- c()
  
  if (experimental_flag == T){
    for (i in 1:length(dataset)){
      if (dataset[[i]]$evidence == 'experimental'){ 
        tss_vector <- c(tss_vector, dataset[[i]]$tss)
        strands_vector <- c(strands_vector, dataset[[i]]$strand)
        ep_vector <- rbind(ep_vector,  dataset[[i]]$mpot)  
        names_vector<-c(names_vector, promoter_names[i])
      }
    }    
    return (list(tss_vector, strands_vector, ep_vector, names_vector))
  }
  else {
    for (i in 1:length(dataset)) {
        tss_vector <- c(tss_vector, dataset[[i]]$tss)
        strands_vector <- c(strands_vector, dataset[[i]]$strand)
        ep_vector <- rbind(ep_vector,  dataset[[i]]$mpot) 
    }
    return(list(tss_vector, strands_vector, ep_vector))
  }
}

#calculation the properties for a given genome with sliding window 200 nt

dynamic_characteristics<-calculate_dynamic_characteristics(ecoli_char, window_size)
E01<-dynamic_characteristics$E01
E02<-dynamic_characteristics$E02

d1<-dynamic_characteristics$d1
d2<-dynamic_characteristics$d2

GC<-dynamic_characteristics$gc #name of the variable is from before

cut_characteristics_dataset <- function(data_set, interval = c(150, 50)) {
  for (i in c('E0Forward', 'E0Reverse', 'dForward', 'dReverse', 'gcForward', 'gcReverse', 'EP')) {
    assign(i, c())
  }
  if (sum(interval) == window_size){
    for (i in 1:length(data_set[[1]])) {
      if (data_set[[2]][i] == 'forward') {
        E0Forward <- rbind(E0Forward, E01[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]][i] + interval[2])])
        dForward <- rbind(dForward, d1[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]] + interval[2])])
        gcForward<-rbind(gcForward, GC[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]] + interval[2])])
      } else {
        E0Reverse<-rbind(E0Reverse, E02[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]][i] + interval[2])])
        dReverse<-rbind(dReverse, d2[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]][i] + interval[2])])
        gcReverse<-rbind(gcReverse, GC[as.numeric(data_set[[1]][i]-interval[1]):(data_set[[1]][i] + interval[2])])
      }
      EP <- rbind(EP, data_set[[3]][i, ])
    }
    #merging matrices for forward and reverse strands  together
    E0Full<-rbind(E0Forward, E0Reverse)
    dFull<-rbind(dForward, dReverse)
    gcFull<-rbind(gcForward, gcReverse)
    return(list(E0Full, dFull, gcFull, EP))
  }
  else {print("The interval is not equal to window size.")}
}

promoter_names <- names(dataset_pro)
gene_names <- names(dataset_gen)
promoter_characteristics <- prepare_sequences(dataset = dataset_pro, experimental_flag = T, promoter_names = promoter_names)
non_promoter_characteristics <- prepare_sequences(dataset = dataset_notpro, experimental_flag = F)
gene_characteristics <- prepare_sequences(dataset = dataset_gen, experimental_flag = F)
island_characteristics <- prepare_sequences(dataset = dataset_isl, experimental_flag = F)
lowscore_characteristics <- prepare_sequences(dataset = dataset_lowscore, experimental_flag = F)

cut_promoters <- cut_characteristics_dataset(promoter_characteristics, interval)
cut_non_promoters <- cut_characteristics_dataset(non_promoter_characteristics, interval)
cut_genes <- cut_characteristics_dataset(gene_characteristics, interval)
cut_islands <- cut_characteristics_dataset(island_characteristics, interval)
cut_lowscores <- cut_characteristics_dataset(lowscore_characteristics, interval)


#setting names
experimental_promoter_names <- sapply(promoter_characteristics[4], as.character)[,1]
make_names <- function(list_of_characteristics, list_of_names) {
  rownames(list_of_characteristics) <- list_of_names
  return(list_of_characteristics)
}

cut_promoters <- lapply(cut_promoters, make_names, list_of_names = experimental_promoter_names)
cut_non_promoters <- lapply(cut_non_promoters, make_names, list_of_names = paste0('Non_promoter_', 1:nrow(cut_non_promoters[[1]])))
cut_genes <- lapply(cut_genes, make_names, list_of_names = paste0('Gene_', 1:nrow(cut_genes[[1]])))
cut_islands <- lapply(cut_islands, make_names, list_of_names = paste0('Island_', 1:nrow(cut_islands[[1]])))
cut_lowscores <- lapply(cut_lowscores, make_names, list_of_names = paste0('Lowscore_', 1:nrow(cut_lowscores[[1]])))

#PCA

pca_set <- rbind(do.call(cbind, cut_promoters),
                 do.call(cbind, cut_non_promoters),
                 do.call(cbind, cut_genes),
                 do.call(cbind, cut_islands),
                 do.call(cbind, cut_lowscores)
)
parameters <- paste0(c(window_size, interval), collapse ="_")
save(pca_set, file=paste0('pre_pca_data_set_', parameters,'.Rdata'))

sc_e0 <- scale(pca_set[,1:(sum(interval)+1)])
sc_d <- scale(pca_set[,(sum(interval)+2) : (2*sum(interval) + 2)])
sc_gc <- scale(pca_set[,(2*sum(interval)+3) : (3*sum(interval) + 3)])
sc_ep <- scale(pca_set[,(3*sum(interval)+4) : (length(pca_set[1,]))])

scaled_pca_set <- cbind(sc_e0, sc_d,sc_ep,sc_gc)
save(scaled_pca_set, file=paste0('scaled_pre_pca_data_set_', parameters,'.Rdata'))

set.seed(999)
prcomp_res <- prcomp(scaled_pca_set)
scaled_pca_set <- as.data.frame(scaled_pca_set)

labels <- c(
  rep('Promoters', length(promoter_characteristics[[1]])),
  rep('Non-promoters', length(non_promoter_characteristics[[1]])),
  rep('Genes', length(gene_characteristics[[1]])),
  rep('Islands', length(island_characteristics[[1]])),
  rep('Lowscore', length(lowscore_characteristics[[1]]))
)

scaled_pca_set$Type <- labels

autoplot(prcomp_res, data=scaled_pca_set, colour='Type')


save(prcomp_res, file=paste0('pca_data_set_', parameters,'.Rdata'))

# PCA Analysis
eigen_values <- (prcomp_res$sdev)^2

# Variances in percentage
variance <- eigen_values*100/sum(eigen_values)

# Cumulative variances
pca_cumulative_varience <- cumsum(variance)
png('cumulative_varience_plot.png', height=1250, width=1250, res=130)
plot(pca_cumulative_varience, type='l', main='Cumulative variance', ylab='Cumulative variance (%)', xlab='Principal components')
vertical_line <- 100 # 65 
horizontal_line <- 99
abline(h=horizontal_line, v=vertical_line, col=12)
text(vertical_line-15, horizontal_line ,paste('Number of \nPCs=', vertical_line,'\n' ,horizontal_line, '% of \nvariance \nretained'), srt = 0.2, pos = 3)

dev.off()

types_to_colors<-c(
  rep('red', length(promoter_characteristics[[1]])),
  rep ('blue', length(non_promoter_characteristics[[1]])),
  rep('green', length(gene_characteristics[[1]])),
  rep('orange', length(dataset_isl)),
  rep('magenta', length(lowscore_characteristics[[1]]))
)

plot3d(
  prcomp_res$x[,1:3],
  col=types_to_colors,
  xlab="PC1", ylab="PC2", zlab="PC3",
  size=3, lwd=5, radius = 4
#  legend("plot")  
  )

# # # SUPERVISED MACHINE LEARNING

registerDoMC(cores = 7)
#load("~/work/2016/iteb/promoters_nonpromoters_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
#load("~/work/2016/iteb/promoters_nonpromoters_initial_ae_data.Rdata")
#load("new_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
#load("new_no_factor_8XII_mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200.Rdata")
# df <- as.data.frame(mixt_5comps_after_pca_on_4_props_ae_size_mpots_gc200)
#load('princ.return.5comps.4props.Rdata')
df <- as.data.frame(prcomp_res$x)
#df<-cbind(habillage_5components_4props, df)
promoters <- df[1:699,]
# non_promoters <- df[700:2579,]
# genes <- df[2580:6006,]
# islands <- df[6007:8234,]
# lowscore <- df[8235:10234,]
# non_promoters <- df[700:(700+699),]
# genes <- df[2580:(2580+699),]
# islands <- df[6007:(6007+699),]
# lowscore <- df[8235:(8235+699),]
non_promoters <- df[sample(700:2579),][1:699,]
genes <- df[sample(2580:6006),][1:699,]
islands <- df[sample(6007:8234),][1:699,]
lowscore <- df[sample(8235:10234),][1:699,]


sequencesTypes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Lowscore', nrow(lowscore))))
promoters_vs_lowscore <- cbind(sequencesTypes, rbind(promoters, lowscore))

sequencesTypes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Nonpromoter', nrow(non_promoters))))
promoters_vs_non_promoters <- cbind(sequencesTypes, rbind(promoters, non_promoters))

sequencesTypes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Island', nrow(islands))))
promoters_vs_islands <- cbind(sequencesTypes, rbind(promoters, islands))

sequencesTypes<-as.factor(c(rep('Promoter', nrow(promoters)), rep('Gene', nrow(genes))))
promoters_vs_genes <- cbind(sequencesTypes, rbind(promoters, genes))

#for (data_set in c(promoters_vs_lowscore, promoters_vs_non_promoters, promoters_vs_islands, promoters_vs_genes)) {

makeClassificationModel <- function(promoter_vs_other, numberOfPCs, trainingSetProportion, MLMethod) {
  promoter_vs_other <- as.data.frame(promoter_vs_other)[, 1:(numberOfPCs + 1)]
  inTraining <- createDataPartition(
    promoter_vs_other$sequencesTypes,
    p = trainingSetProportion,
    list = F)
  training <- promoter_vs_other[inTraining,]
  testing  <- promoter_vs_other[-inTraining,]
  fitControl <- trainControl(method = "repeatedcv", 
                             number = 10, 
                             repeats = 15, 
                             allowParallel = T, 
                             classProbs = T, 
                             summaryFunction = twoClassSummary,
                             savePredictions = T
  )
  fit_model <- train(sequencesTypes ~ .,
                     data = training,
                     method = MLMethod,#nb, lda, nbDiscrete, nbSearch
                     tuneLength = 15,
                     trControl = fitControl,
                     positive = "Promoter"
  )
  predictionClasses <- predict(fit_model, newdata = testing)
  confusionMatrix_promoter_vs_other <- confusionMatrix(data = predictionClasses, testing$sequencesTypes)
  outputFileName <- paste(as.character(promoter_vs_other[,1][700]), numberOfPCs, trainingSetProportion, MLMethod, sep='_')
  confusionName <- paste(outputFileName, '_confusion_matrix', '.Rdata', sep='')
  outputFileName <- paste(outputFileName, '.Rdata', sep='')
  print(outputFileName)
  save(promoter_vs_other, file=outputFileName)
  save(confusionMatrix_promoter_vs_other, file=confusionName)
  
}

listOfSets = list(promoters_vs_genes, promoters_vs_islands, promoters_vs_non_promoters, promoters_vs_lowscore)
# listOfSets = list(promoters_vs_genes)
listOfPCs = list(50, 100, 150)
# listOfPCs = list(100)
listOfLearningSetProportion = list(0.7, 0.8, 0.9)
# listOfLearningSetProportion = list(0.7)
# listOfMethods = list("nb", "rf")
listOfMethods = list("nb")
set.seed(17)
for (s in listOfSets){
  for (p in listOfPCs){
    for (l in listOfLearningSetProportion){
      for (m in listOfMethods){
        makeClassificationModel(s, p, l, m)
      }
    }
  }
}

printConfusionMatrix <- function(dirName) {
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
        spl[[1]][3]
      )
      names(factorRes)[1]  <- 'Type Of Model'
      names(factorRes)[2]  <- 'Learning Method'
      names(factorRes)[3]  <- 'Number Of Components'
      names(factorRes)[4]  <- 'Training Proportion'      
      output <- rbind(output, res)
      outputFactor <- rbind(outputFactor, factorRes)
    }
    
  }
  rownames(output) <- NULL
  rownames(outputFactor) <- NULL
  return(data.frame(output, outputFactor))
}
cfnm <- printConfusionMatrix('/home/artem/work/2017/grant_R_code_and_data/interval_150_50/')
# cfnm07<-cfnm[which(cfnm[,14]==0.7),]


#for promoters_vs_lowscore
set.seed(999)
inTraining <- createDataPartition(promoters_vs_lowscore$factor_to_promoters_vs_lowscore, p = 0.7, list = F)
training <- promoters_vs_lowscore[inTraining,]
testing  <- promoters_vs_lowscore[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 15, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_lowscore <- train(factor_to_promoters_vs_lowscore ~ .,
                                   data = training,
                                   method = "nb",#nb, lda, nbDiscrete, nbSearch
                                   preProcess=c("center", "scale"),
                                   tuneLength = 15,
                                   trControl = fitControl,
                                   metric = "ROC"
)

predictionClasses_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing)
predictionProb_promoters_vs_lowscore <- predict(fit_promoters_vs_lowscore, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_lowscore <- confusionMatrix(data = predictionClasses_promoters_vs_lowscore, testing$factor_to_promoters_vs_lowscore)



#for promoters_vs_non_promoters
set.seed(999)
inTraining <- createDataPartition(promoters_vs_non_promoters$factor_to_promoters_vs_non_promoters, p = 0.7, list = F)
training <- promoters_vs_non_promoters[inTraining,]
testing  <- promoters_vs_non_promoters[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_non_promoters <- train(factor_to_promoters_vs_non_promoters ~ .,
                                        data = training,
                                        method = "nb",#nb, lda, nbDiscrete, nbSearch
                                        preProcess=c("center", "scale"),
                                        tuneLength = 7,
                                        trControl = fitControl,
                                        metric = "ROC"
)
predictionClasses_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing)
predictionProb_promoters_vs_non_promoters <- predict(fit_promoters_vs_non_promoters, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_non_promoters <- confusionMatrix(data = predictionClasses_promoters_vs_non_promoters, testing$factor_to_promoters_vs_non_promoters)

#for promoters_vs_islands
set.seed(999)
inTraining <- createDataPartition(promoters_vs_islands$factor_to_promoters_vs_islands, p = 0.7, list = F)
training <- promoters_vs_islands[inTraining,]
testing  <- promoters_vs_islands[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_islands <- train(factor_to_promoters_vs_islands ~ .,
                                  data = training,
                                  method = "nb",#nb, lda, nbDiscrete, nbSearch
                                  preProcess=c("center", "scale"),
                                  tuneLength = 7,
                                  trControl = fitControl,
                                  metric = "ROC"
)
predictionClasses_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing)
predictionProb_promoters_vs_islands <- predict(fit_promoters_vs_islands, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_islands <- confusionMatrix(data = predictionClasses_promoters_vs_islands, testing$factor_to_promoters_vs_islands)

#for promoters_vs_genes
set.seed(999)
inTraining <- createDataPartition(promoters_vs_genes$factor_to_promoters_vs_genes, p = 0.7, list = F)
training <- promoters_vs_genes[inTraining,]
testing  <- promoters_vs_genes[-inTraining,]
fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, 
                           repeats = 3, 
                           allowParallel = T, 
                           classProbs = T, 
                           summaryFunction = twoClassSummary
)

#train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
fit_promoters_vs_genes <- train(factor_to_promoters_vs_genes ~ .,
                                data = training,
                                method = "nb",#nb, lda, nbDiscrete, nbSearch
                                preProcess=c("center", "scale"),
                                tuneLength = 7,
                                trControl = fitControl,
                                metric = "ROC"
)
predictionClasses_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing)
predictionProb_promoters_vs_genes <- predict(fit_promoters_vs_genes, newdata = testing, type ="prob")
confusionMatrix_promoters_vs_genes <- confusionMatrix(data = predictionClasses_promoters_vs_genes, testing$factor_to_promoters_vs_genes)
