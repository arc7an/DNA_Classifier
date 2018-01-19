#rm(list=ls())
source("/home/artem/work/2018/classifier_on_other_genomes/calculate_sequence_phys_profile.R")
load("pre_pca_data_promoter_set_200_150_50.Rdata")

ecoli <- unlist(read.fasta('/home/artem/work/2016/iteb/e.coli_U00096.2.fasta', as.string = F, seqonly = T))
ecoli_char <- unlist(strsplit(ecoli, ''))


load("/home/artem/work/2018/classifier_on_other_genomes/annotated.clean.promoter.set.Rdata")
string_operations_test <- head(promClean[,1:9][promClean$Sigma =="Sigma70",])[2:3,]
forward_test <- paste(get_forward_substring(ecoli_char, string_operations_test[1,]$TSS, c(60, 20)), collapse="") == toupper(string_operations_test[1,]$seq)
reverse_test <- paste(get_reverse_substring(ecoli_char, string_operations_test[2,]$TSS, c(60, 20)), collapse="") == toupper(string_operations_test[2,]$seq)

#calculated_promoters <- which(rownames(pca_set) %in% c("ynfEp", "aceBp", "aroKp1", "ecpDp1"))
aceBp <- pca_set[which(rownames(pca_set) == "aceBp"),]
ynfEp <- pca_set[which(rownames(pca_set) == "ynfEp"),]
accAp <- pca_set[which(rownames(pca_set) == "accAp"),]

aroKp1 <- pca_set[which(rownames(pca_set) == "aroKp1"),]
ecpDp1 <- pca_set[which(rownames(pca_set) == "ecpDp1"),]
accDp <- pca_set[which(rownames(pca_set) == "accDp"),]

load('/home/artem/work/2016/iteb/spline_dataset_pro.Rdata')
#forward
tested_ynfEp <- dataset_pro$ynfEp
tested_aceBp <- dataset_pro$aceBp
tested_accAp <- dataset_pro$accAp
#reverse
tested_aroKp1 <- dataset_pro$aroKp1
tested_ecpDp1 <- dataset_pro$ecpDp1
tested_accDp <- dataset_pro$accDp

tested_ynfEp_profile <- calculate_profile(ecoli, 200, tested_ynfEp$tss, c(267, 217), -480:239, 1)
tested_aceBp_profile <- calculate_profile(ecoli, 200, tested_aceBp$tss, c(267, 217), -480:239, 1)
tested_accAp_profile <- calculate_profile(ecoli, 200, tested_accAp$tss, c(267, 217), -480:239, 1)

tested_aroKp1_profile <- calculate_profile(ecoli, 200, tested_aroKp1$tss, c(163, 213), -480:239, 2)
tested_ecpDp1_profile <- calculate_profile(ecoli, 200, tested_ecpDp1$tss, c(163, 213), -480:239, 2)
tested_accDp_profile <- calculate_profile(ecoli, 200, tested_accDp$tss, c(163, 213), -480:239, 2)

e0_test_aceBp <- all(aceBp[1:201] == tested_aceBp_profile[1:201])
print(e0_test_aceBp)
e0_test_ynfEp <- all(ynfEp[1:201] == tested_ynfEp_profile[1:201])
print(e0_test_ynfEp)
e0_test_accAp <- all(accAp[1:201] == tested_accAp_profile[1:201])
print(e0_test_accAp)

e0_test_aroKp1 <- all(aroKp1[1:201] == tested_aroKp1_profile[1:201])
print(e0_test_aroKp1)
e0_test_ecpDp1 <- all(ecpDp1[1:201] == tested_ecpDp1_profile[1:201])
print(e0_test_ecpDp1)
e0_test_accDp <- all(accDp[1:201] == tested_accDp_profile[1:201])
print(e0_test_accDp)

ep_test_aceBp <- aceBp[604:1323] - tested_aceBp_profile[403:(403+719)]
ep_test_ynfEp <- ynfEp[604:1323] - tested_ynfEp_profile[403:(403+719)]
ep_test_accAp <- accAp[604:1323] - tested_accAp_profile[403:(403+719)]

ep_test_aroKp1 <- aroKp1[604:1323] - tested_aroKp1_profile[403:(403+719)]
ep_test_ecpDp1 <- ecpDp1[604:1323] - tested_ecpDp1_profile[403:(403+719)]
ep_test_accDp <- accDp[604:1323] - tested_accDp_profile[403:(403+719)]