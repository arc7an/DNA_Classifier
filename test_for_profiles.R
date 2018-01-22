#rm(list=ls())
source("/home/artem/work/2018/classifier_on_other_genomes/calculate_sequence_phys_profile.R")
load("/home/artem/work/2018/classifier_on_other_genomes/pre_pca_data_promoter_set_200_150_50.Rdata")


ecoli <- unlist(read.fasta('/home/artem/work/2016/iteb/e.coli_U00096.2.fasta', as.string = F, seqonly = T))
ecoli_char <- unlist(strsplit(ecoli, ''))
test_txt <- readChar("/home/artem/work/2018/classifier_on_other_genomes/test_sequence_octave_vs_r.txt", 1000)

#octave E0 result for first 1000 bp of ecoli forward strand
octave_E0_forward <- c(229.90, 229.50, 229.24, 229.73, 229.90, 230.30, 230.47, 230.64, 230.87, 230.64, 230.81, 231.21, 231.04, 231.04, 230.78, 230.55, 230.55, 230.95, 230.46, 230.46)
r_E0_forward <- dynchars(test_txt, 200, 201, c(100, 699), "forward")
all((r_E0_forward$E0[1:20] - octave_E0_forward) < 0.01)
octave_E0_reverse <- c(235.61, 236.87, 237.11,  237.27, 237.51, 237.51, 236.87, 236.70, 236.06, 235.90, 235.25, 235.41, 235.66, 235.66, 235.42, 235.42, 235.67, 235.67, 235.42, 235.26)
r_E0_reverse <- dynchars(test_txt, 200, 801, c(99, 700), "reverse")
all((r_E0_reverse$E0[1:20] - octave_E0_reverse) < 0.01)

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

# tested_ynfEp_profile <- calculate_profile(ecoli, 200, tested_ynfEp$tss, c(150, 50), c(267, 217), -480:239, "forward")
# tested_aceBp_profile <- calculate_profile(ecoli, 200, tested_aceBp$tss, c(150, 50), c(267, 217), -480:239, "forward")
# tested_accAp_profile <- calculate_profile(ecoli, 200, tested_accAp$tss, c(150, 50), c(267, 217), -480:239, "forward")
# 
# tested_aroKp1_profile <- calculate_profile(ecoli, 200, tested_aroKp1$tss, c(150, 50), c(163, 213), -480:239, "reverse")
# tested_ecpDp1_profile <- calculate_profile(ecoli, 200, tested_ecpDp1$tss, c(150, 50), c(163, 213), -480:239, "reverse")
# tested_accDp_profile <- calculate_profile(ecoli, 200, tested_accDp$tss, c(150, 50), c(163, 213), -480:239, "reverse")

tested_ynfEp_profile <- calculate_profile(ecoli, 200, tested_ynfEp$tss, c(150, 50), c(250, 150), -480:239, "forward")
tested_aceBp_profile <- calculate_profile(ecoli, 200, tested_aceBp$tss, c(150, 50), c(250, 150), -480:239, "forward")
tested_accAp_profile <- calculate_profile(ecoli, 200, tested_accAp$tss, c(150, 50), c(250, 150), -480:239, "forward")

tested_aroKp1_profile <- calculate_profile(ecoli, 200, tested_aroKp1$tss, c(150, 50), c(250, 150), -480:239, "reverse")
tested_ecpDp1_profile <- calculate_profile(ecoli, 200, tested_ecpDp1$tss, c(150, 50), c(250, 150), -480:239, "reverse")
tested_accDp_profile <- calculate_profile(ecoli, 200, tested_accDp$tss, c(150, 50), c(250, 150), -480:239, "reverse")

#probably these test do not make sense
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
