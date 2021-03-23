install.packages("meta")
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)

#Read in files provided by ISB
demogs <- read.csv(file="demographics_file.csv")
chems <- read.csv(file="chemistries_file.csv")
SNPs <- read.csv(file="SNP_file.csv")
PRS <- read.csv(file="genetics_file.csv")
prots <- read.csv(file="proteins_file.csv")
mets <- read.csv(file="metabolites_file.csv")


#create a file with demographics and polygenic risk scores:
names(demogs)[names(demogs) == 'X'] <- 'public_client_id'
allcovariates <- merge(demogs, PRS)
#chemistries vendor was in the chemistries file, put it in the demographics file
vendor <- subset(chems, select=c(public_client_id, vendor))
chems$vendor<-NULL
allcovariates <- merge(allcovariates, vendor, all=TRUE)
#create numerical variable for sex (male = 0, female = 1)
allcovariates$female <- ifelse(allcovariates$sex=="M",0,1)
allcovariates$vendor_id <- ifelse(allcovariates$vendor=="Quest",0,1)

#create individual snp files, one per snp:
SNPs <- subset(SNPs, select=c(public_client_id, rsid, genotype))
#Kunkle 2019 snps:
rs429358 <- subset(SNPs, SNPs$rsid=='rs429358')
rs429358$rsid<-NULL
names(rs429358)[names(rs429358) == 'genotype'] <- 'rs429358'
rs7412<- subset(SNPs, SNPs$rsid=='rs7412')
rs7412$rsid<-NULL
names(rs7412)[names(rs7412) == 'genotype'] <- 'rs7412'
rs4844610<- subset(SNPs, SNPs$rsid=='rs4844610')
rs4844610$rsid<-NULL
names(rs4844610)[names(rs4844610) == 'genotype'] <- 'rs4844610'
rs6733839<- subset(SNPs, SNPs$rsid=='rs6733839')
rs6733839$rsid<-NULL
names(rs6733839)[names(rs6733839) == 'genotype'] <- 'rs6733839'
rs10933431<- subset(SNPs, SNPs$rsid=='rs10933431')
rs10933431$rsid<-NULL
names(rs10933431)[names(rs10933431) == 'genotype'] <- 'rs10933431'
rs9473117<- subset(SNPs, SNPs$rsid=='rs9473117')
rs9473117$rsid<-NULL
names(rs9473117)[names(rs9473117) == 'genotype'] <- 'rs9473117'
rs9271058<- subset(SNPs, SNPs$rsid=='rs9271058')
rs9271058$rsid<-NULL
names(rs9271058)[names(rs9271058) == 'genotype'] <- 'rs9271058'
rs10808026<- subset(SNPs, SNPs$rsid=='rs10808026')
rs10808026$rsid<-NULL
names(rs10808026)[names(rs10808026) == 'genotype'] <- 'rs10808026'
rs12539172<- subset(SNPs, SNPs$rsid=='rs12539172')
rs12539172$rsid<-NULL
names(rs12539172)[names(rs12539172) == 'genotype'] <- 'rs12539172'
rs9331896<- subset(SNPs, SNPs$rsid=='rs9331896')
rs9331896$rsid<-NULL
names(rs9331896)[names(rs9331896) == 'genotype'] <- 'rs9331896'
rs73223431<- subset(SNPs, SNPs$rsid=='rs73223431')
rs73223431$rsid<-NULL
names(rs73223431)[names(rs73223431) == 'genotype'] <- 'rs73223431'
rs7920721<- subset(SNPs, SNPs$rsid=='rs7920721')
rs7920721$rsid<-NULL
names(rs7920721)[names(rs7920721) == 'genotype'] <- 'rs7920721'
rs7933202<- subset(SNPs, SNPs$rsid=='rs7933202')
rs7933202$rsid<-NULL
names(rs7933202)[names(rs7933202) == 'genotype'] <- 'rs7933202'
rs3851179<- subset(SNPs, SNPs$rsid=='rs3851179')
rs3851179$rsid<-NULL
names(rs3851179)[names(rs3851179) == 'genotype'] <- 'rs3851179'
rs11218343<- subset(SNPs, SNPs$rsid=='rs11218343')
rs11218343$rsid<-NULL
names(rs11218343)[names(rs11218343) == 'genotype'] <- 'rs11218343'
rs3740688<- subset(SNPs, SNPs$rsid=='rs3740688')
rs3740688$rsid<-NULL
names(rs3740688)[names(rs3740688) == 'genotype'] <- 'rs3740688'
rs12881735<- subset(SNPs, SNPs$rsid=='rs12881735')
rs12881735$rsid<-NULL
names(rs12881735)[names(rs12881735) == 'genotype'] <- 'rs12881735'
rs17125924<- subset(SNPs, SNPs$rsid=='rs17125924')
rs17125924$rsid<-NULL
names(rs17125924)[names(rs17125924) == 'genotype'] <- 'rs17125924'
rs593742<- subset(SNPs, SNPs$rsid=='rs593742')
rs593742$rsid<-NULL
names(rs593742)[names(rs593742) == 'genotype'] <- 'rs593742'
rs7185636<- subset(SNPs, SNPs$rsid=='rs7185636')
rs7185636$rsid<-NULL
names(rs7185636)[names(rs7185636) == 'genotype'] <- 'rs7185636'
rs62039712<- subset(SNPs, SNPs$rsid=='rs62039712')
rs62039712$rsid<-NULL
names(rs62039712)[names(rs62039712) == 'genotype'] <- 'rs62039712'
rs138190086<- subset(SNPs, SNPs$rsid=='rs138190086')
rs138190086$rsid<-NULL
names(rs138190086)[names(rs138190086) == 'genotype'] <- 'rs138190086'
rs3752246<- subset(SNPs, SNPs$rsid=='rs3752246')
rs3752246$rsid<-NULL
names(rs3752246)[names(rs3752246) == 'genotype'] <- 'rs3752246'
rs6024870<- subset(SNPs, SNPs$rsid=='rs6024870')
rs6024870$rsid<-NULL
names(rs6024870)[names(rs6024870) == 'genotype'] <- 'rs6024870'
rs2830500 <- subset(SNPs, SNPs$rsid=='rs2830500')
rs2830500 $rsid<-NULL
names(rs2830500 )[names(rs2830500 ) == 'genotype'] <- 'rs2830500'

rs1859788<- subset(SNPs, SNPs$rsid=='rs1859788')
rs1859788$rsid<-NULL
names(rs1859788)[names(rs1859788) == 'genotype'] <- 'rs1859788'


#Recode SNPs to be numerical (homozygous major allele = 0, heterozygous = 1, homozygous minor allele = 2, see supplementary table 1 for major/minor allele designation)
rs429358$rs429358 <- (ifelse(rs429358$rs429358=="T/T", 0,
                             ifelse(rs429358$rs429358=="C/T", 1, 2)))
rs7412$rs7412 <- (ifelse(rs7412$rs7412=="C/C", 0,
                         ifelse(rs7412$rs7412=="C/T", 1, 2)))
rs4844610$rs4844610 <- (ifelse(rs4844610$rs4844610=="C/C", 0,
                               ifelse(rs4844610$rs4844610=="A/C", 1, 2)))
rs6733839$rs6733839 <- (ifelse(rs6733839$rs6733839=="C/C", 0,
                               ifelse(rs6733839$rs6733839=="C/T", 1, 2)))
rs10933431$rs10933431 <- (ifelse(rs10933431$rs10933431=="C/C", 0,
                                 ifelse(rs10933431$rs10933431=="C/G", 1, 2)))
rs9473117$rs9473117 <- (ifelse(rs9473117$rs9473117=="A/A", 0,
                               ifelse(rs9473117$rs9473117=="A/C", 1, 2)))
rs9271058$rs9271058 <- (ifelse(rs9271058$rs9271058=="T/T", 0,
                               ifelse(rs9271058$rs9271058=="A/T", 1, 2)))
rs10808026$rs10808026 <- (ifelse(rs10808026$rs10808026=="C/C", 0,
                                 ifelse(rs10808026$rs10808026=="A/C", 1, 2)))
rs12539172$rs12539172 <- (ifelse(rs12539172$rs12539172=="C/C", 0,
                                 ifelse(rs12539172$rs12539172=="C/T", 1, 2)))
rs9331896$rs9331896 <- (ifelse(rs9331896$rs9331896=="T/T", 0,
                               ifelse(rs9331896$rs9331896=="C/T", 1, 2)))
rs73223431$rs73223431 <- (ifelse(rs73223431$rs73223431=="C/C", 0,
                                 ifelse(rs73223431$rs73223431=="C/T", 1, 2)))
rs7920721$rs7920721 <- (ifelse(rs7920721$rs7920721=="A/A", 0,
                               ifelse(rs7920721$rs7920721=="A/G", 1, 2)))
rs7933202$rs7933202 <- (ifelse(rs7933202$rs7933202=="A/A", 0,
                               ifelse(rs7933202$rs7933202=="A/C", 1, 2)))
rs3851179$rs3851179 <- (ifelse(rs3851179$rs3851179=="C/C", 0,
                               ifelse(rs3851179$rs3851179=="C/T", 1, 2)))
rs11218343$rs11218343 <- (ifelse(rs11218343$rs11218343=="T/T", 0,
                                 ifelse(rs11218343$rs11218343=="C/T", 1, 2)))
rs3740688$rs3740688 <- (ifelse(rs3740688$rs3740688=="T/T", 0,
                               ifelse(rs3740688$rs3740688=="G/T", 1, 2)))
rs12881735$rs12881735 <- (ifelse(rs12881735$rs12881735=="T/T", 0,
                                 ifelse(rs12881735$rs12881735=="C/T", 1, 2)))
rs17125924$rs17125924 <- (ifelse(rs17125924$rs17125924=="A/A", 0,
                                 ifelse(rs17125924$rs17125924=="A/G", 1, 2)))
rs593742$rs593742 <- (ifelse(rs593742$rs593742=="A/A", 0,
                             ifelse(rs593742$rs593742=="A/G", 1, 2)))
rs7185636$rs7185636 <- (ifelse(rs7185636$rs7185636=="T/T", 0,
                               ifelse(rs7185636$rs7185636=="C/T", 1, 2)))
rs62039712$rs62039712 <- (ifelse(rs62039712$rs62039712=="G/G", 0,
                                 ifelse(rs62039712$rs62039712=="A/G", 1, 2)))
rs138190086$rs138190086 <- (ifelse(rs138190086$rs138190086=="G/G", 0,
                                   ifelse(rs138190086$rs138190086=="A/G", 1, 2)))
rs3752246$rs3752246 <- (ifelse(rs3752246$rs3752246=="C/C", 0,
                               ifelse(rs3752246$rs3752246=="C/G", 1, 2)))
rs6024870$rs6024870 <- (ifelse(rs6024870$rs6024870=="G/G", 0,
                               ifelse(rs6024870$rs6024870=="A/G", 1, 2)))
rs2830500$rs2830500 <- (ifelse(rs2830500$rs2830500=="C/C", 0,
                               ifelse(rs2830500$rs2830500=="A/C", 1, 2)))
rs1859788$rs1859788 <- (ifelse(rs1859788$rs1859788=="G/G", 0,
                               ifelse(rs1859788$rs1859788=="A/G", 1, 2)))


#create individual files for the polygenic risk scores (with and without APOE, used the scaled scores)
PRS_apoe <- subset(PRS, select=c(public_client_id, polygenic_alzheimers_APOE_scaled))
PRS_no_apoe <- subset(PRS, select=c(public_client_id, polygenic_alzheimers_scaled))

#create covariate file for chemistries (adjusted for age, sex, vendor, and four PCs)
covariates_chems <- subset(allcovariates, select=c(public_client_id, female, age, vendor_id, PC1, PC2, PC3, PC4))
#create covariate file for proteins & metabolites (adjusted for age, sex, and four PCs)
covariates_protmets <- subset(allcovariates, select=c(public_client_id, female, age, PC1, PC2, PC3, PC4))

#to deal with chemistries needing to be adjusted for vendor (but not proteins or metabolites) we'll separate the phenotypes into chems, prots, and
#mets, run phewas separately on chemistries vs. proteins & mets, then join them back together for the phewas plot (will also need to 
#recalculate FDR adjusted p-values)

#first need to log transform chemistries & metabolites
public_client_id <- chems[,1]
chems2 <- subset(chems, select=c(2:63))
logchems <- log(chems2+1)
logchems <- cbind(public_client_id, logchems)

public_client_id <- mets[,1]
mets2<-mets
mets2$public_client_id<-NULL
logmets <- log(mets2+1)
logmets <- cbind(public_client_id, logmets)

#proteins already log transformed
logprotmets <- merge(prots, logmets, all=TRUE)



###################### PHEWAS ####################################3

#run phewas, first for chemistries then proteins & mets, then join and FDR adjust on a SNP-by-SNP basis. Output full phewas results as csv file

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs4844610, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs4844610, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs4844610_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs6733839, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs6733839, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs6733839_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs10933431, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs10933431, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs10933431_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs9271058, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs9271058, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs9271058_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs9473117, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs9473117, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs9473117_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs12539172, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs12539172, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs12539172_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs10808026, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs10808026, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs10808026_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs73223431, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs73223431, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs73223431_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs9331896, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs9331896, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs9331896_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs7920721, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs7920721, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs7920721_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs3740688, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs3740688, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs3740688_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs7933202, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs7933202, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs7933202_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs3851179, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs3851179, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs3851179_phewas.csv", sep=",", row.names=FALSE)

results<-phewas(phenotypes=logchems, genotypes=rs11218343, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs11218343, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs11218343_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs17125924, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs17125924, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs17125924_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs12881735, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs12881735, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs12881735_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs593742, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs593742, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs593742_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs7185636, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs7185636, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs7185636_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs62039712, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs62039712, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs62039712_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs138190086, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs138190086, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs138190086_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs3752246, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs3752246, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs3752246_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs7412, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs7412, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs7412_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs429358, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs429358, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs429358_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs6024870, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs6024870, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs6024870_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs2830500, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs2830500, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs2830500_phewas.csv", sep=",", row.names=FALSE)

####SNP######
results<-phewas(phenotypes=logchems, genotypes=rs1859788, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=rs1859788, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_rs1859788_phewas.csv", sep=",", row.names=FALSE)




#### POLYGENIC RISK SCORE: WITH APOE ######
results<-phewas(phenotypes=logchems, genotypes=PRS_apoe, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=PRS_apoe, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_PRS_apoe_phewas.csv", sep=",", row.names=FALSE)


#### POLYGENIC RISK SCORE: WITHOUT APOE ######
results<-phewas(phenotypes=logchems, genotypes=PRS_no_apoe, covariates=covariates_chems, cores=1, significance.threshold=c("fdr"))
resultsB<-phewas(phenotypes=logprotmets, genotypes=PRS_no_apoe, covariates=covariates_protmets, cores=1, significance.threshold=c("fdr"))
results_all <- rbind(results, resultsB)
adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)
write.table(adjustresults, file="adjusted_PRS_no_apoe_phewas.csv", sep=",", row.names=FALSE)
