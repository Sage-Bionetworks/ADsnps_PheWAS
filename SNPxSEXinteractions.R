library(plyr)
library(dplyr)

#first run Upload_data_run_phewas.R to upload data, get numeric SNPs, log transform analytes, and create covariate files
#want two data frames: one for chemistries, one for proteins & metabolites (so i can adjust for vendor id for chemistries only)
#each data frame should have the desired covariates and all snps

#create a single file with all the SNPs:
allSNPs <- merge(rs4844610, rs6733839)
allSNPs <- merge(allSNPs, rs10933431)
allSNPs <- merge(allSNPs, rs9271058)
allSNPs <- merge(allSNPs, rs9473117)
allSNPs <- merge(allSNPs, rs12539172)
allSNPs <- merge(allSNPs, rs10808026)
allSNPs <- merge(allSNPs, rs73223431)
allSNPs <- merge(allSNPs, rs9331896)
allSNPs <- merge(allSNPs, rs7920721)
allSNPs <- merge(allSNPs, rs3740688)
allSNPs <- merge(allSNPs, rs7933202)
allSNPs <- merge(allSNPs, rs3851179)
allSNPs <- merge(allSNPs, rs11218343)
allSNPs <- merge(allSNPs, rs17125924)
allSNPs <- merge(allSNPs, rs12881735)
allSNPs <- merge(allSNPs, rs593742)
allSNPs <- merge(allSNPs, rs7185636)
allSNPs <- merge(allSNPs, rs62039712)
allSNPs <- merge(allSNPs, rs138190086)
allSNPs <- merge(allSNPs, rs3752246)
allSNPs <- merge(allSNPs, rs7412)
allSNPs <- merge(allSNPs, rs429358)
allSNPs <- merge(allSNPs, rs6024870)
allSNPs <- merge(allSNPs, rs2830500)
allSNPs <- merge(allSNPs, rs1859788)

#create two data frames: one with chemistries, covariates, and SNPs, one with proteins/mets, covariates, and SNPs
alldatachems <- merge(covariates_chems, logchems, all=TRUE)
alldatachems <- merge(alldatachems, allSNPs, all=TRUE)

alldataprotmets <- merge(covariates_protmets, logprotmets, all=TRUE)
alldataprotmets <- merge(alldataprotmets, allSNPs, all=TRUE)

#a few chemistries had non-variable or missing vendor data (wthey ere all from one vendor), which means the  
#coefficient gets removed from the regression and upsets the ordering of the variables later. Need to 
#take these chemistries out of the chems file and add them to the protmets file, since they do not need vendor adjustment:

movedata <- subset(alldatachems, select=c(public_client_id, DPA,GFR..MDRD..AFRICAN.AM,HDL.PARTICLE.NUMBER,LDL_SIZE,LINOLEIC_ACID,OMEGA_3_TOTAL,OMEGA_6_TOTAL))
alldataprotmets <- merge(alldataprotmets, movedata, all=TRUE)
alldatachems <- subset(alldatachems, select = -c(DPA,GFR..MDRD..AFRICAN.AM,HDL.PARTICLE.NUMBER,LDL_SIZE,LINOLEIC_ACID,OMEGA_3_TOTAL,OMEGA_6_TOTAL))


#SNP by SNP linear regression with genotype x sex as an interaction term, plus covariates (female=1, male=0)

####### rs4844610 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs4844610*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs4844610*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs4844610",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs4844610:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs4844610",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs4844610:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs4844610:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs4844610.csv", sep=",", row.names=FALSE)


####### rs6733839 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs6733839*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs6733839*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs6733839",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs6733839:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs6733839",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs6733839:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs6733839:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs6733839.csv", sep=",", row.names=FALSE)


####### rs10933431 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs10933431*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs10933431*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs10933431",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs10933431:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs10933431",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs10933431:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs10933431:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs10933431.csv", sep=",", row.names=FALSE)


####### rs9473117 ##########
#Run linear regression on each phenotype:

fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs9473117*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs9473117*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs9473117",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs9473117:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs9473117",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs9473117:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs9473117:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs9473117.csv", sep=",", row.names=FALSE)


####### rs9271058 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs9271058*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs9271058*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs9271058",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs9271058:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs9271058",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs9271058:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs9271058:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs9271058.csv", sep=",", row.names=FALSE)


####### rs10808026 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs10808026*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs10808026*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs10808026",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs10808026:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs10808026",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs10808026:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs10808026:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs10808026.csv", sep=",", row.names=FALSE)

####### rs12539172 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs12539172*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs12539172*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs12539172",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs12539172:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs12539172",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs12539172:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs12539172:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs12539172.csv", sep=",", row.names=FALSE)

####### rs9331896 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs9331896*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs9331896*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs9331896",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs9331896:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs9331896",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs9331896:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs9331896:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs9331896.csv", sep=",", row.names=FALSE)


####### rs73223431 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs73223431*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs73223431*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs73223431",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs73223431:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs73223431",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs73223431:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs73223431:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs73223431.csv", sep=",", row.names=FALSE)



####### rs7920721 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs7920721*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs7920721*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs7920721",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs7920721:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs7920721",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs7920721:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs7920721:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs7920721.csv", sep=",", row.names=FALSE)


####### rs7933202 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs7933202*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs7933202*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs7933202",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs7933202:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs7933202",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs7933202:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs7933202:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs7933202.csv", sep=",", row.names=FALSE)


####### rs3851179 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs3851179*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs3851179*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs3851179",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs3851179:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs3851179",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs3851179:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs3851179:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs3851179.csv", sep=",", row.names=FALSE)


####### rs11218343 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs11218343*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs11218343*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs11218343",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs11218343:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs11218343",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs11218343:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs11218343:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs11218343.csv", sep=",", row.names=FALSE)

####### rs3740688 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs3740688*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs3740688*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs3740688",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs3740688:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs3740688",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs3740688:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs3740688:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs3740688.csv", sep=",", row.names=FALSE)

####### rs12881735 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs12881735*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs12881735*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs12881735",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs12881735:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs12881735",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs12881735:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs12881735:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs12881735.csv", sep=",", row.names=FALSE)


####### rs17125924 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs17125924*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs17125924*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs17125924",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs17125924:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs17125924",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs17125924:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs17125924:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs17125924.csv", sep=",", row.names=FALSE)


####### rs593742 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs593742*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs593742*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs593742",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs593742:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs593742",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs593742:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs593742:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs593742.csv", sep=",", row.names=FALSE)


####### rs7185636 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs7185636*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs7185636*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs7185636",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs7185636:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs7185636",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs7185636:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs7185636:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs7185636.csv", sep=",", row.names=FALSE)


####### rs62039712 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs62039712*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs62039712*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs62039712",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs62039712:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs62039712",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs62039712:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs62039712:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs62039712.csv", sep=",", row.names=FALSE)


####### rs138190086 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs138190086*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs138190086*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs138190086",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs138190086:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs138190086",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs138190086:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs138190086:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs138190086.csv", sep=",", row.names=FALSE)


####### rs7412 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs7412*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs7412*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs7412",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs7412:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs7412",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs7412:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs7412:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs7412.csv", sep=",", row.names=FALSE)


####### rs429358 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs429358*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs429358*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs429358",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs429358:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs429358",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs429358:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs429358:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs429358.csv", sep=",", row.names=FALSE)

####### rs3752246 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs3752246*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs3752246*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs3752246",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs3752246:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs3752246",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs3752246:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs3752246:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs3752246.csv", sep=",", row.names=FALSE)


####### rs6024870 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs6024870*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs6024870*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs6024870",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs6024870:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs6024870",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs6024870:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs6024870:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs6024870.csv", sep=",", row.names=FALSE)


####### rs2830500 ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~rs2830500*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~rs2830500*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "rs2830500",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "rs2830500:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "rs2830500",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "rs2830500:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="rs2830500:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_rs2830500.csv", sep=",", row.names=FALSE)

####### PRS_apoe ##########
alldatachems <- merge(alldatachems, PRS_apoe, all=TRUE)
alldatachems <- merge(alldatachems, PRS_no_apoe, all=TRUE)
alldataprotmets <- merge(alldataprotmets, PRS_apoe, all=TRUE)
alldataprotmets <- merge(alldataprotmets, PRS_no_apoe, all=TRUE)
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~polygenic_alzheimers_APOE_scaled*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~polygenic_alzheimers_APOE_scaled*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "PRS_apoe",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "PRS_apoe:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "PRS_apoe",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "PRS_apoe:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="PRS_apoe:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_PRS_apoe.csv", sep=",", row.names=FALSE)




####### PRS_no_apoe ##########
#Run linear regression on each phenotype:
fits1 <- lapply(alldatachems[,9:70], function(x) lm(x~polygenic_alzheimers_scaled*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=alldatachems))
fits2 <- lapply(alldataprotmets[,8:1953], function(x) lm(x~polygenic_alzheimers_scaled*female + age + PC1 + PC2 + PC3 + PC4, data=alldataprotmets))   
summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm1)[1] <- "analyte"
lm2 <- ldply(summaries2, data.frame)
#need to label the rows: first number rows within each group (phenotype), then add labels
colnames(lm2)[1] <- "analyte"
lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                   ifelse(lm1$rownum==2, "PRS_no_apoe",
                          ifelse(lm1$rownum==3, "female",
                                 ifelse(lm1$rownum==4, "age", 
                                        ifelse(lm1$rownum==5, "vendor_id", 
                                               ifelse(lm1$rownum==6, "PC1", 
                                                      ifelse(lm1$rownum==7, "PC2", 
                                                             ifelse(lm1$rownum==8, "PC3",
                                                                    ifelse(lm1$rownum==9, "PC4", "PRS_no_apoe:female")))))))))


lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                   ifelse(lm2$rownum==2, "PRS_no_apoe",
                          ifelse(lm2$rownum==3, "female",
                                 ifelse(lm2$rownum==4, "age", 
                                        ifelse(lm2$rownum==5, "PC1", 
                                               ifelse(lm2$rownum==6, "PC2", 
                                                      ifelse(lm2$rownum==7, "PC3",
                                                             ifelse(lm2$rownum==8, "PC4", "PRS_no_apoe:female"))))))))

#calculate fdr adjusted p-values for the three-way interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
fdr_int <- subset(lm_all, lm_all$coef=="PRS_no_apoe:female")
adjustp <- p.adjust(fdr_int$Pr...t.., "fdr")
fdr2 <- cbind(adjustp, fdr_int)
fdr2$Estimate<-NULL
fdr2$Pr...t..<-NULL
fdr2$rownum<-NULL
lm_adjusted <- merge(lm_all, fdr2, all=TRUE)
write.table(lm_adjusted, file="interaction_sex_PRS_no_apoe.csv", sep=",", row.names=FALSE)
