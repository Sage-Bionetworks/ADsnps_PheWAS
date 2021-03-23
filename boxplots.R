#boxplots for main text and supplementary material

#first run Upload_data_run_phewas.R to upload data, get numeric SNPs, log transform analytes, and create covariate files
#want two data frames: one for chemistries, one for proteins & metabolites (so i can adjust for vendor id for chemistries only)
#each data frame should have the desired covariates and all snps


#Read in files provided by ISB
demogs <- read.csv(file="demographics_file.csv")
chems <- read.csv(file="chemistries_file.csv")
SNPs <- read.csv(file="SNP_file.csv")
PRS <- read.csv(file="genetics_file.csv")
prots <- read.csv(file="proteins_file.csv")
mets <- read.csv(file="metabolites_file.csv")

#replace blank cells with "NA"
demogs[demogs==""]<-NA

names(demogs)[names(demogs) == 'X'] <- 'public_client_id'
vendor <- subset(chems, select=c(public_client_id, vendor))
chems$vendor<-NULL
demogs <- merge(demogs, vendor, all=TRUE)
#create numerical variable for sex (male = 0, female = 1)
demogs$female <- ifelse(demogs$sex=="M",0,1)
demogs$vendor_id <- ifelse(demogs$vendor=="Quest",0,1)
demogs$cholmeds <- ifelse(demogs$cholesterol.enum=="(3) Not at all",0,1)


#need genotypes, not numeric
SNPs <- subset(SNPs, select=c(public_client_id, rsid, genotype))

rs7412<- subset(SNPs, SNPs$rsid=='rs7412')
rs7412$rsid<-NULL
names(rs7412)[names(rs7412) == 'genotype'] <- 'rs7412'
rs429358 <- subset(SNPs, SNPs$rsid=='rs429358')
rs429358$rsid<-NULL
names(rs429358)[names(rs429358) == 'genotype'] <- 'rs429358'
rs10933431<- subset(SNPs, SNPs$rsid=='rs10933431')
rs10933431$rsid<-NULL
names(rs10933431)[names(rs10933431) == 'genotype'] <- 'rs10933431'
rs12539172<- subset(SNPs, SNPs$rsid=='rs12539172')
rs12539172$rsid<-NULL
names(rs12539172)[names(rs12539172) == 'genotype'] <- 'rs12539172'
rs3752246<- subset(SNPs, SNPs$rsid=='rs3752246')
rs3752246$rsid<-NULL
names(rs3752246)[names(rs3752246) == 'genotype'] <- 'rs3752246'
rs1859788<- subset(SNPs, SNPs$rsid=='rs1859788')
rs1859788$rsid<-NULL
names(rs1859788)[names(rs1859788) == 'genotype'] <- 'rs1859788'

#merge into one big file, use dataframes for SNPs and log-transformed chemistries & metabolites generated in Heath_Run_Phewas.R
alldata <- merge(demogs, logchems, all=TRUE)
alldata <- merge(alldata, logprotmets, all=TRUE)
alldata <- merge(alldata, rs7412, all=TRUE)
alldata <- merge(alldata, rs429358, all=TRUE)
alldata <- merge(alldata, rs10933431, all=TRUE)
alldata <- merge(alldata, rs12539172, all=TRUE)
alldata <- merge(alldata, rs3752246, all=TRUE)
alldata <- merge(alldata, rs1859788, all=TRUE)

#create 10-year age categories: 0=18-29, 1=30-39, 2=40-49, 3=50-59, 4=60-69, 5=70+
alldata$age_cat <- ifelse(alldata$age<30, 0,
                          ifelse(alldata$age>29 & alldata$age<40, 1,
                                 ifelse(alldata$age>39 & alldata$age<50, 2,
                                        ifelse(alldata$age>49 & alldata$age<60, 3,
                                               ifelse(alldata$age>69, 5, 4)))))
alldata$age_labels <- ifelse(alldata$age_cat==0, "18-29",
                             ifelse(alldata$age_cat==1, "30-39",
                                    ifelse(alldata$age_cat==2, "40-49",
                                           ifelse(alldata$age_cat==3, "50-59",
                                                  ifelse(alldata$age_cat==4, "60-69", "70+")))))

table(alldata$age_cat)
table(alldata$age_labels)




########NYAP1 plot (rs12539172)
pdf(file="NYAP1_PILRB.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, MET_Q9UKJ0, fill=rs12539172)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("PILRB")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs12539172 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()

pdf(file="NYAP1_PILRA.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, DEV_Q9UKJ1, fill=rs12539172)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("PILRA")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs12539172 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()


########PILRA plot (rs1859788)
alldata$rs1859788<-as.factor(alldata$rs1859788)
alldata$rs1859788 <- factor(alldata$rs1859788, levels=rev(levels(alldata$rs1859788)))
pdf(file="PILRAsnp_PILRAprot.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, DEV_Q9UKJ1, fill=rs1859788)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("PILRA")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("G/G", "A/G", "A/A"), values=c("white", "gray", "black"), name="rs1859788 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()

pdf(file="PILRAsnp_PILRBprot.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, MET_Q9UKJ0, fill=rs1859788)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("PILRB")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("G/G", "A/G", "A/A"), values=c("white", "gray", "black"), name="rs1859788 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()



#########APOE2 plots (rs7412)
pdf(file="APOE2_LDLpartnumber.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LDL.PARTICLE.NUMBER, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL particle number")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(1.8, 4.0)
p
dev.off()

pdf(file="APOE2_LDLcholesterol.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LDL.CHOL.CALCULATION, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(1.8, 4.0)
p
dev.off()

pdf(file="APOE2_Totalcholesterol.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CHOLESTEROL..TOTAL, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("Total cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(1.8, 4.0)
p
dev.off()

pdf(file="APOE2_LDLsmall.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LDL.SMALL, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL small")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(1.8, 4.0)
p
dev.off()

pdf(file="APOE2_LDLR.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CVD3_P01130, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDLR")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(1.8, 4.0)
p
dev.off()

pdf(file="APOE2_HMOX1.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CVD2_P09601, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("HMOX1")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) + ylim(9, 14)
p
dev.off()

pdf(file="APOE2_SLAMF8.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CRE_Q9P0V8, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("SLAMF8")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_RNF31.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, NEX_Q96EP0, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("RNF31")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="Manu_APOE2_CNTNAP2.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, IRE_Q9UHC6, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("CNTNAP2")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_SRP14.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, NEX_P37108, fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("SRP14")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_DGoag18_1_20_4_2.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, oleoyl.arachidonoyl.glycerol..18.1.20.4...2.., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG oag (18:1/20:4) [2]")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_DGoag18_1_20_4_1.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, oleoyl.arachidonoyl.glycerol..18.1.20.4...1.., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG oag (18:1/20:4) [1]")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_DGlag18_2_20_4.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, linoleoyl.arachidonoyl.glycerol..18.2.20.4...1.., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG lag (18:2/20:4)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_DGpag16_0_20_4.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, palmitoyl.arachidonoyl.glycerol..16.0.20.4...2.., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG pag (16:0/20:4)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()


pdf(file="APOE2_DGoog_18_1_18_1_1.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LipidDiacylglycerololeoyl.oleoyl.glycerol..18.1.18.1....1.., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG oog (18:1/18:1)[1]")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE2_MG_1_lg_18_2.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LipidMonoacylglycerol1.linoleoylglycerol..18.2., fill=rs7412)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("MG 1-lg (18:2)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black"), name="rs7412 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()




#########APOE4 plots (rs429358)
#I want to reverse the order of the boxes so minor allele is last in the plot
alldata$rs429358<-as.factor(alldata$rs429358)
alldata$rs429358 <- factor(alldata$rs429358, levels=rev(levels(alldata$rs429358)))

pdf(file="APOE4_LDLpartnumber.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LDL.PARTICLE.NUMBER, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL particle number")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_LDLcholcalc.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LDL.CHOL.CALCULATION, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_Totalcholesterol.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CHOLESTEROL..TOTAL, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("Total cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_trigly_hdl_ratio.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, Triglyceride.HDL.Ratio, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("Triglyceride:HDL ratio")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_HDLchol.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, HDL.CHOL.DIRECT, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("HDL cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_Triglycerides.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, TRIGLYCERIDES, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("Triglycerides")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_PLA2G7.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CAM_Q13093, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("PLA2G7")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()


pdf(file="APOE4_CD28.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, IRE_P10747, fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("CD28")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_DGoag18_1_20_4_1.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, oleoyl.arachidonoyl.glycerol..18.1.20.4...1.., fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG oag (18:1/20:4) [1]")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_DGoag18_1_20_4_2.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, oleoyl.arachidonoyl.glycerol..18.1.20.4...2.., fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG oag (18:1/20:4) [2]")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()

pdf(file="APOE4_DGpag16_1_20_4.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LipidDiacylglycerolpalmitoleoyl.arachidonoyl.glycerol..16.1.20.4...2.., fill=rs429358)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DG pag (16:1/20:4)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()


#########ABCA7 plots (rs3752246)

pdf(file="ABCA7_SL_nervonoyl.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LipidSphingolipid.Metabolismlactosyl.N.nervonoyl.sphingosine..d18.1.24.1.., fill=rs3752246)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("SL (d18:1/24:1)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black"), name="rs3752246 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()

pdf(file="ABCA7_SL_palmitoyl.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, LipidSphingolipid.Metabolismlactosyl.N.palmitoyl.sphingosine..d18.1.16.0., fill=rs3752246)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("SL (d18:1/16:1)")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black"), name="rs3752246 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) 
p
dev.off()


pdf(file="ABCA7_DEFA1.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CAM_P59665, fill=rs3752246)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("DEFA1")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black"), name="rs3752246 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) + ylim(1.8, 4.0)
p
dev.off()



#########INPP5D plot (rs10933431)
pdf(file="INPP5D_IDUA.pdf", width=8, height=7)
p <- ggplot(alldata, aes(x=age_labels, CVD2_P35475, fill=rs10933431)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("IDUA")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black"), name="rs10933431 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) + ylim(3.6, 8.2)
p
dev.off()



#### apoe4 and cholesterol medication use ########

cholmeds <- subset(alldata, !is.na(alldata$cholmeds))
dim(cholmeds)
cholmeds$cholmeds_label <- ifelse(cholmeds$cholmeds==0, "No Cholesterol Meds", "Cholesterol Meds")
#change the levels:
cholmeds$cholmeds_label<-as.factor(cholmeds$cholmeds_label)
cholmeds$cholmeds_label <- factor(cholmeds$cholmeds_label, levels=rev(levels(cholmeds$cholmeds_label)))
pdf(file="APOE4_cholmeds_LDL.pdf", width=10, height=7)
p <- ggplot(cholmeds, aes(x=age_labels, y=LDL.PARTICLE.NUMBER, fill=rs429358)) 
p <- p + geom_boxplot() + theme_linedraw(base_size=20) + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL Particle Number (log)") + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype") + facet_wrap(~cholmeds_label)  
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14))
p <- p + theme(axis.text.x=element_text(angle=90))
p
dev.off()



##### supplementary material plots

apoe4_facets <- subset(alldata, select=c(rs429358, age_labels, LDL.PARTICLE.NUMBER,LDL.CHOL.CALCULATION,CHOLESTEROL..TOTAL,HDL.CHOL.DIRECT,Triglyceride.HDL.Ratio,TRIGLYCERIDES,CAM_Q13093,IRE_P10747,oleoyl.arachidonoyl.glycerol..18.1.20.4...1..,LipidDiacylglycerolpalmitoleoyl.arachidonoyl.glycerol..16.1.20.4...2..,oleoyl.arachidonoyl.glycerol..18.1.20.4...2..))
head(apoe4_facets)

#need to reshape:
long <- reshape2::melt(apoe4_facets, id.vars = c("rs429358", "age_labels"))
head(long)
dim(long)


#change the olink names to common gene names:
long2<-long
long2$variable<-as.character(long2$variable)
long2$variable[long2$variable=="LDL.PARTICLE.NUMBER"]<-"LDL particle number"
long2$variable[long2$variable=="LDL.CHOL.CALCULATION"]<-"LDL cholesterol"
long2$variable[long2$variable=="CHOLESTEROL..TOTAL"]<-"Total cholesterol"
long2$variable[long2$variable=="HDL.CHOL.DIRECT"]<-"HDL cholesterol"
long2$variable[long2$variable=="Triglyceride.HDL.Ratio"]<-"Triglyceride:HDL ratio"
long2$variable[long2$variable=="TRIGLYCERIDES"]<-"Triglycerides"
long2$variable[long2$variable=="CAM_Q13093"]<-"PLA2G7"
long2$variable[long2$variable=="IRE_P10747"]<-"CD28"
long2$variable[long2$variable=="oleoyl.arachidonoyl.glycerol..18.1.20.4...1.."]<-"DG oag (18:1/20:4)[1]"
long2$variable[long2$variable=="LipidDiacylglycerolpalmitoleoyl.arachidonoyl.glycerol..16.1.20.4...2.."]<-"DG pag (16:1/20:4)"
long2$variable[long2$variable=="oleoyl.arachidonoyl.glycerol..18.1.20.4...2.."]<-"DG oag (18:1/20:4)[2]"

long2$rs429358<-as.factor(long2$rs429358)
long2$rs429358 <- factor(long2$rs429358, levels=rev(levels(long2$rs429358)))
long2$variable<-as.factor(long2$variable)
long2$variable <- factor(long2$variable, levels=c("LDL particle number","LDL cholesterol","Total cholesterol","HDL cholesterol","Triglyceride:HDL ratio","Triglycerides","PLA2G7","CD28","DG oag (18:1/20:4)[1]","DG pag (16:1/20:4)","DG oag (18:1/20:4)[2]"))

pdf(file="apoe4_facets.pdf", width=8, height=10)
p <- ggplot(long2, aes(x=age_labels, y=value, fill=rs429358)) 
p <- p + geom_boxplot(outlier.size=0.2) + theme_bw() + xlab("Age group (decade)") + ylab("log values") + guides(fill=guide_legend(title="rs429358")) 
p <- p + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black")) + facet_wrap(~variable, scales="free", ncol=3)
p <- p + theme(legend.position = c(0.9,0.1)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p
dev.off()




######## repeat for apoe2 supplementary data
apoe2_facets <- subset(alldata, select=c(rs7412, age_labels, LDL.PARTICLE.NUMBER,LDL.CHOL.CALCULATION,CHOLESTEROL..TOTAL,LDL.SMALL,oleoyl.arachidonoyl.glycerol..18.1.20.4...2..,oleoyl.arachidonoyl.glycerol..18.1.20.4...1..,linoleoyl.arachidonoyl.glycerol..18.2.20.4...1..,palmitoyl.arachidonoyl.glycerol..16.0.20.4...2..,LipidDiacylglycerololeoyl.oleoyl.glycerol..18.1.18.1....1..,LipidMonoacylglycerol1.linoleoylglycerol..18.2.))
head(apoe2_facets)

#need to reshape:
long <- reshape2::melt(apoe2_facets, id.vars = c("rs7412", "age_labels"))
head(long)
dim(long)

#change the olink names to common gene names:
long2<-long
long2$variable<-as.character(long2$variable)
long2$variable[long2$variable=="LDL.PARTICLE.NUMBER"]<-"LDL particle number"
long2$variable[long2$variable=="LDL.CHOL.CALCULATION"]<-"LDL cholesterol"
long2$variable[long2$variable=="CHOLESTEROL..TOTAL"]<-"Total cholesterol"
long2$variable[long2$variable=="LDL.SMALL"]<-"LDL small"
long2$variable[long2$variable=="oleoyl.arachidonoyl.glycerol..18.1.20.4...2.."]<-"DG oag (18:1/20:4)[2]"
long2$variable[long2$variable=="oleoyl.arachidonoyl.glycerol..18.1.20.4...1.."]<-"DG oag (18:1/20:4)[1]"
long2$variable[long2$variable=="linoleoyl.arachidonoyl.glycerol..18.2.20.4...1.."]<-"DG lag (18:2/20:4) [1]"
long2$variable[long2$variable=="palmitoyl.arachidonoyl.glycerol..16.0.20.4...2.."]<-"DG pag (16:0/20:4)[2]"
long2$variable[long2$variable=="LipidDiacylglycerololeoyl.oleoyl.glycerol..18.1.18.1....1.."]<-"DG oog (18:1/18:1)[1]"
long2$variable[long2$variable=="LipidMonoacylglycerol1.linoleoylglycerol..18.2."]<-"1-Monolinolein"

long2$variable<-as.factor(long2$variable)
long2$variable <- factor(long2$variable, levels=c("LDL particle number","LDL cholesterol","Total cholesterol","LDL small","DG oag (18:1/20:4)[2]","DG oag (18:1/20:4)[1]","DG lag (18:2/20:4) [1]","DG pag (16:0/20:4)[2]","DG oog (18:1/18:1)[1]","1-Monolinolein"))

pdf(file="apoe2_facets.pdf", width=8, height=10)
p <- ggplot(long2, aes(x=age_labels, y=value, fill=rs7412)) 
p <- p + geom_boxplot(outlier.size=0.2) + theme_bw() + xlab("Age group (decade)") + ylab("log values") + guides(fill=guide_legend(title="rs7412")) 
p <- p + scale_fill_manual(breaks = c("C/C", "C/T", "T/T"), values=c("white", "gray", "black")) + facet_wrap(~variable, scales="free", ncol=3)
p <- p + theme(legend.position = c(0.9,0.1)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p
dev.off()




######## repeat for abca7 & INPP5D supplementary data
abca7_facets <- subset(alldata, select=c(rs3752246, age_labels,LipidSphingolipid.Metabolismlactosyl.N.nervonoyl.sphingosine..d18.1.24.1..,LipidSphingolipid.Metabolismlactosyl.N.palmitoyl.sphingosine..d18.1.16.0.,CAM_P59665 ))
head(abca7_facets)

#need to reshape:
long <- reshape2::melt(abca7_facets, id.vars = c("rs3752246", "age_labels"))
head(long)
dim(long)

#change the olink names to common gene names:
long2<-long
long2$variable<-as.character(long2$variable)
long2$variable[long2$variable=="LipidSphingolipid.Metabolismlactosyl.N.nervonoyl.sphingosine..d18.1.24.1.."]<-"LC (d18:1/24:1)"
long2$variable[long2$variable=="LipidSphingolipid.Metabolismlactosyl.N.palmitoyl.sphingosine..d18.1.16.0."]<-"LC (d18:1/16:0)"
long2$variable[long2$variable=="CAM_P59665"]<-"DEFA1"

long2$variable<-as.factor(long2$variable)
long2$variable <- factor(long2$variable, levels=c("LC (d18:1/24:1)","LC (d18:1/16:0)","DEFA1"))

pdf(file="abca7_facets.pdf", width=9, height=3)
p <- ggplot(long2, aes(x=age_labels, y=value, fill=rs3752246)) 
p <- p + geom_boxplot(outlier.size=0.2) + theme_bw() + xlab("Age group (decade)") + ylab("log values") + guides(fill=guide_legend(title="rs3752246")) 
p <- p + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black")) + facet_wrap(~variable, scales="free", ncol=3)
p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p
dev.off()



#####inpp5d
inpp5d_facets <- subset(alldata, select=c(rs10933431, age_labels, CVD2_P35475))
head(inpp5d_facets)

long <- reshape2::melt(inpp5d_facets, id.vars = c("rs10933431", "age_labels"))
head(long)
dim(long)

#change the olink names to common gene names:
long2<-long
long2$variable<-as.character(long2$variable)
long2$variable[long2$variable=="CVD2_P35475"]<-"IDUA"

long2$variable<-as.factor(long2$variable)
long2$variable <- factor(long2$variable, levels="IDUA")

pdf(file="inpp5d_facets.pdf", width=4, height=3)
p <- ggplot(long2, aes(x=age_labels, y=value, fill=rs10933431)) 
p <- p + geom_boxplot(outlier.size=0.2) + theme_bw() + xlab("Age group (decade)") + ylab("log values") + guides(fill=guide_legend(title="rs10933431")) 
p <- p + scale_fill_manual(breaks = c("C/C", "C/G", "G/G"), values=c("white", "gray", "black")) + facet_wrap(~variable, scales="free", ncol=1)
p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p
dev.off()


