# 23-Jan 2025
setwd("/Users/arunkumarramesh/Downloads/")
library(ggplot2)
library(dplyr)
library(pophelper)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(grid)
library(plyr)
library(ggrepel)
library(agricolae)

####  variant site prop (MAF > 0.02) #### 
ceu_cl_smp_prop <- 310295/387615
ibnr_cl_smp_prop <- 623882/769044
ceu_gbm_smp_prop <- 209117/347705
ibnr_gbm_smp_prop <- 418482/688439
ceu_nongbm_smp_prop <- 279348/855999
ibnr_nongbm_smp_prop <- 565856/1881413
ceu_smp_prop <- 545912/1629415
ibnr_smp_prop <- 1135849/3589695
cg_sites <- read.table(file="cg_sites.bed")
387615/sum(cg_sites$V3-cg_sites$V2) ## to get SMP genotype density
769044/sum(cg_sites$V3-cg_sites$V2) ## to get SMP genotype density
ceu_snp_gene_prop <- 264383/17022694
ibnr_snp_gene_prop <- 289088/17022694
ceu_snp_integenic_prop <- 393091/14585335
ibnr_snp_integenic_prop <- 412233/14585335
propsites <- as.data.frame(rbind(ceu_smp_prop,ceu_cl_smp_prop,ceu_gbm_smp_prop,ceu_nongbm_smp_prop,ceu_snp_gene_prop,ceu_snp_integenic_prop,ibnr_smp_prop,ibnr_cl_smp_prop,ibnr_gbm_smp_prop,ibnr_nongbm_smp_prop,ibnr_snp_gene_prop,ibnr_snp_integenic_prop))
colnames(propsites) <- "Proportion"
propsites$Population <- c("CEU","CEU","CEU","CEU","CEU","CEU","IBnr","IBnr","IBnr","IBnr","IBnr","IBnr")
propsites$Type <- c("All methylation sites","Clocklike methylation sites","gbM sites","non-gbM sites","Coding nucleotides","Non-coding nucleotides","All methylation sites","Clocklike methylation sites","gbM sites","non-gbM sites","Coding nucleotides","Non-coding nucleotides")
total <- c("n=17022694","n=14585335","n=1629415","n=387615","n=347705","n=855999","n=17022694","n=14585335","n=3589695","n=769044","n=688439","n=1881413")
propsites$Type <- factor(propsites$Type,levels=c("Coding nucleotides","Non-coding nucleotides","All methylation sites","Clocklike methylation sites","gbM sites","non-gbM sites"))

seg_athaliana <- ggplot(data=propsites,aes(x=Population,y=Proportion,fill=Type)) +
  geom_col(position = "dodge") +
  ylim(c(0,1.05)) +
  scale_fill_manual(values=c("black","grey","brown","darkgreen","lightblue","darkblue")) +
  geom_text(aes(x=c(seq(0.6,1.4,(1.4-0.6)/(6-1)),seq(1.6,2.4,(2.4-1.6)/(6-1))),angle=90,y = 0.95, label = total)) +
  ylab("Proportion of polymorphic sites")

#### pi ####

gbm_genes_poswithnames <- read.table(file="gbm_genes_poswithnames.txt")
gbm_genes_poswithnames$interval <- paste(gbm_genes_poswithnames$V1,":",gbm_genes_poswithnames$V2,"-",gbm_genes_poswithnames$V3,sep="")
non_gbm_genes_poswithnames <- read.table(file="non_gbm_genes_poswithnames.txt")
non_gbm_genes_poswithnames$interval <- paste(non_gbm_genes_poswithnames$V1,":",non_gbm_genes_poswithnames$V2,"-",non_gbm_genes_poswithnames$V3,sep="")
high_sift_genes_poswithnames <- read.table(file="high_sift_genes_poswithnames.txt")
high_sift_genes_poswithnames$interval <- paste(high_sift_genes_poswithnames$V1,":",high_sift_genes_poswithnames$V2,"-",high_sift_genes_poswithnames$V3,sep="")
low_sift_genes_poswithnames <- read.table(file="low_sift_genes_poswithnames.txt")
low_sift_genes_poswithnames$interval <- paste(low_sift_genes_poswithnames$V1,":",low_sift_genes_poswithnames$V2,"-",low_sift_genes_poswithnames$V3,sep="")

ceu_snps_genes <- read.table(file="pi_D_ceu.txt",header=T)
ceu_snps_genes$Pop <- "CEU"
ceu_snps_genes$Type <- "All SNP"
ceu_snps_genes$interval <- gsub(".fa","",ceu_snps_genes$interval)
ceu_snps_genes$pi <- ceu_snps_genes$pi/(as.numeric(gsub(".*-","",gsub(".:","",ceu_snps_genes$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ceu_snps_genes$interval))))
ceu_snps_genes_gbm <- ceu_snps_genes[ceu_snps_genes$interval %in% gbm_genes_poswithnames$interval,]
ceu_snps_genes_non_gbm <- ceu_snps_genes[ceu_snps_genes$interval %in% non_gbm_genes_poswithnames$interval,]
ceu_snps_genes_high_sift <- ceu_snps_genes[ceu_snps_genes$interval %in% high_sift_genes_poswithnames$interval,]
ceu_snps_genes_low_sift <- ceu_snps_genes[ceu_snps_genes$interval %in% low_sift_genes_poswithnames$interval,]
ceu_snps_genes_gbm$Type <- "gbM"
ceu_snps_genes_non_gbm$Type <- "non-gbM"
ceu_snps_genes_high_sift$Type <- "high SIFT-score"
ceu_snps_genes_low_sift$Type <- "low SIFT-score"
ceu_snps_genes$Class <- "SNP"
ceu_snps_genes_gbm$Class <- "SNP"
ceu_snps_genes_non_gbm$Class <- "SNP"
ceu_snps_genes_high_sift$Class <- "SNP"
ceu_snps_genes_low_sift$Class <- "SNP"

ibnr_snps_genes <- read.table(file="pi_D_ibnr.txt",header=T)
ibnr_snps_genes$Pop <- "IBnr"
ibnr_snps_genes$Type <- "All SNP"
ibnr_snps_genes$interval <- gsub(".fa","",ibnr_snps_genes$interval)
ibnr_snps_genes$pi <- ibnr_snps_genes$pi/(as.numeric(gsub(".*-","",gsub(".:","",ibnr_snps_genes$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ibnr_snps_genes$interval))))
ibnr_snps_genes_gbm <- ibnr_snps_genes[ibnr_snps_genes$interval %in% gbm_genes_poswithnames$interval,]
ibnr_snps_genes_non_gbm <- ibnr_snps_genes[ibnr_snps_genes$interval %in% non_gbm_genes_poswithnames$interval,]
ibnr_snps_genes_high_sift <- ibnr_snps_genes[ibnr_snps_genes$interval %in% high_sift_genes_poswithnames$interval,]
ibnr_snps_genes_low_sift <- ibnr_snps_genes[ibnr_snps_genes$interval %in% low_sift_genes_poswithnames$interval,]
ibnr_snps_genes_gbm$Type <- "gbM"
ibnr_snps_genes_non_gbm$Type <- "non-gbM"
ibnr_snps_genes_high_sift$Type <- "high SIFT-score"
ibnr_snps_genes_low_sift$Type <- "low SIFT-score"
ibnr_snps_genes$Class <- "SNP"
ibnr_snps_genes_gbm$Class <- "SNP"
ibnr_snps_genes_non_gbm$Class <- "SNP"
ibnr_snps_genes_high_sift$Class <- "SNP"
ibnr_snps_genes_low_sift$Class <- "SNP"

ceu_snps_intergenic <- read.table(file="pi_D_intergenic_ceu.txt",header=T)
ceu_snps_intergenic$Pop <- "CEU"
ceu_snps_intergenic$Type <- "Non-coding"
ceu_snps_intergenic$Class <- "SNP"
ceu_snps_intergenic$interval <- gsub(".fa","",ceu_snps_intergenic$interval)
ceu_snps_intergenic$pi <- ceu_snps_intergenic$pi/(as.numeric(gsub(".*-","",gsub(".:","",ceu_snps_intergenic$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ceu_snps_intergenic$interval))))

ibnr_snps_intergenic <- read.table(file="pi_D_intergenic_ibnr.txt",header=T)
ibnr_snps_intergenic$Pop <- "IBnr"
ibnr_snps_intergenic$Type <- "Non-coding"
ibnr_snps_intergenic$Class <- "SNP"
ibnr_snps_intergenic$interval <- gsub(".fa","",ibnr_snps_intergenic$interval)
ibnr_snps_intergenic$pi <- ibnr_snps_intergenic$pi/(as.numeric(gsub(".*-","",gsub(".:","",ibnr_snps_intergenic$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ibnr_snps_intergenic$interval))))

ceu_snps_cl <- read.table(file="pi_D_ceu_cg.txt",header=T)
ceu_snps_cl$Pop <- "CEU"
ceu_snps_cl$Type <- "Clocklike"
ceu_snps_cl$Class <- "SNP"
ceu_snps_cl$interval <- gsub(".fa","",ceu_snps_cl$interval)
ceu_snps_cl$pi <- ceu_snps_cl$pi/(as.numeric(gsub(".*-","",gsub(".:","",ceu_snps_cl$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ceu_snps_cl$interval))))

ibnr_snps_cl <- read.table(file="pi_D_ibnr_cg.txt",header=T)
ibnr_snps_cl$Pop <- "IBnr"
ibnr_snps_cl$Type <- "Clocklike"
ibnr_snps_cl$Class <- "SNP"
ibnr_snps_cl$interval <- gsub(".fa","",ibnr_snps_cl$interval)
ibnr_snps_cl$pi <- ibnr_snps_cl$pi/(as.numeric(gsub(".*-","",gsub(".:","",ibnr_snps_cl$interval))) - as.numeric(gsub("-.*","",gsub(".:","",ibnr_snps_cl$interval))))

ceu_smps_genes <- read.table(file="ceu_Dm.txt",header=T)
Dm_filenames_ceu <- read.table(file="Dm_filenames_ceu")
ceu_smps_genes$interval <- Dm_filenames_ceu$V1
ceu_smps_genes$Pop <- "CEU"
ceu_smps_genes$Type <- "All SMP"
colnames(ceu_smps_genes)[6] <- "pi"
colnames(ceu_smps_genes)[4] <- "D"
ceu_smps_genes <- ceu_smps_genes[(ceu_smps_genes$end - ceu_smps_genes$start) > 10,]
ceu_smps_genes$pi <- ceu_smps_genes$pi/(ceu_smps_genes$end)
ceu_smps_genes <- ceu_smps_genes[ceu_smps_genes$segregation_site > 2,]
ceu_smps_genes_gbm <- ceu_smps_genes[ceu_smps_genes$interval %in% gbm_genes_poswithnames$interval,]
ceu_smps_genes_non_gbm <- ceu_smps_genes[ceu_smps_genes$interval %in% non_gbm_genes_poswithnames$interval,]
ceu_smps_genes_high_sift <- ceu_smps_genes[ceu_smps_genes$interval %in% high_sift_genes_poswithnames$interval,]
ceu_smps_genes_low_sift <- ceu_smps_genes[ceu_smps_genes$interval %in% low_sift_genes_poswithnames$interval,]
ceu_smps_genes_gbm$Type <- "gbM"
ceu_smps_genes_non_gbm$Type <- "non-gbM"
ceu_smps_genes_high_sift$Type <- "high SIFT-score"
ceu_smps_genes_low_sift$Type <- "low SIFT-score"
ceu_smps_genes$Class <- "SMP"
ceu_smps_genes_gbm$Class <- "SMP"
ceu_smps_genes_non_gbm$Class <- "SMP"
ceu_smps_genes_high_sift$Class <- "SMP"
ceu_smps_genes_low_sift$Class <- "SMP"

ibnr_smps_genes <- read.table(file="ibnr_Dm.txt",header=T)
Dm_filenames_ibnr <- read.table(file="Dm_filenames_ibnr")
ibnr_smps_genes$interval <- Dm_filenames_ibnr$V1
ibnr_smps_genes$Pop <- "IBnr"
ibnr_smps_genes$Type <- "All SMP"
colnames(ibnr_smps_genes)[6] <- "pi"
colnames(ibnr_smps_genes)[4] <- "D"
ibnr_smps_genes <- ibnr_smps_genes[(ibnr_smps_genes$end - ibnr_smps_genes$start) > 10,]
ibnr_smps_genes$pi <- ibnr_smps_genes$pi/(ibnr_smps_genes$end)
ibnr_smps_genes <- ibnr_smps_genes[ibnr_smps_genes$segregation_site > 2,]
ibnr_smps_genes_gbm <- ibnr_smps_genes[ibnr_smps_genes$interval %in% gbm_genes_poswithnames$interval,]
ibnr_smps_genes_non_gbm <- ibnr_smps_genes[ibnr_smps_genes$interval %in% non_gbm_genes_poswithnames$interval,]
ibnr_smps_genes_high_sift <- ibnr_smps_genes[ibnr_smps_genes$interval %in% high_sift_genes_poswithnames$interval,]
ibnr_smps_genes_low_sift <- ibnr_smps_genes[ibnr_smps_genes$interval %in% low_sift_genes_poswithnames$interval,]
ibnr_smps_genes_gbm$Type <- "gbM"
ibnr_smps_genes_non_gbm$Type <- "non-gbM"
ibnr_smps_genes_high_sift$Type <- "high SIFT-score"
ibnr_smps_genes_low_sift$Type <- "low SIFT-score"
ibnr_smps_genes$Class <- "SMP"
ibnr_smps_genes_gbm$Class <- "SMP"
ibnr_smps_genes_non_gbm$Class <- "SMP"
ibnr_smps_genes_high_sift$Class <- "SMP"
ibnr_smps_genes_low_sift$Class <- "SMP"

ceu_smps_cg <- read.table(file="ceu_cg_Dm.txt",header=T)
Dm_filenames_ceu_cg <- read.table(file="Dm_filenames_ceu_cg")
ceu_smps_cg$interval <- Dm_filenames_ceu_cg$V1
ceu_smps_cg$Pop <- "CEU"
ceu_smps_cg$Type <- "Clocklike"
ceu_smps_cg$Class <- "SMP"
colnames(ceu_smps_cg)[6] <- "pi"
colnames(ceu_smps_cg)[4] <- "D"
ceu_smps_cg <- ceu_smps_cg[(ceu_smps_cg$end - ceu_smps_cg$start) > 10,]
ceu_smps_cg$pi <- ceu_smps_cg$pi/(ceu_smps_cg$end)

ibnr_smps_cg <- read.table(file="ibnr_cg_Dm.txt",header=T)
Dm_filenames_ibnr_cg <- read.table(file="Dm_filenames_ibnr_cg")
ibnr_smps_cg$interval <- Dm_filenames_ibnr_cg$V1
ibnr_smps_cg$Pop <- "IBnr"
ibnr_smps_cg$Type <- "Clocklike"
ibnr_smps_cg$Class <- "SMP"
colnames(ibnr_smps_cg)[6] <- "pi"
colnames(ibnr_smps_cg)[4] <- "D"
ibnr_smps_cg <- ibnr_smps_cg[(ibnr_smps_cg$end - ibnr_smps_cg$start) > 10,]
ibnr_smps_cg$pi <- ibnr_smps_cg$pi/(ibnr_smps_cg$end)

pi_gene <- rbind(ceu_snps_genes_gbm[c(1,4,5,3,6)],ceu_snps_genes_non_gbm[c(1,4,5,3,6)],ceu_snps_genes_high_sift[c(1,4,5,3,6)],ceu_snps_genes_low_sift[c(1,4,5,3,6)],ceu_snps_cl[c(1,4,5,3,6)],ceu_snps_intergenic[c(1,4,5,3,6)],ceu_smps_genes_gbm[c(8,9,10,6,11)],ceu_smps_genes_non_gbm[c(8,9,10,6,11)],ceu_smps_genes_high_sift[c(8,9,10,6,11)],ceu_smps_genes_low_sift[c(8,9,10,6,11)],ceu_smps_cg[c(8,9,10,6,11)],ibnr_snps_genes_gbm[c(1,4,5,3,6)],ibnr_snps_genes_non_gbm[c(1,4,5,3,6)],ibnr_snps_genes_high_sift[c(1,4,5,3,6)],ibnr_snps_genes_low_sift[c(1,4,5,3,6)],ibnr_snps_cl[c(1,4,5,3,6)],ibnr_snps_intergenic[c(1,4,5,3,6)],ibnr_smps_genes_gbm[c(8,9,10,6,11)],ibnr_smps_genes_non_gbm[c(8,9,10,6,11)],ibnr_smps_genes_high_sift[c(8,9,10,6,11)],ibnr_smps_genes_low_sift[c(8,9,10,6,11)],ibnr_smps_cg[c(8,9,10,6,11)])
pi_gene <- pi_gene[!pi_gene$pi == 0,]
pi_gene$Type <- factor(pi_gene$Type,levels = c("Non-coding","Clocklike","gbM","non-gbM","high SIFT-score","low SIFT-score"))
pi_gene_snps <- pi_gene[pi_gene$Class %in% "SNP",]
pi_gene_smps <- pi_gene[pi_gene$Class %in% "SMP",]
pi_gene$pi <- log10(pi_gene$pi)
df_thaliana_pi_gene <- dplyr::count(pi_gene,Type,Class,Pop)
df_thaliana_pi_gene$y <- 0
df_thaliana_pi_gene[df_thaliana_pi_gene$Class %in% "SNP",]$y <- -1

pi_gene_snps$treatment <- paste(pi_gene_snps$Type,pi_gene_snps$Pop,sep="_")
summary(aov(lm(data=pi_gene_snps,pi~treatment)))
snp_test_pvals <- HSD.test(aov(lm(data=pi_gene_snps,pi~treatment)), trt = 'treatment')
snp_test_pvals

pi_gene_smps$treatment <- paste(pi_gene_smps$Type,pi_gene_smps$Pop,sep="_")
summary(aov(lm(data=pi_gene_smps,pi~treatment)))
smp_test_pvals <- HSD.test(aov(lm(data=pi_gene_smps,pi~treatment)), trt = 'treatment')
smp_test_pvals

segments <- data.frame(x=0.7,xend=0.9,y = -3, yend = -3,Class="SNP",Pop="CEU")
segments <- rbind(segments,c(1.1,1.3,-3,-3,"SNP","IBnr"))
segments <- rbind(segments,c(1.7,1.9,-4,-4,"SNP","CEU"))
segments <- rbind(segments,c(2.1,6.3,-4,-4,"SNP","IBnr"))
segments <- rbind(segments,c(0.7,0.9,-1,-1,"SMP","CEU"))
segments <- rbind(segments,c(1.1,1.3,-1,-1,"SMP","IBnr"))
segments <- rbind(segments,c(1.7,5.3,-1.3,-1.3,"SMP","IBnr"))
segments$x <- as.numeric(segments$x)
segments$xend <- as.numeric(segments$xend)
segments$y <- as.numeric(segments$y)
segments$yend <- as.numeric(segments$yend)

grouptext <- data.frame(x=0.8,y=-2.5,label="a",Class="SNP",Pop="CEU")
grouptext <- rbind(grouptext,c(1.2,-2.5,"b", "SNP","IBnr"))
grouptext <- rbind(grouptext,c(1.8,-3.5,"c","SNP","CEU"))
grouptext <- rbind(grouptext,c(4.2,-3.5,"d","SNP","IBnr"))
grouptext <- rbind(grouptext,c(0.8,-0.7,"a","SMP","CEU"))
grouptext <- rbind(grouptext,c(1.2,-0.7,"b","SMP","IBnr"))
grouptext <- rbind(grouptext,c(3.5,-1,"c","SMP","CEU"))
grouptext$x <- as.numeric(grouptext$x)
grouptext$y <- as.numeric(grouptext$y)

pi_gene_plot <- ggplot(data=pi_gene,aes(x=Type,y=pi,fill=Pop))+
  geom_boxplot() +
  ylab(expression('log'[10]*'(Scaled '*pi*')')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.title=element_blank()) + 
  facet_wrap(~Class,scales="free")  +
  xlab("Region") +
  geom_segment(data=segments,aes(x = x, xend = xend, y = y, yend = yend)) +
  geom_text(data=grouptext,aes(x=x,y=y,label=label)) +
  geom_text_repel(data = df_thaliana_pi_gene, aes(y = y, angle=90,label = n),direction="x")

pi_ceu_snp <- ceu_snps_genes[c(1,3)]
colnames(pi_ceu_snp)[2] <- "pi_SNP"
pi_ceu_smp <- ceu_smps_genes[c(8,6)]
colnames(pi_ceu_smp)[2] <- "pi_SMP"
pi_ceu_snp_smp <- inner_join(pi_ceu_snp,pi_ceu_smp,by="interval")
ceu_pi_compare <- ggscatter(data=pi_ceu_snp_smp,x="pi_SNP",y="pi_SMP", shape = 21, cor.coef = T,cor.coef.coord = c(0.000003, 0.01),title = paste("CEU, n=", nrow(pi_ceu_snp_smp),sep=""),size=1) +
  xlab(expression(pi[SNP])) +
  ylab(expression(pi[SMP]))

pi_ibnr_snp <- ibnr_snps_genes[c(1,3)]
colnames(pi_ibnr_snp)[2] <- "pi_SNP"
pi_ibnr_smp <- ibnr_smps_genes[c(8,6)]
colnames(pi_ibnr_smp)[2] <- "pi_SMP"
pi_ibnr_snp_smp <- inner_join(pi_ibnr_snp,pi_ibnr_smp,by="interval")
ibnr_pi_compare <- ggscatter(data=pi_ibnr_snp_smp,x="pi_SNP",y="pi_SMP", shape = 21, cor.coef = T,cor.coef.coord = c(0.000005, 0.02),title = paste("IBnr, n=", nrow(pi_ibnr_snp_smp),sep=""),size=1) +
  xlab(expression(pi[SNP])) +
  ylab(expression(pi[SMP]))

pi_plot <- plot_grid(pi_gene_plot,plot_grid(ceu_pi_compare,ibnr_pi_compare,ncol=1),ncol=2,rel_widths = c(2,1),labels = c("A","C"))

#### TajimasD ####

D_gene <- rbind(ceu_snps_genes_gbm[c(1,4,5,2,6)],ceu_snps_genes_non_gbm[c(1,4,5,2,6)],ceu_snps_genes_high_sift[c(1,4,5,2,6)],ceu_snps_genes_low_sift[c(1,4,5,2,6)],ceu_snps_cl[c(1,4,5,2,6)],ceu_snps_intergenic[c(1,4,5,2,6)],ceu_smps_genes_gbm[c(8,9,10,4,11)],ceu_smps_genes_non_gbm[c(8,9,10,4,11)],ceu_smps_genes_high_sift[c(8,9,10,4,11)],ceu_smps_genes_low_sift[c(8,9,10,4,11)],ceu_smps_cg[c(8,9,10,4,11)],ibnr_snps_genes_gbm[c(1,4,5,2,6)],ibnr_snps_genes_non_gbm[c(1,4,5,2,6)],ibnr_snps_genes_high_sift[c(1,4,5,2,6)],ibnr_snps_genes_low_sift[c(1,4,5,2,6)],ibnr_snps_cl[c(1,4,5,2,6)],ibnr_snps_intergenic[c(1,4,5,2,6)],ibnr_smps_genes_gbm[c(8,9,10,4,11)],ibnr_smps_genes_non_gbm[c(8,9,10,4,11)],ibnr_smps_genes_high_sift[c(8,9,10,4,11)],ibnr_smps_genes_low_sift[c(8,9,10,4,11)],ibnr_smps_cg[c(8,9,10,4,11)])
D_gene$Type <- factor(D_gene$Type,levels = c("Non-coding","Clocklike","gbM","non-gbM","high SIFT-score","low SIFT-score"))
df_thaliana_D_gene <- dplyr::count(D_gene,Type,Class,Pop)

D_gene_snps <- D_gene[D_gene$Class %in% "SNP",]
D_gene_smps <- D_gene[D_gene$Class %in% "SMP",]
D_gene_snps$treatment <- paste(D_gene_snps$Type,D_gene_snps$Pop,sep="_")
summary(aov(lm(data=D_gene_snps,D~treatment)))
snp_Dtest_pvals <- HSD.test(aov(lm(data=D_gene_snps,D~treatment)), trt = 'treatment')
snp_Dtest_pvals

D_gene_smps$treatment <- paste(D_gene_smps$Type,D_gene_smps$Pop,sep="_")
summary(aov(lm(data=D_gene_smps,D~treatment)))
smp_Dtest_pvals <- HSD.test(aov(lm(data=D_gene_smps,D~treatment)), trt = 'treatment')
smp_Dtest_pvals

grouptext_D <- data.frame(x=0.8,y=5.5,label="a",Class="SMP",Pop="CEU")
grouptext_D <- rbind(grouptext_D,c(1.2,5.5,"c","SMP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(1.8,5.5,"b","SMP","CEU"))
grouptext_D <- rbind(grouptext_D,c(2.2,5.5,"de","SMP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(2.8,5.5,"b","SMP","CEU"))
grouptext_D <- rbind(grouptext_D,c(3.2,5.5,"f","SMP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(3.8,5.5,"b","SMP","CEU"))
grouptext_D <- rbind(grouptext_D,c(4.2,5.5,"ef","SMP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(4.8,5.5,"b","SMP","CEU"))
grouptext_D <- rbind(grouptext_D,c(5.2,5.5,"d","SMP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(0.8,5.5,"h","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(1.2,5.5,"ef","SNP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(1.8,5.5,"b","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(2.2,5.5,"a","SNP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(2.8,5.5,"fg","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(3.2,5.5,"cd","SNP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(3.8,5.5,"g","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(4.2,5.5,"d","SNP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(4.8,5.5,"fg","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(5.2,5.5,"c","SNP","IBnr"))
grouptext_D <- rbind(grouptext_D,c(5.8,5.5,"h","SNP","CEU"))
grouptext_D <- rbind(grouptext_D,c(6.2,5.5,"e","SNP","IBnr"))
grouptext_D$x <- as.numeric(grouptext_D$x)
grouptext_D$y <- as.numeric(grouptext_D$y)

D_gene_plot <- ggplot(data=D_gene,aes(x=Type,y=D,fill=Pop))+
  geom_boxplot() +
  ylab(expression("Tajima's "*italic('D'))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),,legend.title=element_blank()) + 
  ylim(-3,7.5) +
  xlab("Region") +
  facet_wrap(~Class,scales="free") +
  geom_text(data=grouptext_D,aes(x=x,y=y,label=label))+
  geom_text_repel(data = df_thaliana_D_gene, aes(y = 7, angle=90,label = n),direction="x")

D_ceu_snp <- ceu_snps_genes[c(1,2)]
colnames(D_ceu_snp)[2] <- "D_SNP"
D_ceu_smp <- ceu_smps_genes[c(8,4)]
colnames(D_ceu_smp)[2] <- "D_SMP"
D_ceu_snp_smp <- inner_join(D_ceu_snp,D_ceu_smp,by="interval")
ceu_D_compare <- ggscatter(data=D_ceu_snp_smp,x="D_SNP",y="D_SMP", shape = 21, cor.coef = T,cor.coef.coord = c(-2, 4),title = paste("CEU, n=", nrow(D_ceu_snp_smp),sep=""),size=1) + 
  xlab(expression(italic('D')[SNP])) +
  ylab(expression(italic('D')[SMP]))

D_ibnr_snp <- ibnr_snps_genes[c(1,2)]
colnames(D_ibnr_snp)[2] <- "D_SNP"
D_ibnr_smp <- ibnr_smps_genes[c(8,4)]
colnames(D_ibnr_smp)[2] <- "D_SMP"
D_ibnr_snp_smp <- inner_join(D_ibnr_snp,D_ibnr_smp,by="interval")
ibnr_D_compare <- ggscatter(data=D_ibnr_snp_smp,x="D_SNP",y="D_SMP", shape = 21, cor.coef = T,cor.coef.coord = c(-2, 4),title = paste("IBnr, n=", nrow(D_ibnr_snp_smp),sep=""),size=1) + 
  xlab(expression(italic('D')[SNP])) +
  ylab(expression(italic('D')[SMP]))

D_plot <- plot_grid(D_gene_plot,plot_grid(ceu_D_compare,ibnr_D_compare,ncol=1),ncol=2,rel_widths = c(2,1),labels = c("B","D"))

pdf(file="pi_D.pdf",height=11,width=12)
plot_grid(pi_plot,D_plot,ncol=1)
dev.off()

smps_cg_all <- inner_join(ceu_smps_cg[c(3:6,8:9)],ibnr_smps_cg[c(3:6,8:9)],by="interval")
smps_genes_all <- inner_join(ceu_smps_genes[c(3:6,8:9)],ibnr_smps_genes[c(3:6,8:9)],by="interval")
pa <- ggscatter(data=smps_genes_all,x="D.x",y="D.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,cor.coef.coord = c(-1, 2.5),title = paste("Coding region, n=", nrow(smps_genes_all),sep=""),size=1) +
  xlab(expression("CEU "*italic('D')[SMP])) +
  ylab(expression("IBnr "*italic('D')[SMP])) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pb <- ggscatter(data=smps_genes_all,x="pi.x",y="pi.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,title = paste("Coding region, n=", nrow(smps_genes_all),sep=""),size=1) +
  xlab(expression(paste("CEU ",pi[SMP]))) +
  ylab(expression(paste("IBnr ",pi[SMP]))) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pc <- ggscatter(data=smps_genes_all,x="segregation_site.x",y="segregation_site.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,title = paste("Coding region, n=", nrow(smps_genes_all),sep=""),size=1) +
  xlab("CEU number of segregating sites") +
  ylab("IBnr number of segregating sites") +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pd <- ggscatter(data=smps_cg_all,x="D.x",y="D.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,cor.coef.coord = c(-1.5, 2.5),title = paste("Clocklike region, n=", nrow(smps_cg_all),sep=""),size=1) +
  xlab(expression("CEU "*italic('D')[SMP])) +
  ylab(expression("IBnr "*italic('D')[SMP])) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pe <- ggscatter(data=smps_cg_all,x="pi.x",y="pi.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,title = paste("Clocklike region, n=", nrow(smps_cg_all),sep=""),size=1) +
  xlab(expression(paste("CEU ",pi[SMP]))) +
  ylab(expression(paste("IBnr ",pi[SMP]))) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pf <- ggscatter(data=smps_cg_all,x="segregation_site.x",y="segregation_site.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,title = paste("Clocklike region, n=", nrow(smps_cg_all),sep=""),size=1) +
  xlab("CEU number of segregating sites") +
  ylab("IBnr number of segregating sites") +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")

pdf(file="D_CEU_IBNR.pdf",height=6,width=10)
plot_grid(pa,pb,pc,pd,pe,pf,ncol=3)
dev.off()

#### SFS  ####
ceu_smps.frq <- read.table(file="ceu_smps.frq",row.names = NULL, header = T)
ceu_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_smps.frq$N_CHR))
ggplot(ceu_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU SMP")


ibnr_smps.frq <- read.table(file="ibnr_smps.frq",row.names = NULL, header = T)
ibnr_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_smps.frq$N_CHR))
ggplot(ibnr_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count)))) +
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites")+
  ggtitle("IBNR SMP")

ceu_dmrs.frq <- read.table(file="ceu_dmrs.frq",row.names = NULL, header = T)
ceu_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs.frq$N_CHR))
ggplot(ceu_dmrs.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU DMR")

ibnr_dmrs.frq <- read.table(file="ibnr_dmrs.frq",row.names = NULL, header = T)
ibnr_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs.frq$N_CHR))
ggplot(ibnr_dmrs.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR DMR")

ceu_gbm.frq <- read.table(file="ceu_gbm.frq",row.names = NULL, header = T)
ceu_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ceu_gbm.frq$N_CHR))
ggplot(ceu_gbm.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU gBM")

ceu_dmrs_gbm.frq <- read.table(file="ceu_dmrs_gbm.frq",row.names = NULL, header = T)
ceu_dmrs_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs_gbm.frq$N_CHR))
ggplot(ceu_dmrs_gbm.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU DMR gBM")


ibnr_gbm.frq <- read.table(file="ibnr_gbm.frq",row.names = NULL, header = T)
ibnr_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_gbm.frq$N_CHR))
ggplot(ibnr_gbm.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR gBM")

ibnr_dmrs_gbm.frq <- read.table(file="ibnr_dmrs_gbm.frq",row.names = NULL, header = T)
ibnr_dmrs_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs_gbm.frq$N_CHR))
ggplot(ibnr_dmrs_gbm.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR DMR gBM")

ceu_sift.frq <- read.table(file="ceu_sift.frq",row.names = NULL, header = T)
ceu_sift.frq$N_CHR <- as.numeric(gsub("A:","",ceu_sift.frq$N_CHR))
ggplot(ceu_sift.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU SIFT")

ceu_dmrs_sift.frq <- read.table(file="ceu_dmrs_sift.frq",row.names = NULL, header = T)
ceu_dmrs_sift.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs_sift.frq$N_CHR))
ggplot(ceu_dmrs_sift.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU DMR SIFT")

ibnr_sift.frq <- read.table(file="ibnr_sift.frq",row.names = NULL, header = T)
ibnr_sift.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_sift.frq$N_CHR))
ggplot(ibnr_sift.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR SIFT")

ibnr_dmrs_sift.frq <- read.table(file="ibnr_dmrs_sift.frq",row.names = NULL, header = T)
ibnr_dmrs_sift.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs_sift.frq$N_CHR))
ggplot(ibnr_dmrs_sift.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR DMR SIFT")


ceu_all_smps.frq <- read.table(file="ceu_all_genes.frq",row.names = NULL, header = T)
ceu_all_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_all_smps.frq$N_CHR))
ggplot(ceu_all_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("CEU all SMP")

ibnr_all_smps.frq <- read.table(file="ibnr_all_genes.frq",row.names = NULL, header = T)
ibnr_all_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_all_smps.frq$N_CHR))
ggplot(ibnr_all_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated alleles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR all SMP")


ceu_non_gbm_smps.frq <- read.table(file="ceu_non_gbm.frq",row.names = NULL, header = T)
ceu_non_gbm_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_non_gbm_smps.frq$N_CHR))
ggplot(ceu_non_gbm_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated non_gbmeles") +
  ylab("Proportion of sites") +
  ggtitle("CEU non gbm SMP")

ibnr_non_gbm_smps.frq <- read.table(file="ibnr_non_gbm.frq",row.names = NULL, header = T)
ibnr_non_gbm_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_non_gbm_smps.frq$N_CHR))
ggplot(ibnr_non_gbm_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated non_gbmeles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR non gbm SMP")


ceu_low_sift_smps.frq <- read.table(file="ceu_low_sift.frq",row.names = NULL, header = T)
ceu_low_sift_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_low_sift_smps.frq$N_CHR))
ggplot(ceu_low_sift_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated low_sifteles") +
  ylab("Proportion of sites") +
  ggtitle("CEU low SIFT SMP")

ibnr_low_sift_smps.frq <- read.table(file="ibnr_low_sift.frq",row.names = NULL, header = T)
ibnr_low_sift_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_low_sift_smps.frq$N_CHR))
ggplot(ibnr_low_sift_smps.frq,aes(x=N_CHR)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))))+
  xlab("Proportion of unmethylated low_sifteles") +
  ylab("Proportion of sites") +
  ggtitle("IBNR low SIFT SMP")



## maf > 0.02

ceu_smps.frq <- read.table(file="ceu_smps.frq",row.names = NULL, header = T)
ceu_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_smps.frq$N_CHR))
ceu_smps.frq <- ceu_smps.frq[(ceu_smps.frq$N_CHR > 0.02) & (ceu_smps.frq$N_CHR < 0.98),]
ceu_smps.frq$Population <- "CEU"
ibnr_smps.frq <- read.table(file="ibnr_smps.frq",row.names = NULL, header = T)
ibnr_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_smps.frq$N_CHR))
ibnr_smps.frq <- ibnr_smps.frq[(ibnr_smps.frq$N_CHR > 0.02) & (ibnr_smps.frq$N_CHR < 0.98),]
ibnr_smps.frq$Population <- "IBnr"
smps.frq <- rbind(ceu_smps.frq,ibnr_smps.frq)
smps.frq <- smps.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
smps.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
smps.frq.plot <- ggplot(smps.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_smps.frq),nrow(ibnr_smps.frq))) +
  ggtitle("Clocklike SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_dmrs.frq <- read.table(file="ceu_dmrs.frq",row.names = NULL, header = T)
ceu_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs.frq$N_CHR))
ceu_dmrs.frq <- ceu_dmrs.frq[(ceu_dmrs.frq$N_CHR > 0.02) & (ceu_dmrs.frq$N_CHR < 0.98),]
ceu_dmrs.frq$Population <- "CEU"
ibnr_dmrs.frq <- read.table(file="ibnr_dmrs.frq",row.names = NULL, header = T)
ibnr_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs.frq$N_CHR))
ibnr_dmrs.frq <- ibnr_dmrs.frq[(ibnr_dmrs.frq$N_CHR > 0.02) & (ibnr_dmrs.frq$N_CHR < 0.98),]
ibnr_dmrs.frq$Population <- "IBnr"
dmrs.frq <- rbind(ceu_dmrs.frq,ibnr_dmrs.frq)
dmrs.frq <- dmrs.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
#dmrs.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
dmrs.frq$N_CHR_binned <- (as.numeric(gsub("\\)|\\]","", gsub(".*,","", dmrs.frq$N_CHR_binned))) + as.numeric(gsub("\\[","", gsub(",.*","", dmrs.frq$N_CHR_binned))) ) / 2
dmrs.frq.plot <- ggplot(dmrs.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="")  +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_dmrs.frq),nrow(ibnr_dmrs.frq))) +
  ggtitle("Clocklike DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_gbm.frq <- read.table(file="ceu_gbm.frq",row.names = NULL, header = T)
ceu_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ceu_gbm.frq$N_CHR))
ceu_gbm.frq <- ceu_gbm.frq[(ceu_gbm.frq$N_CHR > 0.02) & (ceu_gbm.frq$N_CHR < 0.98),]
ceu_gbm.frq$Population <- "CEU"
ibnr_gbm.frq <- read.table(file="ibnr_gbm.frq",row.names = NULL, header = T)
ibnr_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_gbm.frq$N_CHR))
ibnr_gbm.frq <- ibnr_gbm.frq[(ibnr_gbm.frq$N_CHR > 0.02) & (ibnr_gbm.frq$N_CHR < 0.98),]
ibnr_gbm.frq$Population <- "IBnr"
gbm.frq <- rbind(ceu_gbm.frq,ibnr_gbm.frq)
gbm.frq <- gbm.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
gbm.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
gbm.frq.plot <- ggplot(gbm.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_gbm.frq),nrow(ibnr_gbm.frq))) +
  ggtitle("gbM SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_dmrs_gbm.frq <- read.table(file="ceu_dmrs_gbm.frq",row.names = NULL, header = T)
ceu_dmrs_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs_gbm.frq$N_CHR))
ceu_dmrs_gbm.frq <- ceu_dmrs_gbm.frq[(ceu_dmrs_gbm.frq$N_CHR > 0.02) & (ceu_dmrs_gbm.frq$N_CHR < 0.98),]
ceu_dmrs_gbm.frq$Population <- "CEU"
ibnr_dmrs_gbm.frq <- read.table(file="ibnr_dmrs_gbm.frq",row.names = NULL, header = T)
ibnr_dmrs_gbm.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs_gbm.frq$N_CHR))
ibnr_dmrs_gbm.frq <- ibnr_dmrs_gbm.frq[(ibnr_dmrs_gbm.frq$N_CHR > 0.02) & (ibnr_dmrs_gbm.frq$N_CHR < 0.98),]
ibnr_dmrs_gbm.frq$Population <- "IBnr"
dmrs_gbm.frq <- rbind(ceu_dmrs_gbm.frq,ibnr_dmrs_gbm.frq)
dmrs_gbm.frq <- dmrs_gbm.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
dmrs_gbm.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
dmrs_gbm.frq.plot <- ggplot(dmrs_gbm.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_dmrs_gbm.frq),nrow(ibnr_dmrs_gbm.frq))) +
  ggtitle("gbM DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_sift.frq <- read.table(file="ceu_sift.frq",row.names = NULL, header = T)
ceu_sift.frq$N_CHR <- as.numeric(gsub("A:","",ceu_sift.frq$N_CHR))
ceu_sift.frq <- ceu_sift.frq[(ceu_sift.frq$N_CHR > 0.02) & (ceu_sift.frq$N_CHR < 0.98),]
ceu_sift.frq$Population <- "CEU"
ibnr_sift.frq <- read.table(file="ibnr_sift.frq",row.names = NULL, header = T)
ibnr_sift.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_sift.frq$N_CHR))
ibnr_sift.frq <- ibnr_sift.frq[(ibnr_sift.frq$N_CHR > 0.02) & (ibnr_sift.frq$N_CHR < 0.98),]
ibnr_sift.frq$Population <- "IBnr"
sift.frq <- rbind(ceu_sift.frq,ibnr_sift.frq)
sift.frq <- sift.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
sift.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
sift.frq.plot <- ggplot(sift.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_sift.frq),nrow(ibnr_sift.frq))) +
  ggtitle("High SIFT SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_dmrs_sift.frq <- read.table(file="ceu_dmrs_sift.frq",row.names = NULL, header = T)
ceu_dmrs_sift.frq$N_CHR <- as.numeric(gsub("A:","",ceu_dmrs_sift.frq$N_CHR))
ceu_dmrs_sift.frq <- ceu_dmrs_sift.frq[(ceu_dmrs_sift.frq$N_CHR > 0.02) & (ceu_dmrs_sift.frq$N_CHR < 0.98),]
ceu_dmrs_sift.frq$Population <- "CEU"
ibnr_dmrs_sift.frq <- read.table(file="ibnr_dmrs_sift.frq",row.names = NULL, header = T)
ibnr_dmrs_sift.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_dmrs_sift.frq$N_CHR))
ibnr_dmrs_sift.frq <- ibnr_dmrs_sift.frq[(ibnr_dmrs_sift.frq$N_CHR > 0.02) & (ibnr_dmrs_sift.frq$N_CHR < 0.98),]
ibnr_dmrs_sift.frq$Population <- "IBnr"
dmrs_sift.frq <- rbind(ceu_dmrs_sift.frq,ibnr_dmrs_sift.frq)
dmrs_sift.frq <- dmrs_sift.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
dmrs_sift.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
dmrs_sift.frq.plot <- ggplot(dmrs_sift.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_dmrs_sift.frq),nrow(ibnr_dmrs_sift.frq))) +
  ggtitle("High SIFT DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_all_smps.frq <- read.table(file="ceu_all_genes.frq",row.names = NULL, header = T)
ceu_all_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_all_smps.frq$N_CHR))
ceu_all_smps.frq <- ceu_all_smps.frq[(ceu_all_smps.frq$N_CHR > 0.02) & (ceu_all_smps.frq$N_CHR < 0.98),]
ceu_all_smps.frq$Population <- "CEU"
ibnr_all_smps.frq <- read.table(file="ibnr_all_genes.frq",row.names = NULL, header = T)
ibnr_all_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_all_smps.frq$N_CHR))
ibnr_all_smps.frq <- ibnr_all_smps.frq[(ibnr_all_smps.frq$N_CHR > 0.02) & (ibnr_all_smps.frq$N_CHR < 0.98),]
ibnr_all_smps.frq$Population <- "IBnr"
all_smps.frq <- rbind(ceu_all_smps.frq,ibnr_all_smps.frq)
all_smps.frq <- all_smps.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
all_smps.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
all_smps.frq.plot <- ggplot(all_smps.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_all_smps.frq),nrow(ibnr_all_smps.frq))) +
  ggtitle("All SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))


ceu_non_gbm_smps.frq <- read.table(file="ceu_non_gbm.frq",row.names = NULL, header = T)
ceu_non_gbm_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_non_gbm_smps.frq$N_CHR))
ceu_non_gbm_smps.frq <- ceu_non_gbm_smps.frq[(ceu_non_gbm_smps.frq$N_CHR > 0.02) & (ceu_non_gbm_smps.frq$N_CHR < 0.98),]
ceu_non_gbm_smps.frq$Population <- "CEU"
ibnr_non_gbm_smps.frq <- read.table(file="ibnr_non_gbm.frq",row.names = NULL, header = T)
ibnr_non_gbm_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_non_gbm_smps.frq$N_CHR))
ibnr_non_gbm_smps.frq <- ibnr_non_gbm_smps.frq[(ibnr_non_gbm_smps.frq$N_CHR > 0.02) & (ibnr_non_gbm_smps.frq$N_CHR < 0.98),]
ibnr_non_gbm_smps.frq$Population <- "IBnr"
non_gbm_smps.frq <- rbind(ceu_non_gbm_smps.frq,ibnr_non_gbm_smps.frq)
non_gbm_smps.frq <- non_gbm_smps.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
non_gbm_smps.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
non_gbm_smps.frq.plot <- ggplot(non_gbm_smps.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_non_gbm_smps.frq),nrow(ibnr_non_gbm_smps.frq))) +
  ggtitle("non-gbM SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_low_sift_smps.frq <- read.table(file="ceu_low_sift.frq",row.names = NULL, header = T)
ceu_low_sift_smps.frq$N_CHR <- as.numeric(gsub("A:","",ceu_low_sift_smps.frq$N_CHR))
ceu_low_sift_smps.frq <- ceu_low_sift_smps.frq[(ceu_low_sift_smps.frq$N_CHR > 0.02) & (ceu_low_sift_smps.frq$N_CHR < 0.98),]
ceu_low_sift_smps.frq$Population <- "CEU"
ibnr_low_sift_smps.frq <- read.table(file="ibnr_low_sift.frq",row.names = NULL, header = T)
ibnr_low_sift_smps.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_low_sift_smps.frq$N_CHR))
ibnr_low_sift_smps.frq <- ibnr_low_sift_smps.frq[(ibnr_low_sift_smps.frq$N_CHR > 0.02) & (ibnr_low_sift_smps.frq$N_CHR < 0.98),]
ibnr_low_sift_smps.frq$Population <- "IBnr"
low_sift_smps.frq <- rbind(ceu_non_gbm_smps.frq,ibnr_non_gbm_smps.frq)
low_sift_smps.frq <- low_sift_smps.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
low_sift_smps.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
low_sift_smps.frq.plot <- ggplot(low_sift_smps.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_low_sift_smps.frq),nrow(ibnr_low_sift_smps.frq))) +
  ggtitle("Low SIFT SMP") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))


ceu_all_dmrs.frq <- read.table(file="ceu_dmrs_all.frq",row.names = NULL, header = T)
ceu_all_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ceu_all_dmrs.frq$N_CHR))
ceu_all_dmrs.frq <- ceu_all_dmrs.frq[(ceu_all_dmrs.frq$N_CHR > 0.02) & (ceu_all_dmrs.frq$N_CHR < 0.98),]
ceu_all_dmrs.frq$Population <- "CEU"
ibnr_all_dmrs.frq <- read.table(file="ibnr_dmrs_all.frq",row.names = NULL, header = T)
ibnr_all_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_all_dmrs.frq$N_CHR))
ibnr_all_dmrs.frq <- ibnr_all_dmrs.frq[(ibnr_all_dmrs.frq$N_CHR > 0.02) & (ibnr_all_dmrs.frq$N_CHR < 0.98),]
ibnr_all_dmrs.frq$Population <- "IBnr"
all_dmrs.frq <- rbind(ceu_all_dmrs.frq,ibnr_all_dmrs.frq)
all_dmrs.frq <- all_dmrs.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
all_dmrs.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
all_dmrs.frq.plot <- ggplot(all_dmrs.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_all_dmrs.frq),nrow(ibnr_all_dmrs.frq))) +
  ggtitle("All DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_non_gbm_dmrs.frq <- read.table(file="ceu_dmrs_non_gbm.frq",row.names = NULL, header = T)
ceu_non_gbm_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ceu_non_gbm_dmrs.frq$N_CHR))
ceu_non_gbm_dmrs.frq <- ceu_non_gbm_dmrs.frq[(ceu_non_gbm_dmrs.frq$N_CHR > 0.02) & (ceu_non_gbm_dmrs.frq$N_CHR < 0.98),]
ceu_non_gbm_dmrs.frq$Population <- "CEU"
ibnr_non_gbm_dmrs.frq <- read.table(file="ibnr_dmrs_non_gbm.frq",row.names = NULL, header = T)
ibnr_non_gbm_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_non_gbm_dmrs.frq$N_CHR))
ibnr_non_gbm_dmrs.frq <- ibnr_non_gbm_dmrs.frq[(ibnr_non_gbm_dmrs.frq$N_CHR > 0.02) & (ibnr_non_gbm_dmrs.frq$N_CHR < 0.98),]
ibnr_non_gbm_dmrs.frq$Population <- "IBnr"
non_gbm_dmrs.frq <- rbind(ceu_non_gbm_dmrs.frq,ibnr_non_gbm_dmrs.frq)
non_gbm_dmrs.frq <- non_gbm_dmrs.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
non_gbm_dmrs.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
non_gbm_dmrs.frq.plot <- ggplot(non_gbm_dmrs.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_non_gbm_dmrs.frq),nrow(ibnr_non_gbm_dmrs.frq))) +
  ggtitle("non-gbM DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

ceu_low_sift_dmrs.frq <- read.table(file="ceu_dmrs_low_sift.frq",row.names = NULL, header = T)
ceu_low_sift_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ceu_low_sift_dmrs.frq$N_CHR))
ceu_low_sift_dmrs.frq <- ceu_low_sift_dmrs.frq[(ceu_low_sift_dmrs.frq$N_CHR > 0.02) & (ceu_low_sift_dmrs.frq$N_CHR < 0.98),]
ceu_low_sift_dmrs.frq$Population <- "CEU"
ibnr_low_sift_dmrs.frq <- read.table(file="ibnr_dmrs_low_sift.frq",row.names = NULL, header = T)
ibnr_low_sift_dmrs.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_low_sift_dmrs.frq$N_CHR))
ibnr_low_sift_dmrs.frq <- ibnr_low_sift_dmrs.frq[(ibnr_low_sift_dmrs.frq$N_CHR > 0.02) & (ibnr_low_sift_dmrs.frq$N_CHR < 0.98),]
ibnr_low_sift_dmrs.frq$Population <- "IBnr"
low_sift_dmrs.frq <- rbind(ceu_low_sift_dmrs.frq,ibnr_low_sift_dmrs.frq)
forleg <- ggplot(low_sift_dmrs.frq,aes(x=N_CHR,fill=Population)) +
  geom_histogram() +
  theme(legend.position = "top")
forleg <- get_legend(forleg)
low_sift_dmrs.frq <- low_sift_dmrs.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
low_sift_dmrs.frq <- rbind(as.data.frame(low_sift_dmrs.frq),c("CEU","[0.22,0.26)",0,0))
low_sift_dmrs.frq <- rbind(as.data.frame(low_sift_dmrs.frq),c("CEU","[0.38,0.42)",0,0))
low_sift_dmrs.frq <- rbind(as.data.frame(low_sift_dmrs.frq),c("CEU","[0.54,0.58)",0,0))
low_sift_dmrs.frq <- rbind(as.data.frame(low_sift_dmrs.frq),c("CEU","[0.74,0.78)",0,0))
low_sift_dmrs.frq <- rbind(as.data.frame(low_sift_dmrs.frq),c("IBnr","[0.54,0.58)",0,0))
low_sift_dmrs.frq <- low_sift_dmrs.frq[order(low_sift_dmrs.frq$Population,low_sift_dmrs.frq$N_CHR_binned),]
low_sift_dmrs.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
low_sift_dmrs.frq$proportion <- as.numeric(low_sift_dmrs.frq$proportion)
low_sift_dmrs.frq$count <- as.numeric(low_sift_dmrs.frq$count)
low_sift_dmrs.frq.plot <- ggplot(low_sift_dmrs.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="", y="") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_low_sift_dmrs.frq),nrow(ibnr_low_sift_dmrs.frq))) +
  ggtitle("Low SIFT DMR") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))

y.grob <- textGrob("Proportion of sites", gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Proportion of unmethylated alleles", gp=gpar(fontface="bold", col="black", fontsize=15))

sfs_plots <- plot_grid(smps.frq.plot,gbm.frq.plot,sift.frq.plot,all_smps.frq.plot,non_gbm_smps.frq.plot,low_sift_smps.frq.plot,dmrs.frq.plot,dmrs_gbm.frq.plot,dmrs_sift.frq.plot,all_dmrs.frq.plot,non_gbm_dmrs.frq.plot,low_sift_dmrs.frq.plot,ncol=3)
pdf(file="sfs.pdf",height=7.5,width=8.5)
grid.arrange(arrangeGrob(sfs_plots, left = y.grob, bottom = x.grob, top=forleg))
dev.off()

#### Methylation rate estimation #### 

sitemodel_ceu <- read.table(file="model1_meth_rate_ceu.txt")
sitemodel_ceu$Population <- "CEU"
sitemodel_ibnr <- read.table(file="model1_meth_rate_ibnr.txt")
sitemodel_ibnr$Population <- "IBnr"
sitemodel <- rbind(sitemodel_ceu,sitemodel_ibnr)
colnames(sitemodel)[1:2] <- c("Forward","Backward")
sitemodel$Chromosome <- c(seq(1,5),seq(1,5))

a <- ggplot(data = sitemodel,aes(x=Chromosome,y=Forward,fill=Population)) +
  geom_col(position = "dodge") +
  ylab("Forward rate") +
  theme(legend.position = "none")
b <- ggplot(data = sitemodel,aes(x=Chromosome,y=Backward,fill=Population)) +
  geom_col(position = "dodge") +
  ylab("Backward rate") +
  theme(legend.position = "none")
c <- ggplot(data = sitemodel,aes(x=Chromosome,y=Backward/Forward,fill=Population)) +
  geom_col(position = "dodge") +
  ylab("Ratio of backward to forward rate")

pdf("Methylation_rate.pdf",height=3,width=8)
plot_grid(a,b,c,ncol=3,rel_widths = c(1,1,1.3))
dev.off()

#### LD ####
ceu_snps_decay <- read.table(file="ceu_snps_ld_decay.stat",header = T)
ceu_snps_decay$Type <- "SNP"
ceu_snps_decay$Class <- "All sites"
ceu_snps_cg_genes_decay <- read.table(file="ceu_snps_5mb_cg_genes_decay.stat",header = T)
ceu_snps_cg_genes_decay$Type <- "Coding SNP"
ceu_snps_cg_genes_decay$Class <- "Intervals"
ceu_snps_cg_intergenic_decay <- read.table(file="ceu_snps_5mb_cg_intergenic_decay.stat",header = T)
ceu_snps_cg_intergenic_decay$Type <- "Non-coding SNP"
ceu_snps_cg_intergenic_decay$Class <- "Intervals"
ceu_meth_cg_ld_decay.stat <- read.table(file="ceu_meth_cg_ld_decay.stat",header = T)
ceu_meth_cg_ld_decay.stat$Type <- "clocklike SMP"
ceu_meth_cg_ld_decay.stat$Class <- "Intervals"
ceu_meth_ld_decay.stat <- read.table(file="ceu_smps_wholegenome_ld_decay.stat",header = T)
ceu_meth_ld_decay.stat$Type <- "SMP"
ceu_meth_ld_decay.stat$Class <- "All sites"
ceu_meth_gbm_ld_decay.stat <- read.table(file="ceu_meth_gbm_ld_decay.stat",header = T)
ceu_meth_gbm_ld_decay.stat$Type <- "gbM SMP"
ceu_meth_gbm_ld_decay.stat$Class <- "Intervals"
ceu_ld_decay <- rbind(ceu_snps_decay,ceu_snps_cg_genes_decay,ceu_snps_cg_intergenic_decay,ceu_meth_ld_decay.stat,ceu_meth_cg_ld_decay.stat,ceu_meth_gbm_ld_decay.stat)
colnames(ceu_ld_decay) <- c("Dist","Mean_r2","Type","Class")
ceu_ld_decay$Pop <- "CEU"
ceu_ld_decay$Type <- factor(ceu_ld_decay$Type,levels=c("SNP","SMP","Coding SNP", "Non-coding SNP","clocklike SMP","gbM SMP"))
ceu_ld_decay_plot <- ggplot(data=ceu_ld_decay,aes(x=Dist,y=Mean_r2,color=Type)) +
  geom_line(linewidth=1) +
  facet_wrap(~Class) +
  xlab("Base pairs") +
  ylab(expression(paste("r"^"2"))) +
  scale_color_manual(labels = c("SNP","SMP","Coding SNP", "Non-coding SNP","clocklike SMP","gbM SMP"),values=c("purple","brown","black","grey","darkgreen","lightblue")) +
  ggtitle("CEU") +
  xlim(0,100) +
  theme_classic2() +
  theme(legend.position = "none")

for_ld_legend <- ggplot(data=ceu_ld_decay,aes(x=Dist,y=Mean_r2,color=Type)) +
  geom_line(linewidth=1) +
  facet_wrap(~Class) +
  xlab("Base pairs") +
  ylab(expression(paste("r"^"2"))) +
  scale_color_manual(labels = c("SNP","SMP","Coding SNP", "Non-coding SNP","clocklike SMP","gbM SMP"),values=c("purple","brown","black","grey","darkgreen","lightblue")) +
  ggtitle("CEU") +
  xlim(0,100) +
  theme_classic2() +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8))

ibnr_snps_decay <- read.table(file="ibnr_snps_ld_decay.stat",header = T)
ibnr_snps_decay$Type <- "SNP"
ibnr_snps_decay$Class <- "All sites"
ibnr_snps_cg_genes_decay <- read.table(file="ibnr_snps_5mb_cg_genes_decay.stat",header = T)
ibnr_snps_cg_genes_decay$Type <- "Coding SNP"
ibnr_snps_cg_genes_decay$Class <- "Intervals"
ibnr_snps_cg_intergenic_decay <- read.table(file="ibnr_snps_5mb_cg_intergenic_decay.stat",header = T)
ibnr_snps_cg_intergenic_decay$Type <- "Non-coding SNP"
ibnr_snps_cg_intergenic_decay$Class <- "Intervals"
ibnr_meth_cg_ld_decay.stat <- read.table(file="ibnr_meth_cg_ld_decay.stat",header = T)
ibnr_meth_cg_ld_decay.stat$Type <- "clocklike SMP"
ibnr_meth_cg_ld_decay.stat$Class <- "Intervals"
ibnr_meth_ld_decay.stat <- read.table(file="ibnr_smps_wholegenome_ld_decay.stat",header = T)
ibnr_meth_ld_decay.stat$Type <- "SMP"
ibnr_meth_ld_decay.stat$Class <- "All sites"
ibnr_meth_gbm_ld_decay.stat <- read.table(file="ibnr_meth_gbm_ld_decay.stat",header = T)
ibnr_meth_gbm_ld_decay.stat$Type <- "gbM SMP"
ibnr_meth_gbm_ld_decay.stat$Class <- "Intervals"
ibnr_ld_decay <- rbind(ibnr_snps_decay,ibnr_snps_cg_genes_decay,ibnr_snps_cg_intergenic_decay,ibnr_meth_ld_decay.stat,ibnr_meth_cg_ld_decay.stat,ibnr_meth_gbm_ld_decay.stat)
colnames(ibnr_ld_decay) <- c("Dist","Mean_r2","Type","Class")
ibnr_ld_decay$Pop <- "ibnr"
ibnr_ld_decay$Type <- factor(ibnr_ld_decay$Type,levels=c("SNP","SMP","Coding SNP", "Non-coding SNP","clocklike SMP","gbM SMP"))
ibnr_ld_decay_plot <- ggplot(data=ibnr_ld_decay,aes(x=Dist,y=Mean_r2,color=Type)) +
  geom_line(linewidth=1) +
  facet_wrap(~Class) +
  xlab("Base pairs") +
  ylab(expression(paste("r"^"2"))) +
  scale_color_manual(labels = c("SNP","SMP","Coding SNP", "Non-coding SNP","clocklike SMP","gbM SMP"),values=c("purple","brown","black","grey","darkgreen","lightblue")) +
  ggtitle("IBnr") +
  xlim(0,100) +
  theme_classic2() +
  theme(legend.position = "none")

pdf("LD_leg.pdf",height=2.5,width = 8.5)
for_ld_legend
dev.off()

#### IBD blocks  ####

ceu_ld.blocks <- read.table(file="ceu_snps_hmmIBD_output.hmm.txt",header=T)
ibnr_ld.blocks <- read.table(file="ibnr_snps_hmmIBD_output.hmm.txt",header=T)
ceu_ld.blocks <- ceu_ld.blocks[ceu_ld.blocks$different %in% "0",]
ibnr_ld.blocks <- ibnr_ld.blocks[ibnr_ld.blocks$different %in% "0",]
ceu_ld.blocks$Size <- ceu_ld.blocks$end - ceu_ld.blocks$start 
ibnr_ld.blocks$Size <- ibnr_ld.blocks$end - ibnr_ld.blocks$start 

ceu_ld <- ceu_ld.blocks[8]
ceu_ld$Population <- "CEU"
ceu_ld$Type <- "IBD"

ibnr_ld <- ibnr_ld.blocks[8]
ibnr_ld$Population <- "IBnr"
ibnr_ld$Type <- "IBD"

ceu_ld_mid <- ceu_ld.blocks[c(3:5,8)]
colnames(ceu_ld_mid)[1] <- "Chr"
ceu_ld_mid$start <- (ceu_ld_mid$start + ceu_ld_mid$end)
ceu_ld_mid <- ceu_ld_mid[-c(3)]
colnames(ceu_ld_mid)[2] <- "Pos"
ceu_ld_mid$Population <- "CEU"
ceu_ld_mid$Type <- "IBD"

ibnr_ld_mid <- ibnr_ld.blocks[c(3:5,8)]
colnames(ibnr_ld_mid)[1] <- "Chr"
ibnr_ld_mid$start <- (ibnr_ld_mid$start + ibnr_ld_mid$end)
ibnr_ld_mid <- ibnr_ld_mid[-c(3)]
colnames(ibnr_ld_mid)[2] <- "Pos"
ibnr_ld_mid$Population <- "IBnr"
ibnr_ld_mid$Type <- "IBD"

#### dmr size  ####
regions <- read.table(file="popgen5mb.bed")

dmrs_ceu <- read.table(file="10020_methylome_CG.txt",header=T)
dmrs_ceu <- dmrs_ceu[1:3]
dmrs_ceu$size <- dmrs_ceu$end - dmrs_ceu$start
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "1") & (dmrs_ceu$start < regions$V2[1])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "1") & (dmrs_ceu$end > regions$V3[1])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "2") & (dmrs_ceu$start < regions$V2[2])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "2") & (dmrs_ceu$end > regions$V3[2])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "3") & (dmrs_ceu$start < regions$V2[3])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "3") & (dmrs_ceu$end > regions$V3[3])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "4") & (dmrs_ceu$start < regions$V2[4])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "4") & (dmrs_ceu$end > regions$V3[4])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "5") & (dmrs_ceu$start < regions$V2[5])),]
dmrs_ceu <- dmrs_ceu[!((dmrs_ceu$seqnames %in% "5") & (dmrs_ceu$end > regions$V3[5])),]

ggplot(dmrs_ceu,aes(x=size)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))),bins=50)+
  xlab("Size of DMRs (bp)") +
  ylab("Proportion of DMRs") +
  ggtitle("CEU DMR blocks")

dmrs_ibnr <- read.table(file="9950_methylome_CG.txt",header=T)
dmrs_ibnr <- dmrs_ibnr[1:3]
dmrs_ibnr$size <- dmrs_ibnr$end - dmrs_ibnr$start

dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "1") & (dmrs_ibnr$start < regions$V2[1])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "1") & (dmrs_ibnr$end > regions$V3[1])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "2") & (dmrs_ibnr$start < regions$V2[2])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "2") & (dmrs_ibnr$end > regions$V3[2])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "3") & (dmrs_ibnr$start < regions$V2[3])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "3") & (dmrs_ibnr$end > regions$V3[3])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "4") & (dmrs_ibnr$start < regions$V2[4])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "4") & (dmrs_ibnr$end > regions$V3[4])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "5") & (dmrs_ibnr$start < regions$V2[5])),]
dmrs_ibnr <- dmrs_ibnr[!((dmrs_ibnr$seqnames %in% "5") & (dmrs_ibnr$end > regions$V3[5])),]

ggplot(dmrs_ibnr,aes(x=size)) +
  geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))),bins=50)+
  xlab("Size of DMRs (bp)") +
  ylab("Proportion of DMRs") +
  ggtitle("IBNR DMR blocks")

#### ld dmr compare  ####

ceu_dmr <- dmrs_ceu[4]
colnames(ceu_dmr) <- "Size"
ceu_dmr$Population <- "CEU"
ceu_dmr$Type <- "DMR"

ibnr_dmr <- dmrs_ibnr[4]
colnames(ibnr_dmr) <- "Size"
ibnr_dmr$Population <- "IBnr"
ibnr_dmr$Type <- "DMR"

ld_dmr <- rbind(ceu_ld,ibnr_ld,ceu_dmr,ibnr_dmr)
df_thaliana_ld_dmr <- dplyr::count(ld_dmr, Population, Type)
df_thaliana_ld_dmr$n <- paste("n=",df_thaliana_ld_dmr$n,sep="")

LD_DMR_block_plot <- ggplot(data=ld_dmr,aes(x=Population,y=log10(Size),fill=Type))+
  geom_boxplot() +
  theme(legend.title=element_blank()) +
  ylab(expression('Log'[10]*' base pairs'))+
  scale_fill_manual(values=c("purple","gold"))
  #geom_text(data = df_thaliana_ld_dmr,aes(y = log10(160000)+0.5,x=c(0.8,1.2,1.8,2.2),label = n))

ld_decay_plot <- plot_grid(ceu_ld_decay_plot, ibnr_ld_decay_plot)

#### tree  ####

library(ggtree)
library(treeio)
library(ape)
tip.group_1 <- read.csv('ceu_invcf_shuf',header=F,col.names = c("label", "Population"))
tip.group_2 <- read.csv('ibnr_invcf',header=F,col.names = c("label", "Population"))
tip.group_1$Population <- "CEU"
tip.group_2$Population <- "IBNR"
tip.group <- rbind(tip.group_1,tip.group_2)
tip.group$Color <- tip.group$Population
tip.group$Color <- gsub("CEU","blue",tip.group$Color)
tip.group$Color <- gsub("IBNR","orange",tip.group$Color)
tip.group_dmr <- tip.group
tip.group_dmr$label <- gsub("X","",tip.group_dmr$label)

snp_gene_tree<-read.tree("ceu_ibnr_snps_genes_mat_fastme-tree.nwk")
snp_gene_tree_pop <- full_join(snp_gene_tree,tip.group_dmr,by="label")
leg <- ggtree(snp_gene_tree_pop, 
       layout = "ape",
       size=0.5,
       aes(color=Population))+
  geom_tiplab(size=4,face="italic")+
  xlim(-0.03,0.08)+
  ylim(-0.03,0.08)+
  theme(legend.text = element_text(size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.01, y = -0.03, offset = 0.001, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                         breaks=c("CEU", "IBNR"),
                         labels=c("CEU", "IBnr")) +
  ggtitle("Coding")

leg <- get_legend(leg)

pdf(file="tree_leg.pdf",height=1,width=1)
as_ggplot(leg)
dev.off()

snp_gene_tree<-read.tree("ceu_ibnr_snps_genes_mat_fastme-tree.nwk")
snp_gene_tree_pop <- full_join(snp_gene_tree,tip.group_dmr,by="label")
snp_gene_tree_pop_plot <- ggtree(snp_gene_tree_pop, 
              layout = "ape",
              size=0.5,
              aes(color=Population))+
  geom_tiplab(size=4,face="italic")+
    xlim(-0.1,0.5)+
    ylim(-0.1,0.5)+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.01, y = -0.02, offset = 0.01, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBnr")) +
  #ggtitle("Genic")+
  guides(color = "none")

snp_intergenic_tree<-read.tree("ceu_ibnr_snps_intergenic_mat_fastme-tree.nwk")
snp_intergenic_pop_tree <- full_join(snp_intergenic_tree,tip.group_dmr,by="label")
snp_intergenic_pop_tree_plot <- ggtree(snp_intergenic_pop_tree, 
       layout = "ape",
       size=0.5,
       aes(color=Population))+
    geom_tiplab(size=4,face="italic")+
    xlim(-0.1,0.5)+
    ylim(-0.1,0.5)+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.01, y = -0.05, offset = 0.015, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBNR")) +
  #ggtitle("Intergenic")+
  guides(color = "none")

smp_cg_tree<-read.tree("ceu_ibnr_meth_cg_dis_mat_fastme-tree.nwk")
smp_cg_tree_pop <- full_join(smp_cg_tree,tip.group,by="label")
smp_cg_tree_pop_plot <- ggtree(smp_cg_tree_pop, 
       layout = "ape",
       size=0.5,
       aes(color=Population))+
    xlim(-0.3,0.7)+
    ylim(-0.3,0.7)+
    geom_tiplab(size=4,face="italic")+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = -10),
        plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.2, y = -0.2, offset = 0.025, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBNR"))+
  #ggtitle("SMP")+
  guides(color = "none")

dmr_tree<-read.tree("ceu_ibnr_dmrs_cg_mat_fastme-tree.nwk")
dmr_tree_pop <- full_join(dmr_tree,tip.group_dmr,by="label")
dmr_tree_pop_plot <- ggtree(dmr_tree_pop, 
       layout = "ape",
       size=0.5,
       aes(color=Population))+
    xlim(-0.2,0.6)+
    ylim(-0.2,0.6)+
    geom_tiplab(size=4,face="italic")+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.15, y = -0.1, offset = 0.02, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBNR"))+
  #ggtitle("DMR") +
  guides(color = "none")

pdf(file="trees2.pdf",height=5,width=6)
plot_grid(snp_intergenic_pop_tree_plot,snp_gene_tree_pop_plot,smp_cg_tree_pop_plot,dmr_tree_pop_plot,labels = c('Non-coding SNP, n=884064', 'Coding SNP, n=1022482','SMP, n=648376','DMR, n=2142'),hjust=-0.1)
dev.off()

smp_gbm_tree<-read.tree("ceu_ibnr_meth_gbm_dis_mat_fastme-tree.nwk")
smp_gbm_tree_pop <- full_join(smp_gbm_tree,tip.group,by="label")
smp_gbm_tree_pop_plot <- ggtree(smp_gbm_tree_pop, 
                                layout = "ape",
                                size=0.5,
                                aes(color=Population))+
    xlim(-0.3,0.7)+
    ylim(-0.3,0.7)+
    geom_tiplab(size=4,face="italic")+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = -10),
        plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.2, y = -0.2, offset = 0.025, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBNR"))+
  #ggtitle("gBM SMP")+
  guides(color = "none")



dmr_gbm_tree<-read.tree("ceu_ibnr_dmrs_gbm_mat_fastme-tree.nwk")
dmr_gbm_tree_pop <- full_join(dmr_gbm_tree,tip.group_dmr,by="label")
dmr_gbm_tree_pop_plot <- ggtree(dmr_gbm_tree_pop, 
                                layout = "ape",
                                size=0.5,
                                aes(color=Population))+
    xlim(-0.2,0.6)+
    ylim(-0.2,0.6)+
    geom_tiplab(size=4,face="italic")+
    theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "cm"))+
  guides(color = guide_legend(override.aes = list(label = "")))+
  geom_treescale(x = -0.15, y = -0.1, offset = 0.02, fontsize = 3) +
  scale_color_discrete(name  ="Population",
                       breaks=c("CEU", "IBNR"),
                       labels=c("CEU", "IBNR"))+
  #ggtitle("gBM DMR") +
  guides(color = "none")

pdf("trees_gBM.pdf",height = 3, width=6)
plot_grid(smp_gbm_tree_pop_plot,dmr_gbm_tree_pop_plot,labels = c('gbM SMP, n=346521','gbM DMR, n=1200'))
dev.off()

#### DAPC ####

library(ggplot2)
library(gtools)
library(adegenet)
library(dplyr)
library(radiator)
library(scales)

ceu_invcf <- read.table(file="ceu_invcf_shuf2")
colnames(ceu_invcf) <- "INDIVIDUALS"
ceu_invcf$STRATA <- "CEU"
ibnr_invcf <- read.table(file="ibnr_invcf2")
colnames(ibnr_invcf) <- "INDIVIDUALS"
ibnr_invcf$STRATA <- "IBnr"
stratafile <- rbind(ceu_invcf,ibnr_invcf)

smpfilesampleorder <- as.data.frame(c(gsub("^","X", ceu_invcf[order(as.character(ceu_invcf$INDIVIDUALS)),]$INDIVIDUALS),gsub("^","X", ibnr_invcf[order(as.character(ibnr_invcf$INDIVIDUALS)),]$INDIVIDUALS)))
colnames(smpfilesampleorder) <- "INDIVIDUALS"
stratafile2 <- stratafile
stratafile2$INDIVIDUALS <- gsub("^","X",stratafile2$INDIVIDUALS)
stratafile2 <- left_join(smpfilesampleorder,stratafile2,by="INDIVIDUALS")
smp_genlight <- genomic_converter(data = "ceu_ibnr_meth_cg_maf_thin.recode.vcf", strata = stratafile2,output = c("genlight"))
par(mfrow = c(1, 3))
smp_grp <- find.clusters(smp_genlight$genlight,max.n.clust=20)
60
2
smp_mat <- as.matrix(smp_genlight$genlight)
smp_genotype_matrix <- apply(smp_mat,2,na.replace,mean, na.rm = TRUE)
colnames(smp_genotype_matrix) <- seq(1:ncol(smp_genotype_matrix))
smp_xval <- xvalDapc(x=smp_genotype_matrix, grp=smp_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
dev.off()
smp_dapc <- dapc(smp_genlight$genlight, n.pca=5,n.da=2)

smpfilesampleorder <- as.data.frame(c("X10020","X403","X428","X430","X5874","X5907","X6008","X6296","X6396","X6445","X6903","X6951","X6956","X6957","X6976","X6984","X7103","X7120","X7177","X7203","X7207","X7296","X7350","X7424","X7520","X7521","X8284","X8386","X9665","X9668","X9671","X9676","X9679","X9684","X9689","X9690","X9693","X9696","X9727","X9728","X9729","X9730","X9731","X9732","X9733","X9735","X9973","X410","X424","X5890","X5893","X5921","X5950","X5993","X6390","X7067","X8236","X8285","X8290","X8365","X9694","X9914","X9915","X6933","X6970","X7328","X9507","X9510","X9511","X9514","X9515","X9521","X9522","X9524","X9525","X9534","X9535","X9537","X9540","X9541","X9544","X9547","X9556","X9557","X9560","X9562","X9564","X9567","X9568","X9577","X9582","X9588","X9594","X9821","X9822","X9825","X9833","X9834","X9836","X9840","X9841","X9843","X9844","X9845","X9852","X9856","X9857","X9859","X9861","X9864","X9867","X9868","X9870","X9873","X9876","X9878","X9886","X9888","X9899","X9900","X9902","X9904","X9946","X9950","X6971","X7327"))
colnames(smpfilesampleorder) <- "INDIVIDUALS"
stratafile2 <- stratafile
stratafile2$INDIVIDUALS <- gsub("^","X",stratafile2$INDIVIDUALS)
stratafile2 <- left_join(smpfilesampleorder,stratafile2,by="INDIVIDUALS")
smp_gbm_genlight <- genomic_converter(data = "ceu_ibnr_meth_gbm_thin.recode.vcf", strata = stratafile2,output = c("genlight"))
par(mfrow = c(1, 3))
smp_gbm_grp <- find.clusters(smp_gbm_genlight$genlight,max.n.clust=20)
60
2
smp_gbm_mat <- as.matrix(smp_gbm_genlight$genlight)
smp_gbm_genotype_matrix <- apply(smp_gbm_mat,2,na.replace,mean, na.rm = TRUE)
colnames(smp_gbm_genotype_matrix) <- seq(1:ncol(smp_gbm_genotype_matrix))
smp_gbm_xval <- xvalDapc(x=smp_gbm_genotype_matrix, grp=smp_gbm_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
dev.off()
smp_gbm_dapc <- dapc(smp_gbm_genlight$genlight, n.pca=5,n.da=2)

dmrfilesampleorder <-  as.data.frame(c("10020","403","410","424","428","430","5874","5890","5893","5907","5921","5950","5993","6008","6296","6390","6396","6445","6903","6951","6956","6957","6976","6984","7067","7103","7120","7177","7203","7207","7296","7350","7424","7520","7521","8236","8284","8285","8290","8365","8386","9665","9668","9671","9676","9679","9684","9689","9690","9693","9694","9696","9727","9728","9729","9730","9731","9732","9733","9735","9914","9915","9973","6933","6970","6971","7327","7328","9507","9510","9511","9514","9515","9521","9522","9524","9525","9534","9535","9537","9540","9541","9544","9547","9556","9557","9560","9562","9564","9567","9568","9577","9582","9588","9594","9821","9822","9825","9833","9834","9836","9840","9841","9843","9844","9845","9852","9856","9857","9859","9861","9864","9867","9868","9870","9873","9876","9878","9886","9888","9899","9900","9902","9904","9946","9950"))
colnames(dmrfilesampleorder) <- "INDIVIDUALS"
dmrfilesampleorder$INDIVIDUALS <- as.integer(dmrfilesampleorder$INDIVIDUALS)
stratafile3 <- left_join(dmrfilesampleorder,stratafile,by="INDIVIDUALS")
dmr_genlight <- genomic_converter(data = "ceu_ibnr_dmrs_cg_maf.recode.vcf", strata = stratafile3,output = c("genlight"))
par(mfrow = c(1, 3))
dmr_grp <- find.clusters(dmr_genlight$genlight,max.n.clust=20)
60
2
dmr_mat <- as.matrix(dmr_genlight$genlight)
dmr_genotype_matrix <- apply(dmr_mat,2,na.replace,mean, na.rm = TRUE)
colnames(dmr_genotype_matrix) <- seq(1:ncol(dmr_genotype_matrix))
dmr_xval <- xvalDapc(x=dmr_genotype_matrix, grp=dmr_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
dev.off()
dmr_dapc <- dapc(dmr_genlight$genlight, n.pca=5,n.da=2)

ceu_invcf <- read.table(file="ceu_invcf_shuf2")
colnames(ceu_invcf) <- "INDIVIDUALS"
ceu_invcf$STRATA <- "CEU"
ibnr_invcf <- read.table(file="ibnr_invcf2")
colnames(ibnr_invcf) <- "INDIVIDUALS"
ibnr_invcf$STRATA <- "IBnr"
stratafile <- rbind(ceu_invcf,ibnr_invcf)
dmr_gbm_genlight <- genomic_converter(data = "ceu_ibnr_dmrs_gbm.recode.vcf", strata = stratafile,output = c("genlight"))
par(mfrow = c(1, 3))
dmr_gbm_grp <- find.clusters(dmr_gbm_genlight$genlight,max.n.clust=20)
60
2
dmr_gbm_mat <- as.matrix(dmr_gbm_genlight$genlight)
dmr_gbm_genotype_matrix <- apply(dmr_gbm_mat,2,na.replace,mean, na.rm = TRUE)
colnames(dmr_gbm_genotype_matrix) <- seq(1:ncol(dmr_gbm_genotype_matrix))
dmr_gbm_xval <- xvalDapc(x=dmr_gbm_genotype_matrix, grp=dmr_gbm_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
dev.off()
dmr_gbm_dapc <- dapc(dmr_gbm_genlight$genlight, n.pca=5,n.da=2)

ceu_invcf <- read.table(file="ceu_invcf_shuf2")
colnames(ceu_invcf) <- "INDIVIDUALS"
ceu_invcf$STRATA <- "CEU"
ibnr_invcf <- read.table(file="ibnr_invcf2")
colnames(ibnr_invcf) <- "INDIVIDUALS"
ibnr_invcf$STRATA <- "IBnr"
stratafile <- rbind(ceu_invcf,ibnr_invcf)
stratafile <- stratafile[order(stratafile$STRATA,stratafile$INDIVIDUALS),]
genic_genlight <- genomic_converter(data = "ceu_ibnr_snps_genes_maf_thin.recode.vcf", strata = stratafile,output = c("genlight"))
genic_grp <- find.clusters(genic_genlight$genlight,max.n.clust=30) ## select 8 PC for grouping in contrast to methylation
genic_mat <- as.matrix(genic_genlight$genlight)
genic_genotype_matrix <- apply(genic_mat,2,na.replace,mean, na.rm = TRUE)
colnames(genic_genotype_matrix) <- seq(1:ncol(genic_genotype_matrix))
genic_xval <- xvalDapc(x=genic_genotype_matrix, grp=genic_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
genic_dapc <- dapc(genic_genlight$genlight, n.pca=5,n.da=2)

intergenic_genlight <- genomic_converter(data = "ceu_ibnr_snps_intergenic_maf_thin.recode.vcf", strata = stratafile,output = c("genlight"))
intergenic_grp <- find.clusters(intergenic_genlight$genlight,max.n.clust=30)
intergenic_mat <- as.matrix(intergenic_genlight$genlight)
intergenic_genotype_matrix <- apply(intergenic_mat,2,na.replace,mean, na.rm = TRUE)
colnames(intergenic_genotype_matrix) <- seq(1:ncol(intergenic_genotype_matrix))
intergenic_xval <- xvalDapc(x=intergenic_genotype_matrix, grp=intergenic_genlight$genlight$pop,n.pca.max=15,training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.rep = 30, xval.plot = TRUE)
intergenic_dapc <- dapc(intergenic_genlight$genlight, n.pca=5,n.da=2)

pdf("DAPC.pdf",height=5,width=5)
par(mfrow = c(2, 2))
scatter(intergenic_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE, col = hue_pal()(2),posi.leg ="top",leg=T, txt.leg=c("CEU","IBnr"))
title(main = "Non-coding SNP, n=4015")
scatter(genic_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE,  col = hue_pal()(2))
title(main = "Coding SNP, n=5077")
scatter(smp_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE, col=hue_pal()(2))
title(main = "SMP, n=10731")
scatter(dmr_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE,  col = hue_pal()(2))
title(main = "DMR, n=1157")
dev.off()


pdf("DAPC_gBM.pdf",height=3.5,width=7)
par(mfrow = c(1, 2))
scatter(smp_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE, posi.leg ="topleft",leg=T, txt.leg=c("CEU","IBnr"),col=hue_pal()(2))
title(main = "gbM SMP, n=5639")
scatter(dmr_gbm_dapc, ratio.pca=0.2, bg="white", pch=20, cell=0,cstar=0, solid=.5, cex=5, clab=0,mstree=TRUE, col = hue_pal()(2))
title(main = "gbM DMR, n=664")
dev.off()


#### TE frequencies ####

ceu_te.frq <- read.table(file="ceu_te.frq",row.names = NULL, header = T)
ceu_te.frq$N_CHR <- as.numeric(gsub("A:","",ceu_te.frq$N_CHR))
ceu_te.frq <- ceu_te.frq[(ceu_te.frq$N_CHR > 0.02) & (ceu_te.frq$N_CHR < 0.98),]
ceu_te.frq$Population <- "CEU"
ibnr_te.frq <- read.table(file="ibnr_te.frq",row.names = NULL, header = T)
ibnr_te.frq$N_CHR <- as.numeric(gsub("A:","",ibnr_te.frq$N_CHR))
ibnr_te.frq <- ibnr_te.frq[(ibnr_te.frq$N_CHR > 0.02) & (ibnr_te.frq$N_CHR < 0.98),]
ibnr_te.frq$Population <- "IBnr"
te.frq <- rbind(ceu_te.frq,ibnr_te.frq)
te.frq <- te.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0.02, 0.98, by = 0.04), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Population, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Population) %>%
  dplyr::mutate(proportion = count / sum(count))
te.frq$N_CHR_binned <- c(seq(0.04,0.96,0.04),seq(0.04,0.96,0.04))
te.frq.plot <- ggplot(te.frq, aes(x=N_CHR_binned, y=proportion, fill=Population)) +
  geom_bar(stat="identity", position="dodge") +
  labs(y="Proportion of TEs", x="TE absence frequency") +
  theme_minimal()+
  scale_fill_discrete(labels = c(nrow(ceu_te.frq),nrow(ibnr_te.frq))) +
  ggtitle("Transposable elements (TE)") +
  theme(legend.title = element_blank(),legend.position = c(0.5, 0.8))


#### STR analysis   ####

library(GenomicRanges)
library(dplyr)

### get STR data
str <- read.table(file="210217.SuppDataSet2.DiploidUnitNumberCalls.tsv",sep="\t",row.names = 1,header=T)

### only keep STRs that occur in intergenic regions
positions <- colnames(str)
position_data <- data.frame(chr = sub("_.*", "", positions), pos = as.numeric(sub(".*_", "", positions)))
position_data$chr <- sub("chr", "", position_data$chr)
bed_data <- read.table("non_gene_at.bed", col.names = c("chr", "start", "end"))
query_ranges <- GRanges(seqnames = position_data$chr, ranges = IRanges(start = position_data$pos, end = position_data$pos))
bed_ranges <- GRanges(seqnames = as.character(bed_data$chr),  ranges = IRanges(start = bed_data$start, end = bed_data$end))
overlaps <- findOverlaps(query_ranges, bed_ranges)
matched_positions <- positions[queryHits(overlaps)]
str <- str[colnames(str) %in% matched_positions]

### extract polymorphic STRs for CEU
ceusamples <- read.table(file="ceu_invcf_shuf2")
ceu_str <- str[rownames(str) %in% ceusamples$V1,]
ceu_str <- ceu_str[sapply(apply(ceu_str,2,table),length) > 1]

### make all values 2 digits
ceu_str <- as.data.frame(lapply(ceu_str, function(x) sprintf("%02d", x)))

#  function to create diploid genotypes by combining samples
combine_random_rows <- function(df) {
  # Randomly shuffle the row indices
  shuffled_indices <- sample(nrow(df))
  
  # Rearrange the dataframe with shuffled rows
  df <- df[shuffled_indices, , drop = FALSE]
  
  # Create an empty data frame with the same column names
  combined_df <- data.frame(matrix(ncol = ncol(df), nrow = nrow(df) / 2))
  colnames(combined_df) <- colnames(df)
  
  # Loop over each column
  for (i in 1:ncol(df)) {
    # Combine every two randomly paired rows with a '/' for each column
    combined_df[, i] <- sapply(seq(1, nrow(df), by = 2), function(j) {
      return(paste(df[j, i], df[j+1, i], sep = ""))
    })
  }
  
  # Replace all cells that contain "NA" (even partial matches) with actual NA
  combined_df[] <- lapply(combined_df, function(x) {
    x[grepl("NA", x)] <- NA_character_
    return(x)
  })
  
  # Return the combined data frame
  return(combined_df)
}
ceu_str_dip <- combine_random_rows(ceu_str)
ceu_str_dip[is.na(ceu_str_dip)] <- "0000"

# Function to convert haploid data to Genepop format
convert_to_genepop <- function(diploid_data, output_file, title = "Haploid dataset converted for Genepop") {
  # Extract locus names
  loci_names <- colnames(diploid_data)
  sample_ids <- rownames(diploid_data)
  # Open file for writing
  fileConn <- file(output_file, "w")  # "w" ensures a fresh write
  writeLines(title, fileConn)  # Title line
  close(fileConn)
  # Append locus names
  for (l in loci_names) {
    cat(l, "\n", file = output_file, append = TRUE)
  }
  # Append first population marker
  cat("Pop\n", file = output_file, append = TRUE)
  # Append each individual's genotype data
  for (i in seq_along(sample_ids)) {
    line <- paste0(sample_ids[i], " , ", paste(diploid_data[i, ], collapse = "\t"))
    cat(line, "\n", file = output_file, append = TRUE)
  }
  message("Genepop file saved as: ", output_file)
}

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep1_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep2_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep3_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep4_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep5_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep6_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep7_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep8_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep9_100.txt")

ceu_str_100 <- ceu_str_dip[sample(ncol(ceu_str_dip), 100, replace = FALSE)]
convert_to_genepop(ceu_str_100, "ceu_str_rep10_100.txt")

### extract polymorphic STRs for IBnr
ibnrsamples <- read.table(file="ibnr_invcf2")
ibnr_str <- str[rownames(str) %in% ibnrsamples$V1,]
ibnr_str <- ibnr_str[sapply(apply(ibnr_str,2,table),length) > 1]

### make all values 2 digits
ibnr_str <- as.data.frame(lapply(ibnr_str, function(x) sprintf("%02d", x)))

#  function to create diploid genotypes by combining samples
combine_random_rows <- function(df) {
  # Randomly shuffle the row indices
  shuffled_indices <- sample(nrow(df))
  
  # Rearrange the dataframe with shuffled rows
  df <- df[shuffled_indices, , drop = FALSE]
  
  # Create an empty data frame with the same column names
  combined_df <- data.frame(matrix(ncol = ncol(df), nrow = nrow(df) / 2))
  colnames(combined_df) <- colnames(df)
  
  # Loop over each column
  for (i in 1:ncol(df)) {
    # Combine every two randomly paired rows with a '/' for each column
    combined_df[, i] <- sapply(seq(1, nrow(df), by = 2), function(j) {
      return(paste(df[j, i], df[j+1, i], sep = ""))
    })
  }
  
  # Replace all cells that contain "NA" (even partial matches) with actual NA
  combined_df[] <- lapply(combined_df, function(x) {
    x[grepl("NA", x)] <- NA_character_
    return(x)
  })
  
  # Return the combined data frame
  return(combined_df)
}
ibnr_str_dip <- combine_random_rows(ibnr_str)
ibnr_str_dip[is.na(ibnr_str_dip)] <- "0000"

# Function to convert haploid data to Genepop format
convert_to_genepop <- function(diploid_data, output_file, title = "Haploid dataset converted for Genepop") {
  # Extract locus names
  loci_names <- colnames(diploid_data)
  sample_ids <- rownames(diploid_data)
  # Open file for writing
  fileConn <- file(output_file, "w")  # "w" ensures a fresh write
  writeLines(title, fileConn)  # Title line
  close(fileConn)
  # Append locus names
  for (l in loci_names) {
    cat(l, "\n", file = output_file, append = TRUE)
  }
  # Append first population marker
  cat("Pop\n", file = output_file, append = TRUE)
  # Append each individual's genotype data
  for (i in seq_along(sample_ids)) {
    line <- paste0(sample_ids[i], " , ", paste(diploid_data[i, ], collapse = "\t"))
    cat(line, "\n", file = output_file, append = TRUE)
  }
  message("Genepop file saved as: ", output_file)
}

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep1_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep2_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep3_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep4_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep5_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep6_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep7_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep8_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep9_100.txt")

ibnr_str_100 <- ibnr_str_dip[sample(ncol(ibnr_str_dip), 100, replace = FALSE)]
convert_to_genepop(ibnr_str_100, "ibnr_str_rep10_100.txt")
