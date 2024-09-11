library(crayon)
library(dartR)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(adegenet)
library(gtools)
library(vcfR)
library(dplyr)
library(cowplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggrepel)
library(agricolae)

####  variant site prop (MAF > 0.04) #### 

group1.1_smp_prop <- 124137/10062950
group1.2_smp_prop <- 76079/9640492
group6_smp_prop <- 94712/8730867
group1.1_smp_gbm_prop <- 60113/905826
group1.2_smp_gbm_prop <- 34609/879235
group6_smp_gbm_prop <- 38260/835966
gbm_genes_bd <- read.table(file="gbm_genes_bd.bed")
1154347/sum(gbm_genes_bd$V3-gbm_genes_bd$V2) ## to get SMP genotype density
1154347/sum(gbm_genes_bd$V3-gbm_genes_bd$V2) ## to get SMP genotype density

propsites_bd <- as.data.frame(rbind(group1.1_smp_prop,group1.2_smp_prop,group6_smp_prop,group1.1_smp_gbm_prop,group1.2_smp_gbm_prop,group6_smp_gbm_prop))
colnames(propsites_bd) <- "Proportion"
propsites_bd$Group <- c("1.1","1.2","6","1.1","1.2","6")
propsites_bd$Type <- c("All sites","All sites","All sites","gbM sites","gbM sites","gbM sites")
total_bd <- c("n=10062950","n=9640492","n=8730867","n=905826","n=879235","n=835966")

pdf(file="seg_sites_bd.pdf",height=3,width=3.5)
ggplot(data=propsites_bd,aes(x=Type,y=Proportion,fill=Group)) +
  geom_col(position = "dodge") +
  ylim(c(0,0.1)) +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  geom_text(aes(x=c(0.7,1,1.3,1.7,2,2.3),angle=90,y = 0.085, label = total_bd)) +
  ylab("Proportion of polymorphic sites") +
  xlab("")
dev.off()

#### DAPC #### 

smp_vcf <- read.vcfR("smp29_gbm_maf04_thin.recode.vcf")
sample_data <- data.frame(sample_names = c('BdTR11i','BdTR10n','BdTR11c','BdTR11d','BdTR11f','BdTR11g','BdTR11h','BdTR1a',
                                           'BdTR1b','BdTR1e','BdTR1f','BdTR1g','BdTR1h','BdTR1j','BdTR1k','BdTR1m',
                                           'BdTR1n','BdTR2b','BdTR2c','BdTR2d','BdTR2g','BdTR2h','BdTR2j','BdTR2k',
                                           'BdTR2m','BdTR2n','BdTR2p','BdTR2r','BdTR2s'))

group_data <- data.frame(sample_names = c('BdTR11i','BdTR10n','BdTR11c','BdTR11d','BdTR11f','BdTR11g','BdTR11h','BdTR1a',
                                          'BdTR1b','BdTR1e','BdTR1f','BdTR1g','BdTR1h','BdTR1j','BdTR1k','BdTR1m',
                                          'BdTR1n','BdTR2b','BdTR2c','BdTR2d','BdTR2g','BdTR2h','BdTR2j','BdTR2k',
                                          'BdTR2m','BdTR2n','BdTR2p','BdTR2r','BdTR2s'),
                         Group = c("group_6","group_6","group_6","group_6","group_6","group_6","group_6",
                                   "group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1",
                                   "group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2"))
joined_data <- full_join(sample_data, group_data, by = "sample_names")
group_names <- c("group_1.1","group_1.2","group_6")
group_colors <- c( "#FDAE61","#ABDDA4","#2B83BA")
myCols <- setNames(group_colors,group_names)
smp_genlight <- vcfR2genlight(smp_vcf)
smp_grp <- find.clusters(smp_genlight,max.n.clust=20)
smp_1 <- gl.define.pop(smp_genlight, ind.list=c("BdTR1a","BdTR1b","BdTR1e","BdTR1f","BdTR1g","BdTR1h","BdTR1j","BdTR1k","BdTR1m","BdTR1n"), new="group_1.1", verbose = NULL)
smp_2 <- gl.define.pop(smp_1, ind.list=c("BdTR2b","BdTR2c","BdTR2d","BdTR2g","BdTR2h","BdTR2j","BdTR2k","BdTR2m","BdTR2n","BdTR2p","BdTR2r","BdTR2s"), new="group_1.2", verbose = NULL)
smp <- gl.define.pop(smp_2, ind.list=c("BdTR11i","BdTR10n","BdTR11c","BdTR11d","BdTR11f","BdTR11g","BdTR11h"), new="group_6", verbose = NULL)
smp_nafill<-apply(smp,2,na.replace,mean, na.rm = TRUE)
mat <- as.matrix(smp)
genotype_matrix <- apply(mat,2,na.replace,mean, na.rm = TRUE)

xval <- xvalDapc(genotype_matrix, smp$pop, n.pca.max = 15, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)
xval
dapc <- dapc(smp, n.da=2, n.pca=2)

pdf("DAPC_SMP_PC2_Bd.pdf", width = 2.7, height = 2.7)
par(cex=0.6)
scatter(dapc, bg = "white", pch = 20, cstar = 0, col = myCols, solid = 0.4, cex = 4, clab = 0,
        mstree = F, scree.da = FALSE, posi.leg = "topright",
        leg = TRUE, txt.leg = group_names)
title("SMP",font.main=2)
par(xpd=TRUE,new=TRUE)
points(dapc$grp.coord[,1], 	dapc$grp.coord[,2], pch=4,
       cex=2, lwd=2, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,2], pch=4,
       cex=1, lwd=2, col=myCols)
myInset <- function(){
  temp <- dapc$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated\n variance (%)",cex.lab=0.8,
       cex=1, pch=20, type="h", lwd=2)}
add.scatter(myInset(), posi="bottomleft",
            inset=c(0.35,0.0001), ratio=0.2,
            bg=transp("white"))
dev.off()

dmr_vcf <- read.vcfR("dmr29_gbm_maf04.vcf")
sample_data <- data.frame(sample_names = c('BdTR11i','BdTR10n','BdTR11c','BdTR11d','BdTR11f','BdTR11g','BdTR11h','BdTR1a',
                                           'BdTR1b','BdTR1e','BdTR1f','BdTR1g','BdTR1h','BdTR1j','BdTR1k','BdTR1m',
                                           'BdTR1n','BdTR2b','BdTR2c','BdTR2d','BdTR2g','BdTR2h','BdTR2j','BdTR2k',
                                           'BdTR2m','BdTR2n','BdTR2p','BdTR2r','BdTR2s'))

group_data <- data.frame(sample_names = c('BdTR11i','BdTR10n','BdTR11c','BdTR11d','BdTR11f','BdTR11g','BdTR11h','BdTR1a',
                                          'BdTR1b','BdTR1e','BdTR1f','BdTR1g','BdTR1h','BdTR1j','BdTR1k','BdTR1m',
                                          'BdTR1n','BdTR2b','BdTR2c','BdTR2d','BdTR2g','BdTR2h','BdTR2j','BdTR2k',
                                          'BdTR2m','BdTR2n','BdTR2p','BdTR2r','BdTR2s'),
                         Group = c("group_6","group_6","group_6","group_6","group_6","group_6","group_6",
                                   "group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1","group_1.1",
                                   "group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2","group_1.2"))
joined_data <- full_join(sample_data, group_data, by = "sample_names")
group_names <- c("group_1.1","group_1.2","group_6")
group_colors <- c( "#FDAE61","#ABDDA4","#2B83BA")
myCols <- setNames(group_colors,group_names)
dmr_genlight <- vcfR2genlight(dmr_vcf)
dmr_grp <- find.clusters(dmr_genlight,max.n.clust=20)
dmr_1 <- gl.define.pop(dmr_genlight, ind.list=c("BdTR1a","BdTR1b","BdTR1e","BdTR1f","BdTR1g","BdTR1h","BdTR1j","BdTR1k","BdTR1m","BdTR1n"), new="group_1.1", verbose = NULL)
dmr_2 <- gl.define.pop(dmr_1, ind.list=c("BdTR2b","BdTR2c","BdTR2d","BdTR2g","BdTR2h","BdTR2j","BdTR2k","BdTR2m","BdTR2n","BdTR2p","BdTR2r","BdTR2s"), new="group_1.2", verbose = NULL)
dmr <- gl.define.pop(dmr_2, ind.list=c("BdTR11i","BdTR10n","BdTR11c","BdTR11d","BdTR11f","BdTR11g","BdTR11h"), new="group_6", verbose = NULL)
dmr_nafill<-apply(dmr,2,na.replace,mean, na.rm = TRUE)
mat <- as.matrix(dmr)
genotype_matrix <- apply(mat,2,na.replace,mean, na.rm = TRUE)

xval <- xvalDapc(genotype_matrix, dmr$pop, n.pca.max = 15, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)
xval
dapc <- dapc(dmr, n.da=2, n.pca=2)

pdf("DAPC_DMR_PC2_Bd.pdf", width = 2.7, height = 2.7)
par(cex=0.6)
scatter(dapc, bg = "white", pch = 20, cstar = 0, col = myCols, solid = 0.4, cex = 4, clab = 0,
        mstree = F, scree.da = FALSE, posi.leg = "topright",
        leg = TRUE, txt.leg = group_names)
title("DMR",font.main=2)
par(xpd=TRUE,new=TRUE)
points(dapc$grp.coord[,1], 	dapc$grp.coord[,2], pch=4,
       cex=2, lwd=2, col="black")
points(dapc$grp.coord[,1], dapc$grp.coord[,2], pch=4,
       cex=1, lwd=2, col=myCols)
myInset <- function(){
  temp <- dapc$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated\n variance (%)",cex.lab=0.8,
       cex=1, pch=20, type="h", lwd=2)}
add.scatter(myInset(), posi="bottomleft",
            inset=c(0.35,0.0001), ratio=0.2,
            bg=transp("white"))
dev.off()

#### SFS #### 

## gbm regions

group1.1_smps.frq <- read.table(file="smp1.1_gbm_maf10.frq",row.names = NULL, header = T)
group1.1_smps.frq$N_CHR <- as.numeric(gsub("A:","",group1.1_smps.frq$N_CHR))
group1.1_smps.frq <- group1.1_smps.frq[(group1.1_smps.frq$N_CHR > 0.02) & (group1.1_smps.frq$N_CHR < 0.98),]
group1.1_smps.frq$Group <- "1.1"
group1.1_smps.frq$Type <- "SMP"
group1.1_dmr.frq <- read.table(file="dmr1.1_gbm_maf10.frq",row.names = NULL, header = T)
group1.1_dmr.frq$N_CHR <- as.numeric(gsub("A:","",group1.1_dmr.frq$N_CHR))
group1.1_dmr.frq <- group1.1_dmr.frq[(group1.1_dmr.frq$N_CHR > 0.02) & (group1.1_dmr.frq$N_CHR < 0.98),]
group1.1_dmr.frq$Group <- "1.1"
group1.1_dmr.frq$Type <- "DMR"
group1.2_smps.frq <- read.table(file="smp1.2_gbm_maf10.frq",row.names = NULL, header = T)
group1.2_smps.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_smps.frq$N_CHR))
group1.2_smps.frq <- group1.2_smps.frq[(group1.2_smps.frq$N_CHR > 0.02) & (group1.2_smps.frq$N_CHR < 0.98),]
group1.2_smps.frq$Group <- "1.2"
group1.2_smps.frq$Type <- "SMP"
group1.2_dmr.frq <- read.table(file="dmr1.2_gbm_maf10.frq",row.names = NULL, header = T)
group1.2_dmr.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_dmr.frq$N_CHR))
group1.2_dmr.frq <- group1.2_dmr.frq[(group1.2_dmr.frq$N_CHR > 0.02) & (group1.2_dmr.frq$N_CHR < 0.98),]
group1.2_dmr.frq$Group <- "1.2"
group1.2_dmr.frq$Type <- "DMR"
group6_smps.frq <- read.table(file="smp6_gbm_maf04.frq",row.names = NULL, header = T)
group6_smps.frq$N_CHR <- as.numeric(gsub("A:","",group6_smps.frq$N_CHR))
group6_smps.frq <- group6_smps.frq[(group6_smps.frq$N_CHR > 0.02) & (group6_smps.frq$N_CHR < 0.98),]
group6_smps.frq$Group <- "6"
group6_smps.frq$Type <- "SMP"
group6_dmr.frq <- read.table(file="dmr6_gbm_maf04.frq",row.names = NULL, header = T)
group6_dmr.frq$N_CHR <- as.numeric(gsub("A:","",group6_dmr.frq$N_CHR))
group6_dmr.frq <- group6_dmr.frq[(group6_dmr.frq$N_CHR > 0.02) & (group6_dmr.frq$N_CHR < 0.98),]
group6_dmr.frq$Group <- "6"
group6_dmr.frq$Type <- "DMR"
bd_frq <- rbind(group1.1_smps.frq,group1.2_smps.frq,group6_smps.frq,group1.1_dmr.frq,group1.2_dmr.frq,group6_dmr.frq)
bd_frq <- bd_frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0, 1, by = 0.25), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Group,Type, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Group,Type) %>%
  dplyr::mutate(proportion = count / sum(count))
bd_frq$N_CHR_binned <- rep(seq(0.125,0.975,0.25),6)
bd_frq_smp <- bd_frq[bd_frq$Type %in% "SMP",]
bd_frq_smp$Group <- gsub("6",paste("6, n=",nrow(group6_smps.frq),sep=""), bd_frq_smp$Group)
bd_frq_smp$Group <- gsub("1.1",paste("1.1, n=",nrow(group1.1_smps.frq),sep=""), bd_frq_smp$Group)
bd_frq_smp$Group <- gsub("1.2",paste("1.2, n=",nrow(group1.2_smps.frq),sep=""), bd_frq_smp$Group)
bd_frq_dmr <- bd_frq[bd_frq$Type %in% "DMR",]
bd_frq_dmr$Group <- gsub("6",paste("6, n=",nrow(group6_dmr.frq),sep=""), bd_frq_dmr$Group)
bd_frq_dmr$Group <- gsub("1.1",paste("1.1, n=",nrow(group1.1_dmr.frq),sep=""), bd_frq_dmr$Group)
bd_frq_dmr$Group <- gsub("1.2",paste("1.2, n=",nrow(group1.2_dmr.frq),sep=""), bd_frq_dmr$Group)
bd_smp_plot <- ggplot(bd_frq_smp, aes(x=N_CHR_binned, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Proportion of unmethylated alleles", y="Proportion of sites") +
  theme_minimal()+
  ggtitle("SMP") +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  theme(legend.title = element_blank(),legend.position = c(0.8, 0.95))
bd_dmr_plot <- ggplot(bd_frq_dmr, aes(x=N_CHR_binned, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Proportion of unmethylated alleles", y="Proportion of sites") +
  theme_minimal()+
  ggtitle("DMR") +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  theme(legend.title = element_blank(),legend.position = c(0.8, 0.95))

pdf(file="brachypodium_mSFS.pdf",width=6,height=3.5)
plot_grid(bd_smp_plot,bd_dmr_plot,ncol=2)
dev.off()

## non gbm regions
group1.1_smps_nongbm.frq <- read.table(file="smp1.1_exon_nongbm_maf10.frq",row.names = NULL, header = T)
group1.1_smps_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group1.1_smps_nongbm.frq$N_CHR))
group1.1_smps_nongbm.frq <- group1.1_smps_nongbm.frq[(group1.1_smps_nongbm.frq$N_CHR > 0.02) & (group1.1_smps_nongbm.frq$N_CHR < 0.98),]
group1.1_smps_nongbm.frq$Group <- "1.1"
group1.1_smps_nongbm.frq$Type <- "SMP"
group1.1_dmr_nongbm.frq <- read.table(file="dmr1.1_exon_nongbm_maf10.frq",row.names = NULL, header = T)
group1.1_dmr_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group1.1_dmr_nongbm.frq$N_CHR))
group1.1_dmr_nongbm.frq <- group1.1_dmr_nongbm.frq[(group1.1_dmr_nongbm.frq$N_CHR > 0.02) & (group1.1_dmr_nongbm.frq$N_CHR < 0.98),]
group1.1_dmr_nongbm.frq$Group <- "1.1"
group1.1_dmr_nongbm.frq$Type <- "DMR"
group1.2_smps_nongbm.frq <- read.table(file="smp1.2_exon_nongbm_maf10.frq",row.names = NULL, header = T)
group1.2_smps_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_smps_nongbm.frq$N_CHR))
group1.2_smps_nongbm.frq <- group1.2_smps_nongbm.frq[(group1.2_smps_nongbm.frq$N_CHR > 0.02) & (group1.2_smps_nongbm.frq$N_CHR < 0.98),]
group1.2_smps_nongbm.frq$Group <- "1.2"
group1.2_smps_nongbm.frq$Type <- "SMP"
group1.2_dmr_nongbm.frq <- read.table(file="dmr1.2_exon_nongbm_maf10.frq",row.names = NULL, header = T)
group1.2_dmr_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_dmr_nongbm.frq$N_CHR))
group1.2_dmr_nongbm.frq <- group1.2_dmr_nongbm.frq[(group1.2_dmr_nongbm.frq$N_CHR > 0.02) & (group1.2_dmr_nongbm.frq$N_CHR < 0.98),]
group1.2_dmr_nongbm.frq$Group <- "1.2"
group1.2_dmr_nongbm.frq$Type <- "DMR"
group6_smps_nongbm.frq <- read.table(file="smp6_exon_nongbm_maf04.frq",row.names = NULL, header = T)
group6_smps_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group6_smps_nongbm.frq$N_CHR))
group6_smps_nongbm.frq <- group6_smps_nongbm.frq[(group6_smps_nongbm.frq$N_CHR > 0.02) & (group6_smps_nongbm.frq$N_CHR < 0.98),]
group6_smps_nongbm.frq$Group <- "6"
group6_smps_nongbm.frq$Type <- "SMP"
group6_dmr_nongbm.frq <- read.table(file="dmr6_exon_nongbm_maf04.frq",row.names = NULL, header = T)
group6_dmr_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group6_dmr_nongbm.frq$N_CHR))
group6_dmr_nongbm.frq <- group6_dmr_nongbm.frq[(group6_dmr_nongbm.frq$N_CHR > 0.02) & (group6_dmr_nongbm.frq$N_CHR < 0.98),]
group6_dmr_nongbm.frq$Group <- "6"
group6_dmr_nongbm.frq$Type <- "DMR"
bd_nongbm_frq <- rbind(group1.1_smps_nongbm.frq,group1.2_smps_nongbm.frq,group6_smps_nongbm.frq,group1.1_dmr_nongbm.frq,group1.2_dmr_nongbm.frq,group6_dmr_nongbm.frq)
bd_nongbm_frq <- bd_nongbm_frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0, 1, by = 0.25), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Group,Type, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Group,Type) %>%
  dplyr::mutate(proportion = count / sum(count))
bd_nongbm_frq$N_CHR_binned <- rep(seq(0.125,0.975,0.25),6)
bd_nongbm_frq_smp <- bd_nongbm_frq[bd_nongbm_frq$Type %in% "SMP",]
bd_nongbm_frq_smp$Group <- gsub("6",paste("6, n=",nrow(group6_smps_nongbm.frq),sep=""), bd_nongbm_frq_smp$Group)
bd_nongbm_frq_smp$Group <- gsub("1.1",paste("1.1, n=",nrow(group1.1_smps_nongbm.frq),sep=""), bd_nongbm_frq_smp$Group)
bd_nongbm_frq_smp$Group <- gsub("1.2",paste("1.2, n=",nrow(group1.2_smps_nongbm.frq),sep=""), bd_nongbm_frq_smp$Group)
bd_nongbm_frq_dmr <- bd_nongbm_frq[bd_nongbm_frq$Type %in% "DMR",]
bd_nongbm_frq_dmr$Group <- gsub("6",paste("6, n=",nrow(group6_dmr_nongbm.frq),sep=""), bd_nongbm_frq_dmr$Group)
bd_nongbm_frq_dmr$Group <- gsub("1.1",paste("1.1, n=",nrow(group1.1_dmr_nongbm.frq),sep=""), bd_nongbm_frq_dmr$Group)
bd_nongbm_frq_dmr$Group <- gsub("1.2",paste("1.2, n=",nrow(group1.2_dmr_nongbm.frq),sep=""), bd_nongbm_frq_dmr$Group)
bd_smp_nongbm_plot <- ggplot(bd_nongbm_frq_smp, aes(x=N_CHR_binned, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Proportion of unmethylated alleles", y="Proportion of sites") +
  theme_minimal()+
  ggtitle("SMP") +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  theme(legend.title = element_blank(),legend.position = c(0.8, 0.95))
bd_dmr_nongbm_plot <- ggplot(bd_nongbm_frq_dmr, aes(x=N_CHR_binned, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Proportion of unmethylated alleles", y="Proportion of sites") +
  theme_minimal()+
  ggtitle("DMR") +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  theme(legend.title = element_blank(),legend.position = c(0.4, 0.95))

pdf(file="brachypodium_mSFS_nongbm.pdf",width=6,height=3.5)
plot_grid(bd_smp_nongbm_plot,bd_dmr_nongbm_plot,ncol=2)
dev.off()

#### LD ####

group1.1_meth_ld_decay <- read.table(file="smp1.1_gbm_maf10_ld_decay.stat")
group1.1_meth_ld_decay$Group <- "Group 1.1"
group1.2_meth_ld_decay <- read.table(file="smp1.2_gbm_maf10_ld_decay.stat")
group1.2_meth_ld_decay$Group <- "Group 1.2"
group6_meth_ld_decay <- read.table(file="smp6_gbm_maf04_ld_decay.stat")
group6_meth_ld_decay$Group <- "Group 6"
bd_ld_decay <- rbind(group1.1_meth_ld_decay,group1.2_meth_ld_decay,group6_meth_ld_decay)
colnames(bd_ld_decay) <- c("Dist","Mean_r2","Mean_D","Sum_r", "Sum_D","NumberPairs","Group")
bd_ld_decay_plot <- ggplot(data=bd_ld_decay,aes(x=Dist,y=Mean_r2,color=Group)) +
  geom_line() +
  xlab("Base pairs") +
  ylab(expression(paste("r"^"2"))) +
  scale_color_manual(labels = c("Group 1.1, n=97453","Group 1.2, n=34609","Group 6, n=38260"),values=c( "#FDAE61","#ABDDA4","#2B83BA")) +
  xlim(0,300) +
  ylim(0,0.6) +
  theme_classic2() +
  theme(legend.position = c(0.55,0.84),legend.title = element_blank(),legend.text = element_text(size = 8))

pdf("LD_brachypodium.pdf",height=2.5,width = 3.5)
bd_ld_decay_plot
dev.off()


#### DMR sizes ####

## brachypodium

bd_dmrs_1 <- read.table(file="dmr_length_grp1.txt",header=T)
bd_dmrs_6 <- read.table(file="dmr_length_grp6.txt",header=T)
bd_dmrs <- rbind(bd_dmrs_1,bd_dmrs_6)
bd_dmrs <- bd_dmrs[!duplicated(bd_dmrs),]
group1.1_dmr <- read.table(file="dmr1.1_gbm_maf10.frq",row.names = NULL, header = T)
group1.2_dmr <- read.table(file="dmr1.2_gbm_maf10.frq",row.names = NULL, header = T)
group6_dmr <- read.table(file="dmr6_gbm_maf10.frq",row.names = NULL, header = T)
dmrs_all <- rbind(group1.1_dmr,group1.2_dmr,group6_dmr)
bd_dmrs$ID <- paste(bd_dmrs$seqnames,bd_dmrs$middle,sep="_")
bd_dmrs <- bd_dmrs[bd_dmrs$ID %in% paste(dmrs_all$row.names,dmrs_all$CHROM,sep="_"),]
bd_dmrs <- bd_dmrs[c(1:3,5)]
bd_dmrs$Species <- "B. distachyon"

## arabidopsis thaliana
regions <- read.table(file="popgen5mb.bed")
dmrs_ceu <- read.table(file="10020.trim_methylome_CG.txt",header=T)
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
dmrs_ibnr <- read.table(file="9950.trim_methylome_CG.txt",header=T)
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
at_dmrs <- rbind(dmrs_ceu,dmrs_ibnr)
at_dmrs <- at_dmrs[!duplicated(at_dmrs),]
colnames(at_dmrs)[4] <- "length"
at_dmrs$Species <- "A. thaliana"
at_bd_dmrs <- rbind(at_dmrs,bd_dmrs)
mean(bd_dmrs$length)
mean(at_dmrs$length)

pdf(file="DMR_lengths_at_bd.pdf",height = 3, width= 2.5)
ggplot(data=at_bd_dmrs,aes(x=Species,y=log10(length)))+
  geom_boxplot() +
  xlab("") +
  ylab("log10(DMR length in bp)") +
  theme(axis.text.x=element_text(face="italic"))
dev.off()

#### Tajima's D and pi ####

gbm_genes_bd <- read.table(file="gbm_genes_bd.bed")
gbm_genes_bd$V1 <- gsub("NC_016131.3",1,gbm_genes_bd$V1)
gbm_genes_bd$V1 <- gsub("NC_016132.3",2,gbm_genes_bd$V1)
gbm_genes_bd$V1 <- gsub("NC_016133.3",3,gbm_genes_bd$V1)
gbm_genes_bd$V1 <- gsub("NC_016134.3",4,gbm_genes_bd$V1)
gbm_genes_bd$V1 <- gsub("NC_016135.3",5,gbm_genes_bd$V1)
gbm_genes_bd$region <- paste(gbm_genes_bd$V1,":",gbm_genes_bd$V2,"-",gbm_genes_bd$V3,sep="")

Dm_1.1 <- read.csv(file="Dm_theta_exons_group1.1.txt",sep = "\t")
Dm_1.1$region <- gsub("",":",Dm_1.1$region)
Dm_1.1 <- Dm_1.1[Dm_1.1$stop > 10,]
Dm_1.1 <- Dm_1.1[Dm_1.1$segregation_sites > 2,]
Dm_1.1$group <- "1.1"
Dm_1.1$Type <- "non-gbM"
Dm_1.1[Dm_1.1$region %in% gbm_genes_bd$region,]$Type <- "gbM"
Dm_1.1$pi <- Dm_1.1$theta_pi/Dm_1.1$stop

Dm_1.2 <- read.csv(file="Dm_theta_exons_group1.2.txt",sep = "\t")
Dm_1.2$region <- gsub("",":",Dm_1.2$region)
Dm_1.2 <- Dm_1.2[Dm_1.2$stop > 10,]
Dm_1.2 <- Dm_1.2[Dm_1.2$segregation_sites > 2,]
Dm_1.2$group <- "1.2"
Dm_1.2$Type <- "non-gbM"
Dm_1.2[Dm_1.2$region %in% gbm_genes_bd$region,]$Type <- "gbM"
Dm_1.2$pi <- Dm_1.2$theta_pi/Dm_1.2$stop

Dm_6 <- read.csv(file="Dm_theta_exons_group6.txt",sep = "\t")
Dm_6$region <- gsub("",":",Dm_6$region)
Dm_6 <- Dm_6[Dm_6$stop > 10,]
Dm_6 <- Dm_6[Dm_6$segregation_sites > 2,]
Dm_6$group <- "6"
Dm_6$Type <- "non-gbM"
Dm_6[Dm_6$region %in% gbm_genes_bd$region,]$Type <- "gbM"
Dm_6$pi <- Dm_6$theta_pi/Dm_6$stop

Dm_bd <- rbind(Dm_1.1,Dm_1.2,Dm_6)
df_bd_pi_gene <- dplyr::count(Dm_bd,Type,group)

Dm_bd$treatment <- paste(Dm_bd$group,Dm_bd$Type,sep="_")
summary(aov(lm(data=Dm_bd,pi~treatment)))
smp_test_pvals_bd <- HSD.test(aov(lm(data=Dm_bd,pi~treatment)), trt = 'treatment')
smp_test_pvals_bd

pdf("pi_bd.pdf",height=3,width=3)
ggplot(data=Dm_bd,aes(x=Type,y=log10(pi),fill=group)) +
  geom_boxplot() +
  xlab("") +
  ylab(expression('log'[10]*'(Scaled '*pi*')')) +
  ylim(-5,0.5) +
  theme(legend.position="top") +
  scale_fill_manual("Clonal group",values=c("#FDAE61","#ABDDA4","#2B83BA"))+
  geom_text(data = df_bd_pi_gene, aes(x=c(0.75,1,1.25,1.75,2,2.25),y = -0.7, label = c("d","c","ab","c","bc","a"))) +
  geom_text(data = df_bd_pi_gene, aes(x=c(0.75,1,1.25,1.75,2,2.25),y = 0.2, angle=90,label = n))
dev.off()

summary(aov(lm(data=Dm_bd,Dm~treatment)))
smp_test_pvals_D_bd <- HSD.test(aov(lm(data=Dm_bd,Dm~treatment)), trt = 'treatment')
smp_test_pvals_D_bd

pdf("D_bd.pdf",height=3,width=3)
ggplot(data=Dm_bd,aes(x=Type,y=Dm,fill=group)) +
  geom_boxplot() +
  xlab("") +
  ylab(expression("Tajima's "*italic('D'))) +
  theme(legend.position="none") +
  scale_fill_manual("Clonal group",values=c("#FDAE61","#ABDDA4","#2B83BA")) +
  geom_text(data = df_bd_pi_gene, aes(x=c(0.75,1,1.25,1.75,2,2.25),y = 5, label = c("e","bc","c","d","a","ab")))
dev.off()

D_bd_group1 <- inner_join(Dm_1.1[c(3:6,8:10)],Dm_1.2[c(3:6,8:10)],by="region")
pa_gbm_bd <- ggscatter(data=D_bd_group1[D_bd_group1$Type.x %in% "gbM",],x="Dm.x",y="Dm.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,cor.coef.coord = c(-3, 4),size=1) +
  xlab(expression("Group 1.1 Tajima's "*italic('D'))) +
  ylab(expression("Group 1.2 Tajima's "*italic('D'))) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pb_gbm_bd <- ggscatter(data=D_bd_group1[D_bd_group1$Type.x %in% "gbM",],x="theta_pi.x",y="theta_pi.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,size=1) +
  xlab(expression(paste("Group 1.1 ",pi))) +
  ylab(expression(paste("Group 1.2 ",pi))) +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")
pc_gbm_bd <- ggscatter(data=D_bd_group1[D_bd_group1$Type.x %in% "gbM",],x="segregation_sites.x",y="segregation_sites.y", shape = 21, add = "reg.line",  conf.int = TRUE, cor.coef = T,size=1) +
  xlab("Group 1.1 number of segregating sites") +
  ylab("Group 1.2 number of segregating sites") +
  geom_abline(intercept = 0, color="blue",linetype = "dashed")

pdf(file="D_group1.pdf",height=3,width=10)
plot_grid(pa_gbm_bd,pb_gbm_bd,pc_gbm_bd,ncol=3)
dev.off()

