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

group1.1_smp_prop <- 389774/13537770
group1.2_smp_prop <- 76079/9640492
group1.1_smp_gbm_prop <- 40187/1154347
group1.2_smp_gbm_prop <- 34609/879235
gbm_genes_bd <- read.table(file="gbm_genes_bd.bed")
1154347/sum(gbm_genes_bd$V3-gbm_genes_bd$V2) ## to get SMP genotype density
1154347/sum(gbm_genes_bd$V3-gbm_genes_bd$V2) ## to get SMP genotype density

propsites_bd <- as.data.frame(rbind(group1.1_smp_prop,group1.2_smp_prop,group1.1_smp_gbm_prop,group1.2_smp_gbm_prop))
colnames(propsites_bd) <- "Proportion"
propsites_bd$Group <- c("1.1","1.2","1.1","1.2")
propsites_bd$Type <- c("All sites","All sites","gbM sites","gbM sites")
total_bd <- c("n=10062950","n=9640492","n=905826","n=879235")

seg_sites_bd <- ggplot(data=propsites_bd,aes(x=Type,y=Proportion,fill=Group)) +
  geom_col(position = "dodge") +
  ylim(c(0,0.1)) +
  scale_fill_manual(values=c("#FDAE61","#ABDDA4")) +
  geom_text(aes(x=c(0.8,1.2,1.8,2.2),angle=90,y = 0.085, label = total_bd)) +
  ylab("Proportion of polymorphic sites") +
  xlab("")

pdf(file="seg_sites.pdf",height=4,width=9)
plot_grid(seg_athaliana,seg_sites_bd,ncol=2,rel_widths = c(2,1),labels = "AUTO")
dev.off()

#### DAPC #### 

smp_vcf <- read.vcfR("smp1_central_gbm_maf10.vcf")
sample_data <- data.frame(sample_names = colnames(smp_vcf@gt)[-c(1)])
group_data <- data.frame(sample_names = colnames(smp_vcf@gt)[-c(1)],Group = c(rep("group_1.1",6),rep("group_1.2",12)))
joined_data <- full_join(sample_data, group_data, by = "sample_names")
group_names <- c("group_1.1","group_1.2")
group_colors <- c( "#FDAE61","#ABDDA4")
myCols <- setNames(group_colors,group_names)
smp_genlight <- vcfR2genlight(smp_vcf)
smp_grp <- find.clusters(smp_genlight,max.n.clust=15)

smp_1 <- gl.define.pop(smp_genlight, ind.list=c("BdTR1e","BdTR1h","BdTR1j","BdTR1k","BdTR1m","BdTR1n"), new="group_1.1", verbose = NULL)
smp <- gl.define.pop(smp_1, ind.list=c("BdTR2b","BdTR2c","BdTR2d","BdTR2g","BdTR2h","BdTR2j","BdTR2k","BdTR2m","BdTR2n","BdTR2p","BdTR2r","BdTR2s"), new="group_1.2", verbose = NULL)
smp_nafill<-apply(smp,2,na.replace,mean, na.rm = TRUE)
mat <- as.matrix(smp)
genotype_matrix <- apply(mat,2,na.replace,mean, na.rm = TRUE)
xval <- xvalDapc(genotype_matrix, smp$pop, n.pca.max = 5, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)
dapc <- dapc(smp, n.da=1, n.pca=1)

pdf("DAPC_SMP_PC2_Bd.pdf", width = 2.7, height = 2.7)
par(cex=0.6)
scatter(dapc, bg = "white", pch = 20, cstar = 0, col = myCols, solid = 0.4, cex = 4, clab = 0,
        mstree = F, scree.da = FALSE, posi.leg = "topright",
        leg = TRUE, txt.leg = group_names)
title("SMP",font.main=2)
par(xpd=TRUE,new=TRUE)
dev.off()

dmr_vcf <- read.vcfR("dmr1_central_gbm_maf10.vcf")
sample_data <- data.frame(sample_names = colnames(dmr_vcf@gt)[-c(1)])

group_data <- data.frame(sample_names = colnames(dmr_vcf@gt)[-c(1)],Group = c(rep("group_1.1",6),rep("group_1.2",12)))
joined_data <- full_join(sample_data, group_data, by = "sample_names")
group_names <- c("group_1.1","group_1.2")
group_colors <- c( "#FDAE61","#ABDDA4")
myCols <- setNames(group_colors,group_names)
dmr_genlight <- vcfR2genlight(dmr_vcf)
dmr_grp <- find.clusters(dmr_genlight,max.n.clust=15)
10
2
dmr_1 <- gl.define.pop(dmr_genlight, ind.list=c("BdTR1e","BdTR1h","BdTR1j","BdTR1k","BdTR1m","BdTR1n"), new="group_1.1", verbose = NULL)
dmr <- gl.define.pop(dmr_1, ind.list=c("BdTR2b","BdTR2c","BdTR2d","BdTR2g","BdTR2h","BdTR2j","BdTR2k","BdTR2m","BdTR2n","BdTR2p","BdTR2r","BdTR2s"), new="group_1.2", verbose = NULL)
dmr_nafill<-apply(dmr,2,na.replace,mean, na.rm = TRUE)
mat <- as.matrix(dmr)
genotype_matrix <- apply(mat,2,na.replace,mean, na.rm = TRUE)

xval <- xvalDapc(genotype_matrix, dmr$pop, n.pca.max = 5, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)
xval
dapc <- dapc(dmr, n.da=1, n.pca=1)

pdf("DAPC_DMR_PC2_Bd.pdf", width = 2.7, height = 2.7)
par(cex=0.6)
scatter(dapc, bg = "white", pch = 20, cstar = 0, col = myCols, solid = 0.4, cex = 4, clab = 0,
        mstree = F, scree.da = FALSE,txt.leg = group_names)
title("DMR",font.main=2)
par(xpd=TRUE,new=TRUE)
dev.off()

#### SFS #### 

## gbm regions
group1.2_smps.frq <- read.table(file="smp1.2_gbm_maf10.frq",row.names = NULL, header = T)
group1.2_smps.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_smps.frq$N_CHR))
group1.2_smps.frq <- group1.2_smps.frq[(group1.2_smps.frq$N_CHR > 0.02) & (group1.2_smps.frq$N_CHR < 0.98),]
group1.2_smps.frq$Group <- "1.2"
group1.2_smps.frq$Type <- "SMP"
bd_frq <- group1.2_smps.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Group,Type, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Group,Type) %>%
  dplyr::mutate(proportion = count / sum(count))
bd_frq$N_CHR_binned <- seq(0.125,0.975,0.1)
bd_frq_smp <- bd_frq[bd_frq$Type %in% "SMP",]
bd_frq_smp$Group <- gsub("1.2",paste("gbM, n=",nrow(group1.2_smps.frq),sep=""), bd_frq_smp$Group)

## non gbm regions
group1.2_smps_nongbm.frq <- read.table(file="smp1.2_exon_nongbm_maf10.frq",row.names = NULL, header = T)
group1.2_smps_nongbm.frq$N_CHR <- as.numeric(gsub("A:","",group1.2_smps_nongbm.frq$N_CHR))
group1.2_smps_nongbm.frq <- group1.2_smps_nongbm.frq[(group1.2_smps_nongbm.frq$N_CHR > 0.02) & (group1.2_smps_nongbm.frq$N_CHR < 0.98),]
group1.2_smps_nongbm.frq$Group <- "1.2"
group1.2_smps_nongbm.frq$Type <- "SMP"
bd_nongbm_frq <- group1.2_smps_nongbm.frq %>%
  dplyr::mutate(N_CHR_binned = cut(N_CHR, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(Group,Type, N_CHR_binned) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(Group,Type) %>%
  dplyr::mutate(proportion = count / sum(count))
bd_nongbm_frq$N_CHR_binned <- seq(0.125,0.975,0.1)
bd_nongbm_frq_smp <- bd_nongbm_frq[bd_nongbm_frq$Type %in% "SMP",]
bd_nongbm_frq_smp$Group <- gsub("1.2",paste("non-gbM, n=",nrow(group1.2_smps_nongbm.frq),sep=""), bd_nongbm_frq_smp$Group)

bd_combined_frq <- rbind(bd_frq_smp,bd_nongbm_frq_smp)

bd_frq_plot <- ggplot(bd_combined_frq, aes(x=N_CHR_binned, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Proportion of unmethylated alleles", y="Proportion of sites") +
  theme_minimal() +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme(legend.title = element_blank(),legend.position = c(0.7, 0.9))

pdf(file="bd_sfs.pdf",height=3,width=3.5)
bd_frq_plot
dev.off()


#### LD ####

bd_ld_decay <- read.table(file="smp1.2_gbm_maf10_ld_decay.stat",header = T)
bd_ld_decay$Group <- "gbM SMP"
colnames(bd_ld_decay) <- c("Dist","Mean_r2","Group")

bd_all_ld_decay <- read.table(file="smp1.2_maf10_ld_decay.stat",header = T)
bd_all_ld_decay$Group <- "All SMP"
colnames(bd_all_ld_decay) <- c("Dist","Mean_r2","Group")

bd_ld_decay_both <- rbind(bd_ld_decay,bd_all_ld_decay)

bd_ld_decay_plot <- ggplot(data=bd_ld_decay_both,aes(x=Dist,y=Mean_r2,color=Group)) +
  geom_line(linewidth=1) +
  xlab("Base pairs") +
  ylab(expression(paste("r"^"2"))) +
  scale_color_manual(labels = c("All SMP","gbM SMP"),values=c("brown","lightblue")) +
  ggtitle("Group 1.2") +
  xlim(0,100) +
  ylim(0,0.6) +
  theme_classic2() +
  theme(legend.position = c(0.55,0.84),legend.title = element_blank(),legend.text = element_text(size = 8))

pdf("LD.pdf",height=2.5,width = 8.5)
plot_grid(ld_decay_plot,bd_ld_decay_plot,labels=c("A","B"),ncol=2,rel_widths = c(2.8,1))
dev.off()

#### DMR lengths ####

## brachypodium

bd_dmrs <- read.table(file="bd_dmrs_lengths.txt",header=T)
bd_dmrs <- bd_dmrs[!duplicated(bd_dmrs),]
bd_dmrs <- bd_dmrs[!is.na(bd_dmrs$lth),]
bd_dmrs <- bd_dmrs[2]
colnames(bd_dmrs)[1] <- "length"
bd_dmrs$Species <- "B. distachyon"

## arabidopsis thaliana
dmrs_ceu <- read.table(file="ceu_dmrs_lengths.txt",header=T)
#dmrs_ceu$Pop <- "CEU"
dmrs_ibnr <- read.table(file="ibnr_dmrs_lengths.txt",header=T)
#dmrs_ibnr$Pop <- "IBnr"
at_dmrs <- rbind(dmrs_ceu,dmrs_ibnr)
at_dmrs <- at_dmrs[!duplicated(at_dmrs),]
colnames(at_dmrs)[2] <- "length"
at_dmrs$Species <- "A. thaliana"
at_dmrs <- at_dmrs[2:3]
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

Dm_1.2 <- read.csv(file="Dm_theta_exons_group1.2.txt",sep = "\t")
Dm_1.2$region <- gsub("",":",Dm_1.2$region)
Dm_1.2 <- Dm_1.2[Dm_1.2$stop > 10,]
Dm_1.2 <- Dm_1.2[Dm_1.2$segregation_sites > 2,]
Dm_1.2$group <- "1.2"
Dm_1.2$Type <- "non-gbM"
Dm_1.2[Dm_1.2$region %in% gbm_genes_bd$region,]$Type <- "gbM"
Dm_1.2$pi <- Dm_1.2$theta_pi/Dm_1.2$stop
Dm_1.2$Type <- factor(Dm_1.2$Type,levels=c("gbM","non-gbM"))
df_bd_pi_gene <- dplyr::count(Dm_1.2,Type,group)

pdf("pi_bd.pdf",height=3.1,width=2.5)
ggboxplot(Dm_1.2,x="Type",y="log10(pi)",fill="Type") +
  stat_compare_means(comparisons = list(c(1,2)),label.y = -0.9) +
  geom_text(data = df_bd_pi_gene, aes(x=c(1,2),y = 0.2, angle=90,label = n)) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme(legend.position = "none") +
  xlab("")+
  ylab(expression('log'[10]*'(Scaled '*pi[SMP]*')')) +
  ylim(-4,0.3)
dev.off()

pdf("D_bd.pdf",height=3.1,width=2.5)
ggboxplot(Dm_1.2,x="Type",y="Dm",fill="Type") +
  stat_compare_means(comparisons = list(c(1,2)),label.y = 4.7) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme(legend.position = "none") +
  xlab("")+
  ylim(-3,5.2) +
  ylab(expression(italic('D')[SMP]))
dev.off()
