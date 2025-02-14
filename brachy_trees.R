setwd("C:/Users/Andreas/Desktop/brachy_vcf_redo_jan25/")
getwd()

#BiocManager::install("ggtree")
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(ape)
library(RColorBrewer)
theme_set(theme_bw(base_size = 15))

display.brewer.pal(n = 4, name = 'Spectral')
brewer.pal(n = 4, name = "Spectral")
# "#D7191C" "#FDAE61" "#ABDDA4" "#2B83BA" ## groups

display.brewer.pal(n = 3, name = 'RdBu')
brewer.pal(n = 3, name = "RdBu")
# "#EF8A62" "#F7F7F7" "#67A9CF" ## context


tip.group <- read.csv("samplenames_group.txt",sep="\t",header=TRUE)
tip.group
tip.group <- tip.group[1:22,] # only group 1
colnames(tip.group)[1] <- "label"
tip.subgroup <- tip.group
tip.subgroup <- arrange(tip.subgroup, desc(label))
tip.subgroup$subgroup <- c(rep("group_1.2",12),rep("group_1.1",10))
tip.subgroup <- tip.subgroup[,-2]


tip.subgroup2 <- tip.subgroup[-c(18,19,21,22),] # only central region samples BdTR1e, BdTR1h, BdTR1j, BdTR1k, BdTR1m, BdTR1n of group 1.1

#### gBM SMP tree central only
smp_gbm_tree<-read.tree("./for_filtering/distance_matrix/smp1_central_gbm_dis_mat_fastme-tree.nwk")
tree.a13 <- full_join(smp_gbm_tree,tip.subgroup2,by="label")
ggtree(tree.a13, 
       layout = "ape",
       size=0.6,
       aes(color=subgroup), show.legend=F)+
  geom_tiplab(size=5,face="italic")+
  geom_treescale(x = -0.06, y = -0.07, width = NULL, offset = 0.02, color = "black", linesize = 0.5, fontsize = 3.88, family = "sans")+
  labs(title = "Neighbour-joining tree gBM SMPs")+
  coord_cartesian(clip = "on", expand = TRUE, xlim = c(-0.1, 0.6), ylim = c(-0.05, 0.4))+
  theme(legend.position = "bottom")+
  #  geom_tippoint()+
  #  theme(plot.title = element_text(hjust = 1))+
  theme(legend.text = element_text(face="italic",size=12))+
  theme(legend.title = element_blank())+
  scale_colour_manual(labels = c("Group 1.1","Group 1.2"), values = c("#FDAE61","#ABDDA4"), breaks = c("group_1.1","group_1.2"))+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))+
  theme(plot.title = element_text(size=18))+
  theme(plot.subtitle = element_text(size=14))

#### gBM DMR tree central only
dmr_gbm_tree<-read.tree("./for_filtering/distance_matrix/dmr1_central_gbm_dis_mat_fastme-tree.nwk")
tree.a14 <- full_join(dmr_gbm_tree,tip.subgroup2,by="label")
ggtree(tree.a14, 
       layout = "ape",
       size=0.5,
       aes(color=subgroup), show.legend=F)+
  geom_tiplab(size=5,face="italic")+
  geom_treescale(x = 0.45, y = -0.1, width = NULL, offset = 0.02, color = "black", linesize = 0.5, fontsize = 3.88, family = "sans")+
  labs(title = "Neighbour-joining tree gBM DMRs")+
  coord_cartesian(clip = "on", expand = TRUE, xlim = c(-0.1, 0.55), ylim = c(-0.1, 0.45))+
  theme(legend.position = "bottom")+
  #  geom_tippoint()+
  #  theme(plot.title = element_text(hjust = 1))+
  theme(legend.text = element_text(face="italic",size=12))+
  theme(legend.title = element_blank())+
  scale_colour_manual(labels = c("Group 1.1","Group 1.2"), values = c("#FDAE61","#ABDDA4"), breaks = c("group_1.1","group_1.2"))+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))+
  theme(plot.title = element_text(size=18))+
  theme(plot.subtitle = element_text(size=14))

#





###############################################
#trees all of group 1

#### gBM SMP tree
smp_gbm_tree<-read.tree("smp1_gbm_dis_mat_fastme-tree.nwk")
tree.a9 <- full_join(smp_gbm_tree,tip.subgroup,by="label")
ggtree(tree.a9, 
       layout = "ape",
       size=0.6,
       aes(color=subgroup), show.legend=F)+
  geom_tiplab(size=5,face="italic")+
  geom_treescale(x = -0.11, y = -0.15, width = NULL, offset = 0.02, color = "black", linesize = 0.5, fontsize = 3.88, family = "sans")+
  labs(title = "Neighbour-joining tree gBM SMPs")+
  coord_cartesian(clip = "on", expand = TRUE, xlim = c(-0.2, 0.5), ylim = c(-0.15, 0.4))+
  theme(legend.position = "bottom")+
  #  geom_tippoint()+
  #  theme(plot.title = element_text(hjust = 1))+
  theme(legend.text = element_text(face="italic",size=12))+
  theme(legend.title = element_blank())+
  scale_colour_manual(labels = c("Group 1.1","Group 1.2"), values = c("#FDAE61","#ABDDA4"), breaks = c("group_1.1","group_1.2"))+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))+
  theme(plot.title = element_text(size=18))+
  theme(plot.subtitle = element_text(size=14))



#### gBM DMR tree
dmr_gbm_tree<-read.tree("dmr1_gbm_dis_mat_fastme-tree.nwk")
tree.a10 <- full_join(dmr_gbm_tree,tip.subgroup,by="label")
ggtree(tree.a10, 
       layout = "ape",
       size=0.5,
       aes(color=subgroup), show.legend=F)+
  geom_tiplab(size=5,face="italic")+
  geom_treescale(x = -0.11, y = -0.15, width = NULL, offset = 0.02, color = "black", linesize = 0.5, fontsize = 3.88, family = "sans")+
  labs(title = "Neighbour-joining tree gBM DMRs")+
  coord_cartesian(clip = "on", expand = TRUE, xlim = c(-0.2, 0.5), ylim = c(-0.15, 0.4))+
  theme(legend.position = "bottom")+
  #  geom_tippoint()+
  #  theme(plot.title = element_text(hjust = 1))+
  theme(legend.text = element_text(face="italic",size=12))+
  theme(legend.title = element_blank())+
  scale_colour_manual(labels = c("Group 1.1","Group 1.2"), values = c("#FDAE61","#ABDDA4"), breaks = c("group_1.1","group_1.2"))+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7)))+
  theme(plot.title = element_text(size=18))+
  theme(plot.subtitle = element_text(size=14))
