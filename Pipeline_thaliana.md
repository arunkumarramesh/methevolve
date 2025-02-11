# Pipeline for getting SMPs and SNPs for Arabidopsis thaliana

1. Obtain A.thaliana vcf from https://1001genomes.org/data/GMI-MPI/releases/v3.1/ and filter

```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps
gunzip 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip 1001genomes_snp-short-indel_only_ACGTN.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -v snps --max-alleles 2  -O z -o 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz   1001genomes_snp-short-indel_only_ACGTN.vcf.gz
```
2. Create pseudogenomes
```
#invcf is a file with sample names that are found in the WGBS and genomic vcf file

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 1 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr1.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr2.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 3 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr3.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 4 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr4.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa 5 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_chr5.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Mt | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_Mt.fa; done

cat invcf | while read line; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Pt | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $line -p ${line/$/_} 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz >>1001genomes_snp-short-indel_only_ACGTN_snps_Pt.fa; done

mkdir /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas

/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr1.fa chr1/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr2.fa chr2/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr3.fa chr3/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr4.fa chr4/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_chr5.fa chr5/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_Pt.fa Pt/
/proj/popgen/a.ramesh/software/faSplit byname 1001genomes_snp-short-indel_only_ACGTN_snps_Mt.fa Mt/

mkdir /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/snps/fastas/merged
cat ../../invcf | while read line ; do cat $line*[0-9].fa >$line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do java -jar /proj/popgen/a.ramesh/software/picard.jar CreateSequenceDictionary R=$file O=${file/.fa/.dict}; done

```

3. Download methyaltion data from NCBI
```
cd  /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/
/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file PRJNA187927_adabidopsis_single_Acc_List.txt  -O data
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done

/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file PRJNA187927_adabidopsis_paired_Acc_List.txt  -O data/paired
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/paired
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done
```
4. Trim data and concatenate or rename files. Renaming in A.thaliana rename scripts folder
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/
for file in *.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 5 $file ${file/.fastq.gz/_trim.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/paired
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

5. Map methylation reads
```
sed 's/^/\/proj\/popgen\/a.ramesh\/projects\/methylomes\/arabidopsis\/pseudogenomes\//' samplenames5 | paste samplenames5 >samplenames6
sed 's/^/\/proj\/popgen\/a.ramesh\/projects\/methylomes\/arabidopsis\/pseudogenomes\//' samplenames2 | paste samplenames2 >samplenames3

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
cat samplenames6 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 $value1.trim.fq.gz  ; done

cat samplenames3 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
```

6. Deduplicate methylation reads
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
for file in *_bismark_hisat2*.bam ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

7. Call methylation variants
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
cat samplenames3 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam  ; done

cat samplenames6 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.trim_bismark_hisat2.deduplicated.bam  ; done

gunzip *.bismark.cov.gz
gunzip *.CpG_report.txt.gz
```

8. Do binomial test
```
setwd("/proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data")
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")

cov <- read.table(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "_")

context <- read.table(file=contextfiles[1,],header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "_")

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub(gsub("\\..*","",covfile$V1)[1],"",cov_context$chromosome)
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov","",covfile[1,])

b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]

colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  print(i)
  
  cov <- read.table(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "_")
  
  context <- read.table(file=contextfiles[i,],header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "_")
  
  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub(gsub("\\..*","",covfile$V1)[i],"",cov_context$chromosome)
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
  if (is.na(noncoversionrate)) {
    noncoversionrate <- 0
  }
  samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov","",covfile[i,])
  
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  

  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
}
cov_context2$chromosome <- gsub("_.*","",cov_context2$ID)
cov_context2$position   <- as.numeric(gsub(".*_","",cov_context2$ID))

write.table(cov_context2,file="cov_context3.txt",row.names = F)

```
9. Rscript to create methylation vcfs

```
library(dplyr)

cov_context3 <- read.table(file="cov_context3.txt",header=T)
colnames(cov_context3) <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",colnames(cov_context3))
dupsamples <- colnames(cov_context3)[duplicated(colnames(cov_context3))]
alldedupdata <- ""
for (i in 1:length(dupsamples)){
  dupdata <- cov_context3[colnames(cov_context3) %in% (dupsamples[i])]
  dupdata[which(dupdata[,1] != dupdata[,2]),1] <- NA
  dupdata[which(dupdata[,1] != dupdata[,2]),1] <- NA
  dedupdata <- dupdata %>% 
    mutate(new = coalesce(dupdata[,1],dupdata[,2]))
  dedupdata <- dedupdata[3]
  colnames(dedupdata) <- colnames(dupdata)[1]
  alldedupdata <- cbind(alldedupdata,dedupdata)
}
cov_context3 <- cov_context3[!(colnames(cov_context3) %in% dupsamples)]
cov_context3 <- cbind(cov_context3,alldedupdata)

na_count <- apply(cov_context3[6:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[6:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[6:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[1:3]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"

cov_context3 <- cov_context3[-c(1:5)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
cov_context3 <- cbind(meta,cov_context3)
cov_context3$`#CHROM` <- as.numeric(cov_context3$`#CHROM`)
cov_context3$POS <- as.numeric(cov_context3$POS)
cov_context3 <- cov_context3[order(cov_context3$`#CHROM`,cov_context3$POS),]
write.table(cov_context3,file="arabidopsis_meth_cg_combine.vcf",quote = F, row.names = F,sep="\t")

meta2 <- cov_context4[1:3]
colnames(meta2) <- c("#CHROM","POS","ID")
meta2$REF <- "A"
meta2$ALT <- "T"
meta2$QUAL <- 4000
meta2$FILTER <- "PASS"
meta2$INFO <- "DP=1000"
meta2$FORMAT <- "GT"
cov_context4 <- cov_context4[-c(1:5)]
cov_context4[cov_context4 == "U"] <- "0/0"
cov_context4[cov_context4 == "M"] <- "1/1"
cov_context4[is.na(cov_context4)] <- "./."
cov_context4 <- cbind(meta2,cov_context4)
cov_context4$`#CHROM` <- as.numeric(cov_context4$`#CHROM`)
cov_context4$POS <- as.numeric(cov_context4$POS)
cov_context4 <- cov_context4[order(cov_context4$`#CHROM`,cov_context4$POS),]
write.table(cov_context4,file="arabidopsis_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")

```

10. Split reference genome into genes and other regions

```
sed -e 's/\t/:/' -e  's/\t/-/' gene_pos.bed >gene_pos.list
grep -F -f gbm_genes Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'exon' | cut -f 1,4,5 > gbm_genes.bed
grep -v -F -f gbm_genes Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'exon' | cut -f 1,4,5 > non_gbm_genes.bed
# high sift is sift score > 0.75 from File S1 from Hämälä and Tiffin, 2020
grep -F -f high_sift Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'exon' | cut -f 1,4,5 > high_sift.bed
cut -f 1-2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai >chr_lengths
/data/proj2/popgen/a.ramesh/software/bedtools2/bin/bedtools complement -i gene_pos.bed -g chr_lengths >non_gene.bed
sed -e 's/\t/:/' -e  's/\t/-/' non_gene.bed >non_gene.list
grep -F -f low_sift Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'ID=gene' | cut -f 1,4,5,9 | sed -e 's/ID=gene://' -e 's/;Name=.*//' >low_sift_genes_poswithnames.txt
grep -F -f high_sift Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'ID=gene' | cut -f 1,4,5,9 | sed -e 's/ID=gene://' -e 's/;Name=.*//' >high_sift_genes_poswithnames.txt
grep -F -f gbm_genes Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'ID=gene' | cut -f 1,4,5,9 | sed -e 's/ID=gene://' -e 's/;Name=.*//' >gbm_genes_poswithnames.txt
grep -v -F -f gbm_genes Arabidopsis_thaliana.TAIR10.55.gff3 | grep 'ID=gene' | cut -f 1,4,5,9 | sed -e 's/ID=gene://' -e 's/;Name=.*//'  >non_gbm_genes_poswithnames.txt


cat /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $line >>genes.fasta; done
mkdir genes_fasta/
/data/proj2/popgen/a.ramesh/software/faSplit byname genes.fasta genes_fasta/
cd genes_fasta/
ls *fa >filenames

cat /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/cg_sites.list | while read -r line ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $line >>cg.fasta; done
mkdir cg_fasta/
/data/proj2/popgen/a.ramesh/software/faSplit byname cg.fasta cg_fasta/
#Rscript  count_c.R

```

11. Rscript to count number of cytosines

```
library("methimpute",lib.loc="/data/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")

## Only CG context
files <- read.table(file="filenames")
files <- as.character(files$V1)

cytsosine_count <- ""
for (f in 1:length(files)){
  print(f)
  cytosines <- c()
  try(cytosines <- extractCytosinesFromFASTA(files[f], contexts = 'CG'),silent=T)
  if (length(cytosines) > 0){
    cytsosine_count <- rbind(cytsosine_count,(c(files[f],table(cytosines$context))))
  } else {
    cytsosine_count <- rbind(cytsosine_count,(c(files[f],0)))
  }
  #print(cytsosine_count)
}
cytsosine_count <- cytsosine_count[-c(1),]
write.table(cytsosine_count,file="cytsosine_count.txt",row.names=F, col.names=F,quote=F,sep="\t")
```

12. Get list of bismark files for jDMR

```
sed 's/ dna.*//' /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa > Arabidopsis_thaliana.TAIR10.dna.jDMR.fa
sed 's/_bismark_hisat2.*.deduplicated.CpG_report_mod.txt/_methylome.txt/' ceu_reportfiles >ceu_dmrsamples
sed 's/\..*//' ceu_dmrsamples | paste ceu_dmrsamples - >ceu_dmrsamples2
sed 's/_bismark_hisat2.*.deduplicated.CpG_report_mod.txt/_methylome.txt/' ibnr_reportfiles >ibnr_dmrsamples
sed 's/\..*//' ibnr_dmrsamples | paste ibnr_dmrsamples - >ibnr_dmrsamples2
```

13. Rscript to call differentially methylated regions using jDMR (jdmr.R)

```
## modify chr names by removing sample names from them. Was a feature of creating pseudogenomes
contextfiles <- read.table(file="ceu_contextfiles",header = F)

for(i in 1:nrow(contextfiles)){
  contextdata <- read.table(file=contextfiles[i,],header=F)
  contextdata$V1 <- gsub(gsub("\\..*","",contextfiles$V1)[i],"",contextdata$V1)
  write.table(contextdata,file=gsub(".txt","_mod.txt",contextfiles$V1)[i],row.names = F, col.names = F, quote = F)
}

contextfiles <- read.table(file="ibnr_contextfiles",header = F)

for(i in 1:nrow(contextfiles)){
  contextdata <- read.table(file=contextfiles[i,],header=F)
  contextdata$V1 <- gsub(gsub("\\..*","",contextfiles$V1)[i],"",contextdata$V1)
  write.table(contextdata,file=gsub(".txt","_mod.txt",contextfiles$V1)[i],row.names = F, col.names = F, quote = F)
}

library("methimpute",lib.loc="/data/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")
library("jDMR",lib.loc="/data/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")

reportfiles <- read.table(file = "ceu_reportfiles")
for(i in 1:nrow(reportfiles)){
  bismark.data <- importBismark(reportfiles[i,])
  distcor <- distanceCorrelation(bismark.data, separate.contexts = TRUE)
  fit <- estimateTransDist(distcor)
  model <- callMethylationSeparate(data = bismark.data, transDist = fit$transDist, verbosity = 0)
  filename <- as.character(gsub("_bismark_hisat2.*.deduplicated.CpG_report_mod","_methylome",reportfiles[i,]))
  exportMethylome(model, filename = filename)
}

reportfiles <- read.table(file = "ibnr_reportfiles")
for(i in 1:nrow(reportfiles)){
  bismark.data <- importBismark(reportfiles[i,])
  distcor <- distanceCorrelation(bismark.data, separate.contexts = TRUE)
  fit <- estimateTransDist(distcor)
  model <- callMethylationSeparate(data = bismark.data, transDist = fit$transDist, verbosity = 0)
  filename <- as.character(gsub("_bismark_hisat2.*.deduplicated.CpG_report_mod","_methylome",reportfiles[i,]))
  exportMethylome(model, filename = filename)
}


out.dir <- "jDMRresults/"

myfasta <- readDNAStringSet("Arabidopsis_thaliana.TAIR10.dna.jDMR.fa")
CfromFASTAv4(fasta = myfasta, chr = 1, out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", 1, ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = 1, min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("Arabidopsis_thaliana.TAIR10.dna.jDMR.fa")
CfromFASTAv4(fasta = myfasta, chr = 2, out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", 2, ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = 2, min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("Arabidopsis_thaliana.TAIR10.dna.jDMR.fa")
CfromFASTAv4(fasta = myfasta, chr = 3, out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", 3, ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = 3, min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("Arabidopsis_thaliana.TAIR10.dna.jDMR.fa")
CfromFASTAv4(fasta = myfasta, chr = 4, out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", 4, ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = 4, min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("Arabidopsis_thaliana.TAIR10.dna.jDMR.fa")
CfromFASTAv4(fasta = myfasta, chr = 5, out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", 5, ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = 5, min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")


runMethimputeRegions(out.dir = out.dir, samplefiles = "ceu_dmrsamples2", genome = "A_thaliana", context = "CG", nCytosines=8, mincov=3, Regionfiles = out.dir)
runMethimputeRegions(out.dir = out.dir, samplefiles = "ibnr_dmrsamples2", genome = "A_thaliana", context = "CG", nCytosines=8, mincov=3, Regionfiles = out.dir)
```

14. Run jDMR script and generate vcf file

```
Rscript jdmr.R
cd jDMRresults/
Rscript dmr_vcf_create.R
cat vcfheader ceu_dmrs_dip.txt >ceu_dmrs_dip.vcf
cat vcfheader ibnr_dmrs_dip.txt >ibnr_dmrs_dip.vcf
cat vcfheader ceu_dmrs.txt >ceu_dmrs.vcf
cat vcfheader ibnr_dmrs.txt >ibnr_dmrs.vcf
```

15. Rscript to create vcfs for DMRs (dmr_vcf_create.R)

```
library(dplyr)

filenames <- read.table(file="../ceu_dmrsamples2",header = T)
filenames$file <- gsub("_methylome.txt","_methylome.txt_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
dmrs <- dmrs[!duplicated(dmrs),]
dmrs[dmrs$posteriorMax < 0.99,]$status <- NA
gt <- dmrs$status
gt <- gsub("U",0,gt)
gt <- gsub("M",1,gt)
gt[is.na(gt)] <- "."
gt <- as.data.frame(gt)
gt2 <- gt
gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
colnames(gt)[1] <- filenames[1,2]
colnames(gt2)[1] <- filenames[1,2]
lth <- dmrs$end - dmrs$start
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs <- dmrs[3]
colnames(dmrs) <- "ID"
lengths <- as.data.frame(cbind(dmrs,lth))
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)

for (i in 2:nrow(filenames)){
  print(i)
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  dmrs_tmp <- dmrs_tmp[!duplicated(dmrs_tmp),]
  lth <- dmrs_tmp$end - dmrs_tmp$start
  dmrs_tmp$start <- round((dmrs_tmp$start + dmrs_tmp$end)/2)
  dmrs_tmp$end <- paste(dmrs_tmp$seqnames,dmrs_tmp$start,sep="_")
  ID <- dmrs_tmp$end
  lengths <- rbind(lengths,as.data.frame(cbind(ID,lth)))
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  gt <- cbind(gt,dmrs_tmp$end)
  colnames(gt)[2] <- "ID"
  gt2 <- cbind(gt2,dmrs_tmp$end)
  colnames(gt2)[2] <- "ID"
  gt <- gt[!duplicated(gt),]
  gt2 <- gt2[!duplicated(gt2),]
  dmrs <- full_join(dmrs,gt,by="ID")
  dmrs_dip <- full_join(dmrs_dip,gt2,by="ID")
  print(table(duplicated(dmrs_dip$ID)))
}

dmrs[is.na(dmrs)] <- "."
dmrs_dip[is.na(dmrs_dip)] <- "./."

CHROM <- as.integer(gsub("_.*","",dmrs$ID))
POS <- as.integer(gsub(".*_","",dmrs$ID))
ID <- dmrs$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "D"
meta$ALT <- "M"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs <- cbind(meta,dmrs[2:ncol(dmrs)])
dmrs$POS <- as.integer(dmrs$POS)
dmrs <- dmrs[order(dmrs$CHROM,dmrs$POS),]
colnames(dmrs)[1] <- "#CHROM"
write.table(dmrs,file="ceu_dmrs.txt",row.names = F, quote = F,sep="\t")

CHROM <- as.integer(gsub("_.*","",dmrs_dip$ID))
POS <- as.integer(gsub(".*_","",dmrs_dip$ID))
ID <- dmrs_dip$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs_dip <- cbind(meta,dmrs_dip[2:ncol(dmrs_dip)])
dmrs_dip$POS <- as.integer(dmrs_dip$POS)
dmrs_dip <- dmrs_dip[order(dmrs_dip$CHROM,dmrs_dip$POS),]
colnames(dmrs_dip)[1] <- "#CHROM"
write.table(dmrs_dip,file="ceu_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

lengths <- lengths[!duplicated(lengths),]
write.table(lengths,file="ceu_dmrs_lengths.txt",row.names = F, quote = F,sep="\t")

filenames <- read.table(file="../ibnr_dmrsamples2",header = T)
filenames$file <- gsub("_methylome.txt","_methylome.txt_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
dmrs <- dmrs[!duplicated(dmrs),]
dmrs[dmrs$posteriorMax < 0.99,]$status <- NA
gt <- dmrs$status
gt <- gsub("U",0,gt)
gt <- gsub("M",1,gt)
gt[is.na(gt)] <- "."
gt <- as.data.frame(gt)
gt2 <- gt
gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
colnames(gt)[1] <- filenames[1,2]
colnames(gt2)[1] <- filenames[1,2]
lth <- dmrs$end - dmrs$start
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs <- dmrs[3]
colnames(dmrs) <- "ID"
lengths <- as.data.frame(cbind(dmrs,lth))
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)

for (i in 2:nrow(filenames)){
  print(i)
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  dmrs_tmp <- dmrs_tmp[!duplicated(dmrs_tmp),]
  lth <- dmrs_tmp$end - dmrs_tmp$start
  dmrs_tmp$start <- round((dmrs_tmp$start + dmrs_tmp$end)/2)
  dmrs_tmp$end <- paste(dmrs_tmp$seqnames,dmrs_tmp$start,sep="_")
  ID <- dmrs_tmp$end
  lengths <- rbind(lengths,as.data.frame(cbind(ID,lth)))
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  gt <- cbind(gt,dmrs_tmp$end)
  colnames(gt)[2] <- "ID"
  gt2 <- cbind(gt2,dmrs_tmp$end)
  colnames(gt2)[2] <- "ID"
  gt <- gt[!duplicated(gt),]
  gt2 <- gt2[!duplicated(gt2),]
  dmrs <- full_join(dmrs,gt,by="ID")
  dmrs_dip <- full_join(dmrs_dip,gt2,by="ID")
  print(table(duplicated(dmrs_dip$ID)))
}

dmrs[is.na(dmrs)] <- "."
dmrs_dip[is.na(dmrs_dip)] <- "./."

CHROM <- as.integer(gsub("_.*","",dmrs$ID))
POS <- as.integer(gsub(".*_","",dmrs$ID))
ID <- dmrs$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "D"
meta$ALT <- "M"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs <- cbind(meta,dmrs[2:ncol(dmrs)])
dmrs$POS <- as.integer(dmrs$POS)
dmrs <- dmrs[order(dmrs$CHROM,dmrs$POS),]
colnames(dmrs)[1] <- "#CHROM"
write.table(dmrs,file="ibnr_dmrs.txt",row.names = F, quote = F,sep="\t")

CHROM <- as.integer(gsub("_.*","",dmrs_dip$ID))
POS <- as.integer(gsub(".*_","",dmrs_dip$ID))
ID <- dmrs_dip$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs_dip <- cbind(meta,dmrs_dip[2:ncol(dmrs_dip)])
dmrs_dip$POS <- as.integer(dmrs_dip$POS)
dmrs_dip <- dmrs_dip[order(dmrs_dip$CHROM,dmrs_dip$POS),]
colnames(dmrs_dip)[1] <- "#CHROM"
write.table(dmrs_dip,file="ibnr_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

lengths <- lengths[!duplicated(lengths),]
write.table(lengths,file="ibnr_dmrs_lengths.txt",row.names = F, quote = F,sep="\t")
```

16. Get the CEU and IBNR samples that overlap between SNP and methylation vcfs

```
head -n 1 arabidopsis_meth.vcf | cut -f 10- | tr '\t' '\n' >samples_invcf
cat vcfheader arabidopsis_meth.vcf >arabidopsis_meth_all.vcf
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out arabidopsis_meth_cg  --bed  cg_sites.bed --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_cg.recode.vcf  --out arabidopsis_meth_cg_5mb  --bed  popgen5mb.bed --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf arabidopsis_meth_cg_5mb.recode.vcf --keep ceu --maf 0.05 --max-missing 0.5 --out arabidopsis_meth_cg_ceu --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf arabidopsis_meth_cg_5mb.recode.vcf --keep ibnr --maf 0.05 --max-missing 0.2 --out arabidopsis_meth_cg_ibnr --recode
grep 'CHROM' arabidopsis_meth_cg_ceu.recode.vcf | cut -f 10- | tr '\t' '\n' >ceu_invcf
grep 'CHROM' arabidopsis_meth_cg_ibnr.recode.vcf | cut -f 10- | tr '\t' '\n' >ibnr_invcf
shuf -n 63 ceu_invcf >ceu_invcf_shuf # to keep the number of samples analysed for CEU and IBNR the same (n=63)
```

17. Number of variant and invariant SNPs for CEU and IBNR

```
bcftools view  -i "QUAL > 20 || INFO/DP > 20"  -O z -o 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz  1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz
tabix  -f 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz

bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed  -S ceu_invcf_shuf -i "F_MISSING < 0.8 || MAF < 0.0001" -O v -m 1 -M 1 -o ceu_invariant_gene.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz 
wc -l ceu_invariant_gene.vcf
bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed  -S ibnr_invcf -i "F_MISSING < 0.8  || MAF < 0.0001" -O v -m 1 -M 1 -o ibnr_invariant_gene.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ibnr_invariant_gene.vcf

bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed -S ceu_invcf_shuf -i "F_MISSING < 0.8 ||  MAF > 0.02" -O v -m 2 -M 2 -q 0.02:minor -o ceu_variant_gene.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ceu_variant_gene.vcf
bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed -S ibnr_invcf -i "F_MISSING < 0.8 ||  MAF > 0.02" -O v -m 2 -M 2 -q 0.02:minor -o ibnr_variant_gene.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ibnr_variant_gene.vcf

bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed   -S ceu_invcf_shuf -i "F_MISSING < 0.8 || MAF < 0.0001" -O v -m 1 -M 1 -o ceu_invariant_intergenic.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz 
wc -l ceu_invariant_intergenic.vcf
bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed   -S ibnr_invcf -i "F_MISSING < 0.8  || MAF < 0.0001" -O v -m 1 -M 1 -o ibnr_invariant_intergenic.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ibnr_invariant_intergenic.vcf

bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed -S ceu_invcf_shuf -i "F_MISSING < 0.8 ||  MAF > 0.02" -O v -m 2 -M 2 -q 0.02:minor -o ceu_variant_intergenic.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ceu_variant_intergenic.vcf
bcftools view -R /data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed -S ibnr_invcf -i "F_MISSING < 0.8 ||  MAF > 0.02" -O v -m 2 -M 2 -q 0.02:minor -o ibnr_variant_intergenic.vcf 1001genomes_snp-short-indel_with_tair10_only_ACGTN_var_invar.vcf.gz
wc -l ibnr_variant_intergenic.vcf

```

18. Number of variant and invariant SMPs for CEU and IBNR

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_var_invar.vcf --bed  cg_sites.bed --recode --keep ceu_invcf_shuf --max-missing 0.8 --out ceu_meth_var_invar
cat vcfheader ceu_meth_var_invar_5mb.recode.vcf >ceu_meth_var_invar_5mb.vcf
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_invar --maf 0
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_var_invar.vcf  --bed  cg_sites.bed --recode --keep ibnr_invcf --max-missing 0.8 --out ibnr_meth_var_invar
cat vcfheader ibnr_meth_var_invar.recode.vcf >ibnr_meth_var_invar.vcf
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_invar --maf 0
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf arabidopsis_meth_var_invar.vcf --recode  --out test_b  --keep ibnr_invcf --max-missing 0.8
/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools  --vcf test_a.vcf  --out test_a_again --maf 0
/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf arabidopsis_meth_var_invar.vcf --recode  --out test_a  --keep ceu_invcf_shuf --max-missing 0.8
/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf test_b.recode.vcf --out test_b_again --maf 0
```

19. Subset SNP vcfs for CEU and IBNR
```
sed 's/X//' ceu_invcf_shuf > ceu_invcf_shuf2
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz  --keep ceu_invcf_shuf2 --recode  --max-missing 0.8  --out ceu_snps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_snps_wholegenome.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_snps_wholegenome.recode.vcf.gz

sed 's/X//' ibnr_invcf > ibnr_invcf2
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz  --keep ibnr_invcf2 --recode  --max-missing 0.8  --out ibnr_snps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_snps_wholegenome.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_snps_wholegenome.recode.vcf.gz


/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_snps_wholegenome.vcf ceu_snps_wholegenome.recode.vcf.gz ibnr_snps_wholegenome.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_ibnr_snps_wholegenome.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_ibnr_snps_wholegenome.vcf.gz
```

20. Calculate pi and Tajima's D per region for SNPs where region is either the cg interval, gene, or intergenic region

```

mkdir genes_fasta_ceu/
mkdir genes_fasta_ibnr/
mkdir genes_fasta_ceu_ibnr/
mkdir cg_ceu
mkdir cg_ibnr
mkdir intergenic_fasta_ceu/
mkdir intergenic_fasta_ibnr/

cd cg_ceu/

cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_snps_wholegenome.recode.vcf.gz  $line >$line.ceu_snps.vcf; done
wc -l *vcf >vcflengths_ceu_snps
zcat ../ceu_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
for file in *ceu_snps.vcf ; do cat vcfheader_ceu $file >${file/ceu_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ceu_invcf_shuf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript  combine2.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ceu | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a.
ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
ls *fa >fastafiles
Rscript pi_D.R

cd ../cg_ibnr/
cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_snps_wholegenome.recode.vcf.gz  $line >$line.ib
nr_snps.vcf; done
wc -l *vcf >vcflengths_ibnr_snps
zcat ../ibnr_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
for file in *ibnr_snps.vcf ; do cat vcfheader_ibnr $file >${file/ibnr_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ibnr_invcf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript  combine2.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ibnr | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a
.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
ls *fa >fastafiles
Rscript pi_D.R
cd ../

cd genes_fasta_ceu/

cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_snps_wholegenome.recode.vcf.gz  $line >$line.ceu_snps.vcf; done
wc -l *vcf >vcflengths_ceu_snps
zcat ../ceu_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
for file in *ceu_snps.vcf ; do cat vcfheader_ceu $file >${file/ceu_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ceu_invcf_shuf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript  combine2.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ceu | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
ls *fa >fastafiles
Rscript pi_D.R


cd ../genes_fasta_ibnr/
cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_snps_wholegenome.recode.vcf.gz  $line >$line.ibnr_snps.vcf; done
wc -l *vcf >vcflengths_ibnr_snps
zcat ../ibnr_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
for file in *ibnr_snps.vcf ; do cat vcfheader_ibnr $file >${file/ibnr_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ibnr_invcf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript  combine2.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ibnr | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
ls *fa >fastafiles
Rscript pi_D.R
cd ../

cd intergenic_fasta_ceu/

cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_snps_wholegenome.recode.vcf.gz  $line >$line.ceu
_snps.vcf; done
rm Pt*
rm Mt*
wc -l *vcf >vcflengths_ceu_snps
zcat ../ceu_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
for file in *ceu_snps.vcf ; do cat vcfheader_ceu $file >${file/ceu_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ceu_invcf_shuf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript combine.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ceu | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
rm *\:0-*.fa
ls *fa >fastafiles
Rscript pi_D.R

cd ../intergenic_fasta_ibnr/
cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_snps_wholegenome.recode.vcf.gz  $line >$line.ibnr_snps.vcf; done
rm Pt*
rm Mt*
wc -l *vcf >vcflengths_ibnr_snps
zcat ../ibnr_snps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
for file in *ibnr_snps.vcf ; do cat vcfheader_ibnr $file >${file/ibnr_snps.vcf/all.vcf}; done
Rscript vcfstats.R
chmod 777 bad_intervals.sh
./bad_intervals.sh
cp ../ibnr_invcf2 .
ls *.all.vcf | sed 's/.all.vcf//' >intervals
Rscript combine.R
for file in *all.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file ; done
for file in *all.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file ; done
cat intervals_ibnr | while read -r value1 value2 remainder  ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa $value2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p $value1.$value2 -H 1pIu $value2.all.vcf.gz >> $value2.fa ; done
rm *\:0-*.fa
ls *fa >fastafiles
Rscript pi_D.R
cd ../

```

21. Example vcfstats.R
```
vcflengths_var_invar <- read.table(file="vcflengths_ceu_snps")
vcflengths_var_invar <- vcflengths_var_invar[-c(nrow(vcflengths_var_invar)),]
vcflengths_var_invar$V2 <- gsub(".ceu_snps.vcf","",vcflengths_var_invar$V2)
colnames(vcflengths_var_invar) <- c("numvar","interval")
vcflengths_var_invar$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_var_invar$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_var_invar$interval)))
vcflengths_var_invar$prop <- vcflengths_var_invar$numvar/vcflengths_var_invar$length
vcflengths_var_invar <- vcflengths_var_invar[vcflengths_var_invar$numvar < 3,]
write.table(c(paste("rm ",vcflengths_var_invar$interval,".all.vcf",sep=""),paste("rm ",vcflengths_var_invar$interval,".ceu_snps.vcf",sep="")),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)

vcflengths_var_invar <- read.table(file="vcflengths_ceu_snps")
vcflengths_var_invar <- vcflengths_var_invar[-c(nrow(vcflengths_var_invar)),]
vcflengths_var_invar$V2 <- gsub(".ceu_snps.vcf","",vcflengths_var_invar$V2)
colnames(vcflengths_var_invar) <- c("numvar","interval")
vcflengths_var_invar$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_var_invar$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_var_invar$interval)))
vcflengths_var_invar$prop <- vcflengths_var_invar$numvar/vcflengths_var_invar$length
vcflengths_var_invar <- vcflengths_var_invar[vcflengths_var_invar$numvar > 2,]
write.table(cbind(paste(vcflengths_var_invar$interval,".all.vcf",sep=""),vcflengths_var_invar$length),file="goodfiles",sep="\t",quote=F,row.names = F, col.names = F)
```

22. Example combine2.R

```
file1 <- read.table(file="ceu_invcf_shuf2")
file2 <- read.table(file="intervals")
fileall <- expand.grid(file1$V1, file2$V1)  
write.table(fileall, file = "intervals_ceu", quote = F, sep = "\t", row.names = F, col.names = F)

```

23. Example pi_D.R

```
library(ape)
library(adegenet)
require(ape)
library(pegas)
library(dplyr)
library(ggplot2)
library(gridExtra)

fastafiles <- as.character(read.table(file="fastafiles")$V1)
pi_D <- ""

for (f in 1:length(fastafiles)){
  x <- as.DNAbin(read.dna(file=fastafiles[f],format="fasta",as.character = T))
  D_result <- tajima.test(x)
  Pi_result <-nuc.div(x,TRUE)
  pi_D <- rbind(pi_D,c(fastafiles[f],D_result$D,Pi_result[1]))
}
pi_D <- as.data.frame(pi_D[-c(1),])
colnames(pi_D) <- c("interval","D","pi")
write.table(pi_D,file="pi_D_ceu_cg.txt",quote=FALSE,row.names=FALSE,sep = "\t")
```

24. Subset CEU and IBNR SNP vcfs (repeated)

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz --keep ibnr_invcf2 --recode  --max-missing 0.8  --out ibnr_snps
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN_snps.vcf.gz --keep ceu_invcf_shuf2 --recode  --max-missing 0.8  --out ceu_snps


```

25. Calculate IBD block sizes using SNP data

```
cp ceu_snps.recode.vcf.gz ceu_snps.recode2.vcf.gz
gunzip ceu_snps.recode2.vcf.gz
cp ibnr_snps.recode.vcf.gz ibnr_snps.recode2.vcf.gz
gunzip ibnr_snps.recode2.vcf.gz

python3 /proj/popgen/a.ramesh/software/hmmIBD/vcf2hmm.py ceu_snps.recode2.vcf ceu_snps_hmmIBD_input
python3 /proj/popgen/a.ramesh/software/hmmIBD/vcf2hmm.py ibnr_snps.recode2.vcf ibnr_snps_hmmIBD_input
/proj/popgen/a.ramesh/software/hmmIBD/hmmIBD -i ceu_snps_hmmIBD_input_seq.txt -f ceu_snps_hmmIBD_input_freq.txt -o ceu_snps_hmmIBD_output
/proj/popgen/a.ramesh/software/hmmIBD/hmmIBD -i ibnr_snps_hmmIBD_input_seq.txt -f ibnr_snps_hmmIBD_input_freq.txt -o ibnr_snps_hmmIBD_output

```

26. Calculate allele frequencies and linkage disequlibrium using genic SNPs

```
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_snps.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_snps.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_snps.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_snps.recode.vcf.gz

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed --recode  --out ceu_snps_genes
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.bed --recode  --out ibnr_snps_genes

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_snps_genes.recode.vcf --out ceu_snps_genes --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_snps_genes.recode.vcf --out ibnr_snps_genes --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq

/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_snps_genes.recode.vcf -OutStat ceu_snps_5mb_cg_genes_decay  -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_snps_genes.recode.vcf -OutStat ibnr_snps_5mb_cg_genes_decay  -MaxDist 10

python3 ld_bin.py ceu_snps_5mb_cg_genes_decay.LD.gz ceu_snps_5mb_cg_genes_decay.stat 3
python3 ld_bin.py ibnr_snps_5mb_cg_genes_decay.LD.gz ibnr_snps_5mb_cg_genes_decay.stat 3

## also calculate LD for all SNPs

/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_snps.recode.vcf.gz -OutStat ceu_snps_ld_decay -MaxDist 10 -OutType 3
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_snps.recode.vcf.gz -OutStat ibnr_snps_ld_decay -MaxDist 10 -OutType 3

python3 ld_bin.py ceu_snps_ld_decay.LD.gz ceu_snps_ld_decay.stat 3
python3 ld_bin.py ibnr_snps_ld_decay.LD.gz ibnr_snps_ld_decay.stat 3

```

27. Script to calculate average LD every three bases
```
import sys
import pandas as pd

def process_file(input_file, output_file, bin_size):
    # Read input file (gzipped)
    df = pd.read_csv(input_file, sep="\t", compression='gzip')

    # Filter for Dist values between 1 and 2000
    filtered_df = df.loc[(df["Dist"] >= 1) & (df["Dist"] <= 2000)].copy()

    # Bin Dist by the specified bin size and calculate the average r^2 for each bin
    filtered_df["Dist_bin"] = (filtered_df["Dist"] // bin_size) * bin_size
    binned_avg = filtered_df.groupby("Dist_bin")["r^2"].mean().reset_index()

    # Write the result to an output file
    binned_avg.to_csv(output_file, sep="\t", index=False, header=True)

if __name__ == "__main__":
    # Ensure correct number of arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <output_file> <bin_size>")
        sys.exit(1)

    # Get input and output file paths and bin size from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bin_size = int(sys.argv[3])

    # Process the file
    process_file(input_file, output_file, bin_size)
```

28. Estimate genetic distances using genic SNPs

```
grep 'exon' Arabidopsis_thaliana.TAIR10.55.gff3 >gene_pos.bed
sed -e 's/\t/:/' -e  's/\t/-/' gene_pos.bed >gene_pos.list

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps_genes.recode.vcf.gz --out ceu_snps_genes_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps_genes.recode.vcf.gz --out ibnr_snps_genes_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_snps_genes_maf.recode.vcf 
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_snps_genes_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_snps_genes_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_snps_genes_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_snps_genes_maf.recode.vcf ceu_snps_genes_maf.recode.vcf.gz ibnr_snps_genes_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_snps_genes_maf.recode.vcf -OutPut ceu_ibnr_snps_genes.mat

```

29. Calculate allele frequencies and linkage disequlibrium using intergenic SNPs

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed --recode  --out ceu_snps_intergenic
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed --recode  --out ibnr_snps_intergenic

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_snps_intergenic.recode.vcf --out ceu_snps_intergenic --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_snps_intergenic.recode.vcf --out ibnr_snps_intergenic --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq

/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_snps_intergenic.recode.vcf -OutStat ceu_snps_5mb_cg_intergenic_decay  -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_snps_intergenic.recode.vcf -OutStat ibnr_snps_5mb_cg_intergenic_decay  -MaxDist 10

python3 ld_bin.py ceu_snps_5mb_cg_intergenic_decay.LD.gz ceu_snps_5mb_cg_intergenic_decay.stat 3
python3 ld_bin.py ibnr_snps_5mb_cg_intergenic_decay.LD.gz ibnr_snps_5mb_cg_intergenic_decay.stat 3
```

30. Estimate genetic distances using intergenic SNPs

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps_intergenic.recode.vcf.gz --out ceu_snps_intergenic_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps_intergenic.recode.vcf.gz --out ibnr_snps_intergenic_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_snps_intergenic_maf.recode.vcf 
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_snps_intergenic_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_snps_intergenic_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_snps_intergenic_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_snps_intergenic_maf.recode.vcf ceu_snps_intergenic_maf.recode.vcf.gz ibnr_snps_intergenic_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_snps_intergenic_maf.recode.vcf -OutPut ceu_ibnr_snps_intergenic.mat
```

31. Thin SNP VCFs for DAPC inference

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_snps_intergenic_maf_thin --vcf ceu_ibnr_snps_intergenic_maf.recode.vcf --thin 5000
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_snps_genes_maf_thin --vcf ceu_ibnr_snps_genes_maf.recode.vcf --thin 5000

```
32. Filter SMP vcfs for CEU and IBNR

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf arabidopsis_meth_all.vcf --keep ceu_invcf_shuf --recode  --max-missing 0.8  --out ceu_smps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf arabidopsis_meth_all.vcf --keep ibnr_invcf --recode  --max-missing 0.8  --out ibnr_smps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02
sed -i -e '/Pt_/d' -e '/Mt_/d' -e '/5_6000000/d' -e '/5_7000000/d' ceu_smps_wholegenome.recode.vcf
sed -i -e '/Pt_/d' -e '/Mt_/d' -e '/5_6000000/d' -e '/5_7000000/d' ibnr_smps_wholegenome.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_wholegenome.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_wholegenome.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_wholegenome.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_wholegenome.recode.vcf.gz
```
33. Generate individual SMP vcfs for each interval.
```
mkdir genes_fasta_ceu/
mkdir genes_fasta_ibnr/
mkdir cg_ceu
mkdir cg_ibnr

cd cg_ceu/
cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_smps_wholegenome.recode.vcf.gz  $line >$line.ceu_smps.vcf; done
wc -l *vcf >vcflengths_ceu_smps
zcat ../ceu_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
for file in *ceu_smps.vcf; do cat vcfheader_ceu $file >${file/ceu_smps.vcf/all.vcf}; done
wc -l *_smps.vcf >vcflengths_var_invar
Rscript good_intervals.R
grep 'CHROM' 1:5400-5600.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done 
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
tar -cvzf input_files.tar.gz *input.txt

cd ../cg_ibnr/
cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_smps_wholegenome.recode.vcf.gz  line >$line.ibnr_smps.vcf; done
wc -l *vcf >vcflengths_ibnr_smps
zcat ../ibnr_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
for file in *ibnr_smps.vcf; do cat vcfheader_ibnr $file >${file/ibnr_smps.vcf/all.vcf}; done
wc -l *_smps.vcf >vcflengths_var_invar
Rscript good_intervals.R
grep 'CHROM' 1:5400-5600.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
tar -cvzf input_files.tar.gz *input.txt

cd ../genes_fasta_ceu/
cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_smps_wholegenome.recode.vcf.gz  $line >$line.ceu_smps.vcf; done
wc -l *vcf >vcflengths_ceu_smps
zcat ../ceu_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
for file in *_smps.vcf; do cat vcfheader_ceu $file >${file/ceu_smps.vcf/all.vcf}; done
wc -l *all.vcf >vcflengths_var_invar
Rscript good_intervals.R
grep 'CHROM' 1:10003880-10004581.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done 
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
tar -cvzf input_files.tar.gz *input.txt

cd ../genes_fasta_ibnr/
cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_smps_wholegenome.recode.vcf.gz  $line >$line.ibnr_smps.vcf; done
wc -l *vcf >vcflengths_ibnr_smps
zcat ../ibnr_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
for file in *_smps.vcf; do cat vcfheader_ibnr $file >${file/ibnr_smps.vcf/all.vcf}; done
wc -l *all.vcf >vcflengths_var_invar
Rscript good_intervals.R
grep 'CHROM' 1:10003880-10004581.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
cd ../
```

34. Example good_intervals.R
```
```\vcflengths_var_invar <- read.table(file="vcflengths_var_invar")
vcflengths_var_invar <- vcflengths_var_invar[-c(nrow(vcflengths_var_invar)),]
vcflengths_var_invar$V2 <- gsub(".ceu_smps.vcf","",vcflengths_var_invar$V2)
colnames(vcflengths_var_invar) <- c("numvar","interval")
vcflengths_var_invar$interval <- gsub(".all.vcf","",vcflengths_var_invar$interval)
cytosine_count <- read.table(file="../cytosine_count_cg.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

library(dplyr)

merged <- inner_join(cytosine_count,vcflengths_var_invar,by="interval")
merged$prop <- merged$numvar/merged$numc
badintervals <- merged[merged$numvar < 3,]
merged <- merged[merged$numvar > 2,]

write.table(merged,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)
write.table(paste("rm ",badintervals$interval,".all.vcf"),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)

```

35. Estimate per region (gene region or clock like cg) Pi and Tajima's D for SMP.
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/cg_ceu
mkdir input/
for file in *input.txt; do sed -i '2d' $file; done
mv *input.txt input/
cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
Rscript good_intervals.R
mv input/* .
cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ceu.txt -length $value2 -alpha $value5  ; done
head -n 1 5:9986600-9986800.Dm_ceu.txt >Dm_results_header
cat *.Dm_ceu.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ceu_cg_Dm.txt
ls *.Dm_ceu.txt | sed 's/.Dm_ceu.txt//' >Dm_filenames_ceu_cg

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/cg_ibnr
mkdir input/
for file in *input.txt; do sed -i '2d' $file; done
mv *input.txt input/
cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
Rscript good_intervals.R
mv input/* .
cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ibnr.txt -length $value2 -alpha $value5  ; done
head -n 1 5:9986600-9986800.Dm_ibnr.txt >Dm_results_header
cat *.Dm_ibnr.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ibnr_cg_Dm.txt
ls *.Dm_ibnr.txt | sed 's/.Dm_ibnr.txt//' >Dm_filenames_ibnr_cg

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/genes_fasta_ceu
mkdir input/
for file in *input.txt; do sed -i '2d' $file; done
mv *input.txt input/
cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
Rscript good_intervals.R
mv input/* .
cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ceu.txt -length $value2 -alpha $value5  ; done
head -n 1 5:9991685-9992770.Dm_ceu.txt >Dm_results_header
cat *.Dm_ceu.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ceu_Dm.txt
ls *.Dm_ceu.txt | sed 's/.Dm_ceu.txt//' >Dm_filenames_ceu


cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/genes_fasta_ibnr
#mkdir input/
#for file in *input.txt; do sed -i '2d' $file; done
#mv *input.txt input/
#cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
#perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
#Rscript good_intervals.R
#mv input/* .
#cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ibnr.txt -length $value2 -alpha $value5  ; done
#head -n 1 5:9991685-9992770.Dm_ibnr.txt >Dm_results_header
#cat *.Dm_ibnr.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ibnr_Dm.txt
#ls *.Dm_ibnr.txt | sed 's/.Dm_ibnr.txt//' >Dm_filenames_ibnr

```

36. Example good_intervals.R from point 31.
```
library(dplyr)
good_intervals <- read.table(file="good_intervals")
alpha_dm <- read.table(file="alpha_Dm")
alpha_dm$V1 <- gsub(".input.txt","",alpha_dm$V1)
alpha_dm <- alpha_dm[1:2]
good_intervals <- inner_join(good_intervals,alpha_dm,by="V1")
write.table(good_intervals,file="good_intervals_alpha",sep="\t",quote=F,row.names = F, col.names = F)
```

37. Allele frequencies for SMPs and DMRs
```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ceu_meth_cg  --bed  cg_sites.bed  --keep ceu_invcf_shuf  --max-missing 0.8 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ibnr_meth_cg  --bed  cg_sites.bed   --keep ibnr_invcf  --max-missing 0.8 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ceu_all_genes  --keep ceu_invcf_shuf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ibnr_all_genes  --keep ibnr_invcf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ceu_non_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed --keep ceu_invcf_shuf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ibnr_non_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed --keep ibnr_invcf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_cg --recode --bed cg_sites.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_cg --recode --bed cg_sites.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_gbm --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_gbm.recode.vcf --out ceu_dmrs_gbm --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_gbm --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_gbm.recode.vcf --out ibnr_dmrs_gbm --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ceu_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed  --keep ceu_invcf_shuf  --max-missing 0.8 --recode 
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_sift.recode.vcf --out ceu_sift --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ibnr_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed  --keep ibnr_invcf  --max-missing 0.8 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_sift.recode.vcf --out ibnr_sift --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ceu_low_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed --keep ceu_invcf_shuf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf arabidopsis_meth_all.vcf  --out ibnr_low_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed --keep ibnr_invcf  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_sift --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_sift.recode.vcf --out ceu_dmrs_sift --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_sift --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_sift.recode.vcf --out ibnr_dmrs_sift --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_meth_cg.recode.vcf --out ceu_smps  --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_cg.recode.vcf --out ceu_dmrs --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_meth_cg.recode.vcf --out ibnr_smps --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_cg.recode.vcf --out ibnr_dmrs --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_all  --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_all  --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_non_gbm  --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_non_gbm --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_low_sift --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_low_sift --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed
```

38. Calculate LD using SMP and DMR
```
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_meth_cg.recode.vcf -OutStat ceu_meth_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_meth_cg.recode.vcf -OutStat ibnr_meth_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_dmrs_cg.recode.vcf -OutStat ceu_dmrs_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_dmrs_cg.recode.vcf -OutStat ibnr_dmrs_cg_ld_decay -MaxDist 10
```

39. Calculate distance matrices for SMPs and DMRs
```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_meth_cg.recode.vcf --out ceu_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_meth_cg.recode.vcf --out ibnr_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_meth_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_meth_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_meth_cg_maf.recode.vcf ceu_meth_cg_maf.recode.vcf.gz ibnr_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_meth_cg_maf.recode.vcf -OutPut ceu_ibnr_meth_cg_dis.mat

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_cg.recode.vcf --out ceu_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_cg.recode.vcf --out ibnr_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_dmrs_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_dmrs_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_dmrs_cg_maf.recode.vcf ceu_dmrs_cg_maf.recode.vcf.gz ibnr_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_dmrs_cg_maf.recode.vcf -OutPut ceu_ibnr_dmrs_cg.mat

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_ceu --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --gzvcf ceu_smps_wholegenome.recode.vcf.gz --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_ibnr --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --gzvcf ibnr_smps_wholegenome.recode.vcf.gz --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_ceu.recode.vcf --out gbm_ceu_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_ibnr.recode.vcf --out gbm_ibnr_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_ceu_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_ceu_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_ibnr_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_ibnr_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_meth_gbm.recode.vcf gbm_ceu_maf.recode.vcf.gz gbm_ibnr_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_meth_gbm.recode.vcf -OutPut ceu_ibnr_meth_gbm_dis.mat

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_dmr_ceu --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --vcf ceu_dmrs_dip.vcf --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_dmr_ibnr --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --vcf ibnr_dmrs_dip.vcf  --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_dmr_ceu.recode.vcf --out gbm_dmr_ceu_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_dmr_ibnr.recode.vcf --out gbm_dmr_ibnr_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_dmr_ceu_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_dmr_ceu_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_dmr_ibnr_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_dmr_ibnr_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_dmrs_gbm.recode.vcf gbm_dmr_ceu_maf.recode.vcf.gz gbm_dmr_ibnr_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_dmrs_gbm.recode.vcf -OutPut ceu_ibnr_dmrs_gbm.mat
```

40. Get multihetsep files for demographic inference for SNPs
```
cd snps_ceu/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps_wholegenome.recode.vcf.gz  --recode --out ceu_snps_chr1_smc --chr 1
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_snps_chr1_smc.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_snps_chr1_smc.recode.vcf.gz
cat ceu_invcf_shuf2 | while read line ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $line -o $line.filtered.vcf ceu_snps_chr1_smc.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat ceu_invcf_shuf2 | while read line ; do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $line.annotated.vcf.gz --out $line.1.snps --recode --recode-INFO-all  --remove-indels --max-missing 1 --chr 1 ; done
for file in *.snps.recode.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
/proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --chr 1 *.1.snps.recode.vcf.gz >ceu_multihetsep_1
Rscript multihetsep_combined.R
cd ../

cd snps_ibnr/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps_wholegenome.recode.vcf.gz  --recode --out ibnr_snps_chr1_smc --chr 1
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_snps_chr1_smc.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_snps_chr1_smc.recode.vcf.gz
cat ibnr_invcf2 | while read line ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $line -o $line.filtered.vcf ibnr_snps_chr1_smc.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat ibnr_invcf2 | while read line; do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $line.annotated.vcf.gz --out $line.1.snps --recode --recode-INFO-all  --remove-indels --max-missing 1 --chr 1 ; done
for file in *.snps.recode.vcf; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
/proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --chr 1 *.1.snps.recode.vcf.gz >ibnr_multihetsep_1
Rscript multihetsep_combined.R
cd ../
```
41. To create the file combining sample names and chr for the inference above. In R.

```
file1 <- read.table(file="ceu_invcf_shuf2")
file2 <- read.table(file="chrlist")

fileall <- expand.grid(file1$V1, file2$V1)  
write.table(fileall, file = "sample_chr_ceu", quote = F, sep = "\t", row.names = F, col.names = F)

file1 <- read.table(file="ibnr_invcf2")
file2 <- read.table(file="chrlist")

fileall <- expand.grid(file1$V1, file2$V1)
write.table(fileall, file = "sample_chr_ibnr", quote = F, sep = "\t", row.names = F, col.names = F)
```

42. Make the multihetsep files haploid. No hets in file. In R (multihetsep_combined.R).
```
snps <- read.table(file="ceu_multihetsep_1")
geno <- substring(snps$V4,1,1)
for (j in (seq(from = 3, to = 126, by = 2))){
  geno <- paste(geno,substring(snps$V4,j,j),sep="")
} 
snps$V4 <- geno
write.table(snps,file="ceu_multihetsep_snps_1",sep="\t",col.names = F, row.names = F, quote = F)

snps <- read.table(file="ibnr_multihetsep_1")
geno <- substring(snps$V4,1,1)
for (j in (seq(from = 3, to = 126, by = 2))){
  geno <- paste(geno,substring(snps$V4,j,j),sep="")
} 
snps$V4 <- geno
write.table(snps,file="ibnr_multihetsep_snps_1",sep="\t",col.names = F, row.names = F, quote = F)
```


43. Get list of cytosine sites. In R.

```
library("methimpute",lib.loc="/data/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")
cytosines <- extractCytosinesFromFASTA("Arabidopsis_thaliana.TAIR10.dna.toplevel.fa", contexts = 'CG')
cytosines <- as.data.frame(cytosines)
cytosines <- cytosines[1:3]
cytosines$seqnames <- gsub(" dna.*","",cytosines$seqnames)
cytosines$start <- cytosines$start - 1
options(scipen=10)
write.table(cytosines[cytosines$seqnames %in% "1",],file="1.cytosines_tair10.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cytosines[cytosines$seqnames %in% "2",],file="2.cytosines_tair10.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cytosines[cytosines$seqnames %in% "3",],file="3.cytosines_tair10.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cytosines[cytosines$seqnames %in% "4",],file="4.cytosines_tair10.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cytosines[cytosines$seqnames %in% "5",],file="5.cytosines_tair10.bed",sep="\t",quote=F,row.names=F,col.names=F)
```


44. Get multihetsep files for demographic inference for SMPs
```
cd smp_ceu/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_smps_wholegenome.recode.vcf.gz --recode  --out ceu_smps_5mb --bed popgen5mb2.bed
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_5mb.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_5mb.recode.vcf.gz
cat ceu_invcf_shuf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ceu_smps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat sample_chr_ceu |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ceu_multihetsep_meth_$line ; done
Rscript multihetsep_combined.R
cd ../

cd smp_ibnr/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_smps_wholegenome.recode.vcf.gz --recode  --out ibnr_smps_5mb --bed popgen5mb2.bed
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_5mb.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_5mb.recode.vcf.gz
cat ibnr_invcf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ibnr_smps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat sample_chr_ibnr |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ibnr_multihetsep_meth_$line ; done
Rscript multihetsep_combined.R
cd ../
```

45. Make the multihetsep files haploid. No hets in file. In R. (multihetsep_combined.R).
```
popgen5mb <- read.table(file="../popgen5mb2.bed",header = T)
for (i in 1:5){
  smps <- read.table(paste("ceu_multihetsep_meth_",i,sep=""))
  smps$V2 <- smps$V2 - popgen5mb$start[i]
  smps$V3[1] <- 1
  geno <- substring(smps$V4,1,1)
  for (j in (seq(from = 3, to = 126, by = 2))){
    geno <- paste(geno,substring(smps$V4,j,j),sep="")
  } 
  smps$V4 <- geno
  write.table(smps,paste("ceu_multihetsep_smps_5mb_",i,sep=""),sep="\t",col.names = F, row.names = F, quote = F)
}

popgen5mb <- read.table(file="../popgen5mb2.bed",header = T)
for (i in 1:5){
  smps <- read.table(paste("ibnr_multihetsep_meth_",i,sep=""))
  smps$V2 <- smps$V2 - popgen5mb$start[i]
  smps$V3[1] <- 1
  geno <- substring(smps$V4,1,1)
  for (j in (seq(from = 3, to = 126, by = 2))){
    geno <- paste(geno,substring(smps$V4,j,j),sep="")
  } 
  smps$V4 <- geno
  write.table(smps,paste("ibnr_multihetsep_smps_5mb_",i,sep=""),sep="\t",col.names = F, row.names = F, quote = F)
}
```

46. Get multihetsep files for demographic inference for 5mb regions for SNPs

```
cd ceu_smcm/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps_wholegenome.recode.vcf.gz  --recode  --out ceu_snps_5mb --bed popgen5mb2.bed
/data/proj2/popgen/a.ramesh/software/htslib/bgzip -f ceu_snps_5mb.recode.vcf 
/data/proj2/popgen/a.ramesh/software/htslib/tabix -f ceu_snps_5mb.recode.vcf.gz
cat ceu_invcf_shuf |  while read -r  sample remainder; do /data/proj2/popgen/a.ramesh/software/bcftools/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ceu_snps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /data/proj2/popgen/a.ramesh/software/bcftools/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /data/proj2/popgen/a.ramesh/software/htslib/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /data/proj2/popgen/a.ramesh/software/htslib/tabix -f $file; done
cat sample_chr_ceu |  while read -r value1 value2 remainder ;  do /data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /data/proj2/popgen/a.ramesh/software/htslib/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /data/proj2/popgen/a.ramesh/software/htslib/tabix -f $file; done
cat chrlist | while read line; do /data/proj2/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --chr $line *.$line.snps.recode.vcf.gz >ceu_multihetsep_snp_$line ; done
Rscript multihetsep_combined.R
cd ../

cd ibnr_smcm/
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps_wholegenome.recode.vcf.gz --recode  --out ibnr_snps_5mb --bed popgen5mb2.bed
/data/proj2/popgen/a.ramesh/software/htslib/bgzip -f ibnr_snps_5mb.recode.vcf 
/data/proj2/popgen/a.ramesh/software/htslib/tabix -f ibnr_snps_5mb.recode.vcf.gz
cat ibnr_invcf |  while read -r  sample remainder; do /data/proj2/popgen/a.ramesh/software/bcftools/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ibnr_snps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /data/proj2/popgen/a.ramesh/software/bcftools/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /data/proj2/popgen/a.ramesh/software/htslib/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /data/proj2/popgen/a.ramesh/software/htslib/tabix -f $file; done
cat sample_chr_ibnr |  while read -r value1 value2 remainder ;  do /data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /data/proj2/popgen/a.ramesh/software/htslib/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /data/proj2/popgen/a.ramesh/software/htslib/tabix -f $file; done
cat chrlist | while read line; do /data/proj2/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --chr $line *.$line.snps.recode.vcf.gz >ibnr_multihetsep_snp_$line ; done
Rscript multihetsep_combined.R
cd ../
```

47. Make the multihetsep files haploid. For SNPs in 5mb No hets in file. In R. (multihetsep_combined.R). Then combine with methylation 5Mb files.

```
popgen5mb <- read.table(file="../popgen5mb2.bed",header = T)
for (i in 1:5){
  snps <- read.table(paste("ceu_multihetsep_snp_",i,sep=""))
  snps$V2 <- snps$V2 - popgen5mb$start[i]
  snps$V3[1] <- 1
  geno <- substring(snps$V4,1,1)
  for (j in (seq(from = 3, to = 126, by = 2))){
    geno <- paste(geno,substring(snps$V4,j,j),sep="")
  } 
  snps$V4 <- geno
  write.table(snps,paste("ceu_multihetsep_snps_5mb_",i,sep=""),sep="\t",col.names = F, row.names = F, quote = F)
}

samples <- c(56,52,55,8,29,57,46,4,23,32,61,40,49,3,16,26,28,53,5,15)
extract_letters <- function(string, positions) {
  paste0(sapply(positions, function(pos) substr(string, pos, pos)), collapse = "")
}
for (i in 1:5){
  snps <- read.table(paste("ceu_multihetsep_snps_5mb_",i,sep=""))
  smps <- read.table(paste("ceu_multihetsep_smps_5mb_",i,sep=""))
  smps$V4 <- gsub("A","D",smps$V4)
  smps$V4 <- gsub("T","M",smps$V4)
  both <- rbind(snps,smps)
  both <- both[order(both$V1,both$V2),]
  both$V3 <- c(both$V2[1]-1,diff(both$V2))
  both$V4 <- sapply(both$V4, extract_letters, positions = samples)
  write.table(both,paste("ceu_multihetsep_snp_meth_20_",i,sep=""),sep=" ",col.names = F, row.names = F, quote = F)
}

popgen5mb <- read.table(file="../popgen5mb2.bed",header = T)
for (i in 1:5){
  snps <- read.table(paste("ibnr_multihetsep_snp_",i,sep=""))
  snps$V2 <- snps$V2 - popgen5mb$start[i]
  snps$V3[1] <- 1
  geno <- substring(snps$V4,1,1)
  for (j in (seq(from = 3, to = 126, by = 2))){
    geno <- paste(geno,substring(snps$V4,j,j),sep="")
  } 
  snps$V4 <- geno
  write.table(snps,paste("ibnr_multihetsep_snps_5mb_",i,sep=""),sep="\t",col.names = F, row.names = F, quote = F)
}

samples <- c(56,52,55,8,29,57,46,4,23,32,61,40,49,3,16,26,28,53,5,15)
extract_letters <- function(string, positions) {
  paste0(sapply(positions, function(pos) substr(string, pos, pos)), collapse = "")
}
for (i in 1:5){
  snps <- read.table(paste("ibnr_multihetsep_snps_5mb_",i,sep=""))
  smps <- read.table(paste("ibnr_multihetsep_smps_5mb_",i,sep=""))
  smps$V4 <- gsub("A","D",smps$V4)
  smps$V4 <- gsub("T","M",smps$V4)
  both <- rbind(snps,smps)
  both <- both[order(both$V1,both$V2),]
  both$V3 <- c(both$V2[1]-1,diff(both$V2))
  both$V4 <- sapply(both$V4, extract_letters, positions = samples)
  write.table(both,paste("ibnr_multihetsep_snp_meth_20_",i,sep=""),sep=" ",col.names = F, row.names = F, quote = F)
}

```
48. Estimate methylation rates using SMCm 

```
###################
# Source function #
###################
library(eSMC2,lib.loc="/data/proj2/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")
################################

########
#Script#
########
nsim=5
M=20
M_s=M
mu <- 7*10^-9
r=3.5*10^-8
rho=r/mu
SB=F
SF=F
ER=T
rate_m=2.56*10^(-4)
rate_d=6.30*10^(-4)
model=1
results_site=list()
for(x in 1:5){
  M=M_s
  O=as.matrix(Get_real_data("/data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/snps",M=M,paste("ceu_multihetsep_snp_meth_20_",x,sep ="" )))
  M=dim(O)[1]-2
# Estimating rates
  Estimation_results=Methylation_rate_estimation(n=40,rho= rho,methylation=c(abs(log10(rate_m)),abs(log10(rate_d))),O,mu_r=(mu),NC=1,mu_b=1,Region=F,Free = T)
  print(c("Meth_rate",Estimation_results$mu_m,x,model))
  write.table(paste(Estimation_results$mu_m,collapse = " "),file=paste(x,model,"meth_rate_ceu.txt",sep="_"),row.names = F,col.names = F,sep="\t",quote = F)
}


###################
# Source function #
###################
library(eSMC2,lib.loc="/data/proj2/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.1/")
################################

########
#Script#
########
nsim=5
M=20
M_s=M
mu <- 7*10^-9
r=3.5*10^-8
rho=r/mu
SB=F
SF=F
ER=T
rate_m=2.56*10^(-4)
rate_d=6.30*10^(-4)
model=1
results_site=list()
for(x in 1:5){
  M=M_s
  O=as.matrix(Get_real_data("/data/proj2/popgen/a.ramesh/projects/methylomes/arabidopsis/snps",M=M,paste("ibnr_multihetsep_snp_meth_20_",x,sep ="" )))
  M=dim(O)[1]-2
# Estimating rates
  Estimation_results=Methylation_rate_estimation(n=40,rho= rho,methylation=c(abs(log10(rate_m)),abs(log10(rate_d))),O,mu_r=(mu),NC=1,mu_b=1,Region=F,Free = T)
  print(c("Meth_rate",Estimation_results$mu_m,x,model))
  write.table(paste(Estimation_results$mu_m,collapse = " "),file=paste(x,model,"meth_rate_ibnr.txt",sep="_"),row.names = F,col.names = F,sep="\t",quote = F)
}

## in bash, combine results
#cat *_meth_rate_ceu.txt >model1_meth_rate_ceu.txt
#cat *_meth_rate_ibnr.txt >model1_meth_rate_ibnr.txt

```

49. Calculate site frequency spectra for transposable elements.

```
Rscript tevcf.R
Rscript te_create_bed.R
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip arabidopsis_te.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix arabidopsis_te.vcf.gz

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_te.vcf.gz --keep ceu_invcf_shuf --recode  --max-missing 0.8  --out ceu_te --min-alleles 2 --max-alleles 2 --maf 0.02
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_te.vcf.gz --keep ibnr_invcf --recode  --max-missing 0.8  --out ibnr_te --min-alleles 2 --max-alleles 2 --maf 0.02

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_te.recode.vcf --freq --out ceu_te
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_te.recode.vcf --freq --out ibnr_te
```

50. tevcf.R script.

```
cov_context3 <- read.table(file="13059_2021_2348_MOESM2_ESM.bed",header=T)
meta <- cov_context3[1:3]
meta$start <- round(((meta$start+meta$stop)/2),digits=0)
meta$stop <- paste(meta$Chr,meta$start,sep="_")
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"
cov_context3 <- cov_context3[-c(1:6)]
cov_context3[is.na(cov_context3)] <- "."
cov_context3 <- cbind(meta,cov_context3)
cov_context3$`#CHROM` <- as.numeric(cov_context3$`#CHROM`)
cov_context3$POS <- as.numeric(cov_context3$POS)
cov_context3 <- cov_context3[order(cov_context3$`#CHROM`,cov_context3$POS),]
write.table(cov_context3,file="arabidopsis_te.vcf",quote = F, row.names = F,sep="\t")

```

51. te_create_bed.R script.

```
tefile <- read.table(file="13059_2021_2348_MOESM2_ESM.bed",header=T)

ceu_invcf_shuf <- read.table(file="ceu_invcf_shuf")
ibnr_invcf <- read.table(file="ibnr_invcf")
meta <- tefile[1:6]

ceu_dat <- tefile[7:ncol(tefile)]
ceu_dat <- ceu_dat[colnames(ceu_dat) %in% ceu_invcf_shuf$V1]
ceu_na <- rowSums(is.na(ceu_dat))/ncol(ceu_dat)
ceu_maf <- rowSums(ceu_dat == 1, na.rm = T) / ((rowSums(ceu_dat == 1, na.rm = T)) + rowSums(ceu_dat == 0, na.rm = T))
ceu_dat <- cbind(meta,ceu_dat)
ceu_dat <- ceu_dat[which(ceu_maf >= 0.02 & ceu_maf <= 0.98 & ceu_na <= 0.2),]
write.table(ceu_dat[1:3],file="ceu_te.bed",row.names = F, col.names = F, quote = F, sep="\t")


ibnr_dat <- tefile[7:ncol(tefile)]
ibnr_dat <- ibnr_dat[colnames(ibnr_dat) %in% ibnr_invcf$V1]
ibnr_na <- rowSums(is.na(ibnr_dat))/ncol(ibnr_dat)
ibnr_maf <- rowSums(ibnr_dat == 1, na.rm = T) / ((rowSums(ibnr_dat == 1, na.rm = T)) + rowSums(ibnr_dat == 0, na.rm = T))
ibnr_dat <- cbind(meta,ibnr_dat)
ibnr_dat <- ibnr_dat[which(ibnr_maf >= 0.02 & ibnr_maf <= 0.98 & ibnr_na <= 0.2),]
write.table(ibnr_dat[1:3],file="ibnr_te.bed",row.names = F, col.names = F, quote = F, sep="\t")
```

Appendix 1a. Repeating entire analysis and intermediate files were deleted. Only focusing on CEU and IBnr. Logic same as before but these were the results that were eventually used.
```
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
sed -e 's/X/ /' -e 's/$/_/' ceu | grep -F -f - rename_se2.sh | cut -d ' ' -f 2 | sed 's/_trim.fq.gz//' >>ceu_sra_files
sed -e 's/X/ /' -e 's/$/_/' ceu | grep -F -f - rename_pe2.sh | cut -d ' ' -f 2 | sed 's/_.*//' >>ceu_sra_files
sed -e 's/X/ /' -e 's/$/_/' ceu | grep -F -f - /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/scripts/rename_se.sh | sed -e 's/cat //' -e 's/>.*//' -e 's/_trim.fq.gz//g' | tr -s ' ' '\n' >>ceu_sra_files
sed -e 's/X/ /' -e 's/$/_/' ceu | grep -F -f - /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/scripts/rename_pe.sh | sed -e 's/cat //' -e 's/>.*//' -e 's/_..paired.fq.gz//g' | tr -s ' ' '\n' >>ceu_sra_files


sed -e 's/X/ /' -e 's/$/_/' ibnr | grep -F -f - rename_se2.sh | cut -d ' ' -f 2 | sed 's/_trim.fq.gz//' >>ibnr_sra_files
sed -e 's/X/ /' -e 's/$/_/' ibnr | grep -F -f - rename_pe2.sh | cut -d ' ' -f 2 | sed 's/_.*//' >>ibnr_sra_files
sed -e 's/X/ /' -e 's/$/_/' ibnr | grep -F -f - /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/scripts/rename_se.sh | sed -e 's/cat //' -e 's/>.*//' -e 's/_trim.fq.gz//g' | tr -s ' ' '\n' >>ibnr_sra_files
sed -e 's/X/ /' -e 's/$/_/' ibnr | grep -F -f - /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/scripts/rename_pe.sh | sed -e 's/cat //' -e 's/>.*//' -e 's/_..paired.fq.gz//g' | tr -s ' ' '\n' >>ibnr_sra_files


/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file ceu_sra_files -O ceu_data
/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file ibnr_sra_files -O ibnr_data
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done
for file in *.sralite; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done
cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --gzip --split-3  $file; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
mkdir paired
mv *_1.fastq.gz paired/
mv *_2.fastq.gz paired/
cd paired/
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
cd ../
for file in *.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 20 $file ${file/.fastq.gz/_trim.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done


cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
mkdir paired
mv *_1.fastq.gz paired/
mv *_2.fastq.gz paired/
cd paired/
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
cd ../
for file in *.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 20 $file ${file/.fastq.gz/_trim.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done


mv paired/*_1.paired.fq.gz .
mv paired/*_2.paired.fq.gz .
for file in *_trim.fq.gz; do mv $file ${file/_trim.fq.gz/.trim.fq.gz}; done
for file in *_1.paired.fq.gz; do mv $file ${file/_1.paired.fq.gz/.1.paired.fq.gz}; done
for file in *_2.paired.fq.gz; do mv $file ${file/_2.paired.fq.gz/.2.paired.fq.gz}; done



cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
mv SRR771710_trim.fq.gz 7120_trim.fq.gz
mv SRR771735_trim.fq.gz 7177_trim.fq.gz
mv SRR771745_trim.fq.gz 7203_trim.fq.gz
mv SRR771746_trim.fq.gz 7207_trim.fq.gz
mv SRR771755_trim.fq.gz 7520_trim.fq.gz
mv SRR771822_trim.fq.gz 6984_trim.fq.gz
mv SRR4295448_trim.fq.gz 9731_trim.fq.gz
mv SRR4295449_trim.fq.gz 9735_trim.fq.gz
mv SRR3299958_trim.fq.gz 5921_trim.fq.gz
mv SRR3299964_trim.fq.gz 5984_trim.fq.gz
mv SRR3300341_trim.fq.gz 8236_trim.fq.gz
mv SRR3300365_trim.fq.gz 8285_trim.fq.gz
mv SRR3300383_trim.fq.gz 8365_trim.fq.gz
mv SRR3300957_trim.fq.gz 9694_trim.fq.gz
mv SRR3301449_trim.fq.gz 9914_trim.fq.gz
mv SRR4295311_1.paired.fq.gz 5874_1.paired.fq.gz
mv SRR4295312_1.paired.fq.gz 6008_1.paired.fq.gz
mv SRR4295321_1.paired.fq.gz 6396_1.paired.fq.gz
mv SRR4295327_1.paired.fq.gz 6957_1.paired.fq.gz
mv SRR4295335_1.paired.fq.gz 6976_1.paired.fq.gz
mv SRR4295350_1.paired.fq.gz 7296_1.paired.fq.gz
mv SRR4295377_1.paired.fq.gz 8235_1.paired.fq.gz
mv SRR4295380_1.paired.fq.gz 8284_1.paired.fq.gz
mv SRR3299872_1.paired.fq.gz 410_1.paired.fq.gz
mv SRR3299876_1.paired.fq.gz 424_1.paired.fq.gz
mv SRR3299947_1.paired.fq.gz 5890_1.paired.fq.gz
mv SRR3299951_1.paired.fq.gz 5893_1.paired.fq.gz
mv SRR3299957_1.paired.fq.gz 5921_1.paired.fq.gz
mv SRR3299959_1.paired.fq.gz 5950_1.paired.fq.gz
mv SRR3299963_1.paired.fq.gz 5984_1.paired.fq.gz
mv SRR3299965_1.paired.fq.gz 5993_1.paired.fq.gz
mv SRR3300062_1.paired.fq.gz 6390_1.paired.fq.gz
mv SRR3300340_1.paired.fq.gz 8236_1.paired.fq.gz
mv SRR3300364_1.paired.fq.gz 8285_1.paired.fq.gz
mv SRR3300366_1.paired.fq.gz 8290_1.paired.fq.gz
mv SRR3300382_1.paired.fq.gz 8365_1.paired.fq.gz
mv SRR3300956_1.paired.fq.gz 9694_1.paired.fq.gz
mv SRR3301448_1.paired.fq.gz 9914_1.paired.fq.gz
mv SRR3301450_1.paired.fq.gz 9915_1.paired.fq.gz
mv SRR4295311_2.paired.fq.gz 5874_2.paired.fq.gz
mv SRR4295312_2.paired.fq.gz 6008_2.paired.fq.gz
mv SRR4295321_2.paired.fq.gz 6396_2.paired.fq.gz
mv SRR4295327_2.paired.fq.gz 6957_2.paired.fq.gz
mv SRR4295335_2.paired.fq.gz 6976_2.paired.fq.gz
mv SRR4295350_2.paired.fq.gz 7296_2.paired.fq.gz
mv SRR4295377_2.paired.fq.gz 8235_2.paired.fq.gz
mv SRR4295380_2.paired.fq.gz 8284_2.paired.fq.gz
mv SRR3299872_2.paired.fq.gz 410_2.paired.fq.gz
mv SRR3299876_2.paired.fq.gz 424_2.paired.fq.gz
mv SRR3299947_2.paired.fq.gz 5890_2.paired.fq.gz
mv SRR3299951_2.paired.fq.gz 5893_2.paired.fq.gz
mv SRR3299957_2.paired.fq.gz 5921_2.paired.fq.gz
mv SRR3299959_2.paired.fq.gz 5950_2.paired.fq.gz
mv SRR3299963_2.paired.fq.gz 5984_2.paired.fq.gz
mv SRR3299965_2.paired.fq.gz 5993_2.paired.fq.gz
mv SRR3300062_2.paired.fq.gz 6390_2.paired.fq.gz
mv SRR3300340_2.paired.fq.gz 8236_2.paired.fq.gz
mv SRR3300364_2.paired.fq.gz 8285_2.paired.fq.gz
mv SRR3300366_2.paired.fq.gz 8290_2.paired.fq.gz
mv SRR3300382_2.paired.fq.gz 8365_2.paired.fq.gz
mv SRR3300956_2.paired.fq.gz 9694_2.paired.fq.gz
mv SRR3301448_2.paired.fq.gz 9914_2.paired.fq.gz
mv SRR3301450_2.paired.fq.gz 9915_2.paired.fq.gz

cat SRR771683_trim.fq.gz SRR771684_trim.fq.gz > 6903_trim.fq.gz
cat SRR771722_trim.fq.gz SRR3299880_trim.fq.gz > 430_trim.fq.gz
cat SRR771776_trim.fq.gz SRR771777_trim.fq.gz > 6951_trim.fq.gz
cat SRR4295282_trim.fq.gz SRR4295283_trim.fq.gz > 428_trim.fq.gz
cat SRR4295319_trim.fq.gz SRR4295320_trim.fq.gz > 6296_trim.fq.gz
cat SRR4295433_trim.fq.gz SRR4295434_trim.fq.gz > 9665_trim.fq.gz
cat SRR4295435_trim.fq.gz SRR4295436_trim.fq.gz > 9689_trim.fq.gz
cat SRR3299808_trim.fq.gz SRR3299809_trim.fq.gz SRR3299810_trim.fq.gz > 10020_trim.fq.gz
cat SRR3299870_trim.fq.gz SRR3299871_trim.fq.gz > 403_trim.fq.gz
cat SRR3299873_trim.fq.gz SRR3299874_trim.fq.gz SRR3299875_trim.fq.gz > 410_trim.fq.gz
cat SRR3299877_trim.fq.gz SRR3299878_trim.fq.gz SRR3299879_trim.fq.gz > 424_trim.fq.gz
cat SRR3299938_trim.fq.gz SRR3299939_trim.fq.gz SRR3299940_trim.fq.gz SRR3299941_trim.fq.gz SRR3299942_trim.fq.gz SRR3299943_trim.fq.gz > 5837_trim.fq.gz
cat SRR3299948_trim.fq.gz SRR3299949_trim.fq.gz SRR3299950_trim.fq.gz > 5890_trim.fq.gz
cat SRR3299952_trim.fq.gz SRR3299953_trim.fq.gz SRR3299954_trim.fq.gz > 5893_trim.fq.gz
cat SRR3299955_trim.fq.gz SRR3299956_trim.fq.gz > 5907_trim.fq.gz
cat SRR3299960_trim.fq.gz SRR3299961_trim.fq.gz SRR3299962_trim.fq.gz > 5950_trim.fq.gz
cat SRR3299966_trim.fq.gz SRR3299967_trim.fq.gz > 5993_trim.fq.gz
cat SRR3300063_trim.fq.gz SRR3300064_trim.fq.gz SRR3300065_trim.fq.gz > 6390_trim.fq.gz
cat SRR3300068_trim.fq.gz SRR3300069_trim.fq.gz SRR3300070_trim.fq.gz SRR3300071_trim.fq.gz > 6445_trim.fq.gz
cat SRR3300105_trim.fq.gz SRR3300106_trim.fq.gz > 6956_trim.fq.gz
cat SRR3300148_trim.fq.gz SRR3300149_trim.fq.gz > 7067_trim.fq.gz
cat SRR3300159_trim.fq.gz SRR3300160_trim.fq.gz SRR3300161_trim.fq.gz > 7103_trim.fq.gz
cat SRR3300253_trim.fq.gz SRR3300254_trim.fq.gz > 7350_trim.fq.gz
cat SRR3300291_trim.fq.gz SRR3300292_trim.fq.gz > 7424_trim.fq.gz
cat SRR3300367_trim.fq.gz SRR3300368_trim.fq.gz > 8290_trim.fq.gz
cat SRR3300386_trim.fq.gz SRR3300387_trim.fq.gz > 8386_trim.fq.gz
cat SRR3300934_trim.fq.gz SRR3300935_trim.fq.gz SRR3300936_trim.fq.gz > 9668_trim.fq.gz
cat SRR3300937_trim.fq.gz SRR3300938_trim.fq.gz SRR3300939_trim.fq.gz > 9671_trim.fq.gz
cat SRR3300940_trim.fq.gz SRR3300941_trim.fq.gz SRR3300942_trim.fq.gz > 9676_trim.fq.gz
cat SRR3300943_trim.fq.gz SRR3300944_trim.fq.gz > 9678_trim.fq.gz
cat SRR3300945_trim.fq.gz SRR3300946_trim.fq.gz SRR3300947_trim.fq.gz > 9679_trim.fq.gz
cat SRR3300948_trim.fq.gz SRR3300949_trim.fq.gz > 9684_trim.fq.gz
cat SRR3300950_trim.fq.gz SRR3300951_trim.fq.gz SRR3300952_trim.fq.gz > 9690_trim.fq.gz
cat SRR3300953_trim.fq.gz SRR3300954_trim.fq.gz SRR3300955_trim.fq.gz > 9693_trim.fq.gz
cat SRR3300958_trim.fq.gz SRR3300959_trim.fq.gz SRR3300960_trim.fq.gz > 9696_trim.fq.gz
cat SRR3301041_trim.fq.gz SRR3301042_trim.fq.gz > 9727_trim.fq.gz
cat SRR3301043_trim.fq.gz SRR3301044_trim.fq.gz > 9728_trim.fq.gz
cat SRR3301045_trim.fq.gz SRR3301046_trim.fq.gz SRR3301047_trim.fq.gz SRR3301048_trim.fq.gz > 9729_trim.fq.gz
cat SRR3301049_trim.fq.gz SRR3301050_trim.fq.gz > 9730_trim.fq.gz
cat SRR3301051_trim.fq.gz SRR3301052_trim.fq.gz SRR3301053_trim.fq.gz SRR3301054_trim.fq.gz > 9732_trim.fq.gz
cat SRR3301055_trim.fq.gz SRR3301056_trim.fq.gz > 9733_trim.fq.gz
cat SRR3301451_trim.fq.gz SRR3301452_trim.fq.gz > 9915_trim.fq.gz
cat SRR3301564_trim.fq.gz SRR3301565_trim.fq.gz SRR3301566_trim.fq.gz > 9973_trim.fq.gz

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data

mv SRR4295325_1.paired.fq.gz 6933_1.paired.fq.gz
mv SRR3300113_1.paired.fq.gz 6971_1.paired.fq.gz
mv SRR3300236_1.paired.fq.gz 7327_1.paired.fq.gz
mv SRR4295325_2.paired.fq.gz 6933_2.paired.fq.gz
mv SRR3300113_2.paired.fq.gz 6971_2.paired.fq.gz
mv SRR3300236_2.paired.fq.gz 7327_2.paired.fq.gz

#cat SRR3300110_trim.fq.gz SRR3300111_trim.fq.gz SRR3300112_trim.fq.gz > 6970_trim.fq.gz
#cat SRR3300114_trim.fq.gz SRR3300115_trim.fq.gz SRR3300116_trim.fq.gz > 6971_trim.fq.gz
#cat SRR3300237_trim.fq.gz SRR3300238_trim.fq.gz SRR3300239_trim.fq.gz > 7327_trim.fq.gz
#cat SRR3300240_trim.fq.gz SRR3300241_trim.fq.gz SRR3300242_trim.fq.gz SRR3300243_trim.fq.gz > 7328_trim.fq.gz
#cat SRR3300512_trim.fq.gz SRR3300513_trim.fq.gz SRR3300514_trim.fq.gz > 9507_trim.fq.gz
#cat SRR3300521_trim.fq.gz SRR3300522_trim.fq.gz SRR3300523_trim.fq.gz > 9510_trim.fq.gz
#cat SRR3300524_trim.fq.gz SRR3300525_trim.fq.gz SRR3300526_trim.fq.gz > 9511_trim.fq.gz
#cat SRR3300533_trim.fq.gz SRR3300534_trim.fq.gz SRR3300535_trim.fq.gz > 9514_trim.fq.gz
#cat SRR3300536_trim.fq.gz SRR3300537_trim.fq.gz SRR3300538_trim.fq.gz > 9515_trim.fq.gz
#cat SRR3300553_trim.fq.gz SRR3300554_trim.fq.gz SRR3300555_trim.fq.gz > 9521_trim.fq.gz
#cat SRR3300556_trim.fq.gz SRR3300557_trim.fq.gz SRR3300558_trim.fq.gz > 9522_trim.fq.gz
#cat SRR3300562_trim.fq.gz SRR3300563_trim.fq.gz SRR3300564_trim.fq.gz > 9524_trim.fq.gz
#cat SRR3300565_trim.fq.gz SRR3300566_trim.fq.gz SRR3300567_trim.fq.gz > 9525_trim.fq.gz
#cat SRR3300589_trim.fq.gz SRR3300590_trim.fq.gz SRR3300591_trim.fq.gz > 9534_trim.fq.gz
#cat SRR3300592_trim.fq.gz SRR3300593_trim.fq.gz SRR3300594_trim.fq.gz > 9535_trim.fq.gz
#cat SRR3300595_trim.fq.gz SRR3300596_trim.fq.gz SRR3300597_trim.fq.gz > 9537_trim.fq.gz
#cat SRR3300601_trim.fq.gz SRR3300602_trim.fq.gz > 9540_trim.fq.gz
#cat SRR3300603_trim.fq.gz SRR3300604_trim.fq.gz SRR3300605_trim.fq.gz > 9541_trim.fq.gz
#cat SRR3300611_trim.fq.gz SRR3300612_trim.fq.gz SRR3300613_trim.fq.gz > 9544_trim.fq.gz
#cat SRR3300620_trim.fq.gz SRR3300621_trim.fq.gz SRR3300622_trim.fq.gz > 9547_trim.fq.gz
#cat SRR3300644_trim.fq.gz SRR3300645_trim.fq.gz SRR3300646_trim.fq.gz > 9556_trim.fq.gz
#cat SRR3300647_trim.fq.gz SRR3300648_trim.fq.gz SRR3300649_trim.fq.gz > 9557_trim.fq.gz
#cat SRR3300656_trim.fq.gz SRR3300657_trim.fq.gz SRR3300658_trim.fq.gz > 9560_trim.fq.gz
#cat SRR3300662_trim.fq.gz SRR3300663_trim.fq.gz SRR3300664_trim.fq.gz > 9562_trim.fq.gz
#cat SRR3300668_trim.fq.gz SRR3300669_trim.fq.gz SRR3300670_trim.fq.gz > 9564_trim.fq.gz
#cat SRR3300677_trim.fq.gz SRR3300678_trim.fq.gz SRR3300679_trim.fq.gz > 9567_trim.fq.gz
#cat SRR3300680_trim.fq.gz SRR3300681_trim.fq.gz SRR3300682_trim.fq.gz > 9568_trim.fq.gz
#cat SRR3300703_trim.fq.gz SRR3300704_trim.fq.gz SRR3300705_trim.fq.gz > 9577_trim.fq.gz
#cat SRR3300717_trim.fq.gz SRR3300718_trim.fq.gz SRR3300719_trim.fq.gz > 9582_trim.fq.gz
#cat SRR3300734_trim.fq.gz SRR3300735_trim.fq.gz SRR3300736_trim.fq.gz > 9588_trim.fq.gz
#cat SRR3300752_trim.fq.gz SRR3300753_trim.fq.gz SRR3300754_trim.fq.gz > 9594_trim.fq.gz
#cat SRR3301234_trim.fq.gz SRR3301235_trim.fq.gz > 9821_trim.fq.gz
#cat SRR3301236_trim.fq.gz SRR3301237_trim.fq.gz SRR3301238_trim.fq.gz > 9822_trim.fq.gz
#cat SRR3301243_trim.fq.gz SRR3301244_trim.fq.gz > 9825_trim.fq.gz
#cat SRR3301259_trim.fq.gz SRR3301260_trim.fq.gz > 9833_trim.fq.gz
#cat SRR3301261_trim.fq.gz SRR3301262_trim.fq.gz SRR3301263_trim.fq.gz > 9834_trim.fq.gz
#cat SRR3301267_trim.fq.gz SRR3301268_trim.fq.gz > 9836_trim.fq.gz
#cat SRR3301276_trim.fq.gz SRR3301277_trim.fq.gz SRR3301278_trim.fq.gz > 9840_trim.fq.gz
#cat SRR3301279_trim.fq.gz SRR3301280_trim.fq.gz > 9841_trim.fq.gz
#cat SRR3301281_trim.fq.gz SRR3301282_trim.fq.gz SRR3301283_trim.fq.gz SRR3301284_trim.fq.gz SRR3301285_trim.fq.gz SRR3301286_trim.fq.gz > 9843_trim.fq.gz
#cat SRR3301287_trim.fq.gz SRR3301288_trim.fq.gz SRR3301289_trim.fq.gz SRR3301290_trim.fq.gz SRR3301291_trim.fq.gz SRR3301292_trim.fq.gz > 9844_trim.fq.gz
#cat SRR3301293_trim.fq.gz SRR3301294_trim.fq.gz > 9845_trim.fq.gz
#cat SRR3301313_trim.fq.gz SRR3301314_trim.fq.gz > 9852_trim.fq.gz
#cat SRR3301323_trim.fq.gz SRR3301324_trim.fq.gz SRR3301325_trim.fq.gz > 9856_trim.fq.gz
#cat SRR3301326_trim.fq.gz SRR3301327_trim.fq.gz SRR3301328_trim.fq.gz > 9857_trim.fq.gz
#cat SRR3301331_trim.fq.gz SRR3301332_trim.fq.gz SRR3301333_trim.fq.gz > 9859_trim.fq.gz
#cat SRR3301336_trim.fq.gz SRR3301337_trim.fq.gz > 9861_trim.fq.gz
#cat SRR3301340_trim.fq.gz SRR3301341_trim.fq.gz > 9864_trim.fq.gz
#cat SRR3301345_trim.fq.gz SRR3301346_trim.fq.gz > 9867_trim.fq.gz
#cat SRR3301347_trim.fq.gz SRR3301348_trim.fq.gz > 9868_trim.fq.gz
#cat SRR3301351_trim.fq.gz SRR3301352_trim.fq.gz SRR3301353_trim.fq.gz SRR3301354_trim.fq.gz > 9870_trim.fq.gz
#cat SRR3301359_trim.fq.gz SRR3301360_trim.fq.gz SRR3301361_trim.fq.gz > 9873_trim.fq.gz
#cat SRR3301366_trim.fq.gz SRR3301367_trim.fq.gz SRR3301368_trim.fq.gz > 9876_trim.fq.gz
#cat SRR3301371_trim.fq.gz SRR3301372_trim.fq.gz SRR3301373_trim.fq.gz > 9878_trim.fq.gz
#cat SRR3301389_trim.fq.gz SRR3301390_trim.fq.gz > 9886_trim.fq.gz
#cat SRR3301391_trim.fq.gz SRR3301392_trim.fq.gz SRR3301393_trim.fq.gz SRR3301394_trim.fq.gz > 9888_trim.fq.gz
#cat SRR3301414_trim.fq.gz SRR3301415_trim.fq.gz SRR3301416_trim.fq.gz SRR3301417_trim.fq.gz SRR3301418_trim.fq.gz > 9899_trim.fq.gz
#cat SRR3301419_trim.fq.gz SRR3301420_trim.fq.gz SRR3301421_trim.fq.gz > 9900_trim.fq.gz
#cat SRR3301425_trim.fq.gz SRR3301426_trim.fq.gz SRR3301427_trim.fq.gz SRR3301428_trim.fq.gz > 9902_trim.fq.gz
#cat SRR3301432_trim.fq.gz SRR3301433_trim.fq.gz > 9904_trim.fq.gz
#cat SRR3301516_trim.fq.gz SRR3301517_trim.fq.gz SRR3301518_trim.fq.gz > 9946_trim.fq.gz
#cat SRR3301527_trim.fq.gz SRR3301528_trim.fq.gz > 9950_trim.fq.gz


cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data
sed -e 's/X//' ceu | grep -F -f - samplenames3 >../ceu_data/samplenames_ceu_pe
sed -e 's/X//' ceu | grep -F -f - samplenames6 >../ceu_data/samplenames_ceu_se
sed -e 's/X//' ibnr | grep -F -f - samplenames3 >../ibnr_data/samplenames_ibnr_pe
sed -e 's/X//' ibnr | grep -F -f - samplenames6 >../ibnr_data/samplenames_ibnr_se

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
for file in *_trim.fq.gz; do mv $file ${file/_trim.fq.gz/.trim.fq.gz}; done
cat samplenames_ceu_se |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 $value1.trim.fq.gz  ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
for file in *_trim.fq.gz; do mv $file ${file/_trim.fq.gz/.trim.fq.gz}; done
cat samplenames_ibnr_se |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 $value1.trim.fq.gz  ; done
bismark_se.sh

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
for file in *_1.paired.fq.gz; do mv $file ${file/_1.paired.fq.gz/.1.paired.fq.gz}; done
for file in *_2.paired.fq.gz; do mv $file ${file/_2.paired.fq.gz/.2.paired.fq.gz}; done
cat samplenames_ceu_pe |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
for file in *_1.paired.fq.gz; do mv $file ${file/_1.paired.fq.gz/.1.paired.fq.gz}; done
for file in *_2.paired.fq.gz; do mv $file ${file/_2.paired.fq.gz/.2.paired.fq.gz}; done
cat samplenames_ibnr_pe |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/  --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
bismark_pe.sh

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
for file in *_bismark_hisat2*.bam ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
for file in *_bismark_hisat2*.bam ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
cat samplenames_ceu_se |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor  --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.trim_bismark_hisat2.deduplicated.bam  ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
cat samplenames_ibnr_se |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor  --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.trim_bismark_hisat2.deduplicated.bam  ; done


cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
cat samplenames_ceu_pe |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam  ; done

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
cat samplenames_ibnr_pe |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam  ; done


cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_data
#mkdir ../ceu_ibnr_data
#for file in *.bismark.cov.gz; do gunzip $file; done
#for file in *.CpG_report.txt.gz ; do gunzip $file; done
#for file in *.CpG_report.txt; do sed -e '/Mt/d' -e '/Pt/d' $file | awk -F'\t' '{$1 = substr($1, length($1), 1); print $0}' OFS='\t'  > ${file/_report.txt/_report2.txt} ; done
#python3 merge_files_byname_report.py . .

#for file in *.cov; do  cut -f 1,2,5,6 $file | sed '/Mt/d' | awk -F'\t' '{$1 = substr($1, length($1), 1); print $0}' OFS='\t' | sed 's/t/Pt/' > ${file/cov/cov2} ; done
#python3 merge_files_byname.py . .
#for file in *.cov2; do  awk -F'\t' '($3 + $4) > 9' $file > ${file/cov2/cov4} ; done
#for file in *.cov4; do python3 combine_cgs_mp.py $file >${file/.cov4/.cov5} ; done
#mv *cov4 ../ceu_ibnr_data/
#mv *cov5 ../ceu_ibnr_data/

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ibnr_data
#for file in *.bismark.cov.gz; do gunzip $file; done
#for file in *.CpG_report.txt.gz ; do gunzip $file; done
#for file in *.CpG_report.txt; do sed -e 's/^....//' -e '/Mt/d' -e '/Pt/d' $file  > ${file/_report.txt/_report2.txt} ; done
#python3 merge_files_byname_report.py . .

#for file in *.cov; do  cut -f 1,2,5,6 $file | sed -e 's/^....//'  > ${file/cov/cov2} ; done
#python3 merge_files_byname.py . .
#for file in *.cov2; do  awk -F'\t' '($3 + $4) > 4' $file > ${file/cov2/cov4} ; done
#for file in *.cov4; do python3 combine_cgs_mp.py $file >${file/.cov4/.cov5} ; done
#mv *cov4 ../ceu_ibnr_data/
#mv *cov5 ../ceu_ibnr_data/

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/ceu_ibnr_data
#ls *.cov4 >covfiles4 
#ls *.cov5 >covfiles5
#for file in  *combined.cov4; do sed -i 's/\.0//' $file ; done 
#Rscript binomtest_mC2_mod2025.R
#Rscript binomtest_mC2_cg_combine.R


############ next 7 lines run in 10.152.152.200!!!
############ Rscript jdmr.R >>jdmr_out 2>>jdmr_err 
############ for file in *_combined_methylome; do mv $file ${file/_combined_methylome/_methylome} ; done
############ sed -e 's/.trim//' -e 's/.1.paired//' -e s'/.txt//' ../data/ceu_dmrsamples2 >ceu_dmrsamples2 
############ sed -e 's/.trim//' -e 's/.1.paired//' -e s'/.txt//' ../data/ibnr_dmrsamples2 >ibnr_dmrsamples2
############ Rscript dmr_vcf_create.R
############ cat vcfheader ceu_dmrs_dip.txt >ceu_dmrs_dip.vcf
############ cat vcfheader ibnr_dmrs_dip.txt >ibnr_dmrs_dip.vcf

#Rscript methvcf.R
#sed -i '/^NA/d' arabidopsis_meth.vcf
#sed -i '/^NA/d' arabidopsis_meth_var_invar.vcf 
#sed -i 's/4e+05/400000/' arabidopsis_meth.vcf
#sed -i 's/3e+06/3000000/' arabidopsis_meth.vcf
#sed -i 's/1e+07/10000000/' arabidopsis_meth.vcf
#sed -i 's/7e+06/7000000/' arabidopsis_meth.vcf
#sed -i 's/5e+06/5000000/' arabidopsis_meth.vcf
#cat vcfheader arabidopsis_meth.vcf >arabidopsis_meth_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f arabidopsis_meth_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f arabidopsis_meth_all.vcf.gz

#sed -i -e 's/1e+07/10000000/' -e 's/1.2e+07/12000000/' -e 's/e+06/000000/' -e 's/e+05/00000/'  arabidopsis_meth_var_invar.vcf
#cat vcfheader arabidopsis_meth_var_invar.vcf >arabidopsis_meth_var_invar_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f arabidopsis_meth_var_invar_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f arabidopsis_meth_var_invar_all.vcf.gz

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz --bed  cg_sites.bed --recode --keep ceu_invcf_shuf --max-missing 0.8 --out ceu_meth_var_invar
#cat vcfheader ceu_meth_var_invar.recode.vcf >ceu_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz  --bed  cg_sites.bed --recode --keep ibnr_invcf --max-missing 0.8 --out ibnr_meth_var_invar
#cat vcfheader ibnr_meth_var_invar.recode.vcf >ibnr_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz --recode --keep ceu_invcf_shuf --max-missing 0.8 --out ceu_meth_var_invar
#cat vcfheader ceu_meth_var_invar.recode.vcf >ceu_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz  --recode --keep ibnr_invcf --max-missing 0.8 --out ibnr_meth_var_invar
#cat vcfheader ibnr_meth_var_invar.recode.vcf >ibnr_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz --bed  gbm_genes_poswithnames.txt --recode --keep ceu_invcf_shuf --max-missing 0.8 --out ceu_meth_var_invar
#cat vcfheader ceu_meth_var_invar.recode.vcf >ceu_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz  --bed  gbm_genes_poswithnames.txt --recode --keep ibnr_invcf --max-missing 0.8 --out ibnr_meth_var_invar
#cat vcfheader ibnr_meth_var_invar.recode.vcf >ibnr_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz --bed  non_gbm_genes_poswithnames.txt --recode --keep ceu_invcf_shuf --max-missing 0.8 --out ceu_meth_var_invar
#cat vcfheader ceu_meth_var_invar.recode.vcf >ceu_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ceu_meth_var_invar.vcf   --out ceu_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_var_invar_all.vcf.gz  --bed  non_gbm_genes_poswithnames.txt --recode --keep ibnr_invcf --max-missing 0.8 --out ibnr_meth_var_invar
#cat vcfheader ibnr_meth_var_invar.recode.vcf >ibnr_meth_var_invar.vcf
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_invar --maf 0
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf ibnr_meth_var_invar.vcf  --out ibnr_meth_var --maf 0.02 --min-alleles 2 --max-alleles 2

#Rscript methvcf_cg_combine.R
#sed -i '/^NA/d' arabidopsis_meth_cg_combine.vcf
#sed -i '/^NA/d' arabidopsis_meth_var_invar_cg_combine.vcf
#sed -i 's/4e+05/400000/' arabidopsis_meth_cg_combine.vcf
#sed -i 's/3e+06/3000000/' arabidopsis_meth_cg_combine.vcf
#sed -i 's/7e+06/7000000/' arabidopsis_meth_cg_combine.vcf
#sed -i 's/5e+06/5000000/' arabidopsis_meth_cg_combine.vcf
#cat vcfheader arabidopsis_meth_cg_combine.vcf >arabidopsis_meth_cg_combine_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f arabidopsis_meth_cg_combine_all.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f arabidopsis_meth_cg_combine_all.vcf.gz

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_meth_all.vcf.gz --keep ceu_invcf_shuf --recode  --max-missing 0.8  --out ceu_smps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_meth_all.vcf.gz --keep ibnr_invcf --recode  --max-missing 0.8  --out ibnr_smps_wholegenome --min-alleles 2 --max-alleles 2 --maf 0.02

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_meth_cg_combine_all.vcf.gz --keep ceu_invcf_shuf --recode  --max-missing 0.8  --out ceu_smps_wholegenome_cg_combine --min-alleles 2 --max-alleles 2 --maf 0.02
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf arabidopsis_meth_cg_combine_all.vcf.gz --keep ibnr_invcf --recode  --max-missing 0.8  --out ibnr_smps_wholegenome_cg_combine --min-alleles 2 --max-alleles 2 --maf 0.02

#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_wholegenome.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_wholegenome.recode.vcf.gz

#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_wholegenome.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_wholegenome.recode.vcf.gz

#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_wholegenome_cg_combine.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_wholegenome_cg_combine.recode.vcf.gz

#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_wholegenome_cg_combine.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_wholegenome_cg_combine.recode.vcf.gz

#mkdir smp_ceu_combine
#mkdir smp_ibnr_combine

#cp ceu_invcf_shuf smp_ceu/
#cp popgen5mb2.bed smp_ceu/
#cp chrlist smp_ceu/
#cp sample_chr_ceu smp_ceu/
#cp *.cytosines_tair10.bed smp_ceu/

#cp ibnr_invcf smp_ibnr/
#cp popgen5mb2.bed smp_ibnr/
#cp chrlist smp_ibnr/
#cp sample_chr_ibnr smp_ibnr/
#cp *.cytosines_tair10.bed smp_ibnr/

#cp ceu_invcf_shuf smp_ceu_combine/
#cp popgen5mb2.bed smp_ceu_combine/
#cp chrlist smp_ceu_combine/
#cp sample_chr_ceu smp_ceu_combine/
#cp *.cytosines_tair10.bed smp_ceu_combine/

#cp ibnr_invcf smp_ibnr_combine/
#cp popgen5mb2.bed smp_ibnr_combine/
#cp chrlist smp_ibnr_combine/
#cp sample_chr_ibnr smp_ibnr_combine/
#cp *.cytosines_tair10.bed smp_ibnr_combine/

cd smp_ceu/
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ../ceu_smps_wholegenome.recode.vcf.gz --recode  --out ceu_smps_5mb --bed popgen5mb2.bed
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_5mb.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_5mb.recode.vcf.gz
#cat ceu_invcf_shuf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ceu_smps_5mb.recode.vcf.gz ; done
#for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
#for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat sample_chr_ceu |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
#for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ceu_multihetsep_meth_$line ; done
#Rscript multihetsep_combined.R
cd ../

cd smp_ibnr/
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ../ibnr_smps_wholegenome.recode.vcf.gz --recode  --out ibnr_smps_5mb --bed popgen5mb2.bed
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_5mb.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_5mb.recode.vcf.gz
#cat ibnr_invcf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ibnr_smps_5mb.recode.vcf.gz ; done
#for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
#for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat sample_chr_ibnr |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
#for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ibnr_multihetsep_meth_$line ; done
#Rscript multihetsep_combined.R
cd ../

cd smp_ceu_combine/
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ../ceu_smps_wholegenome_cg_combine.recode.vcf.gz --recode  --out ceu_smps_5mb --bed popgen5mb2.bed
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_smps_5mb.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_smps_5mb.recode.vcf.gz
#cat ceu_invcf_shuf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ceu_smps_5mb.recode.vcf.gz ; done
#for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
#for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat sample_chr_ceu |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
#for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ceu_multihetsep_meth_$line ; done
#Rscript multihetsep_combined.R
cd ../

cd smp_ibnr_combine/
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ../ibnr_smps_wholegenome_cg_combine.recode.vcf.gz --recode  --out ibnr_smps_5mb --bed popgen5mb2.bed
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_smps_5mb.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_smps_5mb.recode.vcf.gz
#cat ibnr_invcf |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c 1 -O v -s $sample -o $sample.filtered.vcf ibnr_smps_5mb.recode.vcf.gz ; done
#for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
#for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat sample_chr_ibnr |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --max-missing 1.0 --recode --recode-INFO-all  --chr $value2 ; done
#for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
#for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
#cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >ibnr_multihetsep_meth_$line ; done
#Rscript multihetsep_combined.R
cd ../


#mkdir genes_fasta_ceu/
#mkdir genes_fasta_ibnr/
#mkdir cg_ceu
#mkdir cg_ibnr

cd cg_ceu/
#cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_smps_wholegenome.recode.vcf.gz  $line >$line.ceu_smps.vcf; done
#wc -l *vcf >vcflengths_ceu_smps
#zcat ../ceu_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
#for file in *ceu_smps.vcf; do cat vcfheader_ceu $file >${file/ceu_smps.vcf/all.vcf}; done
#wc -l *_smps.vcf >vcflengths_var_invar
#Rscript good_intervals_part1.R
#chmod 777 bad_intervals.sh
#./bad_intervals.sh
#grep 'CHROM' 1:5400-5600.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
#for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done 
#cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
#tar -cvzf input_files.tar.gz *input.txt
#mkdir -p input/
#for file in *input.txt; do sed -i '2d' $file; done
#mv *input.txt input/
#cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
#perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
##Rscript good_intervals_part2.R 
#mv input/* .
#cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ceu.txt -length $value2 -alpha $value5  ; done
#head -n 1 5:9942800-9943000.Dm_ceu.txt >Dm_results_header
#cat *.Dm_ceu.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ceu_cg_Dm.txt
#ls *.Dm_ceu.txt | sed 's/.Dm_ceu.txt//' >Dm_filenames_ceu_cg

cd ../cg_ibnr/
#cat ../cg_sites.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_smps_wholegenome.recode.vcf.gz  $line >$line.ibnr_smps.vcf; done
#wc -l *vcf >vcflengths_ibnr_smps
#zcat ../ibnr_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
#for file in *ibnr_smps.vcf; do cat vcfheader_ibnr $file >${file/ibnr_smps.vcf/all.vcf}; done
#wc -l *_smps.vcf >vcflengths_var_invar
#Rscript good_intervals_part1.R
#chmod 777 bad_intervals.sh
#./bad_intervals.sh
#grep 'CHROM' 1:5400-5600.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
#for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done
#cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
#tar -cvzf input_files.tar.gz *input.txt
#mkdir -p input/
#for file in *input.txt; do sed -i '2d' $file; done
#mv *input.txt input/
#cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
#perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
#Rscript good_intervals_part2.R 
#mv input/* .
#cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ibnr.txt -length $value2 -alpha $value5  ; done
#head -n 1 5:9959000-9959400.Dm_ibnr.txt >Dm_results_header
#cat *.Dm_ibnr.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ibnr_cg_Dm.txt
#ls *.Dm_ibnr.txt | sed 's/.Dm_ibnr.txt//' >Dm_filenames_ibnr_cg

cd ../genes_fasta_ceu/
#cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ceu_smps_wholegenome.recode.vcf.gz  $line >$line.ceu_smps.vcf; done
#wc -l *vcf >vcflengths_ceu_smps
#zcat ../ceu_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ceu
#for file in *_smps.vcf; do cat vcfheader_ceu $file >${file/ceu_smps.vcf/all.vcf}; done
#wc -l *all.vcf >vcflengths_var_invar
#Rscript good_intervals_part1.R
#chmod 777 bad_intervals.sh
#./bad_intervals.sh
#grep 'CHROM' 1:10003880-10004581.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
#for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done 
#cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
#tar -cvzf input_files.tar.gz *input.txt
#mkdir -p input/
#for file in *input.txt; do sed -i '2d' $file; done
#mv *input.txt input/
#cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
#perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
#Rscript good_intervals_part2.R 
#mv input/* .
#cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ceu.txt -length $value2 -alpha $value5  ; done
#head -n 1 5:9991685-9992770.Dm_ceu.txt >Dm_results_header
#cat *.Dm_ceu.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ceu_Dm.txt
#ls *.Dm_ceu.txt | sed 's/.Dm_ceu.txt//' >Dm_filenames_ceu

cd ../genes_fasta_ibnr/
#cat /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gene_pos.list | while read -r line ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix ../ibnr_smps_wholegenome.recode.vcf.gz  $line >$line.ibnr_smps.vcf; done
#wc -l *vcf >vcflengths_ibnr_smps
#zcat ../ibnr_smps_wholegenome.recode.vcf.gz | grep '#' >vcfheader_ibnr
#for file in *_smps.vcf; do cat vcfheader_ibnr $file >${file/ibnr_smps.vcf/all.vcf}; done
#wc -l *all.vcf >vcflengths_var_invar
#Rscript good_intervals_part1.R
#chmod 777 bad_intervals.sh
#./bad_intervals.sh
#grep 'CHROM' 1:10003880-10004581.all.vcf | cut --complement  -f 3-9 | sed -e 's/CHROM/chr/' -e 's/POS/position/' >Dm_header
#for file in *.all.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' | cat Dm_header -  >${file/.all.vcf/.input.txt} ; done
#cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
#mkdir -p input/
#for file in *input.txt; do sed -i '2d' $file; done
#mv *input.txt input/
#cp /proj/popgen/a.ramesh/software/diff_two_seq.pl .
#perl /proj/popgen/a.ramesh/software/alpha_estimation.pl -dir input -output  alpha_Dm -length_list length_list
#Rscript good_intervals_part2.R 
#mv input/* .
#cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl /proj/popgen/a.ramesh/software/Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_ibnr.txt -length $value2 -alpha $value5  ; done
#head -n 1 5:9991685-9992770.Dm_ibnr.txt >Dm_results_header
#cat *.Dm_ibnr.txt | sed '/chr/d' | cat Dm_results_header - | sed 's/#//' >ibnr_Dm.txt
#ls *.Dm_ibnr.txt | sed 's/.Dm_ibnr.txt//' >Dm_filenames_ibnr

cd ../


#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_ceu --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --gzvcf ceu_smps_wholegenome.recode.vcf.gz --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_ibnr --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --gzvcf ibnr_smps_wholegenome.recode.vcf.gz --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_ceu.recode.vcf --out gbm_ceu_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_ibnr.recode.vcf --out gbm_ibnr_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_ceu_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_ceu_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_ibnr_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_ibnr_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_meth_gbm.recode.vcf gbm_ceu_maf.recode.vcf.gz gbm_ibnr_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_meth_gbm.recode.vcf -OutPut ceu_ibnr_meth_gbm_dis.mat

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_dmr_ceu --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --vcf ceu_dmrs_dip.vcf --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --out gbm_dmr_ibnr --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed --vcf ibnr_dmrs_dip.vcf  --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_dmr_ceu.recode.vcf --out gbm_dmr_ceu_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf gbm_dmr_ibnr.recode.vcf --out gbm_dmr_ibnr_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_dmr_ceu_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_dmr_ceu_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f gbm_dmr_ibnr_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f gbm_dmr_ibnr_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_dmrs_gbm.recode.vcf gbm_dmr_ceu_maf.recode.vcf.gz gbm_dmr_ibnr_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_dmrs_gbm.recode.vcf -OutPut ceu_ibnr_dmrs_gbm.mat



#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ceu_meth_cg  --bed  cg_sites.bed  --keep ceu_invcf_shuf  --max-missing 0.8 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ibnr_meth_cg  --bed  cg_sites.bed   --keep ibnr_invcf  --max-missing 0.8 --recode

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ceu_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed  --keep ceu_invcf_shuf  --max-missing 0.8 --recode 
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_gbm.recode.vcf --out ceu_gbm --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ibnr_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed  --keep ibnr_invcf  --max-missing 0.8 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_gbm.recode.vcf --out ibnr_gbm --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ceu_all_genes  --keep ceu_invcf_shuf  --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ibnr_all_genes  --keep ibnr_invcf  --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ceu_non_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed --keep ceu_invcf_shuf  --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ibnr_non_gbm  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed --keep ibnr_invcf  --max-missing 0.8 --freq


#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_cg --recode --bed cg_sites.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_cg --recode --bed cg_sites.bed

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_gbm --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_gbm.recode.vcf --out ceu_dmrs_gbm --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_gbm --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/gbm_genes.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_gbm.recode.vcf --out ibnr_dmrs_gbm --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ceu_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed  --keep ceu_invcf_shuf  --max-missing 0.8 --recode 
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_sift.recode.vcf --out ceu_sift --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz  --out ibnr_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed  --keep ibnr_invcf  --max-missing 0.8 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_sift.recode.vcf --out ibnr_sift --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ceu_low_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed --keep ceu_invcf_shuf  --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --gzvcf arabidopsis_meth_all.vcf.gz   --out ibnr_low_sift  --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed --keep ibnr_invcf  --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_sift --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_sift.recode.vcf --out ceu_dmrs_sift --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_sift --recode --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/high_sift.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_sift.recode.vcf --out ibnr_dmrs_sift --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_meth_cg.recode.vcf --out ceu_smps  --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_cg.recode.vcf --out ceu_dmrs --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_meth_cg.recode.vcf --out ibnr_smps --max-missing 0.8 --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_cg.recode.vcf --out ibnr_dmrs --max-missing 0.8 --freq

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_all  --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_all  --freq
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_non_gbm  --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_non_gbm --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gbm_genes.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_dip.vcf --out ceu_dmrs_low_sift --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_dip.vcf --out ibnr_dmrs_low_sift --freq --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/low_sift.bed

#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_meth_cg.recode.vcf -OutStat ceu_meth_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_meth_cg.recode.vcf -OutStat ibnr_meth_cg_ld_decay -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF gbm_ceu.recode.vcf -OutStat ceu_meth_gbm_ld_decay -MaxDist 10 -OutType 3
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF gbm_ibnr.recode.vcf -OutStat ibnr_meth_gbm_ld_decay -MaxDist 10 -OutType 3
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_smps_wholegenome.recode.vcf.gz -OutStat ceu_smps_wholegenome_ld_decay -MaxDist 10 -OutType 3
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_smps_wholegenome.recode.vcf.gz -OutStat ibnr_smps_wholegenome_ld_decay -MaxDist 10 -OutType 3

python3 ld_bin.py ceu_meth_cg_ld_decay.LD.gz ceu_meth_cg_ld_decay.stat 3
python3 ld_bin.py ibnr_meth_cg_ld_decay.LD.gz ibnr_meth_cg_ld_decay.stat 3
python3 ld_bin.py ceu_meth_gbm_ld_decay.LD.gz ceu_meth_gbm_ld_decay.stat 3
python3 ld_bin.py ibnr_meth_gbm_ld_decay.LD.gz ibnr_meth_gbm_ld_decay.stat 3
python3 ld_bin.py ceu_smps_wholegenome_ld_decay.LD.gz ceu_smps_wholegenome_ld_decay.stat 3
python3 ld_bin.py ibnr_smps_wholegenome_ld_decay.LD.gz ibnr_smps_wholegenome_ld_decay.stat 3

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_meth_cg.recode.vcf --out ceu_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_meth_cg.recode.vcf --out ibnr_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_meth_cg_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_meth_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_meth_cg_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_meth_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_meth_cg_maf.recode.vcf ceu_meth_cg_maf.recode.vcf.gz ibnr_meth_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_meth_cg_maf.recode.vcf -OutPut ceu_ibnr_meth_cg_dis.mat

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_cg.recode.vcf --out ceu_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_cg.recode.vcf --out ibnr_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_dmrs_cg_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_dmrs_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_dmrs_cg_maf.recode.vcf
#/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_dmrs_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_dmrs_cg_maf.recode.vcf ceu_dmrs_cg_maf.recode.vcf.gz ibnr_dmrs_cg_maf.recode.vcf.gz
#/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_dmrs_cg_maf.recode.vcf -OutPut ceu_ibnr_dmrs_cg.mat

#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_meth_cg_maf_thin --vcf ceu_ibnr_meth_cg_maf.recode.vcf --thin 1000
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_meth_gbm_thin --vcf ceu_ibnr_meth_gbm.recode.vcf --thin 1000
#sed -i 's/-/_/'  ceu_ibnr_meth_cg_maf_thin.recode.vcf
#sed -i 's/-/_/'  ceu_ibnr_meth_gbm_thin.recode.vcf

```

Appendix 1b. One difference from before is that when there are paired and single end methylation reads for the same sample, the counts for methyated and unmethylated reads were summed prior to the binomial test. First for the coverage files. Implemented in merge_files_byname.py.
```
import os
import argparse
import pandas as pd
from collections import defaultdict
import shutil

def get_file_pairs(file_list):
    """Group files by the prefix before the first period and return pairs."""
    file_dict = defaultdict(list)

    # Group files based on the prefix before the first period
    for file in file_list:
        prefix = file.split('.')[0]  # Get the part before the first period
        file_dict[prefix].append(file)
    
    # Only return pairs of files
    file_pairs = [(files[0], files[1]) for files in file_dict.values() if len(files) == 2]
    return file_pairs

def read_file(file_path):
    """Read the input file as a tab-separated file."""
    return pd.read_csv(file_path, sep="\t", header=None, names=["Column1", "Column2", "Column3", "Column4"], dtype={"Column1": str})

def merge_files(file1, file2):
    """Merge the two files based on the unique ID and sum the columns."""
    # Read files
    df1 = read_file(file1)
    df2 = read_file(file2)

    # Convert "Column1" to string (if it makes sense for your data)
    df1["Column1"] = df1["Column1"].astype(str)
    df2["Column1"] = df2["Column1"].astype(str)

    # Convert "Column2" to integer (if it makes sense for your data)
    df1["Column2"] = df1["Column2"].astype(int)
    df2["Column2"] = df2["Column2"].astype(int)

    # Add unique ID
    df1["ID"] = df1["Column1"] + "_" + df1["Column2"].astype(str)
    df2["ID"] = df2["Column1"] + "_" + df2["Column2"].astype(str)

    # Merge dataframes
    merged_df = pd.merge(df1, df2, on="ID", how="inner", suffixes=("_file1", "_file2"))
    merged_df.fillna(0, inplace=True)

    # Summing columns
    merged_df["Column3_sum"] = merged_df["Column3_file1"] + merged_df["Column3_file2"]
    merged_df["Column4_sum"] = merged_df["Column4_file1"] + merged_df["Column4_file2"]

    unique_df1 = df1[~df1["ID"].isin(df2["ID"])]
    unique_df2 = df2[~df2["ID"].isin(df1["ID"])]

    unique_df = pd.concat([unique_df1, unique_df2], ignore_index=True)
    unique_df = unique_df.iloc[:, :4]
    unique_df.columns = ['Column1', 'Column2', 'Column3_sum', 'Column4_sum']

    # Prepare the final output
    result_df = merged_df[["Column1_file1", "Column2_file1", "Column3_sum", "Column4_sum"]].copy()  # Explicitly make a copy
    result_df.rename(columns={"Column1_file1": "Column1", "Column2_file1": "Column2"}, inplace=True)
    

    # Convert "Column1" to string (if it makes sense for your data)
    result_df["Column1"] = result_df["Column1"].astype(str)

    # Convert "Column2" to integer (if it makes sense for your data)
    result_df["Column2"] = result_df["Column2"].astype(int)

    # Convert only numeric columns to int, leave non-numeric columns as they are
    result_df.loc[:, "Column3_sum"] = result_df["Column3_sum"].astype(int)
    result_df.loc[:, "Column4_sum"] = result_df["Column4_sum"].astype(int)

    result_df = pd.concat([result_df, unique_df], ignore_index=True)

    result_df[result_df.columns[0]] = result_df[result_df.columns[0]].astype(str)
    result_df[result_df.columns[1]] = result_df[result_df.columns[1]].astype(int)
    result_df = result_df.sort_values(by=[result_df.columns[0], result_df.columns[1]])

    result_df.fillna(0, inplace=True)
    return result_df


def run_merge(file1, file2, output_dir):
    """Run the merge for a pair of files and save the result."""
    # Get the sample name (prefix before the first period) from file1 or file2
    sample_name = file1.split('.1')[0]  # Extract prefix before the first period
    
    # Output file paths
    output_file = os.path.join(output_dir, f"{sample_name}_combined.cov2")
    
    # Merge files
    merged_df = merge_files(file1, file2)
    
    # Save the merged result to the output file
    merged_df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Merged files {file1} and {file2} -> {output_file}")
    
    
    # Move the original files to the 'old' folder
    old_folder = os.path.join(os.path.dirname(file1), 'old')
    if not os.path.exists(old_folder):
        os.makedirs(old_folder)
    
    # Move files
    shutil.move(file1, os.path.join(old_folder, os.path.basename(file1)))
    shutil.move(file2, os.path.join(old_folder, os.path.basename(file2)))
    print(f"Moved {file1} and {file2} to {old_folder}")

def main():
    parser = argparse.ArgumentParser(description="Process file pairs for merging.")
    parser.add_argument("input_dir", help="Directory containing the files")
    parser.add_argument("output_dir", help="Directory to save the merged files")
    args = parser.parse_args()

    # Get the list of files from the input directory
    file_list = [f for f in os.listdir(args.input_dir) if f.endswith("bismark.cov2")]
    
    # Get the pairs of files to merge
    file_pairs = get_file_pairs(file_list)

    # Run the merge for each pair
    for file1, file2 in file_pairs:
        run_merge(os.path.join(args.input_dir, file1), os.path.join(args.input_dir, file2), args.output_dir)

if __name__ == "__main__":
    main()

```

Appendix 1c. Also for the report files. Implemented in merge_files_byname_report.py 
```
import os
import argparse
import pandas as pd
from collections import defaultdict
import shutil

def get_file_pairs(file_list):
    """Group files by the prefix before the first period and return pairs."""
    file_dict = defaultdict(list)

    # Group files based on the prefix before the first period
    for file in file_list:
        prefix = file.split('.')[0]  # Get the part before the first period
        file_dict[prefix].append(file)
    
    # Only return pairs of files
    file_pairs = [(files[0], files[1]) for files in file_dict.values() if len(files) == 2]
    return file_pairs

def read_file(file_path):
    """Read the input file as a tab-separated file."""
    return pd.read_csv(file_path, sep="\t", header=None, names=["Column1", "Column2", "Column3", "Column4", "Column5", "Column6", "Column7"], dtype={"Column1": str})

def merge_files(file1, file2):
    """Merge the two files based on the unique ID and sum the columns."""
    # Read files
    df1 = read_file(file1)
    df2 = read_file(file2)

    # Convert "Column1" to string (if it makes sense for your data)
    df1["Column1"] = df1["Column1"].astype(str)
    df2["Column1"] = df2["Column1"].astype(str)

    # Convert "Column2" to integer (if it makes sense for your data)
    df1["Column2"] = df1["Column2"].astype(int)
    df2["Column2"] = df2["Column2"].astype(int)

    # Add unique ID
    df1["ID"] = df1["Column1"] + "_" + df1["Column2"].astype(str)
    df2["ID"] = df2["Column1"] + "_" + df2["Column2"].astype(str)

    df1 = df1.add_suffix("_file1")
    df2 = df2.add_suffix("_file2")

    # Merge dataframes
    merged_df = pd.concat([df1, df2], axis=1)

    # Summing columns
    merged_df["Column4_sum"] = merged_df["Column4_file1"] + merged_df["Column4_file2"]
    merged_df["Column5_sum"] = merged_df["Column5_file1"] + merged_df["Column5_file2"]

    # Prepare the final output
    result_df = merged_df[["Column1_file1", "Column2_file1", "Column3_file1", "Column4_sum", "Column5_sum", "Column6_file1", "Column7_file1",]].copy()  # Explicitly make a copy
    result_df.rename(columns={"Column1_file1": "Column1", "Column2_file1": "Column2", "Column3_file1": "Column4", "Column6_file1": "Column6", "Column7_file1": "Column7"}, inplace=True)
    

    # Convert "Column1" to string (if it makes sense for your data)
    result_df["Column1"] = result_df["Column1"].astype(str)

    # Convert "Column2" to integer (if it makes sense for your data)
    result_df["Column2"] = result_df["Column2"].astype(int)

    # Convert only numeric columns to int, leave non-numeric columns as they are
    result_df.loc[:, "Column4_sum"] = result_df["Column4_sum"].astype(int)
    result_df.loc[:, "Column5_sum"] = result_df["Column5_sum"].astype(int)


    result_df[result_df.columns[0]] = result_df[result_df.columns[0]].astype(str)
    result_df[result_df.columns[1]] = result_df[result_df.columns[1]].astype(int)

    result_df.fillna(0, inplace=True)
    return result_df


def run_merge(file1, file2, output_dir):
    """Run the merge for a pair of files and save the result."""
    # Get the sample name (prefix before the first period) from file1 or file2
    sample_name = file1.split('.1')[0]  # Extract prefix before the first period
    
    # Output file paths
    output_file = os.path.join(output_dir, f"{sample_name}_combined.CpG_report2.txt")
    
    # Merge files
    merged_df = merge_files(file1, file2)
    
    # Save the merged result to the output file
    merged_df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Merged files {file1} and {file2} -> {output_file}")
    
    
    # Move the original files to the 'old' folder
    old_folder = os.path.join(os.path.dirname(file1), 'old2')
    if not os.path.exists(old_folder):
        os.makedirs(old_folder)
    
    # Move files
    shutil.move(file1, os.path.join(old_folder, os.path.basename(file1)))
    shutil.move(file2, os.path.join(old_folder, os.path.basename(file2)))
    print(f"Moved {file1} and {file2} to {old_folder}")

def main():
    parser = argparse.ArgumentParser(description="Process file pairs for merging.")
    parser.add_argument("input_dir", help="Directory containing the files")
    parser.add_argument("output_dir", help="Directory to save the merged files")
    args = parser.parse_args()

    # Get the list of files from the input directory
    file_list = [f for f in os.listdir(args.input_dir) if f.endswith(".CpG_report2.txt")]
    
    # Get the pairs of files to merge
    file_pairs = get_file_pairs(file_list)

    # Run the merge for each pair
    for file1, file2 in file_pairs:
        run_merge(os.path.join(args.input_dir, file1), os.path.join(args.input_dir, file2), args.output_dir)

if __name__ == "__main__":
    main()

```
Appendix 1d. For one set of coverage files, the CG counts for the + and - strands were summed. Implemented in  These results are not included in the paper.

```
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import os

# Function to process a chunk of data
def process_chunk(chunk):
    results = []

    # Ensure columns are numeric, handling any non-numeric values
    chunk["Column 2"] = pd.to_numeric(chunk["Column 2"], errors='coerce').fillna(0).astype(int)
    chunk["Column 3"] = pd.to_numeric(chunk["Column 3"], errors='coerce').fillna(0).astype(int)
    chunk["Column 4"] = pd.to_numeric(chunk["Column 4"], errors='coerce').fillna(0).astype(int)

    # Convert chunk to numpy array for processing
    data = chunk.to_numpy()

    i = 0
    while i < len(data) - 1:
        current_row = data[i]
        next_row = data[i + 1]

        # Check if the difference in Column 2 is 1 and the first row is odd, the second is even
        if (next_row[1] - current_row[1] == 1 and
            current_row[1] % 2 != 0 and
            next_row[1] % 2 == 0):
           
            # Sum Column 3 and Column 4 values
            sum_column_3 = current_row[2] + next_row[2]
            sum_column_4 = current_row[3] + next_row[3]
           
            # Append formatted result
            results.append(f"{current_row[0]}     {current_row[1]}     {sum_column_3}     {sum_column_4}")
           
            # Skip the next row
            i += 2
        else:
            # Append the current row as is
            results.append(f"{current_row[0]}     {current_row[1]}     {current_row[2]}     {current_row[3]}")
            i += 1

    # Add any remaining unprocessed row
    if i < len(data):
        last_row = data[i]
        results.append(f"{last_row[0]}     {last_row[1]}     {last_row[2]}     {last_row[3]}")

    return results

# Function to read the file in chunks and process them in parallel
def process_file_in_parallel(file_path, chunk_size=100000):
    # Pool for parallel processing
    pool = Pool(cpu_count())

    # Read the file in chunks and apply the process function in parallel
    results = []
    for chunk in pd.read_csv(file_path, header=None, delimiter='\s+',
                             names=["ID", "Column 2", "Column 3", "Column 4"],
                             chunksize=chunk_size):
        # Submit each chunk to the pool
        chunk_results = pool.apply_async(process_chunk, args=(chunk,))
        results.append(chunk_results)

    # Close and join the pool
    pool.close()
    pool.join()

    # Collect results
    all_results = []
    for result in results:
        all_results.extend(result.get())

    return all_results

# Main function
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process large files in parallel and find row pairs.")
    parser.add_argument('file', type=str, help="Path to the input data file")
    parser.add_argument('--chunk-size', type=int, default=100000, help="Number of rows per chunk (default: 100,000)")

    args = parser.parse_args()

    # Check if file exists
    if not os.path.exists(args.file):
        print(f"Error: File {args.file} not found.")
        exit(1)

    # Process the file
    processed_results = process_file_in_parallel(args.file, args.chunk_size)

    # Print the results
    for line in processed_results:
        print(line)

```

Appendix 1e. binomtest_mC2_mod2025 has implements the binomial test for the new set of files

```

library(dplyr)
library(ggplot2)
library(data.table)

covfile <- read.table(file="covfiles4")

cov <- fread(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "count.methylated", "count.unmethylated")
cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
cov$count.total <- cov$count.methylated + cov$count.unmethylated
cov$pval <- 0
noncoversionrate <- sum(cov[cov$chromosome %in% "Pt",]$count.methylated)/sum(cov[cov$chromosome %in% "Pt",]$count.total)
b <- apply(cov[,c(3,6)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
cov <- cov[,c(1,2,5,7)]
cov$fdr <- p.adjust(cov$pval,method = "fdr")
cov$call <- "U"
cov[cov$fdr < 0.01,]$call <- "M"
cov <- cov[,-c(4,5)]
samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov4","",covfile[1,])
samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov4","",covfile[1,])
colnames(cov)[ncol(cov)] <- samplename
cov2 <- cov

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles4")
  print(i)
  cov <- fread(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "count.methylated", "count.unmethylated")
  cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
  cov$count.total <- cov$count.methylated + cov$count.unmethylated
  cov$pval <- 0
  noncoversionrate <- sum(cov[cov$chromosome %in% "Pt",]$count.methylated)/sum(cov[cov$chromosome %in% "Pt",]$count.total)
  b <- apply(cov[,c(3,6)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
  cov <- cov[,c(1,2,5,7)]
  cov$fdr <- p.adjust(cov$pval,method = "fdr")
  cov$call <- "U"
  cov[cov$fdr < 0.01,]$call <- "M"
  cov <- cov[,-c(4,5)]
  
  samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov4","",covfile[i,])
  samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov4","",covfile[i,])
  colnames(cov)[ncol(cov)] <- samplename
  cov <- cov[,c(3,4)]
  cov2 <- merge(cov2,cov,by="ID",all=TRUE)
  write.table(cov2,file="cov_tmp.txt",row.names = F)
}

cov2$chromosome <- gsub("-.*","",cov2$ID)
cov2$position <- as.numeric(gsub(".*-","",cov2$ID))
write.table(cov2,file="cov_context3_new.txt",row.names = F)

```

Appendix 1f. binomtest_mC2_cg_combine.R has the binomial test for the coverage counts where symmertic CG sites  were combined. Not shown in paper.

```

library(dplyr)
library(ggplot2)
library(data.table)

covfile <- read.table(file="covfiles5")

cov <- fread(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "count.methylated", "count.unmethylated")
cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
cov$count.total <- cov$count.methylated + cov$count.unmethylated
cov$pval <- 0
noncoversionrate <- sum(cov[cov$chromosome %in% "Pt",]$count.methylated)/sum(cov[cov$chromosome %in% "Pt",]$count.total)
b <- apply(cov[,c(3,6)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
cov <- cov[,c(1,2,5,7)]
cov$fdr <- p.adjust(cov$pval,method = "fdr")
cov$call <- "U"
cov[cov$fdr < 0.01,]$call <- "M"
cov <- cov[,-c(4,5)]
samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov5","",covfile[1,])
samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov5","",covfile[1,])
colnames(cov)[ncol(cov)] <- samplename
cov2 <- cov

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles5")
  print(i)
  cov <- fread(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "count.methylated", "count.unmethylated")
  cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
  cov$count.total <- cov$count.methylated + cov$count.unmethylated
  cov$pval <- 0
  noncoversionrate <- sum(cov[cov$chromosome %in% "Pt",]$count.methylated)/sum(cov[cov$chromosome %in% "Pt",]$count.total)
  b <- apply(cov[,c(3,6)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
  cov <- cov[,c(1,2,5,7)]
  cov$fdr <- p.adjust(cov$pval,method = "fdr")
  cov$call <- "U"
  cov[cov$fdr < 0.01,]$call <- "M"
  cov <- cov[,-c(4,5)]
  
  samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov5","",covfile[i,])
  samplename <- gsub(".trim_bismark_hisat2.deduplicated.bismark.cov5","",covfile[i,])
  colnames(cov)[ncol(cov)] <- samplename
  cov <- cov[,c(3,4)]
  cov2 <- merge(cov2,cov,by="ID",all=TRUE)
  write.table(cov2,file="cov_tmp_cg_combine.txt",row.names = F)
}

cov2$chromosome <- gsub("-.*","",cov2$ID)
cov2$position <- as.numeric(gsub(".*-","",cov2$ID))
write.table(cov2,file="cov_context3_cg_combine.txt",row.names = F)

```

Appendix 1g. Updated methvcf.R file to process methylation state calls into VCF.

```
cov_context3 <- read.table(file="cov_context3_new.txt",header=T)
colnames(cov_context3) <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov4","",colnames(cov_context3))
colnames(cov_context3) <- gsub("_combined.cov4","",colnames(cov_context3))
na_count <- apply(cov_context3[4:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[4:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[4:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[c(2,3,1)]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"

cov_context3 <- cov_context3[-c(1:3)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
cov_context3 <- cbind(meta,cov_context3)
cov_context3$`#CHROM` <- as.numeric(cov_context3$`#CHROM`)
cov_context3$POS <- as.numeric(cov_context3$POS)
cov_context3 <- cov_context3[order(cov_context3$`#CHROM`,cov_context3$POS),]
write.table(cov_context3,file="arabidopsis_meth.vcf",quote = F, row.names = F,sep="\t")

#meta2 <- cov_context4[c(2,3,1)]
#colnames(meta2) <- c("#CHROM","POS","ID")
#meta2$REF <- "A"
#meta2$ALT <- "T"
#meta2$QUAL <- 4000
#meta2$FILTER <- "PASS"
#meta2$INFO <- "DP=1000"
#meta2$FORMAT <- "GT"
#cov_context4 <- cov_context4[-c(1:3)]
#cov_context4[cov_context4 == "U"] <- "0/0"
#cov_context4[cov_context4 == "M"] <- "1/1"
#cov_context4[is.na(cov_context4)] <- "./."
#cov_context4 <- cbind(meta2,cov_context4)
#cov_context4$`#CHROM` <- as.numeric(cov_context4$`#CHROM`)
#cov_context4$POS <- as.numeric(cov_context4$POS)
#cov_context4 <- cov_context4[order(cov_context4$`#CHROM`,cov_context4$POS),]
#write.table(cov_context4,file="arabidopsis_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")
```

Appendix 1h. Updated methvcf_cg_combine.R file to process methylation state calls into VCF where symmertic CG sites  were combined. Not shown in paper.

```
cov_context3 <- read.table(file="cov_context3_cg_combine.txt",header=T)
colnames(cov_context3) <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov5","",colnames(cov_context3))
colnames(cov_context3) <- gsub("_combined.cov5","",colnames(cov_context3))
na_count <- apply(cov_context3[4:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[4:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[4:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[c(2,3,1)]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"

cov_context3 <- cov_context3[-c(1:3)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
cov_context3 <- cbind(meta,cov_context3)
cov_context3$`#CHROM` <- as.numeric(cov_context3$`#CHROM`)
cov_context3$POS <- as.numeric(cov_context3$POS)
cov_context3 <- cov_context3[order(cov_context3$`#CHROM`,cov_context3$POS),]
write.table(cov_context3,file="arabidopsis_meth_cg_combine.vcf",quote = F, row.names = F,sep="\t")


meta2 <- cov_context4[c(2,3,1)]
colnames(meta2) <- c("#CHROM","POS","ID")
meta2$REF <- "A"
meta2$ALT <- "T"
meta2$QUAL <- 4000
meta2$FILTER <- "PASS"
meta2$INFO <- "DP=1000"
meta2$FORMAT <- "GT"
cov_context4 <- cov_context4[-c(1:3)]
cov_context4[cov_context4 == "U"] <- "0/0"
cov_context4[cov_context4 == "M"] <- "1/1"
cov_context4[is.na(cov_context4)] <- "./."
cov_context4 <- cbind(meta2,cov_context4)
cov_context4$`#CHROM` <- as.numeric(cov_context4$`#CHROM`)
cov_context4$POS <- as.numeric(cov_context4$POS)
cov_context4 <- cov_context4[order(cov_context4$`#CHROM`,cov_context4$POS),]
write.table(cov_context4,file="arabidopsis_meth_var_invar_cg_combine.vcf",quote = F, row.names = F,sep="\t")

```

Appendix 1i. Updated dmr_vcf_create.R to create vcf from DMR call files. jdmr.R remains the same as before. Only a filename change in this script.

```
library(dplyr)

filenames <- read.table(file="../ceu_dmrsamples2",header = T)
filenames$file <- gsub("_methylome","_methylome_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
dmrs <- dmrs[!duplicated(dmrs),]
dmrs[dmrs$posteriorMax < 0.99,]$status <- NA
gt <- dmrs$status
gt <- gsub("U",0,gt)
gt <- gsub("M",1,gt)
gt[is.na(gt)] <- "."
gt <- as.data.frame(gt)
gt2 <- gt
gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
colnames(gt)[1] <- filenames[1,2]
colnames(gt2)[1] <- filenames[1,2]
lth <- dmrs$end - dmrs$start
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs <- dmrs[3]
colnames(dmrs) <- "ID"
lengths <- as.data.frame(cbind(dmrs,lth))
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)

for (i in 2:nrow(filenames)){
  print(i)
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  dmrs_tmp <- dmrs_tmp[!duplicated(dmrs_tmp),]
  lth <- dmrs_tmp$end - dmrs_tmp$start
  dmrs_tmp$start <- round((dmrs_tmp$start + dmrs_tmp$end)/2)
  dmrs_tmp$end <- paste(dmrs_tmp$seqnames,dmrs_tmp$start,sep="_")
  ID <- dmrs_tmp$end
  lengths <- rbind(lengths,as.data.frame(cbind(ID,lth)))
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  gt <- cbind(gt,dmrs_tmp$end)
  colnames(gt)[2] <- "ID"
  gt2 <- cbind(gt2,dmrs_tmp$end)
  colnames(gt2)[2] <- "ID"
  gt <- gt[!duplicated(gt),]
  gt2 <- gt2[!duplicated(gt2),]
  dmrs <- full_join(dmrs,gt,by="ID")
  dmrs_dip <- full_join(dmrs_dip,gt2,by="ID")
  print(table(duplicated(dmrs_dip$ID)))
}

dmrs[is.na(dmrs)] <- "."
dmrs_dip[is.na(dmrs_dip)] <- "./."

CHROM <- as.integer(gsub("_.*","",dmrs$ID))
POS <- as.integer(gsub(".*_","",dmrs$ID))
ID <- dmrs$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "D"
meta$ALT <- "M"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs <- cbind(meta,dmrs[2:ncol(dmrs)])
dmrs$POS <- as.integer(dmrs$POS)
dmrs <- dmrs[order(dmrs$CHROM,dmrs$POS),]
colnames(dmrs)[1] <- "#CHROM"
write.table(dmrs,file="ceu_dmrs.txt",row.names = F, quote = F,sep="\t")

CHROM <- as.integer(gsub("_.*","",dmrs_dip$ID))
POS <- as.integer(gsub(".*_","",dmrs_dip$ID))
ID <- dmrs_dip$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs_dip <- cbind(meta,dmrs_dip[2:ncol(dmrs_dip)])
dmrs_dip$POS <- as.integer(dmrs_dip$POS)
dmrs_dip <- dmrs_dip[order(dmrs_dip$CHROM,dmrs_dip$POS),]
colnames(dmrs_dip)[1] <- "#CHROM"
write.table(dmrs_dip,file="ceu_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

lengths <- lengths[!duplicated(lengths),]
write.table(lengths,file="ceu_dmrs_lengths.txt",row.names = F, quote = F,sep="\t")

filenames <- read.table(file="../ibnr_dmrsamples2",header = T)
filenames$file <- gsub("_methylome","_methylome_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
dmrs <- dmrs[!duplicated(dmrs),]
dmrs[dmrs$posteriorMax < 0.99,]$status <- NA
gt <- dmrs$status
gt <- gsub("U",0,gt)
gt <- gsub("M",1,gt)
gt[is.na(gt)] <- "."
gt <- as.data.frame(gt)
gt2 <- gt
gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
colnames(gt)[1] <- filenames[1,2]
colnames(gt2)[1] <- filenames[1,2]
lth <- dmrs$end - dmrs$start
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs <- dmrs[3]
colnames(dmrs) <- "ID"
lengths <- as.data.frame(cbind(dmrs,lth))
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)

for (i in 2:nrow(filenames)){
  print(i)
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  dmrs_tmp <- dmrs_tmp[!duplicated(dmrs_tmp),]
  lth <- dmrs_tmp$end - dmrs_tmp$start
  dmrs_tmp$start <- round((dmrs_tmp$start + dmrs_tmp$end)/2)
  dmrs_tmp$end <- paste(dmrs_tmp$seqnames,dmrs_tmp$start,sep="_")
  ID <- dmrs_tmp$end
  lengths <- rbind(lengths,as.data.frame(cbind(ID,lth)))
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  gt <- cbind(gt,dmrs_tmp$end)
  colnames(gt)[2] <- "ID"
  gt2 <- cbind(gt2,dmrs_tmp$end)
  colnames(gt2)[2] <- "ID"
  gt <- gt[!duplicated(gt),]
  gt2 <- gt2[!duplicated(gt2),]
  dmrs <- full_join(dmrs,gt,by="ID")
  dmrs_dip <- full_join(dmrs_dip,gt2,by="ID")
  print(table(duplicated(dmrs_dip$ID)))
}

dmrs[is.na(dmrs)] <- "."
dmrs_dip[is.na(dmrs_dip)] <- "./."

CHROM <- as.integer(gsub("_.*","",dmrs$ID))
POS <- as.integer(gsub(".*_","",dmrs$ID))
ID <- dmrs$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "D"
meta$ALT <- "M"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs <- cbind(meta,dmrs[2:ncol(dmrs)])
dmrs$POS <- as.integer(dmrs$POS)
dmrs <- dmrs[order(dmrs$CHROM,dmrs$POS),]
colnames(dmrs)[1] <- "#CHROM"
write.table(dmrs,file="ibnr_dmrs.txt",row.names = F, quote = F,sep="\t")

CHROM <- as.integer(gsub("_.*","",dmrs_dip$ID))
POS <- as.integer(gsub(".*_","",dmrs_dip$ID))
ID <- dmrs_dip$ID
meta <- as.data.frame(cbind(CHROM,POS,ID))
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "."
meta$FORMAT <- "GT"
dmrs_dip <- cbind(meta,dmrs_dip[2:ncol(dmrs_dip)])
dmrs_dip$POS <- as.integer(dmrs_dip$POS)
dmrs_dip <- dmrs_dip[order(dmrs_dip$CHROM,dmrs_dip$POS),]
colnames(dmrs_dip)[1] <- "#CHROM"
write.table(dmrs_dip,file="ibnr_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

lengths <- lengths[!duplicated(lengths),]
write.table(lengths,file="ibnr_dmrs_lengths.txt",row.names = F, quote = F,sep="\t")

```
