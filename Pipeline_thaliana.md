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
meta$`#CHROM` <- as.numeric(meta$`#CHROM`)
meta$POS <- as.numeric(meta$POS)
meta <- meta[order(meta$`#CHROM`,meta$POS),]
write.table(cbind(meta,cov_context3),file="arabidopsis_meth.vcf",quote = F, row.names = F,sep="\t")

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
meta2$`#CHROM` <- as.numeric(meta2$`#CHROM`)
meta2$POS <- as.numeric(meta2$POS)
meta2 <- meta2[order(meta2$`#CHROM`,meta2$POS),]
write.table(cbind(meta2,cov_context4),file="arabidopsis_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")

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
fasta <- "Arabidopsis_thaliana.TAIR10.dna.jDMR.fa" 
samplefile_ceu <- "ceu_dmrsamples2" #
samplefile_ibnr <- "ibnr_dmrsamples2" 
runjDMRregions(fasta.file = fasta, out.dir = out.dir, samplefiles = samplefile_ceu, genome = "A_thaliana", contexts = "CG")
runjDMRregions(fasta.file = fasta, out.dir = out.dir, samplefiles = samplefile_ibnr, genome = "A_thaliana", contexts = "CG")
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

filenames <- read.table(file="../ceu_dmrsamples2",header = T)
filenames$file <- gsub("_methylome.txt","_methylome_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
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
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs$context <- "D"
dmrs$posteriorMax <- "M"
dmrs <- dmrs[1:5]
colnames(dmrs) <- c("#CHROM","POS","ID","REF","ALT")
dmrs$QUAL <- 4000
dmrs$FILTER <- "PASS"
dmrs$INFO <- "."
dmrs$FORMAT <- "GT"
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)
for (i in 2:nrow(filenames)){
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  dmrs <- cbind(dmrs,gt)
  dmrs_dip <- cbind(dmrs_dip,gt2)
}
dmrs_dip$REF <- "A"
dmrs_dip$ALT <- "T"
write.table(dmrs[dmrs$`#CHROM` %in% seq(1,5),],file="ceu_dmrs.txt",row.names = F, quote = F,sep="\t")
write.table(dmrs_dip[dmrs_dip$`#CHROM` %in% seq(1,5),],file="ceu_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

filenames <- read.table(file="../ibnr_dmrsamples2",header = T)
filenames$file <- gsub("_methylome.txt","_methylome_CG.txt",filenames$file)
dmrs <- read.table(file=filenames[1,1],header=T)
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
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs$context <- "D"
dmrs$posteriorMax <- "M"
dmrs <- dmrs[1:5]
colnames(dmrs) <- c("#CHROM","POS","ID","REF","ALT")
dmrs$QUAL <- 4000
dmrs$FILTER <- "PASS"
dmrs$INFO <- "."
dmrs$FORMAT <- "GT"
dmrs_dip <- cbind(dmrs,gt2)
dmrs <- cbind(dmrs,gt)
for (i in 2:nrow(filenames)){
  dmrs_tmp <- read.table(file=filenames[i,1],header=T)
  gt <- dmrs_tmp$status
  gt <- gsub("U",0,gt)
  gt <- gsub("M",1,gt)
  gt[is.na(gt)] <- "."
  gt <- as.data.frame(gt)
  gt2 <- gt
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "/")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  dmrs <- cbind(dmrs,gt)
  dmrs_dip <- cbind(dmrs_dip,gt2)
}
dmrs_dip$REF <- "A"
dmrs_dip$ALT <- "T"
write.table(dmrs[dmrs$`#CHROM` %in% seq(1,5),],file="ibnr_dmrs.txt",row.names = F, quote = F,sep="\t")
write.table(dmrs_dip[dmrs_dip$`#CHROM` %in% seq(1,5),],file="ibnr_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

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
```

27. Estimate genetic distances using genic SNPs

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

28. Calculate allele frequencies and linkage disequlibrium using intergenic SNPs

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ceu_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed --recode  --out ceu_snps_intergenic
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf ibnr_snps.recode.vcf.gz --bed /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/genomes/non_gene.bed --recode  --out ibnr_snps_intergenic

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_snps_intergenic.recode.vcf --out ceu_snps_intergenic --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_snps_intergenic.recode.vcf --out ibnr_snps_intergenic --min-alleles 2 --max-alleles 2 --max-missing 0.8 --freq

/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_snps_intergenic.recode.vcf -OutStat ceu_snps_5mb_cg_intergenic_decay  -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_snps_intergenic.recode.vcf -OutStat ibnr_snps_5mb_cg_intergenic_decay  -MaxDist 10

```

29. Estimate genetic distances using intergenic SNPs

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

30. Thin SNP VCFs for DAPC inference

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_snps_intergenic_maf_thin --vcf ceu_ibnr_snps_intergenic_maf.recode.vcf --thin 5000
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out ceu_ibnr_snps_genes_maf_thin --vcf ceu_ibnr_snps_genes_maf.recode.vcf --thin 5000

```
31. Filter SMP vcfs for CEU and IBNR

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
32. Generate individual SMP vcfs for each interval.
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

33. Example good_intervals.R
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

34. Estimate per region (gene region or clock like cg) Pi and Tajima's D for SMP.
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

35. Example good_intervals.R from point 31.
```
library(dplyr)
good_intervals <- read.table(file="good_intervals")
alpha_dm <- read.table(file="alpha_Dm")
alpha_dm$V1 <- gsub(".input.txt","",alpha_dm$V1)
alpha_dm <- alpha_dm[1:2]
good_intervals <- inner_join(good_intervals,alpha_dm,by="V1")
write.table(good_intervals,file="good_intervals_alpha",sep="\t",quote=F,row.names = F, col.names = F)
```

36. Allele frequencies for SMPs and DMRs
```
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

37. Calculate LD using SMP and DMR
```
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_meth_cg.recode.vcf -OutStat ceu_meth_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_meth_cg.recode.vcf -OutStat ibnr_meth_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ceu_dmrs_cg.recode.vcf -OutStat ceu_dmrs_cg_ld_decay -MaxDist 10
#/proj/popgen/a.ramesh/software/PopLDdecay/PopLDdecay -InVCF ibnr_dmrs_cg.recode.vcf -OutStat ibnr_dmrs_cg_ld_decay -MaxDist 10
```

38. Calculate distance matrices for SMPs and DMRs
```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_meth_cg.recode.vcf.gz --out ceu_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_meth_cg.recode.vcf.gz --out ibnr_meth_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_meth_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_meth_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_meth_cg_maf.recode.vcf ceu_meth_cg_maf.recode.vcf.gz ibnr_meth_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_meth_cg_maf.recode.vcf -OutPut ceu_ibnr_meth_cg_dis.mat

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ceu_dmrs_cg.recode.vcf.gz --out ceu_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf ibnr_dmrs_cg.recode.vcf.gz --out ibnr_dmrs_cg_maf --max-missing 0.8  --maf 0.02 --recode
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ceu_dmrs_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ceu_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f ibnr_dmrs_cg_maf.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f ibnr_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/bcftools-1.16/bcftools merge -m none -Ov -o ceu_ibnr_dmrs_cg_maf.recode.vcf ceu_dmrs_cg_maf.recode.vcf.gz ibnr_dmrs_cg_maf.recode.vcf.gz
/proj/popgen/a.ramesh/software/VCF2Dis-1.50/bin/VCF2Dis -InPut ceu_ibnr_dmrs_cg_maf.recode.vcf -OutPut ceu_ibnr_dmrs_cg.mat
```

39. Get multihetsep files for demographic inference for SNPs
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
40. To create the file combining sample names and chr for the inference above. In R.

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

41. Make the multihetsep files haploid. No hets in file. In R (multihetsep_combined.R).
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


42. Get list of cytosine sites. In R.

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


43. Get multihetsep files for demographic inference for SMPs
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

44. Make the multihetsep files haploid. No hets in file. In R. (multihetsep_combined.R).
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

45. Get multihetsep files for demographic inference for 5mb regions for SNPs

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

46. Make the multihetsep files haploid. For SNPs in 5mb No hets in file. In R. (multihetsep_combined.R). Then combine with methylation 5Mb files.

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
47. Estimate methylation rates using SMCm 

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

