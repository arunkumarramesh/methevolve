1. Get data

```
/data/proj/teaching/NGS_course/Softwares/sratoolkit.2.11.0-centos_linux64/bin/prefetch --option-file SRRlist_wgbs.txt -O data_wgbs

/data/proj/teaching/NGS_course/Softwares/sratoolkit.2.11.0-centos_linux64/bin/prefetch --option-file acc_list_wgs.txt -O data_wgs
```

2. Fastq-dump
>*for wgs data*
```
for file in SRR*; do 
/data/proj/teaching/NGS_course/Softwares/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done
```

>*for wgbs data*
```
for file in SRR*; do 
/data/proj/teaching/NGS_course/Softwares/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done
```
---
>*--split-3 splits .sra files into forward and reverse reads*\
>*--gzip option outputs files in gzipped .gz format*
---

3. Trim adapters
>*for wgs data (paired end)*
```
for file in *_1.fastq.gz; do java -jar /data/proj/teaching/NGS_course/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```
---
>*PE is paired end mode*\
>*-phred33 speciefies the phred score option (usually 33)*\
>*-threads 20 can be run on mulitple threads*\
>*requires 2 input files and generates 4 output files*\
>*ILLUMINACLIP: TruSeq3-PE.fa is fasta file for adapters in paired end mode*\
---
>*for wgbs data (single end)*
```
for file in *.fastq.gz; do java -jar /data/proj/teaching/NGS_course/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 20 $file ${file/.fastq.gz/.trimmed.fq.gz} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

4. Index genome
>*download from NCBI*
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005505.3/

>*index*
```
picard CreateSequenceDictionary -R GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna -O GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.dict
```
5. Map wgs data to reference genome
```
for file in *_1.paired.fq.gz; do bwa mem -t 24 /data/home/students/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} 2>paired_map_err | samtools view -S -b > ${file/_1.paired.fq.gz/_mapped.bam} 2>paired_map_err; done
```

6. Add readgroups, filter reads, sort reads, merge reads by accession using markduplicates, index bams
```
for file in *_mapped.bam ; do picard AddOrReplaceReadGroups -I $file -O ${file/_mapped.bam/_readgroup.bam} -LB brachypodium_distachyon -PL illumina -PU 1 -SM $file; done
for file in *_readgroup.bam; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b $file >${file/_readgroup.bam/_mq20.bam}; done
for file in *_mq20.bam; do samtools sort -@ 5 $file >${file/_mq20.bam/_sort.bam} ; done
# for each accession provide all SRR.bam files, repeat for all accessions
picard MarkDuplicates -I SRR8190034_sort.bam -I SRR8190670_sort.bam -O BdTR1a_marked.bam -M BdTR1a_metrics.txt
for file in *_marked.bam ; do picard AddOrReplaceReadGroups -I $file -O ${file/_marked.bam/_marked_readgroup.bam} -LB brachypodium_distachyon -PL illumina -PU 1 -SM $file; done
for file in *_marked_readgroup.bam ; do picard BuildBamIndex -I $file; done
```

7. Call variants
>*In R, create submission script for variant calling*
```
#variant_call.R

names <- read.table(file = "samplenames30")

output_file <- "variant_call.txt"
file_conn <- file(output_file, "w")

# add bash stuff
cat("#$-cwd", file = file_conn)
cat("\n", file = file_conn)
cat("#$-pe serial 1", file = file_conn)
cat("\n", file = file_conn)
cat("#$-l vf=1g", file = file_conn)
cat("\n", file = file_conn)
cat("#$-N variant_call", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("source ~/.bashrc", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("cd /data/proj2/popgen/a.zauchner/brachy/data_wgs/bam/", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)


for (i in 1:nrow(names)){
    cat("gatk HaplotypeCaller -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna", 
        paste("-I ", names[1][i,], "_marked_readgroup.bam", sep = ""), 
        paste("-O ", names[1][i,], ".g.vcf.gz -ERC GVCF", sep = ""), 
        "\n", 
        file = file_conn)
}

close(file_conn)

# don't forget to dos2unix using vim combine_gvcfs.txt -c "set ff=unix" -c ":wq"
```
>*command (exammple for first 2 samples)*
```
cd /data/proj2/popgen/a.zauchner/brachy/data_wgs/bam/

gatk HaplotypeCaller -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna -I Bd21_marked_readgroup.bam -O Bd21.g.vcf.gz -ERC GVCF 
gatk HaplotypeCaller -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna -I BdTR10n_marked_readgroup.bam -O BdTR10n.g.vcf.gz -ERC GVCF 
```
8. Combine GVCFs
>*In R, write script to combine samples*
```
# combine.R
names <- read.table(file = "samplenames30")

output_file <- "combine_GVCFs.txt"
file_conn <- file(output_file, "w")

# add bash stuff
cat("#$-cwd", file = file_conn)
cat("\n", file = file_conn)
cat("#$-pe serial 1", file = file_conn)
cat("\n", file = file_conn)
cat("#$-l vf=1g", file = file_conn)
cat("\n", file = file_conn)
cat("#$-N combine_GVCFs", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("source ~/.bashrc", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("cd /data/proj2/popgen/a.zauchner/brachy/data_wgs/bam/", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)

cat("gatk CombineGVCFs -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna",
    file = file_conn)

for (i in 1:nrow(names)){
      cat(" --variant ", names[1][i,], ".g.vcf.gz", sep = "",
      file = file_conn)
}

cat(" -O brachypodium.all.g.vcf.gz", file = file_conn)

close(file_conn)

# don't forget to dos2unix!
```
>*command*
```
cd /data/proj2/popgen/a.zauchner/brachy/data_wgs/bam/

gatk CombineGVCFs -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna --variant Bd21.g.vcf.gz --variant BdTR10n.g.vcf.gz --variant BdTR11c.g.vcf.gz --variant BdTR11d.g.vcf.gz --variant BdTR11f.g.vcf.gz --variant BdTR11g.g.vcf.gz --variant BdTR11h.g.vcf.gz --variant BdTR11i.g.vcf.gz --variant BdTR1a.g.vcf.gz --variant BdTR1b.g.vcf.gz --variant BdTR1e.g.vcf.gz --variant BdTR1f.g.vcf.gz --variant BdTR1g.g.vcf.gz --variant BdTR1h.g.vcf.gz --variant BdTR1j.g.vcf.gz --variant BdTR1k.g.vcf.gz --variant BdTR1m.g.vcf.gz --variant BdTR1n.g.vcf.gz --variant BdTR2b.g.vcf.gz --variant BdTR2c.g.vcf.gz --variant BdTR2d.g.vcf.gz --variant BdTR2g.g.vcf.gz --variant BdTR2h.g.vcf.gz --variant BdTR2j.g.vcf.gz --variant BdTR2k.g.vcf.gz --variant BdTR2m.g.vcf.gz --variant BdTR2n.g.vcf.gz --variant BdTR2p.g.vcf.gz --variant BdTR2r.g.vcf.gz --variant BdTR2s.g.vcf.gz -O brachypodium.all.g.vcf.gz
```

9. Genotype GVCFs 
```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna -V brachypodium.all.g.vcf.gz -O brachypodium.output.vcf.gz
```
10. Filter vcfs and count variants
```
gatk SelectVariants -V brachypodium.output.vcf.gz -select-type SNP -O brachypodium.snps.vcf.gz
bcftools view -i 'FMT/DP>=5 & QUAL>=30' brachypodium.snps.vcf > brachypodium.snps.filtered.vcf
```
```
bcftools stats brachypodium_snps_filtered.recode.vcf > brachypodium_snps_filtered.recode.stats
```
11. Create pseudoreferences
>*In R, create list with sample name and chromosome identifyer*
```
# pseudorefs.R
# make list with sample name and chromosome

#wrote sample names to list
#file1 <- unique(read.table(file="list3.tsv")[1])
#file1 <- rbind(file1, "BdTR1j")
#file1 <- file1[order(file1$V1), , drop = FALSE]
#write.table(file1, file = "samplenames30", quote = F, sep = "\t", row.names = F, col.names = F)

file1 <- read.table(file = "samplenames30")
for (i in 1:nrow(file1)){
  file1[[1]][i] <- paste(file1[[1]][i], "_marked.bam", sep = "")
}

chrlist <- c("NC_016131.3", "NC_016132.3", "NC_016133.3", "NC_016134.3", "NC_016135.3", "NC_011032.1")
file2 <- data.frame(chrlist)

fileall <- expand.grid(file1$V1, file2$chrlist)  # expand.grid creates data frame from all combinations
write.table(fileall, file = "sample_chr", quote = F, sep = "\t", row.names = F, col.names = F)
```

>*Use file to make 5 .fa files for the 5 chromosomes + Plastid containing samples with variants applied*
```
bgzip brachypodium_snps_filtered.recode.vcf
tabix brachypodium_snps_filtered.recode.vcf.gz
cat sample_chr | while read -r value1 value2 remainder ; do samtools faidx /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna $value2 | bcftools consensus -M N -s $value1 -p ${value1%marked.bam} -H 1pIu brachypodium_snps_filtered.recode.vcf.gz >> ${value2}_brachypodium.fa; done
```
>*Use faSplit to split .fa files by sample names*
```
faSplit byname NC_016131.3_brachypodium.fa
```
>*merge and index (samplenames30 is list of samplenames e.g. Bd21, BDTR11f, ...)* 
```
cat samplenames30 | while read line ; do cat $line*.fa > $line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do picard CreateSequenceDictionary -R $file -O ${file/.fa/.dict}; done
```
>*move into new folders*
```
cat samplenames30 | while read line ; do mkdir $line  ; done
cat samplenames30 | while read line ; do mv $line.merged.fa $line  ; done
cat samplenames30 | while read line ; do mv $line.merged.dict $line  ; done
cat samplenames30 | while read line ; do mv $line.merged.fa.fai $line  ; done
```

12. Move to WGBS data and index genomes for methylation calling
```
cat samplenames30 | while read line ; do /data/proj2/popgen/a.zauchner/softwares/Bismark-0.24.1/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /data/proj2/popgen/a.zauchner/softwares/hisat2-2.2.1/ /data/proj2/popgen/a.zauchner/brachy/refs/$line ; done
```

13. Combine SRR fastq files of each accession
>*In R, write file with one command line per accession using "SRR_sample_wgbs" which is a list with two columns containing the SRR number of the read in column 1 and the corresponding sample name in column 2*
```
SRR_sample <- read.table(file="SRR_sample_wgbs")

#write file with commands cat SRR1.trimmed.fq.gz SRR2.trimmed.fq.gz ... > name.wgbs.fq.gz

output_file <- "combine_fq_wgbs.txt"
file_conn <- file(output_file, "w")

#add bash stuff
cat("#$-cwd", file = file_conn)
cat("\n", file = file_conn)
cat("#$-pe serial 8", file = file_conn)
cat("\n", file = file_conn)
cat("#$-l vf=1g", file = file_conn)
cat("\n", file = file_conn)
cat("#$-N combine_fq_wgbs", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("source ~/.bashrc", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)
cat("cd /data/proj2/popgen/a.zauchner/brachy/data_wgbs/", file = file_conn)
cat("\n", file = file_conn)
cat("\n", file = file_conn)


for(i in 1:30){
  cat("cat ", file = file_conn)
  name <- unique(SRR_sample$V2)[i]
  df_name <- SRR_sample$V1[SRR_sample$V2 == name]
  
  for (j in 1:length(df_name)){
    cat(
      paste(df_name[j], ".trimmed.fq.gz ", sep = ""),
      file = file_conn)
  }
  cat(paste("> ", name, ".wgbs.fq.gz", sep = ""),
      "\n",
      file = file_conn)
}

close(file_conn)

# don't forget to dos2unix!
```
>*Command in bash (only 2 of 30 lines are shown)*
```
cd /data/proj2/popgen/a.zauchner/brachy/data_wgbs/

cat SRR5050864.trimmed.fq.gz SRR5032129.trimmed.fq.gz SRR5046527.trimmed.fq.gz SRR5050902.trimmed.fq.gz SRR5032126.trimmed.fq.gz SRR5046521.trimmed.fq.gz SRR5273351.trimmed.fq.gz SRR5032109.trimmed.fq.gz SRR5046553.trimmed.fq.gz SRR5481677.trimmed.fq.gz SRR5481675.trimmed.fq.gz SRR5273493.trimmed.fq.gz > BdTR11h.wgbs.fq.gz 
cat SRR5050907.trimmed.fq.gz SRR5032153.trimmed.fq.gz SRR5046482.trimmed.fq.gz SRR5273422.trimmed.fq.gz SRR5032075.trimmed.fq.gz SRR5046477.trimmed.fq.gz SRR5273355.trimmed.fq.gz SRR5032071.trimmed.fq.gz SRR5046555.trimmed.fq.gz SRR5481676.trimmed.fq.gz SRR5273562.trimmed.fq.gz SRR5273496.trimmed.fq.gz > BdTR2c.wgbs.fq.gz 
```

14. Map single-end reads to pseudoreferences
```
cat samplenames30 |  while read line ; do /data/proj2/popgen/a.zauchner/softwares/Bismark-0.24.1/bismark --multicore 24 --hisat2 --path_to_hisat2 /data/proj2/popgen/a.zauchner/softwares/hisat2-2.2.1/ --genome_folder /data/proj2/popgen/a.zauchner/brachy/refs/$line/ $line.wgbs.fq.gz  ; done
```

15. Deduplicate methylation reads
```
for file in *.wgbs_bismark_hisat2.bam; do /data/proj2/popgen/a.zauchner/softwares/Bismark-0.24.1/deduplicate_bismark --bam $file ; done
```

16. Call methylation variants
```
cat samplenames30 |  while read line ; do /data/proj2/popgen/a.zauchner/softwares/Bismark-0.24.1/bismark_methylation_extractor --multicore 24 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder /data/proj2/popgen/a.zauchner/brachy/refs/$line $line.wgbs_bismark_hisat2.deduplicated.bam ; done
```

17. Binomial test
```
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")


cov <- read.table(file=as.character(covfile$V1[1]), header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
print(cov$ID[1])

context <- read.table(file=as.character(contextfiles$V1[1]),header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "-")
print(context$ID[1])

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub("_N","N", gsub(gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile$V1)[1],"",cov_context$chromosome))
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 9, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.total)

print(noncoversionrate)

rm(cov)
rm(context)


b <- apply(cov_context[c(5,11)], 1, binom.test, p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
dim(cov_context)

rm(b)

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]
dim(cov_context)

samplename <- gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile[1,])
colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

rm(cov_context)

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  cov <- read.table(file=as.character(covfile[i,]), header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
  print(cov$ID[1])

  context <- read.table(file=as.character(contextfiles[i,]),header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "-")
  print(context$ID[1])
  print(i)

  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub("_N","N", gsub(gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile$V1)[i],"",cov_context$chromosome))
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 9, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
  cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.total)
  
  print(noncoversionrate)

  rm(cov)
  rm(context)

  
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
  dim(cov_context)

  rm(b)

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  dim(cov_context)

  samplename <- gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile[i,])
  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
  dim(cov_context2)
  rm(cov_context)
}

#### write to file for further processing

#cov_context2$chromosome <- gsub("-.*","",cov_context2$ID)
#cov_context2$position   <- as.numeric(gsub(".*-","",cov_context2$ID))
write.table(cov_context2,file="cov_context3_all.txt",row.names = F)

```

18. Make vcfs
```
library(dplyr)
library(ggplot2)
library(data.table)

cov_context3 <- fread(file="cov_context3_all.txt",header=T)
na_count <- apply(cov_context3[,6:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[,6:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[,6:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[,1:3]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"
cov_context3 <- cov_context3[,-c(1:5)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
fwrite(cbind(meta,cov_context3),file="brachy_meth.vcf",quote = F, row.names = F,sep="\t")

meta2 <- cov_context4[,1:3]
colnames(meta2) <- c("#CHROM","POS","ID")
meta2$REF <- "A"
meta2$ALT <- "T"
meta2$QUAL <- 4000
meta2$FILTER <- "PASS"
meta2$INFO <- "DP=1000"
meta2$FORMAT <- "GT"
cov_context4 <- cov_context4[,-c(1:5)]
cov_context4[cov_context4 == "U"] <- "0/0"
cov_context4[cov_context4 == "M"] <- "1/1"
cov_context4[is.na(cov_context4)] <- "./."
fwrite(cbind(meta2,cov_context4),file="brachy_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")
```

19. get DMRs with jDMR
```
## Change chromosome names, consistency necessary (namechange_for_jDMR.R)

# this reads in the CpG report files from bismark and changes the chromosome names from "BdTR11i_NC_016132.3" to "NC_016132.3" and writes them as a new file with a different name "BdTR11i_CpG_report.txt" -> input for methimpute

library(data.table)

contextfiles <- read.table(file="contextfiles")

for(i in 1:nrow(contextfiles)){
data <- fread(contextfiles[i,], header = FALSE)
name <- as.character(gsub(".wgbs_bismark_hisat2.deduplicated.CpG_report.txt","_",contextfiles[i,]))
data$V1 <- gsub(name,"", data$V1)

out.filename <- as.character(gsub(".wgbs_bismark_hisat2.deduplicated.CpG_report","_CpG_report",contextfiles[i,]))
fwrite(data, file = out.filename, eol = "\n", sep = "\t", row.names = F, col.names = F)
}

#data_after <- fread("BdTR11i_CpG_report.txt", header = TRUE)

### run jDMR

library(devtools)
library(methimpute) #for cluster, load package from local account with library(methimpute, lib.loc = ...)
library(jDMR)
library(data.table)
  
##read file
contextfiles <- read.table(file="contextfiles")
reportfiles <- read.table(file = "reportfiles") # reportfiles are context files with chromosome names adjusted (see namechange_for_jDMR.R)

## without inflateMethylome

# get methylome with methimpute
# i <-1 # run line by line if memory allocation error occurs

for(i in 1:nrow(reportfiles)){
bismark.data <- importBismark(reportfiles[i,])
#bismark.data <- importBismark("BdTR11i_CpG_report.txt")

distcor <- distanceCorrelation(bismark.data, separate.contexts = TRUE)
fit <- estimateTransDist(distcor)
#print(fit)
model <- callMethylationSeparate(data = bismark.data, transDist = fit$transDist, verbosity = 0)
filename <- as.character(gsub(".wgbs_bismark_hisat2.deduplicated.CpG_report","_methylome",contextfiles[i,]))
exportMethylome(model, filename = filename)
}

#move methylome files to /methylome/ folder

# run jDMR
library("methimpute",lib.loc="/data/proj2/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.3/")
library("jDMR",lib.loc="/data/proj2/home/users/a.ramesh/R/x86_64-redhat-linux-gnu-library/4.3/")
library(Biostrings)
library(data.table)

out.dir <- "jDMRresults/"

myfasta <- readDNAStringSet("GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff.fna")
CfromFASTAv4(fasta = myfasta, chr = "NC_016131.3", out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", "NC_016131.3", ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = "NC_016131.3", min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff.fna")
CfromFASTAv4(fasta = myfasta, chr = "NC_016132.3", out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", "NC_016132.3", ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = "NC_016132.3", min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff.fna")
CfromFASTAv4(fasta = myfasta, chr = "NC_016133.3", out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", "NC_016133.3", ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = "NC_016133.3", min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff.fna")
CfromFASTAv4(fasta = myfasta, chr = "NC_016134.3", out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", "NC_016134.3", ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = "NC_016134.3", min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")

myfasta <- readDNAStringSet("GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff.fna")
CfromFASTAv4(fasta = myfasta, chr = "NC_016135.3", out.dir = out.dir, write.output = TRUE)
ref.genome <- fread(paste0(out.dir, "/cytosine_positions_chr", "NC_016135.3", ".csv", sep = ""))
makeReg(ref.genome = ref.genome, contexts = c("CG"), makeRegnull = c(FALSE), chr = "NC_016135.3", min.C = 8, N.boot = 10^5, N.sim.C = "all", fp.rate = 0.01, set.tol = 0.01, out.dir = out.dir, out.name = "Arabidopsis")


runMethimputeRegions(out.dir = out.dir, samplefiles = "samplefiles", genome = "B_dist", context = "CG", nCytosines=8, mincov=3, Regionfiles = out.dir)
```

20. Make vcf from jDMR results
```
library(dplyr)

filenames <- read.table(file="../samplefiles",header = T)
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
gt2$gt <- paste(gt2$gt,gt2$gt,sep = "|")
colnames(gt)[1] <- filenames[1,2]
colnames(gt2)[1] <- filenames[1,2]
lth <- dmrs$end - dmrs$start
dmrs$start <- round((dmrs$start + dmrs$end)/2)
dmrs$end <- paste(dmrs$seqnames,dmrs$start,sep="_")
dmrs <- dmrs[3]
colnames(dmrs) <- "ID"
dmrs$ID <- gsub("NC_01613","",dmrs$ID)
dmrs$ID <- gsub("\\.3","",dmrs$ID)
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
  gt2$gt <- paste(gt2$gt,gt2$gt,sep = "|")
  colnames(gt)[1] <- filenames[i,2]
  colnames(gt2)[1] <- filenames[i,2]
  gt <- cbind(gt,dmrs_tmp$end)
  colnames(gt)[2] <- "ID"
  gt$ID <- gsub("NC_01613","",gt$ID)
  gt$ID <- gsub("\\.3","",gt$ID)
  gt2 <- cbind(gt2,dmrs_tmp$end)
  colnames(gt2)[2] <- "ID"
  gt2$ID <- gsub("NC_01613","",gt2$ID)
  gt2$ID <- gsub("\\.3","",gt2$ID)
  gt <- gt[!duplicated(gt),]
  gt2 <- gt2[!duplicated(gt2),]
  dmrs <- full_join(dmrs,gt,by="ID")
  dmrs_dip <- full_join(dmrs_dip,gt2,by="ID")
  print(table(duplicated(dmrs_dip$ID)))
}

dmrs[is.na(dmrs)] <- "."
dmrs_dip[is.na(dmrs_dip)] <- ".|."

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
write.table(dmrs,file="bd_dmrs.txt",row.names = F, quote = F,sep="\t")

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
write.table(dmrs_dip,file="bd_dmrs_dip.txt",row.names = F, quote = F,sep="\t")

lengths <- lengths[!duplicated(lengths),]
write.table(lengths,file="bd_dmrs_lengths.txt",row.names = F, quote = F,sep="\t")
```

21. Make vcfs for group 1, group 6 and all 29 (30 minus reference Bd21) for SNPs, SMPs and DMRs each
```
# get list with samplenames and group
setwd("G:/Meine Ablage/Master/RStudio/Thesis/")
df1 <- fread("group1_names", header = F, drop = "V1")
df6 <- fread("group6_names", header = F, drop = "V1")

dmr_30 <- fread("bd_dmrs_dip.vcf", skip = "##")
dmr_29 <- dmr_30[,-10]
colnames(dmr_29) %in% df1$V2
#add header

dmr_grp1 <- dmr_29[,-c(10:16)]
colnames(dmr_grp1) %in% df1$V2
fwrite(dmr_grp1,file="bd_dmrs_dip1.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")
#add header

dmr_grp6 <- dmr_29[,-c(17:38)]
colnames(dmr_grp6) %in% df6$V2
fwrite(dmr_grp6,file="bd_dmrs_dip6.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")
#add header


### snps
setwd("C:/Users/Excellaptop/ubuntu/vcfs/popLDecay")

snp <- fread("brachypodium.snps.filtered.vcf", header = T)
snp$FILTER <- gsub(".", "PASS", snp$FILTER)

snp_29 <- snp[,-10]
colnames(snp_29)[10:length(colnames(snp_29))] <- gsub("_marked.bam", "", colnames(snp_29)[10:length(colnames(snp_29))])
fwrite(snp_29, file = "snp_29.vcf", quote = F, row.names = F,sep ="\t", eol = "\n")

colnames(snp_29) %in% df1$V2

snp_grp1 <- snp_29[,-c(10:16)]
colnames(snp_grp1) %in% df1$V2
fwrite(snp_grp1,file="snp_grp1.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

snp_grp6 <- snp_29[,-c(17:38)]
colnames(snp_grp6) %in% df6$V2
fwrite(snp_grp6,file="snp_grp6.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

setwd("C:/Users/Excellaptop/Desktop/data/vcfs")

### smps

smp <- fread("brachy_meth.vcf", header = T)
smp_29 <- smp[,-11]
fwrite(smp_29, file = "smp_29.vcf", quote = F, row.names = F,sep ="\t", eol = "\n")

colnames(smp_29) %in% df1$V2

smp_grp1 <- smp_29[,-c(10:16)]
colnames(smp_grp1) %in% df1$V2
fwrite(smp_grp1,file="smp_grp1.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

smp_grp6 <- smp_29[,-c(17:38)]
colnames(smp_grp6) %in% df6$V2
fwrite(smp_grp6,file="smp_grp6.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

```
22. Prepare vcfs for filtering
>*vcf_cleanup.R, Change chromosomes names in vcfs to 1,2, etc., add headers from snp vcfs and rewrite to ..._unfiltered.vcf*
>*count variants*
```
grep -v '^#' snp_29_unfiltered.vcf | wc -l
```
>*change chromosome names in bed file to 1,2, etc.*

23. Replace heterozygoes with NA in snp vcf files
```
# snp_nohet.R
library(data.table)
library(dplyr)

snp29.unfiltered <- fread(input = "snp29_unfiltered.vcf")
snp1.unfiltered <- fread(input = "snp1_unfiltered.vcf")
snp6.unfiltered <- fread(input = "snp6_unfiltered.vcf")


snp6 <- snp6.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

snp1 <- snp1.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

snp29 <- snp29.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

fwrite(snp1, file = "snp1_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)
fwrite(snp6, file = "snp6_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)
fwrite(snp29, file = "snp29_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)

```

24. filter with vcftools and plink ---> see step 48 for all final filtering steps
```
### see step 48

# filter max 20% NA, only gene body
vcftools --vcf snp29_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gene2.bed --out snp29_na --recode-INFO-all
vcftools --vcf snp1_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gene2.bed --out snp1_na --recode-INFO-all
vcftools --vcf snp6_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gene2.bed --out snp6_na --recode-INFO-all
 
vcftools --vcf smp29_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out smp29_na --recode-INFO-all
vcftools --vcf smp1_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out smp1_na --recode-INFO-all
vcftools --vcf smp6_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out smp6_na --recode-INFO-all

vcftools --vcf dmr29_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out dmr29_na --recode-INFO-all
vcftools --vcf dmr1_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out dmr1_na --recode-INFO-all
vcftools --vcf dmr6_unfiltered.vcf --max-missing 0.8 --recode --bed gene2.bed --out dmr6_na --recode-INFO-all

# count
grep -v '^#' snp29_na.vcf | wc -l

# filter MAF with plink --> generates binary fileset

plink --vcf snp29_na.recode.vcf --maf 0.04 --make-bed --out ./plink/snp29_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp1_na.recode.vcf --maf 0.04 --make-bed --out ./plink/snp1_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp6_na.recode.vcf --maf 0.04 --make-bed --out ./plink/snp6_maf --set-missing-var-ids @:# --keep-allele-order

plink --vcf smp29_na.recode.vcf --maf 0.04 --make-bed --out ./plink/smp29_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1_na.recode.vcf --maf 0.04 --make-bed --out ./plink/smp1_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_na.recode.vcf --maf 0.04 --make-bed --out ./plink/smp6_maf --set-missing-var-ids @:# --keep-allele-order

plink --vcf dmr29_na.recode.vcf --maf 0.04 --make-bed --out ./plink/dmr29_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1_na.recode.vcf --maf 0.04 --make-bed --out ./plink/dmr1_maf --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_na.recode.vcf --maf 0.04 --make-bed --out ./plink/dmr6_maf --set-missing-var-ids @:# --keep-allele-order


# make vcfs from binary fileset

plink --bfile ./plink/snp29_maf --recode vcf-iid -out ./plink/snp29_maf --keep-allele-order
plink --bfile ./plink/snp1_maf --recode vcf-iid -out ./plink/snp1_maf --keep-allele-order
plink --bfile ./plink/snp6_maf --recode vcf-iid -out ./plink/snp6_maf --keep-allele-order

plink --bfile ./plink/smp29_maf --recode vcf-iid -out ./plink/smp29_maf --keep-allele-order
plink --bfile ./plink/smp1_maf --recode vcf-iid -out ./plink/smp1_maf --keep-allele-order
plink --bfile ./plink/smp6_maf --recode vcf-iid -out ./plink/smp6_maf --keep-allele-order

plink --bfile ./plink/dmr29_maf --recode vcf-iid -out ./plink/dmr29_maf --keep-allele-order
plink --bfile ./plink/dmr1_maf --recode vcf-iid -out ./plink/dmr1_maf --keep-allele-order
plink --bfile ./plink/dmr6_maf --recode vcf-iid -out ./plink/dmr6_maf --keep-allele-order

```
25. calculate distance matrix and plot trees
```
# calculate distance matrix for tree
/mnt/c/Users/Excellaptop/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./plink/snp29_maf.vcf  -OutPut ./tree/snp29_dis.mat
/mnt/c/Users/Excellaptop/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./plink/smp29_maf.vcf  -OutPut ./tree/smp29_dis.mat
/mnt/c/Users/Excellaptop/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./plink/dmr29_maf.vcf  -OutPut ./tree/dmr29_dis.mat


# http://www.atgc-montpellier.fr/fastme/
# make trees in R

tree.R
```
26. Calculate and plot LD decay
```
## popLDdecay

/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/snp29_maf.vcf -OutStat ./LDdecay/snp29_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/snp1_maf.vcf -OutStat ./LDdecay/snp1_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/snp6_maf.vcf -OutStat ./LDdecay/snp6_LDdecay.stat

/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/smp29_maf.vcf -OutStat ./LDdecay/smp29_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/smp1_maf.vcf -OutStat ./LDdecay/smp1_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/smp6_maf.vcf -OutStat ./LDdecay/smp6_LDdecay.stat

/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/dmr29_maf.vcf -OutStat ./LDdecay/dmr29_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/dmr1_maf.vcf -OutStat ./LDdecay/dmr1_LDdecay.stat
/mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/PopLDdecay -InVCF ./plink/dmr6_maf.vcf -OutStat ./LDdecay/dmr6_LDdecay.stat

## plot LDdecay 1kb window

perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/snp29_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_snp29_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/snp1_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_snp1_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/snp6_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_snp6_LDdecay_1kb

perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/smp29_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_smp29_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/smp1_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_smp1_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/smp6_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_smp6_LDdecay_1kb

perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/dmr29_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_dmr29_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/dmr1_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_dmr1_LDdecay_1kb
perl /mnt/c/Users/Excellaptop/ubuntu/softwares/PopLDdecay-3.42/bin/Plot_OnePop.pl -inFile ./LDdecay/dmr6_LDdecay.stat.gz -maxX 1 -output ./LDdecay/fig_dmr6_LDdecay_1kb

```
27. Calculate site-pi using vcftools
```
# calculate site-pi

vcftools --vcf ./plink/snp1_maf.vcf --site-pi --out ./pi/snp1_maf
vcftools --vcf ./plink/snp6_maf.vcf --site-pi --out ./pi/snp6_maf

vcftools --vcf ./plink/smp1_maf.vcf --site-pi --out ./pi/smp1_maf
vcftools --vcf ./plink/smp6_maf.vcf --site-pi --out ./pi/smp6_maf

vcftools --vcf ./plink/dmr1_maf.vcf --site-pi --out ./pi/dmr1_maf
vcftools --vcf ./plink/dmr6_maf.vcf --site-pi --out ./pi/dmr6_maf

```
28. Get dmr_pos.list and gene_pos.list in the format chr:start-end for all dmrs that are in gene regions (from bed file)
```
library(data.table)

dmr.reg <- fread(input = "DMR_all_pre.vcf", header = TRUE)
dmr.reg <- dmr.reg[!dmr.reg$seqnames =="NC_011032.1"] # remove plastid
dmr.reg$seqnames <- gsub("NC_016131.3", "1", dmr.reg$seqnames)
dmr.reg$seqnames <- gsub("NC_016132.3", "2", dmr.reg$seqnames)
dmr.reg$seqnames <- gsub("NC_016133.3", "3", dmr.reg$seqnames)
dmr.reg$seqnames <- gsub("NC_016134.3", "4", dmr.reg$seqnames)
dmr.reg$seqnames <- gsub("NC_016135.3", "5", dmr.reg$seqnames)
unique(dmr.reg$seqnames) # change chromosome names
  
dmr.reg <- dmr.reg[,c(1,2,3)]
dmr.reg$middle <- round(apply(dmr.reg[,c(2,3)],1, mean)) # calculate middle of dmr


# gene body bed file
genebody <- fread(input = "gene.bed")
genebody <- genebody[,-c(4,5,6)]
colnames(genebody) <- c("seqnames", "start", "stop")
genebody$seqnames <- gsub("NC_016131.3", "1", genebody$seqnames)
genebody$seqnames <- gsub("NC_016132.3", "2", genebody$seqnames)
genebody$seqnames <- gsub("NC_016133.3", "3", genebody$seqnames)
genebody$seqnames <- gsub("NC_016134.3", "4", genebody$seqnames)
genebody$seqnames <- gsub("NC_016135.3", "5", genebody$seqnames)
unique(genebody$seqnames)
genebody<- genebody[!genebody$seqnames =="NC_011032.1"] # remove plastid & other contigs
genebody<- genebody[!genebody$seqnames =="NW_019841143.1"]
genebody<- genebody[!genebody$seqnames =="NW_019841144.1"]
genebody<- genebody[!genebody$seqnames =="NW_019841145.1"]
genebody<- genebody[!genebody$seqnames =="NW_019841146.1"]
genebody<- genebody[!genebody$seqnames =="NW_019841147.1"]
unique(genebody$seqnames)

# only keep dmrs with middle point in gene body

get_pos <- function(keep) {
  return(which(keep))
}

genebody1 <- genebody[genebody$seqnames== 1]
genebody2 <- genebody[genebody$seqnames== 2]
genebody3 <- genebody[genebody$seqnames== 3]
genebody4 <- genebody[genebody$seqnames== 4]
genebody5 <- genebody[genebody$seqnames== 5]

dmr.reg1 <- dmr.reg[dmr.reg$seqnames==1]
dmr.reg2 <- dmr.reg[dmr.reg$seqnames==2]
dmr.reg3 <- dmr.reg[dmr.reg$seqnames==3]
dmr.reg4 <- dmr.reg[dmr.reg$seqnames==4]
dmr.reg5 <- dmr.reg[dmr.reg$seqnames==5]

#1
keep1 <- lapply(1:nrow(genebody1), function(i) {
data.table::between(dmr.reg1$middle, genebody1$start[i], genebody1$stop[i])
})
keep <- keep1
pos1 <- lapply(keep, get_pos)
head(pos1)

#2
keep2 <- lapply(1:nrow(genebody2), function(i) {
  data.table::between(dmr.reg2$middle, genebody2$start[i], genebody2$stop[i])
})
keep <- keep2
pos2 <- lapply(keep, get_pos)
head(pos2)


#3
keep3 <- lapply(1:nrow(genebody3), function(i) {
  data.table::between(dmr.reg3$middle, genebody3$start[i], genebody3$stop[i])
})
keep <- keep3
pos3 <- lapply(keep, get_pos)
head(pos3)

#4
keep4 <- lapply(1:nrow(genebody4), function(i) {
  data.table::between(dmr.reg4$middle, genebody4$start[i], genebody4$stop[i])
})
keep <- keep4
pos4 <- lapply(keep, get_pos)
head(pos4)

#5
keep5 <- lapply(1:nrow(genebody5), function(i) {
  data.table::between(dmr.reg5$middle, genebody5$start[i], genebody5$stop[i])
})
keep <- keep5
pos5 <- lapply(keep, get_pos)
head(pos5)


# take indices --> each list corresponds to elements in genebody e.g. first list in pos5 is first gene region in genebody chromosome 5
# numbers in list are indices of dmr.reg starting from each chromosome

# get gene regions where there are elements in the lists
# e.g. pos5 list 1 has elements inside --> use pos 1 in genebody 5


a5 <- vector()
for(i in 1:length(pos5)){
  if (lengths(pos5)[i]>0){
    a5 <- append(a5, pos5[[i]])
    print(i)
  }
}

a4 <- vector()
for(i in 1:length(pos4)){
  if (lengths(pos4)[i]>0){
    a4 <- append(a4, pos4[[i]])
    print(i)
  }
}

a3 <- vector()
for(i in 1:length(pos3)){
  if (lengths(pos3)[i]>0){
    a3 <- append(a3, pos3[[i]])
    print(i)
  }
}

a2 <- vector()
for(i in 1:length(pos2)){
  if (lengths(pos2)[i]>0){
    a2 <- append(a2, pos2[[i]])
    print(i)
  }
}

a1 <- vector()
for(i in 1:length(pos1)){
  if (lengths(pos1)[i]>0){
    a1 <- append(a1, pos1[[i]])
    print(i)
  }
}



c5 <- dmr.reg5[a5]
c4 <- dmr.reg4[a4]
c3 <- dmr.reg3[a3]
c2 <- dmr.reg2[a2]
c1 <- dmr.reg1[a1]



dmr_regions_in_gene <- rbind(c1,c2,c3,c4,c5)
fwrite(dmr_regions_in_gene, file = "dmr_regions_in_gene.txt", quote = F, row.names = F,sep ="\t", eol = "\n")

list <- dmr_regions_in_gene[, c(1,2,3)]
list$start_end <- paste(list$start, list$end, sep = "-")
list$chr_start_end <- paste(list$seqnames, list$start_end, sep = ":")
dmr_pos.list <- list$chr_start_end

write.table(dmr_pos.list, file = "dmr_pos.list", quote = F, row.names = F,sep ="\t", eol = "\n", append = TRUE, )

```
29. Split SNP vcfs by gene and dmr regions
```
bgzip -f snp1_maf.vcf
tabix -f snp1_maf.vcf.gz

bgzip -f snp6_maf.vcf
tabix -f snp6_maf.vcf.gz


# split vcfs by region
dos2unix dmr_pos.list
dos2unix gene_pos.list

### SNPs
# dmr positions
cat dmr_pos.list | while read -r line ; do tabix snp1_maf.vcf.gz $line >./group1_snp/dmr_regions/$line.snp.vcf; done
cat dmr_pos.list | while read -r line ; do tabix snp6_maf.vcf.gz $line >./group6_snp/dmr_regions/$line.snp.vcf; done

# gene positions
cat gene_pos.list | while read -r line ; do tabix snp1_maf.vcf.gz $line >./group1_snp/gene_regions/$line.snp.vcf; done
cat gene_pos.list | while read -r line ; do tabix snp6_maf.vcf.gz $line >./group6_snp/gene_regions/$line.snp.vcf; done


# count lines and write to file
# do for both groups
for file in *vcf ; do wc -l $file >>vcflengths_snp_dmr ; done
for file in *vcf ; do wc -l $file >>vcflengths_snp_gene ; done
```

30. For SNPs: get good intervals and remove bad interval files
```
#############################################################

#####   SNPs

#############################################################

####### GROUP 1

###### DMR regions
setwd("C:/Users/Excellaptop/Desktop/stats/snp/group1_snp/dmr_regions/")

vcflengths_snp_dmr <- read.table(file="vcflengths_snp_dmr")
vcflengths_snp_dmr$V2 <- gsub(".snp.vcf","",vcflengths_snp_dmr$V2)
colnames(vcflengths_snp_dmr) <- c("numvar","interval")
vcflengths_snp_dmr$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_dmr$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_dmr$interval)))
vcflengths_snp_dmr <- vcflengths_snp_dmr[vcflengths_snp_dmr$numvar < 3,]

write.table(paste("rm ",vcflengths_snp_dmr$interval,".snp.vcf",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)


vcflengths_snp_dmr <- read.table(file="vcflengths_snp_dmr")
vcflengths_snp_dmr$V2 <- gsub(".snp.vcf","",vcflengths_snp_dmr$V2)
colnames(vcflengths_snp_dmr) <- c("numvar","interval")
vcflengths_snp_dmr$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_dmr$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_dmr$interval)))
vcflengths_snp_dmr <- vcflengths_snp_dmr[vcflengths_snp_dmr$numvar >= 3,]

write.table(cbind(paste(vcflengths_snp_dmr$interval,".all.vcf",sep=""),vcflengths_snp_dmr$length),file="goodfiles",sep="\t",quote=F,row.names = F, col.names = F)
# got 0 here



###### GENE regions
setwd("C:/Users/Excellaptop/Desktop/stats/snp/group1_snp/gene_regions/")

vcflengths_snp_gene <- read.table(file="vcflengths_snp_gene")
vcflengths_snp_gene$V2 <- gsub(".snp.vcf","",vcflengths_snp_gene$V2)
colnames(vcflengths_snp_gene) <- c("numvar","interval")
vcflengths_snp_gene$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_gene$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_gene$interval)))
vcflengths_snp_gene <- vcflengths_snp_gene[vcflengths_snp_gene$numvar < 3,]

write.table(paste("rm ",vcflengths_snp_gene$interval,".snp.vcf",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)


vcflengths_snp_gene <- read.table(file="vcflengths_snp_gene")
vcflengths_snp_gene$V2 <- gsub(".snp.vcf","",vcflengths_snp_gene$V2)
colnames(vcflengths_snp_gene) <- c("numvar","interval")
vcflengths_snp_gene$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_gene$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_gene$interval)))
vcflengths_snp_gene <- vcflengths_snp_gene[vcflengths_snp_gene$numvar >= 3,]

write.table(cbind(paste(vcflengths_snp_gene$interval,".all.vcf",sep=""),vcflengths_snp_gene$length),file="goodfiles",sep="\t",quote=F,row.names = F, col.names = F)


#############################################################

####### GROUP 6

###### DMR regions
setwd("C:/Users/Excellaptop/Desktop/stats/snp/group6_snp/dmr_regions/")

vcflengths_snp_dmr <- read.table(file="vcflengths_snp_dmr")
vcflengths_snp_dmr$V2 <- gsub(".snp.vcf","",vcflengths_snp_dmr$V2)
colnames(vcflengths_snp_dmr) <- c("numvar","interval")
vcflengths_snp_dmr$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_dmr$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_dmr$interval)))
vcflengths_snp_dmr <- vcflengths_snp_dmr[vcflengths_snp_dmr$numvar < 3,]

write.table(paste("rm ",vcflengths_snp_dmr$interval,".snp.vcf",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)


vcflengths_snp_dmr <- read.table(file="vcflengths_snp_dmr")
vcflengths_snp_dmr$V2 <- gsub(".snp.vcf","",vcflengths_snp_dmr$V2)
colnames(vcflengths_snp_dmr) <- c("numvar","interval")
vcflengths_snp_dmr$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_dmr$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_dmr$interval)))
vcflengths_snp_dmr <- vcflengths_snp_dmr[vcflengths_snp_dmr$numvar >= 3,]

write.table(cbind(paste(vcflengths_snp_dmr$interval,".all.vcf",sep=""),vcflengths_snp_dmr$length),file="goodfiles",sep="\t",quote=F,row.names = F, col.names = F)
# 0 here


###### GENE regions
setwd("C:/Users/Excellaptop/Desktop/stats/snp/group6_snp/gene_regions/")

vcflengths_snp_gene <- read.table(file="vcflengths_snp_gene")
vcflengths_snp_gene$V2 <- gsub(".snp.vcf","",vcflengths_snp_gene$V2)
colnames(vcflengths_snp_gene) <- c("numvar","interval")
vcflengths_snp_gene$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_gene$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_gene$interval)))
vcflengths_snp_gene <- vcflengths_snp_gene[vcflengths_snp_gene$numvar < 3,]

write.table(paste("rm ",vcflengths_snp_gene$interval,".snp.vcf",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)


vcflengths_snp_gene <- read.table(file="vcflengths_snp_gene")
vcflengths_snp_gene$V2 <- gsub(".snp.vcf","",vcflengths_snp_gene$V2)
colnames(vcflengths_snp_gene) <- c("numvar","interval")
vcflengths_snp_gene$length <- as.numeric(gsub(".*-","", gsub(".*:","",vcflengths_snp_gene$interval))) - as.numeric(gsub("-.*","", gsub(".*:","",vcflengths_snp_gene$interval)))
vcflengths_snp_gene <- vcflengths_snp_gene[vcflengths_snp_gene$numvar >= 3,]

write.table(cbind(paste(vcflengths_snp_gene$interval,".all.vcf",sep=""),vcflengths_snp_gene$length),file="goodfiles",sep="\t",quote=F,row.names = F, col.names = F)

# got 2 here


### dos2unix bad_intervals.sh files and run in linux
# ./bad_intervals.sh

```
31. For SNPs: calculate D and Pi for regions
```
####### calculate Tajima's D and pi

# Group 1

setwd("/mnt/c/Users/Excellaptop/Desktop/stats/snp/group1_snp/gene_regions/")

goodfiles1 <- read.table(file="./goodfiles")
for(f in 1:nrow(goodfiles1)){
  i=0
  while(i<200){
    print(i)
    i <- i + 1
    system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles1$V1[f]," --out tmp --TajimaD ",goodfiles1$V2[f]+i,sep=""))
    tmpfile <- read.table(file="tmp.Tajima.D",header = T)
    print(nrow(tmpfile))
    if (nrow(tmpfile) == 1) {
      system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles1$V1[f]," --out ",goodfiles1$V1[f]," --TajimaD ",goodfiles1$V2[f]+i,sep=""))
      break
    }
  }
}

for(f in 1:nrow(goodfiles1)){
  i=0
  while(i<200){
    print(i)
    i <- i + 1
    system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles1$V1[f]," --out tmp --window-pi ",goodfiles1$V2[f]+i,sep=""))
    tmpfile <- read.table(file="tmp.windowed.pi",header = T)
    print(nrow(tmpfile))
    if (nrow(tmpfile) == 1) {
      system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles1$V1[f]," --out ",goodfiles1$V1[f]," --window-pi ",goodfiles1$V2[f]+i,sep=""))
      break
    }
  }
}

# Group 6

setwd("/mnt/c/Users/Excellaptop/Desktop/stats/snp/group6_snp/gene_regions/")

goodfiles6 <- read.table(file="./goodfiles")
for(f in 1:nrow(goodfiles6)){
  i=0
  while(i<200){
    print(i)
    i <- i + 1
    system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles6$V1[f]," --out tmp --TajimaD ",goodfiles6$V2[f]+i,sep=""))
    tmpfile <- read.table(file="tmp.Tajima.D",header = T)
    print(nrow(tmpfile))
    if (nrow(tmpfile) == 1) {
      system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles6$V1[f]," --out ",goodfiles6$V1[f]," --TajimaD ",goodfiles6$V2[f]+i,sep=""))
      break
    }
  }
}

for(f in 1:nrow(goodfiles6)){
  i=0
  while(i<200){
    print(i)
    i <- i + 1
    system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles6$V1[f]," --out tmp --window-pi ",goodfiles6$V2[f]+i,sep=""))
    tmpfile <- read.table(file="tmp.windowed.pi",header = T)
    print(nrow(tmpfile))
    if (nrow(tmpfile) == 1) {
      system(paste("/mnt/c/Users/Excellaptop/ubuntu/softwares/vcftools/src/cpp/vcftools --vcf ",goodfiles6$V1[f]," --out ",goodfiles6$V1[f]," --window-pi ",goodfiles6$V2[f]+i,sep=""))
      break
    }
  }
}
```
32. For SMPs: split reference in gene regions and dmr regions
```
cat ./gene_pos.list | while read -r line ; do samtools faidx GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff_chrnames.fna $line >>genes.fasta; done
# in conda environment
faSplit byname genes.fasta genes_fasta/

cat ./dmr_pos.list | while read -r line ; do samtools faidx GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff_chrnames.fna $line >>dmr.fasta; done
# in conda environment
faSplit byname dmr.fasta dmr_fasta/


cd genes_fasta/
ls *fa >filenames_genes

cd dmr_fasta/
#ls *fa >filenames_dmr
for file in *fa ; do ls $file >>filenames_dmr ; done
```
33. For SMPs: Count cytosines
```
library("methimpute")

## Only CG context
## genes
files <- read.table(file="./filenames_genes")
files$V1 <- gsub(":", "", files$V1)
files <- as.character(files$V1)

cytosine_count <- ""
for (f in 1:length(files)){
  print(f)
  cytosines <- c()
  try(cytosines <- extractCytosinesFromFASTA(files[f], contexts = 'CG'),silent=T)
  if (length(cytosines) > 0){
    cytosine_count <- rbind(cytosine_count,(c(files[f],table(cytosines$context))))
  } else {
    cytosine_count <- rbind(cytosine_count,(c(files[f],0)))
  }
  #print(cytosine_count)
}
cytosine_count <- cytosine_count[-c(1),]
write.table(cytosine_count,file="cytosine_count_genes2.txt",row.names=F, col.names=F,quote=F,sep="\t")


## dmrs
setwd("C:/Users/Excellaptop/Desktop/stats/smp/dmr_fasta/")
files <- read.table(file="./filenames_dmr")
files$V1 <- gsub(":", "", files$V1)
files <- as.character(files$V1)

cytosine_count <- ""
for (f in 1:length(files)){
  print(f)
  cytosines <- c()
  try(cytosines <- extractCytosinesFromFASTA(files[f], contexts = 'CG'),silent=T)
  if (length(cytosines) > 0){
    cytosine_count <- rbind(cytosine_count,(c(files[f],table(cytosines$context))))
  } else {
    cytosine_count <- rbind(cytosine_count,(c(files[f],0)))
  }
  #print(cytosine_count)
}
cytosine_count <- cytosine_count[-c(1),]
write.table(cytosine_count,file="cytosine_count_dmr.txt",row.names=F, col.names=F,quote=F,sep="\t")
```
34. For SMPs: split vcfs
```
bgzip -f smp1_maf.vcf
tabix -f smp1_maf.vcf.gz

bgzip -f smp6_maf.vcf
tabix -f smp6_maf.vcf.gz


# dmr positions
cat dmr_pos.list | while read -r line ; do tabix smp1_maf.vcf.gz $line >./group1_smp/dmr_regions/$line.smp.vcf; done
cat dmr_pos.list | while read -r line ; do tabix smp6_maf.vcf.gz $line >./group6_smp/dmr_regions/$line.smp.vcf; done

# gene positions
cat gene_pos.list | while read -r line ; do tabix smp1_maf.vcf.gz $line >./group1_smp/gene_regions/$line.smp.vcf; done
cat gene_pos.list | while read -r line ; do tabix smp6_maf.vcf.gz $line >./group6_smp/gene_regions/$line.smp.vcf; done

for file in *vcf ; do wc -l $file >>vcflengths_smp_dmr ; done
for file in *vcf ; do wc -l $file >>vcflengths_smp_gene ; done
```
35. For SMPs: Get good intervals and remove other files
```
#############################################################

#####   sMPs

#############################################################

####### GROUP 1

###### DMR regions
library(dplyr)
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group1_smp/dmr_regions/")

vcflengths_smp_dmr <- read.table(file="vcflengths_smp_dmr")
vcflengths_smp_dmr$V2 <- gsub(".smp.vcf","",vcflengths_smp_dmr$V2)
colnames(vcflengths_smp_dmr) <- c("numvar","interval")
cytosine_count <- read.table(file="/Users/Excellaptop/Desktop/stats/smp/cytosine_count_dmr.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
cytosine_count$V1 <- gsub("", ":", cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

merged <- inner_join(cytosine_count,vcflengths_smp_dmr,by="interval")
merged$prop <- merged$numvar/merged$numc
merged1 <- merged[merged$numvar >= 3 & merged$prop > 0.05,]
#merged2 <- merged[merged$prop > 0.05,]

write.table(merged1,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)

# bad_intervals
merged4 <- merged[!(merged$numvar >= 3 & merged$prop > 0.05),]
write.table(paste("rm ",merged4$interval,".input.txt",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)



###### GENE regions
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group1_smp/gene_regions/")

vcflengths_smp_gene <- read.table(file="vcflengths_smp_gene")
vcflengths_smp_gene$V2 <- gsub(".smp.vcf","",vcflengths_smp_gene$V2)
colnames(vcflengths_smp_gene) <- c("numvar","interval")
cytosine_count <- read.table(file="/Users/Excellaptop/Desktop/stats/smp/cytosine_count_genes.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
cytosine_count$V1 <- gsub("", ":", cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

merged <- inner_join(cytosine_count,vcflengths_smp_gene,by="interval")
merged$prop <- merged$numvar/merged$numc
#merged1 <- merged[merged$numvar >= 3,]
merged2 <- merged[merged$numvar >= 3 & merged$prop > 0.05,]

#write.table(merged1,file="good_intervals1",sep="\t",quote=F,row.names = F, col.names = F)
write.table(merged2,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)

# bad_intervals
#merged3 <- merged[merged$numvar < 3,]
merged4 <- merged[!(merged$numvar >= 3 & merged$prop > 0.05),]
#write.table(paste("rm ",merged3$interval,".input.txt",sep=""),file="bad_intervals1.sh",sep="\t",quote=F,row.names = F, col.names = F)
write.table(paste("rm ",merged4$interval,".input.txt",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)


#############################################################

####### GROUP 6

###### DMR regions
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group6_smp/dmr_regions/")

vcflengths_smp_dmr <- read.table(file="vcflengths_smp_dmr")
vcflengths_smp_dmr$V2 <- gsub(".smp.vcf","",vcflengths_smp_dmr$V2)
colnames(vcflengths_smp_dmr) <- c("numvar","interval")
cytosine_count <- read.table(file="/Users/Excellaptop/Desktop/stats/smp/cytosine_count_dmr.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
cytosine_count$V1 <- gsub("", ":", cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

merged <- inner_join(cytosine_count,vcflengths_smp_dmr,by="interval")
merged$prop <- merged$numvar/merged$numc
merged1 <- merged[merged$numvar >= 3 & merged$prop > 0.05,]
#merged2 <- merged[merged$prop > 0.05,]

write.table(merged1,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)

# bad_intervals
merged4 <- merged[!(merged$numvar >= 3 & merged$prop > 0.05),]
write.table(paste("rm ",merged4$interval,".input.txt",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)




###### GENE regions
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group6_smp/gene_regions/")

vcflengths_smp_gene <- read.table(file="vcflengths_smp_gene")
vcflengths_smp_gene$V2 <- gsub(".smp.vcf","",vcflengths_smp_gene$V2)
colnames(vcflengths_smp_gene) <- c("numvar","interval")
cytosine_count <- read.table(file="/Users/Excellaptop/Desktop/stats/smp/cytosine_count_genes.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
cytosine_count$V1 <- gsub("", ":", cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

merged <- inner_join(cytosine_count,vcflengths_smp_gene,by="interval")
merged$prop <- merged$numvar/merged$numc
merged1 <- merged[merged$numvar >= 3 & merged$prop > 0.05,]
#merged2 <- merged[merged$prop > 0.05,]

write.table(merged1,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)

# bad_intervals
merged4 <- merged[!(merged$numvar >= 3 & merged$prop > 0.05),]
write.table(paste("rm ",merged4$interval,".input.txt",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)
```
36. For SMPs: Get alpha
>*repeat for both groups and gene and dmr regions*
```
## group1_smp/gene_regions
cd /mnt/c/Users/Excellaptop/Desktop/stats/smp/group1_smp/gene_regions
for file in *.smp.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' -  >${file/.smp.vcf/.input.txt} ; done 
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
mkdir input
mv *input.txt input/
dos2unix bad_intervals.sh
mv bad_intervals.sh input/
cd input/
./bad_intervals.sh
cd ../
perl alpha_estimation.pl -dir ./input -output  alpha_Dm_brachy -length_list length_list
```

37. For SMPs: impute NA because Dm test script interpretes NA as 0
```
#impute.R
impute <- function(inputfile){
  count_na <- apply(inputfile[, 3:ncol(inputfile)], 1, function(x) sum(x == "."))
  count_1 <- apply(inputfile[, 3:ncol(inputfile)], 1, function(x) sum(x == "1"))
  prob <- count_1/(ncol(inputfile)-2)
  
  for(i in 1:nrow(inputfile)){
    miss <- inputfile[i, 3:ncol(inputfile)] == "."
    replacement <- rbinom(count_na[i], 1, prob[i]) # 1 because one row/trial/sample
    inputfile[i, 3:ncol(inputfile)][miss] <- replacement
  }
  return(inputfile)
}

# loop through all input files and write them to a new file
# copy input/ and rename to impute/

# group 1 dmr
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group1_smp/dmr_regions/impute")
list <- read.table("length_list", header = FALSE)
list$V1 <- gsub(":", "", list$V1)

for(f in 1:nrow(list)){
  file <- read.table(list[f,1], header = FALSE)
  file_imputed <- impute(inputfile = file)
  filename <- gsub("input.txt", "input.imputed.txt", list[f,1])
  write.table(file_imputed, file = filename, sep = "\t", quote = F , col.names = F, row.names = F)
}

# group 1 gene
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group1_smp/gene_regions/impute")
list <- read.table("length_list", header = FALSE)
list$V1 <- gsub(":", "", list$V1)

for(f in 1:nrow(list)){
  file <- read.table(list[f,1], header = FALSE)
  file_imputed <- impute(inputfile = file)
  filename <- gsub("input.txt", "input.imputed.txt", list[f,1])
  write.table(file_imputed, file = filename, sep = "\t", quote = F , col.names = F, row.names = F)
}

# group 6 dmr
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group6_smp/dmr_regions/impute")
list <- read.table("length_list", header = FALSE)
list$V1 <- gsub(":", "", list$V1)

for(f in 1:nrow(list)){
  file <- read.table(list[f,1], header = FALSE)
  file_imputed <- impute(inputfile = file)
  filename <- gsub("input.txt", "input.imputed.txt", list[f,1])
  write.table(file_imputed, file = filename, sep = "\t", quote = F , col.names = F, row.names = F)
}

# group 6 gene
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group6_smp/gene_regions/impute")
list <- read.table("length_list", header = FALSE)
list$V1 <- gsub(":", "", list$V1)

for(f in 1:nrow(list)){
  file <- read.table(list[f,1], header = FALSE)
  file_imputed <- impute(inputfile = file)
  filename <- gsub("input.txt", "input.imputed.txt", list[f,1])
  write.table(file_imputed, file = filename, sep = "\t", quote = F , col.names = F, row.names = F)
}
```
38. For SMPs: get Dm estimates for all groups and regions
```
# in R (good_intervals_Dm.R)
library(dplyr)

# run for each group and region

# group 1 dmr
setwd("C:/Users/Excellaptop/Desktop/stats/smp/group1_smp/dmr_regions")
good_intervals <- read.table(file="good_intervals")
alpha_dm <- read.table(file="alpha_Dm_brachy_impute")
alpha_dm$V1 <- gsub(".input.imputed.txt","",alpha_dm$V1)
alpha_dm <- alpha_dm[1:2]
good_intervals <- inner_join(good_intervals,alpha_dm,by="V1")

write.table(good_intervals,file="good_intervals_alpha",sep="\t",quote=F,row.names = F, col.names = F)

#######
# Dm estimates

cd /mnt/c/Users/Excellaptop/Desktop/stats/smp/group1_smp/dmr_regions
dos2unix good_intervals_alpha
cp good_intervals_alpha impute/
cd impute/
cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl ../Dm_test_new.pl -input $value1.input.imputed.txt -output $value1.Dm_group1_dmr.txt -length $value2 -alpha $value5  ; done
```

39. Calculate admixture plots for snps, smps and dmrs
```
### admixture
cd ./admixture/
for K in 2 3 4 5 6 7 8 9 10 ; do /mnt/c/Users/Excellaptop/ubuntu/softwares/admixture_linux-1.3.0/admixture32 --cv snp29_maf.bed $K | tee log${K}.out; done
grep -h CV log*.out
grep -h CV log*.out|awk -F '[=):]' '{print $2"\t"$4}'|sed '1i\\K\tCV_error' > CV_for_plot_snp.txt
```

40. filtering 10% maf
```
### subgroups
# 10% filtering
plink --vcf snp1_na.recode.vcf --maf 0.1 --make-bed --out ./snp1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp6_na.recode.vcf --maf 0.1 --make-bed --out ./snp6_maf10 --set-missing-var-ids @:# --keep-allele-order

plink --vcf smp1_na.recode.vcf --maf 0.1 --make-bed --out ./smp1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_na.recode.vcf --maf 0.1 --make-bed --out ./smp6_maf10 --set-missing-var-ids @:# --keep-allele-order

plink --vcf dmr1_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_na.recode.vcf --maf 0.1 --make-bed --out ./dmr6_maf10 --set-missing-var-ids @:# --keep-allele-order



plink --bfile ./snp1_maf10 --recode vcf-iid -out ./snp1_maf10 --keep-allele-order
plink --bfile ./snp6_maf10 --recode vcf-iid -out ./snp6_maf10 --keep-allele-order

plink --bfile ./smp1_maf10 --recode vcf-iid -out ./smp1_maf10 --keep-allele-order
plink --bfile ./smp6_maf10 --recode vcf-iid -out ./smp6_maf10 --keep-allele-order

plink --bfile ./dmr1_maf10 --recode vcf-iid -out ./dmr1_maf10 --keep-allele-order
plink --bfile ./dmr6_maf10 --recode vcf-iid -out ./dmr6_maf10 --keep-allele-order

grep -v '^#' ./plink/snp29_maf.vcf | wc -l
```
41. get subgroup vcfs
```
plink --vcf snp1_na.recode.vcf --make-bed --out ./snp1_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1_na.recode.vcf --make-bed --out ./smp1_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1_na.recode.vcf --make-bed --out ./dmr1_maf100 --set-missing-var-ids @:# --keep-allele-order

plink --bfile ./snp1_maf100 --recode vcf-iid -out ./snp1_maf100 --keep-allele-order
plink --bfile ./smp1_maf100 --recode vcf-iid -out ./smp1_maf100 --keep-allele-order
plink --bfile ./dmr1_maf100 --recode vcf-iid -out ./dmr1_maf100 --keep-allele-order


bcftools view -S pop1.1.txt snp1_maf100.vcf > snp1.1_na.recode.vcf
bcftools view -S pop1.2.txt snp1_maf100.vcf > snp1.2_na.recode.vcf

bcftools view -S pop1.1.txt smp1_maf100.vcf > smp1.1_na.recode.vcf
bcftools view -S pop1.2.txt smp1_maf100.vcf > smp1.2_na.recode.vcf

bcftools view -S pop1.1.txt dmr1_maf100.vcf > dmr1.1_na.recode.vcf
bcftools view -S pop1.2.txt dmr1_maf100.vcf > dmr1.2_na.recode.vcf

plink --vcf snp1.1_na.recode.vcf --maf 0.1 --make-bed --out ./snp1.1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp1.2_na.recode.vcf --maf 0.1 --make-bed --out ./snp1.2_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1.1_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1.2_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.2_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1.1_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.1_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1.2_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.2_maf10 --set-missing-var-ids @:# --keep-allele-order

plink --bfile ./snp1.1_maf10 --recode vcf-iid -out ./snp1.1_maf10 --keep-allele-order
plink --bfile ./snp1.2_maf10 --recode vcf-iid -out ./snp1.2_maf10 --keep-allele-order
plink --bfile ./smp1.1_maf10 --recode vcf-iid -out ./smp1.1_maf10 --keep-allele-order
plink --bfile ./smp1.2_maf10 --recode vcf-iid -out ./smp1.2_maf10 --keep-allele-order
plink --bfile ./dmr1.1_maf10 --recode vcf-iid -out ./dmr1.1_maf10 --keep-allele-order
plink --bfile ./dmr1.2_maf10 --recode vcf-iid -out ./dmr1.2_maf10 --keep-allele-order

```

42. sitepi and D with 10% maf
```
# calculate site-pi and D with 10% filtering

vcftools --vcf ./snp1_maf10.vcf --site-pi --out ./pi/snp1_maf10
vcftools --vcf ./snp1.1_maf10.vcf --site-pi --out ./pi/snp1.1_maf10
vcftools --vcf ./snp1.2_maf10.vcf --site-pi --out ./pi/snp1.2_maf10
vcftools --vcf ./snp6_maf10.vcf --site-pi --out ./pi/snp6_maf10

vcftools --vcf ./smp1_maf10.vcf --site-pi --out ./pi/smp1_maf10
vcftools --vcf ./smp1.1_maf10.vcf --site-pi --out ./pi/smp1.1_maf10
vcftools --vcf ./smp1.2_maf10.vcf --site-pi --out ./pi/smp1.2_maf10
vcftools --vcf ./smp6_maf10.vcf --site-pi --out ./pi/smp6_maf10

vcftools --vcf ./dmr1_maf10.vcf --site-pi --out ./pi/dmr1_maf10
vcftools --vcf ./dmr1.1_maf10.vcf --site-pi --out ./pi/dmr1.1_maf10
vcftools --vcf ./dmr1.2_maf10.vcf --site-pi --out ./pi/dmr1.2_maf10
vcftools --vcf ./dmr6_maf10.vcf --site-pi --out ./pi/dmr6_maf10
```
43. Fst
```
vcf-query -l smp1_maf.vcf | grep 'BdTR1' > pop1.1.txt
vcf-query -l smp1_maf.vcf | grep 'BdTR2' > pop1.1.txt

vcftools --vcf smp1_maf.vcf --weir-fst-pop pop1.1.txt --weir-fst-pop pop1.2.txt --out smp1.1_vs_smp1.2
vcftools --vcf dmr1_maf.vcf --weir-fst-pop pop1.1.txt --weir-fst-pop pop1.2.txt --out dmr1.1_vs_dmr1.2
vcftools --vcf snp1_maf.vcf --weir-fst-pop pop1.1.txt --weir-fst-pop pop1.2.txt --out snp1.1_vs_snp1.2
```
44. randomly downsample by 50% of variants
```
(cat smp29_maf.vcf | head -n 10000 | grep ^# ; grep -v ^# smp29_maf.vcf| shuf -n 6947 | LC_ALL=C sort -k1,1V -k2,2n) > smp29_small.vcf
(cat dmr29_maf.vcf | head -n 10000 | grep ^# ; grep -v ^# dmr29_maf.vcf| shuf -n 1161 | LC_ALL=C sort -k1,1V -k2,2n) > dmr29_small.vcf

/mnt/c/Users/Excellaptop/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./smp29_small.vcf  -OutPut ./smp29_small_dis.mat
/mnt/c/Users/Excellaptop/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./dmr29_small.vcf  -OutPut ./dmr29_small_dis.mat
```
45. Repeat steps to generate new SMP vcf
>*Repeat step 17 with lower threshold (>4)*
```
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")


cov <- read.table(file=as.character(covfile$V1[1]), header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
print(cov$ID[1])

context <- read.table(file=as.character(contextfiles$V1[1]),header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "-")
print(context$ID[1])

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub("_N","N", gsub(gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile$V1)[1],"",cov_context$chromosome))
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.total)

print(noncoversionrate)

rm(cov)
rm(context)

#cov_context <- cov_context[cov_context$count.total > 4,]

b <- apply(cov_context[c(5,11)], 1, binom.test, p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
dim(cov_context)

rm(b)

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]
dim(cov_context)

samplename <- gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile[1,])
colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

rm(cov_context)

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  cov <- read.table(file=as.character(covfile[i,]), header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
  print(cov$ID[1])

  context <- read.table(file=as.character(contextfiles[i,]),header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "-")
  print(context$ID[1])
  print(i)

  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub("_N","N", gsub(gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile$V1)[i],"",cov_context$chromosome))
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
  cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_011032.1",]$count.total)
  
  print(noncoversionrate)

  rm(cov)
  rm(context)

  #cov_context <- cov_context[cov_context$count.total > 4,]
  
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))
  dim(cov_context)

  rm(b)

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  dim(cov_context)

  samplename <- gsub(".wgbs_bismark_hisat2.deduplicated.bismark.cov","",covfile[i,])
  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
  dim(cov_context2)
  rm(cov_context)
}

#### write to file for further processing

#cov_context2$chromosome <- gsub("-.*","",cov_context2$ID)
#cov_context2$position   <- as.numeric(gsub(".*-","",cov_context2$ID))
write.table(cov_context2,file="cov_context3_all_thresh4.txt",row.names = F)

```
>*steps 18, 21*
```
library(dplyr)
library(ggplot2)
library(data.table)

cov_context3 <- fread(file="cov_context3_all_thresh4.txt",header=T)
na_count <- apply(cov_context3[,6:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[,6:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[,6:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[,1:3]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"
cov_context3 <- cov_context3[,-c(1:5)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
fwrite(cbind(meta,cov_context3),file="brachy_meth.thresh4.vcf",quote = F, row.names = F,sep="\t")

meta2 <- cov_context4[,1:3]
colnames(meta2) <- c("#CHROM","POS","ID")
meta2$REF <- "A"
meta2$ALT <- "T"
meta2$QUAL <- 4000
meta2$FILTER <- "PASS"
meta2$INFO <- "DP=1000"
meta2$FORMAT <- "GT"
cov_context4 <- cov_context4[,-c(1:5)]
cov_context4[cov_context4 == "U"] <- "0/0"
cov_context4[cov_context4 == "M"] <- "1/1"
cov_context4[is.na(cov_context4)] <- "./."
fwrite(cbind(meta2,cov_context4),file="brachy_meth_var_invar.thresh4.vcf",quote = F, row.names = F,sep="\t")

################
df1 <- fread("group1_names", header = F, drop = "V1")
df6 <- fread("group6_names", header = F, drop = "V1")

smp <- fread("brachy_meth.thresh4.vcf", header = T)
smp_29 <- smp[,-11]
fwrite(smp_29, file = "smp_29_thresh4.vcf", quote = F, row.names = F,sep ="\t", eol = "\n")

colnames(smp_29) %in% df1$V2

smp_grp1 <- smp_29[,-c(10:16)]
colnames(smp_grp1) %in% df1$V2
fwrite(smp_grp1,file="smp_grp1_thresh4.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

smp_grp6 <- smp_29[,-c(17:38)]
colnames(smp_grp6) %in% df6$V2
fwrite(smp_grp6,file="smp_grp6_thresh4.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")
```
>*cleanup vcfs*
```
library(data.table)
library(dplyr)

# read in vcf files

smp29 <- fread("smp_29_thresh4.vcf")
smp1 <- fread("smp_grp1_thresh4.vcf")
smp6 <- fread("smp_grp6_thresh4.vcf")


# remove plastid
smp29 <- dplyr::filter(smp29, `#CHROM` != "NC_011032.1")
smp1 <- dplyr::filter(smp1, `#CHROM` != "NC_011032.1")
smp6 <- dplyr::filter(smp6, `#CHROM` != "NC_011032.1")

# change chromosome names to 1, 2, etc. in #CHROM and ID columns


## smp
smp29$`#CHROM` <- gsub("NC_016131.3", "1", smp29$`#CHROM`)
smp29$`#CHROM` <- gsub("NC_016132.3", "2", smp29$`#CHROM`)
smp29$`#CHROM` <- gsub("NC_016133.3", "3", smp29$`#CHROM`)
smp29$`#CHROM` <- gsub("NC_016134.3", "4", smp29$`#CHROM`)
smp29$`#CHROM` <- gsub("NC_016135.3", "5", smp29$`#CHROM`)
smp29$ID <- gsub("NC_016131.3-", "1:", smp29$ID)
smp29$ID <- gsub("NC_016132.3-", "2:", smp29$ID)
smp29$ID <- gsub("NC_016133.3-", "3:", smp29$ID)
smp29$ID <- gsub("NC_016134.3-", "4:", smp29$ID)
smp29$ID <- gsub("NC_016135.3-", "5:", smp29$ID)

smp1$`#CHROM` <- gsub("NC_016131.3", "1", smp1$`#CHROM`)
smp1$`#CHROM` <- gsub("NC_016132.3", "2", smp1$`#CHROM`)
smp1$`#CHROM` <- gsub("NC_016133.3", "3", smp1$`#CHROM`)
smp1$`#CHROM` <- gsub("NC_016134.3", "4", smp1$`#CHROM`)
smp1$`#CHROM` <- gsub("NC_016135.3", "5", smp1$`#CHROM`)
smp1$ID <- gsub("NC_016131.3-", "1:", smp1$ID)
smp1$ID <- gsub("NC_016132.3-", "2:", smp1$ID)
smp1$ID <- gsub("NC_016133.3-", "3:", smp1$ID)
smp1$ID <- gsub("NC_016134.3-", "4:", smp1$ID)
smp1$ID <- gsub("NC_016135.3-", "5:", smp1$ID)


smp6$`#CHROM` <- gsub("NC_016131.3", "1", smp6$`#CHROM`)
smp6$`#CHROM` <- gsub("NC_016132.3", "2", smp6$`#CHROM`)
smp6$`#CHROM` <- gsub("NC_016133.3", "3", smp6$`#CHROM`)
smp6$`#CHROM` <- gsub("NC_016134.3", "4", smp6$`#CHROM`)
smp6$`#CHROM` <- gsub("NC_016135.3", "5", smp6$`#CHROM`)
smp6$ID <- gsub("NC_016131.3-", "1:", smp6$ID)
smp6$ID <- gsub("NC_016132.3-", "2:", smp6$ID)
smp6$ID <- gsub("NC_016133.3-", "3:", smp6$ID)
smp6$ID <- gsub("NC_016134.3-", "4:", smp6$ID)
smp6$ID <- gsub("NC_016135.3-", "5:", smp6$ID)



# write to file again
fwrite(smp29, file = "smp29_thresh4_unfiltered.vcf", quote = F, sep = "\t", eol = "\n")
fwrite(smp1, file = "smp1_thresh4_unfiltered.vcf", quote = F, sep = "\t", eol = "\n")
fwrite(smp6, file = "smp6_thresh4_unfiltered.vcf", quote = F, sep = "\t", eol = "\n")

test <-fread("smp6_thresh4_unfiltered.vcf")

## add header to vcf files
```

46. get SNP invariants
```
bcftools view -i 'FMT/DP>=5 & QUAL>=30'
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.zauchner/brachy/GCF_000005505.3/GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic.fna -all-sites -V brachypodium.all.g.vcf.gz -O brachypodium.var_invar.vcf.gz
vcftools --gzvcf brachypodium.var_invar.vcf.gz --recode --maf 0 --out brachypodium.var_invar_maf0
```
>*make vcf*
```
snp_invar <- fread("brachypodium.var_invar.filtered.maf0.recode.vcf", header = T)
snp_invar$FILTER <- gsub(".", "PASS", snp_invar$FILTER)

snp_invar_29 <- snp_invar[,-10]
colnames(snp_invar_29)[10:length(colnames(snp_invar_29))] <- gsub("_marked.bam", "", colnames(snp_invar_29)[10:length(colnames(snp_invar_29))])
fwrite(snp_invar_29, file = "snp_invar_29.vcf", quote = F, row.names = F,sep ="\t", eol = "\n")

colnames(snp_invar_29) %in% df1$V2

snp_invar_grp1 <- snp_invar_29[,-c(10:16)]
colnames(snp_invar_grp1) %in% df1$V2
fwrite(snp_invar_grp1,file="snp_invar_1.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")

snp_invar_grp6 <- snp_invar_29[,-c(17:38)]
colnames(snp_invar_grp6) %in% df6$V2
fwrite(snp_invar_grp6,file="snp_invar_6.vcf",quote = F, row.names = F,sep ="\t", eol = "\n")
```
>*remove heterozygoes*
```
library(data.table)
library(dplyr)

snp29.unfiltered <- fread(input = "snp29_invar_unfiltered.vcf")
snp1.unfiltered <- fread(input = "snp1_invar_unfiltered.vcf")
snp6.unfiltered <- fread(input = "snp6_invar_unfiltered.vcf")


snp6 <- snp6.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp6 <- snp6 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

snp1 <- snp1.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp1 <- snp1 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

snp29 <- snp29.unfiltered %>%
  mutate_all(~gsub("0/1", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1/0", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("0\\|1", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1\\|0", "./.", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("0\\|0", "0/0", .))
snp29 <- snp29 %>%
  mutate_all(~gsub("1\\|1", "1/1", .))

fwrite(snp1, file = "snp1_invar_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)
fwrite(snp6, file = "snp6_invar_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)
fwrite(snp29, file = "snp29_invar_unfiltered_nohet.vcf", quote = F, eol = "\n", sep = "\t", row.names = F)
```
47. get .bed files for noncoding regions, gbM genes and exons
>**
```
grep 'exon' genomic.gff | cut -f 1,4,5 > exon_pos.bed
sed -i '/KZ62297/d' exon_pos.bed
sort -k 1,2 -V exon_pos.bed >exon_pos_sorted.bed

# change chromosome names to match vcfs -> exon_pos2.bed (1, 2, 3, etc.) and exon_pos3.bed (NC_016131, NC_016132, etc.)

grep 'exon' GCA_000005505.4_Brachypodium_distachyon_v3.0_genomic.gff | cut -f 1,4,5 >gene_pos.bed
sed -i '/KZ62297/d' gene_pos.bed
sort -k 1,2 -V gene_pos.bed >gene_pos2.bed
/data/proj2/popgen/a.ramesh/software/bedtools2/bin/bedtools complement -i gene_pos2.bed -g chr_lengths >non_gene.bed
sed -e 's/\t/:/' -e  's/\t/-/' non_gene.bed >non_gene.list

# change chromosome names to match vcfs -> non_gene2.bed (1, 2, 3, etc.) and non_gene3.bed (NC_016131, NC_016132, etc.)

grep -F -f gBM_genes GCA_000005505.4_Brachypodium_distachyon_v3.0_genomic.gff  | grep 'exon' | cut -f 1,4,5 > gbm_genes.bed
sed -i '/KZ62297/d' gbm_genes.bed
sed -e 's/\t/:/' -e  's/\t/-/' gbm_genes.bed >gbm_genes.list

# change chromosome names to match vcfs -> gbm_genes2.bed (1, 2, 3, etc.) and gbm_genes3.bed (NC_016131, NC_016132, etc.)
```

48. filtering of VCFs
```
# filtering of vcfs for gBM SMP, gBM DMR and coding and noncoding SNPs, and invariants


#### SNPs ######################################################################################################################################

# filter max 20% NA, only non coding SNP
vcftools --vcf snp29_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene2.bed --out ./snp/snp29_noncoding_na --recode-INFO-all
vcftools --vcf snp1_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene2.bed --out ./snp/snp1_noncoding_na --recode-INFO-all
vcftools --vcf snp6_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene2.bed --out ./snp/snp6_noncoding_na --recode-INFO-all 

grep -v '^#' snp29_noncoding_na.recode.vcf | wc -l


plink --vcf snp29_noncoding_na.recode.vcf --maf 0.04 --make-bed --out ./snp29_noncoding_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp1_noncoding_na.recode.vcf --maf 0.04 --make-bed --out ./snp1_noncoding_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp6_noncoding_na.recode.vcf --maf 0.04 --make-bed --out ./snp6_noncoding_maf04 --set-missing-var-ids @:# --keep-allele-order

plink --bfile ./snp29_noncoding_maf04 --recode vcf-iid -out ./snp29_noncoding_maf04 --keep-allele-order
plink --bfile ./snp1_noncoding_maf04 --recode vcf-iid -out ./snp1_noncoding_maf04 --keep-allele-order
plink --bfile ./snp6_noncoding_maf04 --recode vcf-iid -out ./snp6_noncoding_maf04 --keep-allele-order



### coding and non coding invariants

vcftools --vcf snp29_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene3.bed --out ./snp/snp29_invar_noncoding_na --recode-INFO-all
vcftools --vcf snp1_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene3.bed --out ./snp/snp1_invar_noncoding_na --recode-INFO-all
vcftools --vcf snp6_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed non_gene3.bed --out ./snp/snp6_invar_noncoding_na --recode-INFO-all 


vcftools --vcf snp29_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./snp/snp29_invar_coding_na --recode-INFO-all
vcftools --vcf snp1_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./snp/snp1_invar_coding_na --recode-INFO-all
vcftools --vcf snp6_invar_unfiltered_nohet.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./snp/snp6_invar_coding_na --recode-INFO-all 




#### SMPs ####################################################################################################################################
# filter max 20% NA, gBM genes and all exons

vcftools --vcf smp29_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./smp/smp29_gbm_na --recode-INFO-all
vcftools --vcf smp1_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./smp/smp1_gbm_na --recode-INFO-all
vcftools --vcf smp6_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed gbm_genes3.bed --out ./smp/smp6_gbm_na --recode-INFO-all

vcftools --vcf smp29_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed exon_pos3.bed --out ./smp/smp29_exon_na --recode-INFO-all
vcftools --vcf smp1_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed exon_pos3.bed --out ./smp/smp1_exon_na --recode-INFO-all
vcftools --vcf smp6_thresh4_unfiltered.vcf --max-missing 0.8 --recode --bed exon_pos3.bed --out ./smp/smp6_exon_na --recode-INFO-all



# change chromosome numbers
plink --vcf smp29_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./smp29_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./smp1_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./smp6_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp29_gbm_maf04 --recode vcf-iid -out ./smp29_gbm_maf04 --keep-allele-order
plink --bfile ./smp1_gbm_maf04 --recode vcf-iid -out ./smp1_gbm_maf04 --keep-allele-order
plink --bfile ./smp6_gbm_maf04 --recode vcf-iid -out ./smp6_gbm_maf04 --keep-allele-order


plink --vcf smp29_exon_na.recode.vcf --maf 0.04 --make-bed --out ./smp29_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1_exon_na.recode.vcf --maf 0.04 --make-bed --out ./smp1_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_exon_na.recode.vcf --maf 0.04 --make-bed --out ./smp6_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp29_exon_maf04 --recode vcf-iid -out ./smp29_exon_maf04 --keep-allele-order
plink --bfile ./smp1_exon_maf04 --recode vcf-iid -out ./smp1_exon_maf04 --keep-allele-order
plink --bfile ./smp6_exon_maf04 --recode vcf-iid -out ./smp6_exon_maf04 --keep-allele-order


## get non-gbM exons
from exon filtered vcfs

vcftools --vcf smp29_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/smp29_exon_nongbm --recode-INFO-all
vcftools --vcf smp1_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/smp1_exon_nongbm --recode-INFO-all
vcftools --vcf smp6_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/smp6_exon_nongbm --recode-INFO-all


plink --vcf smp29_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./smp29_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./smp1_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./smp6_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp29_exon_nongbm_maf04 --recode vcf-iid -out ./smp29_exon_nongbm_maf04 --keep-allele-order
plink --bfile ./smp1_exon_nongbm_maf04 --recode vcf-iid -out ./smp1_exon_nongbm_maf04 --keep-allele-order
plink --bfile ./smp6_exon_nongbm_maf04 --recode vcf-iid -out ./smp6_exon_nongbm_maf04 --keep-allele-order




#### DMRs ####################################################################################################################################

vcftools --vcf bd_dmrs_dip29.vcf --max-missing 0.8 --recode --bed gbm_genes2.bed --out ./dmr/dmr29_gbm_na --recode-INFO-all
vcftools --vcf bd_dmrs_dip1.vcf --max-missing 0.8 --recode --bed gbm_genes2.bed --out ./dmr/dmr1_gbm_na --recode-INFO-all
vcftools --vcf bd_dmrs_dip6.vcf --max-missing 0.8 --recode --bed gbm_genes2.bed --out ./dmr/dmr6_gbm_na --recode-INFO-all


plink --vcf dmr29_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./dmr29_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./dmr1_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_gbm_na.recode.vcf --maf 0.04 --make-bed --out ./dmr6_gbm_maf04 --set-missing-var-ids @:# --keep-allele-order

plink --bfile ./dmr29_gbm_maf04 --recode vcf-iid -out ./dmr29_gbm_maf04 --keep-allele-order
plink --bfile ./dmr1_gbm_maf04 --recode vcf-iid -out ./dmr1_gbm_maf04 --keep-allele-order
plink --bfile ./dmr6_gbm_maf04 --recode vcf-iid -out ./dmr6_gbm_maf04 --keep-allele-order

## get non-gbM dmrs
vcftools --vcf bd_dmrs_dip29.vcf --max-missing 0.8 --recode --bed exon_pos2.bed --out ./dmr/dmr29_exon_na --recode-INFO-all
vcftools --vcf bd_dmrs_dip1.vcf --max-missing 0.8 --recode --bed exon_pos2.bed --out ./dmr/dmr1_exon_na --recode-INFO-all
vcftools --vcf bd_dmrs_dip6.vcf --max-missing 0.8 --recode --bed exon_pos2.bed --out ./dmr/dmr6_exon_na --recode-INFO-all

#from exon filtered dmr vcfs
vcftools --vcf dmr29_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/dmr29_exon_nongbm --recode-INFO-all
vcftools --vcf dmr1_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/dmr1_exon_nongbm --recode-INFO-all
vcftools --vcf dmr6_exon_na.recode.vcf --recode --exclude-bed gbm_genes2.bed --out ./non/dmr6_exon_nongbm --recode-INFO-all


plink --vcf dmr29_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./dmr29_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./dmr1_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_exon_nongbm.recode.vcf --maf 0.04 --make-bed --out ./dmr6_exon_nongbm_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr29_exon_nongbm_maf04 --recode vcf-iid -out ./dmr29_exon_nongbm_maf04 --keep-allele-order
plink --bfile ./dmr1_exon_nongbm_maf04 --recode vcf-iid -out ./dmr1_exon_nongbm_maf04 --keep-allele-order
plink --bfile ./dmr6_exon_nongbm_maf04 --recode vcf-iid -out ./dmr6_exon_nongbm_maf04 --keep-allele-order

# DMR exons
plink --vcf dmr29_exon_na.recode.vcf --maf 0.04 --make-bed --out ./dmr29_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1_exon_na.recode.vcf --maf 0.04 --make-bed --out ./dmr1_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_exon_na.recode.vcf --maf 0.04 --make-bed --out ./dmr6_exon_maf04 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr29_exon_maf04 --recode vcf-iid -out ./dmr29_exon_maf04 --keep-allele-order
plink --bfile ./dmr1_exon_maf04 --recode vcf-iid -out ./dmr1_exon_maf04 --keep-allele-order
plink --bfile ./dmr6_exon_maf04 --recode vcf-iid -out ./dmr6_exon_maf04 --keep-allele-order



#### subgroups ####################################################################################################################################

# 10% filtering

plink --vcf snp1_noncoding_na.recode.vcf --maf 0.1 --make-bed --out ./snp1_noncoding_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp6_noncoding_na.recode.vcf --maf 0.1 --make-bed --out ./snp6_noncoding_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./snp1_noncoding_maf10 --recode vcf-iid -out ./snp1_noncoding_maf10 --keep-allele-order
plink --bfile ./snp6_noncoding_maf10 --recode vcf-iid -out ./snp6_noncoding_maf10 --keep-allele-order


plink --vcf smp1_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./smp1_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./smp6_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_gbm_maf10 --recode vcf-iid -out ./smp1_gbm_maf10 --keep-allele-order
plink --bfile ./smp6_gbm_maf10 --recode vcf-iid -out ./smp6_gbm_maf10 --keep-allele-order


plink --vcf smp1_exon_na.recode.vcf --maf 0.1 --make-bed --out ./smp1_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_exon_na.recode.vcf --maf 0.1 --make-bed --out ./smp6_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_exon_maf10 --recode vcf-iid -out ./smp1_exon_maf10 --keep-allele-order
plink --bfile ./smp6_exon_maf10 --recode vcf-iid -out ./smp6_exon_maf10 --keep-allele-order



plink --vcf dmr1_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./dmr6_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_gbm_maf10 --recode vcf-iid -out ./dmr1_gbm_maf10 --keep-allele-order
plink --bfile ./dmr6_gbm_maf10 --recode vcf-iid -out ./dmr6_gbm_maf10 --keep-allele-order


plink --vcf smp1_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./smp1_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp6_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./smp6_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_exon_nongbm_maf10 --recode vcf-iid -out ./smp1_exon_nongbm_maf10 --keep-allele-order
plink --bfile ./smp6_exon_nongbm_maf10 --recode vcf-iid -out ./smp6_exon_nongbm_maf10 --keep-allele-order


plink --vcf dmr1_exon_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_exon_na.recode.vcf --maf 0.1 --make-bed --out ./dmr6_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_exon_maf10 --recode vcf-iid -out ./dmr1_exon_maf10 --keep-allele-order
plink --bfile ./dmr6_exon_maf10 --recode vcf-iid -out ./dmr6_exon_maf10 --keep-allele-order


plink --vcf dmr1_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./dmr1_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr6_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./dmr6_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_exon_nongbm_maf10 --recode vcf-iid -out ./dmr1_exon_nongbm_maf10 --keep-allele-order
plink --bfile ./dmr6_exon_nongbm_maf10 --recode vcf-iid -out ./dmr6_exon_nongbm_maf10 --keep-allele-order


### make subgroups
plink --vcf snp1_noncoding_na.recode.vcf --make-bed --out ./snp1_noncoding_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./snp1_noncoding_maf100 --recode vcf-iid -out ./snp1_noncoding_maf100 --keep-allele-order
bcftools view -S pop1.1.txt snp1_noncoding_maf100.vcf > snp1.1_noncoding_na.recode.vcf
bcftools view -S pop1.2.txt snp1_noncoding_maf100.vcf > snp1.2_noncoding_na.recode.vcf
plink --vcf snp1.1_noncoding_na.recode.vcf --maf 0.1 --make-bed --out ./snp1.1_noncoding_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf snp1.2_noncoding_na.recode.vcf --maf 0.1 --make-bed --out ./snp1.2_noncoding_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./snp1.1_noncoding_maf10 --recode vcf-iid -out ./snp1.1_noncoding_maf10 --keep-allele-order
plink --bfile ./snp1.2_noncoding_maf10 --recode vcf-iid -out ./snp1.2_noncoding_maf10 --keep-allele-order


plink --vcf smp1_gbm_na.recode.vcf --make-bed --out ./smp1_gbm_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_gbm_maf100 --recode vcf-iid -out ./smp1_gbm_maf100 --keep-allele-order
bcftools view -S pop1.1.txt smp1_gbm_maf100.vcf > smp1.1_gbm_na.recode.vcf
bcftools view -S pop1.2.txt smp1_gbm_maf100.vcf > smp1.2_gbm_na.recode.vcf
plink --vcf smp1.1_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.1_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1.2_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.2_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1.1_gbm_maf10 --recode vcf-iid -out ./smp1.1_gbm_maf10 --keep-allele-order
plink --bfile ./smp1.2_gbm_maf10 --recode vcf-iid -out ./smp1.2_gbm_maf10 --keep-allele-order


plink --vcf dmr1_gbm_na.recode.vcf --make-bed --out ./dmr1_gbm_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_gbm_maf100 --recode vcf-iid -out ./dmr1_gbm_maf100 --keep-allele-order
bcftools view -S pop1.1.txt dmr1_gbm_maf100.vcf > dmr1.1_gbm_na.recode.vcf
bcftools view -S pop1.2.txt dmr1_gbm_maf100.vcf > dmr1.2_gbm_na.recode.vcf
plink --vcf dmr1.1_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.1_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1.2_gbm_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.2_gbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1.1_gbm_maf10 --recode vcf-iid -out ./dmr1.1_gbm_maf10 --keep-allele-order
plink --bfile ./dmr1.2_gbm_maf10 --recode vcf-iid -out ./dmr1.2_gbm_maf10 --keep-allele-order


# exons
plink --vcf smp1_exon_na.recode.vcf --make-bed --out ./smp1_exon_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_exon_maf100 --recode vcf-iid -out ./smp1_exon_maf100 --keep-allele-order
bcftools view -S pop1.1.txt smp1_exon_maf100.vcf > smp1.1_exon_na.recode.vcf
bcftools view -S pop1.2.txt smp1_exon_maf100.vcf > smp1.2_exon_na.recode.vcf
plink --vcf smp1.1_exon_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.1_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1.2_exon_na.recode.vcf --maf 0.1 --make-bed --out ./smp1.2_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1.1_exon_maf10 --recode vcf-iid -out ./smp1.1_exon_maf10 --keep-allele-order
plink --bfile ./smp1.2_exon_maf10 --recode vcf-iid -out ./smp1.2_exon_maf10 --keep-allele-order


plink --vcf dmr1_exon_na.recode.vcf --make-bed --out ./dmr1_exon_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_exon_maf100 --recode vcf-iid -out ./dmr1_exon_maf100 --keep-allele-order
bcftools view -S pop1.1.txt dmr1_exon_maf100.vcf > dmr1.1_exon_na.recode.vcf
bcftools view -S pop1.2.txt dmr1_exon_maf100.vcf > dmr1.2_exon_na.recode.vcf
plink --vcf dmr1.1_exon_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.1_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1.2_exon_na.recode.vcf --maf 0.1 --make-bed --out ./dmr1.2_exon_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1.1_exon_maf10 --recode vcf-iid -out ./dmr1.1_exon_maf10 --keep-allele-order
plink --bfile ./dmr1.2_exon_maf10 --recode vcf-iid -out ./dmr1.2_exon_maf10 --keep-allele-order

# exons non gbm

plink --vcf smp1_exon_nongbm.recode.vcf --make-bed --out ./smp1_exon_nongbm_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1_exon_nongbm_maf100 --recode vcf-iid -out ./smp1_exon_nongbm_maf100 --keep-allele-order
bcftools view -S pop1.1.txt smp1_exon_nongbm_maf100.vcf > smp1.1_exon_nongbm.recode.vcf
bcftools view -S pop1.2.txt smp1_exon_nongbm_maf100.vcf > smp1.2_exon_nongbm.recode.vcf
plink --vcf smp1.1_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./smp1.1_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf smp1.2_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./smp1.2_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./smp1.1_exon_nongbm_maf10 --recode vcf-iid -out ./smp1.1_exon_nongbm_maf10 --keep-allele-order
plink --bfile ./smp1.2_exon_nongbm_maf10 --recode vcf-iid -out ./smp1.2_exon_nongbm_maf10 --keep-allele-order

plink --vcf dmr1_exon_nongbm.recode.vcf --make-bed --out ./dmr1_exon_nongbm_maf100 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1_exon_nongbm_maf100 --recode vcf-iid -out ./dmr1_exon_nongbm_maf100 --keep-allele-order
bcftools view -S pop1.1.txt dmr1_exon_nongbm_maf100.vcf > dmr1.1_exon_nongbm.recode.vcf
bcftools view -S pop1.2.txt dmr1_exon_nongbm_maf100.vcf > dmr1.2_exon_nongbm.recode.vcf
plink --vcf dmr1.1_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./dmr1.1_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --vcf dmr1.2_exon_nongbm.recode.vcf --maf 0.1 --make-bed --out ./dmr1.2_exon_nongbm_maf10 --set-missing-var-ids @:# --keep-allele-order
plink --bfile ./dmr1.1_exon_nongbm_maf10 --recode vcf-iid -out ./dmr1.1_exon_nongbm_maf10 --keep-allele-order
plink --bfile ./dmr1.2_exon_nongbm_maf10 --recode vcf-iid -out ./dmr1.2_exon_nongbm_maf10 --keep-allele-order


#############################################################################

# Distance matrix

/mnt/c/Users/Andreas/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./smp/smp29_gbm_maf04.vcf  -OutPut ./distance_matrix/smp29_gbm_dis.mat
/mnt/c/Users/Andreas/ubuntu/softwares/VCF2Dis-1.50/bin/VCF2Dis -InPut  ./dmr/dmr29_gbm_maf04.vcf  -OutPut ./distance_matrix/dmr29_gbm_dis.mat

```

49. Calculate Theta Pi and Dm
>*split file into exon regions*
```
cat ./exons2.list | while read -r line ; do samtools faidx GCF_000005505.3_Brachypodium_distachyon_v3.0_genomic_chromnames_rmscaff_chrnames.fna $line >>exons.fasta; done

# in conda environment
faSplit byname exons.fasta exons_fasta/

cd exons_fasta/
for file in *fa ; do ls $file >>filenames_exons ; done
```
>*Count CG sites per exon region*
```
library("methimpute")
library(data.table)

## Only CG context
## genes
files <- fread(file="./filenames_exons", header = F)
files$V1 <- gsub(":", "", files$V1)
files <- as.character(files$V1)

results <- vector("list", length(files))

for (f in seq_along(files)) {
  print(f)
  cytosines <- NULL
  try(cytosines <- extractCytosinesFromFASTA(files[f], contexts = 'CG'), silent = TRUE)
  if (!is.null(cytosines) && length(cytosines) > 0) {
    context_table <- table(cytosines$context)
    results[[f]] <- c(files[f], as.list(context_table))
  } else {
    results[[f]] <- c(files[f], 0)
  }
}

results <- rbindlist(lapply(results, as.data.table))
fwrite(results, "cytosine_count_exons.txt", sep = "\t", col.names = FALSE, quote = FALSE)
```
>*Split vcfs*
```
bgzip -f smp1_exon_maf04.vcf
tabix -f smp1_exon_maf04.vcf.gz

bgzip -f smp6_exon_maf04.vcf
tabix -f smp6_exon_maf04.vcf.gz

bgzip -f smp1.1_exon_maf10.vcf
tabix -f smp1.1_exon_maf10.vcf.gz

bgzip -f smp1.2_exon_maf10.vcf
tabix -f smp1.2_exon_maf10.vcf.gz



cat exons2.list | while read -r line ; do tabix smp1_exon_maf04.vcf.gz $line >./group1/exons/$line.smp.vcf; done
cat exons2.list | while read -r line ; do tabix smp6_exon_maf04.vcf.gz $line >./group6/exons/$line.smp.vcf; done

cat exons2.list | while read -r line ; do tabix smp1.1_exon_maf10.vcf.gz $line >./group1.1/exons/$line.smp.vcf; done
cat exons2.list | while read -r line ; do tabix smp1.2_exon_maf10.vcf.gz $line >./group1.2/exons/$line.smp.vcf; done

# for all groups
for file in *vcf ; do wc -l $file >>vcflengths_exon ; done
```

>*Get good and bad intervals*
```
library(dplyr)

###### GROUP 1
###### EXON regions
setwd("/group1/exons/")

vcflengths_exon <- read.table(file="vcflengths_exon")
vcflengths_exon$V2 <- gsub(".smp.vcf","",vcflengths_exon$V2)
colnames(vcflengths_exon) <- c("numvar","interval")
cytosine_count <- read.table(file="C:/Users/Andreas/Desktop/desktop_thesis/data/vcf_redo/stats_thresh4/cytosine_count_exons.txt")
cytosine_count$V1 <- gsub(".fa","",cytosine_count$V1)
cytosine_count$V1 <- gsub("", ":", cytosine_count$V1)
colnames(cytosine_count) <- c("interval","numc")

merged <- inner_join(cytosine_count,vcflengths_exon,by="interval")
merged$prop <- merged$numvar/merged$numc
merged2 <- merged[merged$numvar >= 3 & merged$prop > 0.05,]
write.table(merged2,file="good_intervals",sep="\t",quote=F,row.names = F, col.names = F)

merged4 <- merged[!(merged$numvar >= 3 & merged$prop > 0.05),]
write.table(paste("rm ",merged4$interval,".input.txt",sep=""),file="bad_intervals.sh",sep="\t",quote=F,row.names = F, col.names = F)

## Repeat for groups 6, 1.1 and 1.2 with according vcf files
```
>*remove bad intervals and alpha estimation*
```
cd group1/exons
for file in *.smp.vcf; do sed '/##/d' $file | cut -f 1,2,10- | sed 's/\/.//g' -  >${file/.smp.vcf/.input.txt} ; done 
cut -f 1-2 good_intervals | sed 's/\t/.input.txt\t/' >length_list
mkdir input
for file in *.input.txt; do mv $file input/ ; done
dos2unix bad_intervals.sh
mv bad_intervals.sh input/
cd input/
./bad_intervals.sh
cd ../
perl alpha_estimation.pl -dir ./input -output  alpha_Dm_brachy -length_list length_list
# repeat for other groups



library(dplyr)
# run for each group and region
# group 1 exons
setwd("C:/Users/Andreas/Desktop/desktop_thesis/data/vcf_redo/stats_thresh4/group1/exons/")
good_intervals <- read.table(file="good_intervals")
alpha_dm <- read.table(file="alpha_Dm_brachy")
alpha_dm$V1 <- gsub(".input.txt","",alpha_dm$V1)
alpha_dm <- alpha_dm[1:2]
good_intervals <- inner_join(good_intervals,alpha_dm,by="V1")

write.table(good_intervals,file="good_intervals_alpha",sep="\t",quote=F,row.names = F, col.names = F)
```
>*calculate statistics for all groups*
```
#group1
cd /mnt/c/Users/Andreas/Desktop/desktop_thesis/data/vcf_redo/stats_thresh4/group1/exons
dos2unix good_intervals_alpha
cp good_intervals_alpha input/
cd input/
cat good_intervals_alpha |  while read -r value1 value2 value3 value4 value5 remainder ;  do perl Dm_test_new.pl -input $value1.input.txt -output $value1.Dm_group1_exons.txt -length $value2 -alpha $value5  ; done


#in R 
library(dplyr)
# EXONS group 1

setwd("C:/Users/Andreas/Desktop/desktop_thesis/data/vcf_redo/stats_thresh4/group1/exons/input/")
list <- read.table("length_list", header = FALSE, sep = "\t")
list$V1 <- gsub("input.txt", "Dm_group1_exons.txt", list$V1)
list$V1 <- gsub(":", "", list$V1)



e <- read.table(list[1,1], header = F, sep = "\t")
name <- gsub(".Dm_group1_exons.txt", "", list[1,1])
e$region <- name
for (i in 2:nrow(list)){
  region <- gsub(".Dm_group1_exons.txt", "", list[i,1])
  e2 <- cbind(read.table(list[i,1], header = F), region)
  e <- rbind(e, e2)
}
print(e)
colnames(e) <- c("#chr", "start", "stop", "Dm", "segregation_sites", "theta_pi", "theta_s", "region")
e$group <- "Group 1"
print(e)

write.table(e, file = "Dm_theta_exons_group1.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#repeat for other groups
```

50. Thin vcf for DAPC
```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --recode --out smp29_gbm_maf04_thin --vcf smp29_gbm_maf04.vcf --thin 3000
```

51. Calculate LD decay

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp1.1_gbm_maf10.vcf --out smp1.1_gbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp1.2_gbm_maf10.vcf --out smp1.2_gbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp6_gbm_maf04.vcf --out smp6_gbm_maf04 --freq

/proj/popgen/a.ramesh/software/PopLDdecay/bin/PopLDdecay  -InVCF smp1.1_gbm_maf10.vcf -OutStat smp1.1_gbm_maf10_ld_decay -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/bin/PopLDdecay  -InVCF smp1.2_gbm_maf10.vcf -OutStat smp1.2_gbm_maf10_ld_decay -MaxDist 10
/proj/popgen/a.ramesh/software/PopLDdecay/bin/PopLDdecay  -InVCF smp6_gbm_maf04.vcf -OutStat smp6_gbm_maf04_ld_decay -MaxDist 10
```

51. Calculate mSFS.

```

/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf smp1.1_gbm_maf10.vcf --out smp1.2_gbm_maf10 --freq
/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf smp1.2_gbm_maf10.vcf --out smp1.2_gbm_maf10 --freq
/data/proj2/popgen/a.ramesh/software/vcftools/src/cpp/vcftools --vcf smp6_gbm_maf04.vcf --out smp6_gbm_maf04 --freq

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp1.1_exon_nongbm_maf10.vcf --out smp1.1_exon_nongbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp1.2_exon_nongbm_maf10.vcf --out smp1.2_exon_nongbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf smp6_exon_nongbm_maf04.vcf --out smp6_exon_nongbm_maf04 --freq

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr1.1_gbm_maf10.vcf --out dmr1.1_gbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr1.2_gbm_maf10.vcf --out dmr1.2_gbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr6_gbm_maf04.vcf --out dmr6_gbm_maf04 --freq

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr1.1_exon_nongbm_maf10.vcf --out dmr1.1_exon_nongbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr1.2_exon_nongbm_maf10.vcf --out dmr1.2_exon_nongbm_maf10 --freq
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf dmr6_exon_nongbm_maf04.vcf --out dmr6_exon_nongbm_maf04 --freq

```

52. Add header to variant plus invariant vcf
```
cat vcf_header brachy_meth_var_invar.vcf >brachy_meth_var_invar_all.vcf
```

53. get number of invariant sites for all sites and for gbm sites only.

```
cp gbm_genes.bed gbm_genes_bd.bed 

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --bed  gbm_genes.bed --max-missing 0.8 --recode --vcf brachy_meth_var_invar_all.vcf --out brachy_meth_var_invar_all
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group1.1 --vcf brachy_meth_var_invar_all.vcf  --out group1.1_meth_invar --maf 0 --max-missing 0.9
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group1.2 --vcf brachy_meth_var_invar_all.vcf  --out group1.2_meth_invar --maf 0 --max-missing 0.9
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group6 --vcf brachy_meth_var_invar_all.vcf  --out group6_meth_invar --maf 0 --max-missing 0.9

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group1.1 --vcf brachy_meth_var_invar_all.recode.vcf  --out group1.1_meth_invar_gbm --maf 0 --max-missing 0.9
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group1.2 --vcf brachy_meth_var_invar_all.recode.vcf  --out group1.2_meth_invar_gbm  --maf 0 --max-missing 0.9 
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --keep group6 --vcf brachy_meth_var_invar_all.recode.vcf  --out group6_meth_invar_gbm --maf 0 --max-missing 0.9


/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools   --vcf smp1.1_gbm_maf10.vcf   --out test1.1 --maf 0.005
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf smp1.2_gbm_maf10.vcf   --out test1.2 --maf 0.005
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools  --vcf smp6_gbm_maf10.vcf   --out test6 --maf 0.005

```

54. Generate multihetsep files for demographic inference. 

```
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf brachy_meth_var_invar_all.vcf --recode  --out group1.1_smps_5mb --keep group1.1 --max-missing 0.9 --bed popgen5mb2_bd.bed
/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --vcf brachy_meth_var_invar_all.vcf --recode  --out group1.2_smps_5mb --keep group1.2 --max-missing 0.9 --bed popgen5mb2_bd.bed

cd smps_group1.1/
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f group1.1_smps_5mb.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f group1.1_smps_5mb.recode.vcf.gz
cat group1.1 |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c1 -O v -s $sample -o $sample.filtered.vcf group1.1_smps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat sample_chr_group1.1 |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >group1.1_multihetsep_meth_$line ; done
Rscript multihetsep_combined.R
cd ../

cd smps_group1.2/
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f group1.2_smps_5mb.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix -f group1.2_smps_5mb.recode.vcf.gz
cat group1.2 |  while read -r  sample remainder; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools view -c1 -O v -s $sample -o $sample.filtered.vcf group1.2_smps_5mb.recode.vcf.gz ; done
for file in *.filtered.vcf ; do /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools annotate -x INFO,^FORMAT/GT -O v -o ${file/.filtered/.annotated} $file ; done
for file in *.annotated.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.annotated.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat sample_chr_group1.2 |  while read -r value1 value2 remainder ;  do /proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf $value1.annotated.vcf.gz --out $value1.$value2.snps --min-alleles 1 --max-alleles 3 --recode --recode-INFO-all  --chr $value2 ; done
for file in *.snps.recode.vcf  ; do /proj/popgen/a.ramesh/software/htslib-1.16/bgzip -f $file; done
for file in *.snps.recode.vcf.gz  ; do /proj/popgen/a.ramesh/software/htslib-1.16/tabix -f $file; done
cat chrlist | while read line; do /proj/popgen/a.ramesh/software/msmc-tools/generate_multihetsep.py --mask ../$line.cytosines_tair10.bed  --chr $line *.$line.snps.recode.vcf.gz >group1.2_multihetsep_meth_$line ; done
Rscript multihetsep_combined.R
cd ../
```

