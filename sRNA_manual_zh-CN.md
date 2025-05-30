## 1. 测序数据质量评估

This part is comprised of several steps and software to check quality, trim reads, remove adapters, filter length, and trace contamination. It is involved in the following packages: fastqc, fastp, and mirtrace.

### 1.1 安装

``` bash
conda create -n smrna -y fastqc fastp mirtrace
conda activate smrna
```

### 1.2 检查质量

``` bash
fastqc -j $THREADS -o $OUTPUT_DIR ${FASTQ}.fq
multiqc .
```
> 指控部分待完成  

### 1.3 提取UMI

(*Optional*) Please refer to `UMI-tools` and `nf-core/smrnaseq`

### 1.4 接头去除

``` bash
fastp
```

------------------------------------------------------------------------

## 2. UMI去重

(*Optional*) 参考GitHub `nf-core/smrnaseq`.

------------------------------------------------------------------------

## 3. miRNA指控

miRNA quality check, including PHRED score, taxonomic classification, read length, and proportion of different RNAs.

### 3.1 安装

``` bash
conda install mirtrace
```

### 3.2 使用

``` bash
mirtrace qc -t 30 -s hsa ${FASTQ}.fq  ## run qc mode
mirtrace trace qc -t 30 -s hsa ${FASTQ}.fq  # run trace mode
```

------------------------------------------------------------------------

## 4. 去除RNA污染

去除包含tRNA, rRNA, ncRNA, cDNA, piRNA的reads

### 4.1 安装

``` bash
conda install -y bowtie blat seqkit
```

### 4.2 准备参考数据

cDNA: Homo_sapiens.GRCh38.cdna.all.fa.gz <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/>\
ncRNA: Homo_sapiens.GRCh38.ncrna.fa.gz <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/>\
tRNA: hg38-tRNAs.fa <https://gtrnadb.ucsc.edu/>\
piRNA: hsa.v3.0.fa.gz <http://bigdata.ibp.ac.cn/piRBase/download.php>\
rRNA: rRNA_reference.fa <https://github.com/friedlanderlab/mirtrace>\
miRNA: hairpin and mature <https://www.mirbase.org/>

### 4.3 过滤

#### 4.3.1 过滤rRNA

``` bash
# Align raw reads to rRNA reference, and output unaligned reads (without rRNA) for the next step.
# Warning: some bowtie1 versions don't support ".gz" format.
bowtie -v 1 --threads $THREADS --un ${RRNA_UNALIGNED}.fq $RRNA_INDEX ${FASTQ}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${RRNA_FASTA}.fa -o ${RRNA_BAM}.bam -

# After all alignments, you can combine $LOG to a single file showing mapping statistics.
function stats_contaminant(){
  output="$1"
  
  shift
  
  if [ -f "$output" ];then
    rm "$output"
  fi
  
  echo -e "Sample\tReads_Processed\tAligned_Reads\tFailed_Reads" > "$output"
  
  for file in "$@";do
    id=$(basename "$file" .log)
    reads_processed=$(awk 'NR==1{print $4}' "$file")
    aligned=$(awk 'NR==2{print $9,$10}' "$file")
    failed=$(awk 'NR==3{print $7,$8}' "$file")
    echo -e "${id}\t${reads_processed}\t${aligned}\t${failed}" >> "$output"
  done
  
}

# Define the output file at the first position, and input $LOG behind output file name.
stats_contaminant ${LOG_SUM}.txt ${LOG}.log <$LOG2,$LOG3,$LOG4,...>
```

#### 4.3.2 过滤tRNA

``` bash
# Align unmapped reads from the last step, ${RRNA_UNALIGNED}.fq, to tRNA reference
bowtie -v 2 --threads $THREADS --un ${TRNA_UNALIGNED}.fq $TRNA_INDEX ${RRNA_UNALIGNED}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${TRNA_FASTA}.fa -o ${TRNA_BAM}.bam -

# Mapping summary
stats_contaminant ${LOG_SUM}.txt ${LOG}.log <$LOG2,$LOG3,$LOG4,...>
```

#### 4.3.3 过滤cDNA

``` bash
# Search which hairpin miRNA are present in cDNA data

# blat <database> <query>
blat -out=blast8 -q=rna ${HSA_CDNA}.fa ${HSA_HAIRPIN}.fa ${OUTPUT}.blast8

# Extract significant BLAT hits
awk -v FS="\t" '{if($11 < 1e-5) print $2}' ${OUTPUT}.blast8 | sort | uniq > ${HSA_HAIRPIN_UNIQ}.txt  

# Remove the hairpin miRNA from the cDNA data
seqkit grep -v -f ${HSA_HAIRPIN_UNIQ}.txt ${HSA_CDNA}.fa > ${HSA_CDNA_NO_MIRNA}.fa

# Build index for ${HSA_CDNA_NO_MIRNA}.fa
bowtie-build ${HSA_CDNA_NO_MIRNA}.fa ${HSA_CDNA_NO_MIRNA}.ebwt

# Align filtered reads from 4.3.2 to cDNA
bowtie -v 1 --threads $THREADS --un ${CDNA_UNALIGNED}.fq ${HSA_CDNA_NO_MIRNA}.ebwt ${TRNA_UNALIGNED}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${HSA_CDNA_NO_MIRNA}.fa -o ${CDNA_BAM}.bam -

# Mapping summary
stats_contaminant ${LOG_SUM}.txt ${LOG}.log <$LOG2,$LOG3,$LOG4,...>


# ${HSA_CDNA_NO_MIRNA}.fa: ~/zniu_ws/ref/human/hg38_sRNA/GRCh38_cDNA_no_mirna.fa
# ${HSA_CDNA_NO_MIRNA}.ebwt: ~/zniu_ws/ref/human/hg38_sRNA/bowtie_index/GRCh38_cDNA_no_mirna
# ${OUTPUT}.blast8: ~/zniu_ws/ref/human/hg38_sRNA/BLAT/which_hairpin_in_cDNA.blast8
# ${HSA_HAIRPIN_UNIQ}.txt: ~/zniu_ws/ref/human/hg38_sRNA/sig_hsa_mirna_cdna.txt
# ${HSA_CDNA}.fa: ~/zniu_ws/ref/human/hg38_sRNA/Homo_sapiens.GRCh38.cdna.all.fa.gz
# ${HAIRPIN_MIRNA}.fa: ~/zniu_ws/ref/human/hg38_sRNA/hairpin.fa
# ${HSA_HAIRPIN}.fa: ~/zniu_ws/ref/human/hg38_sRNA/hairpin_hsa.fa
```

#### 4.3.4 过滤ncRNA

``` bash
# Search which hairpin miRNAs are present in the ncRNA data
blat -out=blast8 -q=rna ${HSA_NCRNA}.fa ${HSA_HAIRPIN}.fa ${OUTPUT}.blast8

# Extract the significant hits
awk -v FS="\t" '{if($11 < 1e-5) print $2}' ${OUTPUT}.blast8 | sort | uniq > ${HSA_HAIRPIN_UNIQ}.txt 

# Remove the hairpin miRNAs from the ncRNA data
seqkit grep -v -f ${HSA_HAIRPIN_UNIQ}.txt ${HSA_NCRNA}.fa > ${HSA_NCRNA_NO_MIRNA}.fa

# Build bowtie index for ${HSA_NCRNA_NO_MIRNA}.fa
bowtie-build ${HSA_NCRNA_NO_MIRNA}.fa ${HSA_NCRNA_NO_MIRNA}.ebwt

# Map which reads are ncRNA
bowtie -v 1 --threads $THREADS --un ${NCRNA_UNALIGNED}.fq ${HSA_NCRNA_NO_MIRNA}.ebwt ${CDNA_UNALIGNED}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${HSA_NCRNA_NO_MIRNA}.fa -o ${NCRNA_BAM}.bam -

# Mapping summary
stats_contaminant ${LOG_SUM}.txt ${LOG}.log <$LOG2,$LOG3,$LOG4,...>


# ${HSA_NCRNA}.fa: ~/zniu_ws/ref/human/hg38_sRNA/Homo_sapiens.GRCh38.ncrna.fa.gz
# ${HSA_HAIRPIN}.fa: ~/zniu_ws/ref/human/hg38_sRNA/hairpin_hsa.fa
# ${OUTPUT}.blast8: ~/zniu_ws/ref/human/hg38_sRNA/BLAT/which_hairpin_in_ncRNA.blast8
# ${HSA_HAIRPIN_UNIQ}.txt: ~/zniu_ws/ref/human/hg38_sRNA/sig_hsa_mirna_ncrna.txt
# ${HSA_NCRNA_NO_MIRNA}.fa: ~/zniu_ws/ref/human/hg38_sRNA/GRCh38_ncRNA_no_mirna.fa
# ${HSA_NCRNA_NO_MIRNA}.ebwt: ~/zniu_ws/ref/human/hg38_sRNA/bowtie_index/GRCh38_ncRNA_no_mirna
# RRNA_INDEX=/mnt/workspace_lgu/zniu/ref/human/hg38_sRNA/bowtie_index/rRNA_bowtie
# TRNA_INDEX=/mnt/workspace_lgu/zniu/ref/human/hg38_sRNA/bowtie_index/tRNA_bowtie
# CDNA_INDEX=/mnt/workspace_lgu/zniu/ref/human/hg38_sRNA/bowtie_index/GRCh38_cDNA_no_mirna_new
# NCRNA_INDEX=/mnt/workspace_lgu/zniu/ref/human/hg38_sRNA/bowtie_index/GRCh38_ncRNA_no_mirna_new
```

#### 4.3.5 过滤piRNA

``` bash
Not defined
```



## 5. miRNA定量

### 5.1 解析miRBase数据库的mature.fa和hairpin.fa 

#### 5.1.1 安装

``` bash
conda install fastx_toolkit
```

#### 5.1.2 解析miRBase
``` bash
# 第一个参数：miRNA fasta文件；第二个参数：物种名<hsa/mmu/dme>；第三个参数：输出文件名
function parse_mir(){
  FASTA="$1"
  SPS="$2"
  OUT="$3"
  
  if [[ "${FASTA: -3}" == ".gz" ]];then
    gunzip -f "$FASTA"
    FASTA="${FASTA%%.*}"
  fi
  
  sed 's/&gt;/>/g' "$FASTA"|sed 's#<br>#\n#g'|sed 's#</p>##g'|sed 's#<p>##g'|sed -e :a -e '/^\n*$/{$d;N;};/\n$/ba' > tmp_clean_html.fa
  sed '/^[^>]/s/[^AUGCaugc]/N/g' tmp_clean_html.fa > tmp_parsed.fa
  sed -i 's/\s.*//' tmp_parsed.fa
  seqkit grep -r -p ".*${SPS}.*" tmp_parsed.fa > tmp_sps.fa
  seqkit seq --rna2dna tmp_sps.fa > tmp_dna.fa
  fasta_formatter -w 0 -i tmp_dna.fa -o ${OUT} 
  rm tmp_sps.fa tmp_parsed.fa tmp_clean_html.fa tmp_dna.fa
}

parse_mir hairpin.fa hsa hairpin_hsa_parsed.fa
parse_mir mature.fa hsa mature_hsa_parsed.fa
```

#### 5.1.3 比对reads到mature miRNA

``` bash
# 构建maure_hsa_parsed.fa的索引
bowtie-build mature_hsa_parsed.fa ${MATURE_HSA}.ebwt

# 比对到mature miRNA
bowtie -k 50 -e 99999 --best --strata --sam --chunkmbs 2048 --threads 20 --un ${MATURE_UNMAP}.fq ${MATURE_HSA}.ebwt ${NCRNA_UNALIGNED}.fq 2>${LOG}.log |\
samtools view -@ $THREADS -bS -o ${PREFIX}.bam - 

# 构建hairpin_hsa_parsed.fa的索引
bowtie-build hairpin_hsa_parsed.fa ${HAIRPIN_HSA}.ebwt

# 将未比对到mature miRNA的reads再重新比对至hairpin miRNA
bowtie -k 50 -e 99999 --best --strata --sam --chunkmbs 2048 --threads 20 --un ${HAIRPIN_UNMAP}.fq ${HAIRPIN_HSA}.ebwt ${MATURE_UNMAP}.fq 2>${LOG}.log |\
samtools view -@ $THREADS -bS -o ${PREFIX}.bam - 
```
#### 5.1.4 统计定量比对结果
```bash

# 统计比对到mature miRNA
samtools sort -@ 20 --reference mature_hsa_parsed.fa -o ${PREFIX}_sort.bam ${PREFIX}.bam
samtools index -@ 20 ${PREFIX}_sort.bam
samtools stats -@ 20 --reference mature_hsa_parsed.fa ${PREFIX}_sort.bam > ${PREFIX}.stats
samtools flagstat -@ 20 ${PREFIX}_sort.bam > ${PREFIX}.flagstat
samtools idxstats -@ 20 ${PREFIX}_sort.bam > ${PREFIX}.idxstats

function bam_stats(){
  ref="$1"
  shift
  
  for file in "$@"
  do
    prefix=${file%%.*}
    samtools sort -@ 5 --reference $ref -o ${prefix}_sort.bam $file
    samtools index -@ 5 ${prefix}_sort.bam
    samtools stats -@ 5 --reference $ref ${prefix}_sort.bam > ${prefix}.stats
    samtools flagstat -@ 5 ${prefix}_sort.bam > ${prefix}.flagstat
    samtools idxstats -@ 5 ${prefix}_sort.bam > ${prefix}.idxstats
  done
}

bam_stats ${prefix} *.bam
```

在R中输出count table
``` R
# Originally written by Phil Ewels and Chuan Wang and released under the MIT license.
# Contributions by Alexander Peltzer, Anabella Trigila, James Fellows Yates, Sarah Djebali, Kevin Menden, Konrad Stawinski and Lorena Pantano also released under the MIT license. See LICENSE https://github.com/nf-core/smrnaseq/blob/master/LICENSE for details.


library("limma")
library("edgeR")
library("statmod")
library("data.table")
library("gplots")
library("methods")


# 将mature和hairpin的idxstats分别放在一个列表中
filelist <- list()
filelist[[1]] <- list.files('/mnt/workspace_lgu/zniu/project/epi_clock/HL60_RNA/BMK250321-CP281-ZX01-1101/BMK_DATA_20250527124829_1/Data/4_Mapping',pattern = '.*mature.idxstats',full.names = T)  # 把mature的idxstats放在第一个元素中
filelist[[2]] <- list.files('/mnt/workspace_lgu/zniu/project/epi_clock/HL60_RNA/BMK250321-CP281-ZX01-1101/BMK_DATA_20250527124829_1/Data/4_Mapping',pattern = '.*hairpin.idxstats',full.names = T)  # 把hairpin的idxstats放在第一个元素中
names(filelist)<-c("mature","hairpin")

for (i in 1:2){
  header <- names(filelist)[i]
  
  # 准备以geneID为行名，sampleID为列名的表格
  data <- do.call("cbind", lapply(filelist[[i]], fread, header=FALSE, select=c(3)))
  unmapped<-do.call("cbind", lapply(filelist[[i]], fread, header=FALSE, select=c(4)))
  data<-as.data.frame(data)
  unmapped<-as.data.frame(unmapped)
  temp <- fread(filelist[[i]][1],header=FALSE, select=c(1))
  rownames(data)<-temp$V1
  rownames(unmapped)<-temp$V1
  colnames(data)<-gsub("_mature.*","",basename(filelist[[i]]))
  colnames(unmapped)<-gsub("_mature.*","",basename(filelist[[i]]))
  
  data<-data[rownames(data)!="*",,drop=FALSE]
  unmapped<-unmapped[rownames(unmapped)=="*",,drop=FALSE]
  
  # Write the summary table of unmapped reads
  write.table(unmapped,file=paste(header,"_unmapped_read_counts.txt",sep=""),sep='\t',quote=FALSE)
  
  # Remove genes with 0 reads in all samples
  row_sub = apply(data, 1, function(row) all(row ==0 ))
  # Only subset if at least one sample is remaining
  nr_keep <- sum(row_sub)
  if (nr_keep > 0){
    data<-data[!row_sub,, drop=FALSE]
  }
  #Also check for colSums > 0, otherwise DGEList will fail if samples have entirely colSum == 0 #Fixes #134
  drop_colsum_zero <- (colSums(data, na.rm=T) != 0) # T if colSum is not 0, F otherwise
  data <- data[, drop_colsum_zero] # all the non-zero columns
  
  write.csv(t(data),file=paste(header,"_counts.csv",sep=""))
  
  # Normalization
  dataDGE<-DGEList(counts=data,genes=rownames(data))
  o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
  dataDGE <- dataDGE[o,]
  
  # Save log10(TPM)
  tpm = cpm(dataDGE, normalized.lib.sizes=F, log = F, prior.count = 0.001)
  tpm = tpm + 0.001
  tpm = log10(tpm)
  ttpm = t(tpm)
  write.table(ttpm,file=paste(header,"_logtpm.txt",sep=""),sep='\t',quote=FALSE)
  write.csv(ttpm,file=paste(header,"_logtpm.csv",sep=""))
  
  # TMM
  dataNorm <- calcNormFactors(dataDGE)
  
  # Print normalized read counts to file
  dataNorm_df<-as.data.frame(cpm(dataNorm))
  write.table(dataNorm_df,file=paste(header,"_normalized_CPM.txt",sep=""),sep='\t',quote=FALSE)
  
  if (length(filelist[[1]]) > 1){ # with more than 1 sample
    # Print heatmap based on normalized read counts
    pdf(paste(header,"_CPM_heatmap.pdf",sep=""))
    heatmap.2(cpm(dataNorm),col=redgreen(100),key=TRUE,scale="row",density.info="none",trace="none")
    dev.off()
  }
  
  # Make MDS plot (only perform with 3 or more samples)
  if (ncol(dataNorm$counts) > 2){
    pdf(paste(header,"_edgeR_MDS_plot.pdf",sep=""))
    MDSdata <- plotMDS(dataNorm)
    dev.off()
    
    # Print distance matrix to file
    write.table(MDSdata$distance.matrix, paste(header,"_edgeR_MDS_distance_matrix.txt",sep=""), quote=FALSE, sep="\t")
    
    # Print plot x,y co-ordinates to file
    MDSxy = data.frame(x=MDSdata$x, y=MDSdata$y)
    colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))
    
    write.table(MDSxy, paste(header,"_edgeR_MDS_plot_coordinates.txt",sep=""), quote=FALSE, sep="\t")
    
    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)
    
    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))
    
    # Plot a heatmap of correlations
    pdf(paste(header,"_log2CPM_sample_distances_heatmap.pdf",sep=""))
    hmap <- heatmap.2(as.matrix(dists),main="Sample Correlations", key.title="Distance", trace="none",dendrogram="row", margin=c(9, 9))
    dev.off()
    
    # Plot the heatmap dendrogram
    pdf(paste(header,"_log2CPM_sample_distances_dendrogram.pdf",sep=""))
    plot(hmap$rowDendrogram, main="Sample Dendrogram")
    dev.off()
    
    # Write clustered distance values to file
    write.table(hmap$carpet, paste(header,"_log2CPM_sample_distances.txt",sep=""), quote=FALSE, sep="\t")
  } else {
    warning("Not enough samples to create an MDS plot. At least 3 samples are required.")
  }
}


```

















