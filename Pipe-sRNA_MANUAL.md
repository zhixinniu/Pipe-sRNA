## 1. Quality control

This part is comprised of several steps and software to check quality, trim reads, remove adapters, filter length, and trace contamination. It is involved in the following packages: fastqc, fastp, and mirtrace.

### 1.1 Installation

``` bash
conda create -n smrna -y fastqc fastp mirtrace
conda activate smrna
```

### 1.2 Quality check

``` bash
fastqc -j $THREADS -o $OUTPUT_DIR ${FASTQ}.fq
```

### 1.3 UMI extraction and miRNA adapter trimming

(*Optional*) Please refer to `UMI-tools` and `nf-core/smrnaseq`

### 1.4 Adapter trimming

``` bash
fastp
```

------------------------------------------------------------------------

## 2. UMI deduplication

(*Optional*) Please refer to `nf-core/smrnaseq`.

------------------------------------------------------------------------

## 3. miRNA QC

miRNA quality check, including PHRED score, taxonomic classification, read length, and proportion of different RNAs.

### 3.1 Installation

``` bash
conda install mirtrace
```

### 3.2 Usage

``` bash
mirtrace qc -t 30 -s hsa ${FASTQ}.fq  ## run qc mode
mirtrace trace qc -t 30 -s hsa ${FASTQ}.fq  # run trace mode
```

------------------------------------------------------------------------

## 4. Contamination filtering

Remove reads of tRNA, rRNA, ncRNA, cDNA, piRNA

### 4.1 Installation

``` bash
conda install -y bowtie blat seqkit
```

### 4.2 Make bowtie index

``` bash
bowtie-build ${FASTA}.fa ${PREFIX_INDEX}
```

> cDNA: Homo_sapiens.GRCh38.cdna.all.fa.gz <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/>\
> ncRNA: Homo_sapiens.GRCh38.ncrna.fa.gz <https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/>\
> tRNA: hg38-tRNAs.fa <https://gtrnadb.ucsc.edu/>\
> piRNA: hsa.v3.0.fa.gz <http://bigdata.ibp.ac.cn/piRBase/download.php>\
> rRNA: rRNA_reference.fa <https://github.com/friedlanderlab/mirtrace>\
> miRNA: hairpin and mature <https://www.mirbase.org/>

### 4.3 Filtering

#### 4.3.1 Filter rRNA

``` bash
# Align raw reads to rRNA reference, and output unaligned reads (without rRNA) for the next step.
# Warning: some bowtie1 versions don't support ".gz" format.
bowtie -v 1 --threads $THREADS --un ${RRNA_UNALIGNED}.fq $RRNA_INDEX ${FASTQ}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${RRNA_FASTA}.fa -o ${RRNA_BAM}.bam -

# After all alignments, you can combine $LOG to a single file showing mapping statistics.
echo -e "Sample\tReads_Processed\tAligned_Reads\tFailed_Reads" > align_sum.txt
for file in $(ls *.log)
do
  id=$(basename "$file" .log)
  reads_processed=$(awk 'NR==1{print $4}' "$file")
  aligned=$(awk 'NR==2{print $9,$10}' "$file")
  failed=$(awk 'NR==3{print $7,$8}' "$file")
  echo -e "${id}\t${reads_processed}\t${aligned}\t${failed}"
done >> align_sum.txt
```

#### 4.3.2 Filter tRNA

``` bash
# Align unmapped reads from the last step, ${RRNA_UNALIGNED}.fq, to tRNA reference
bowtie -v 1 --threads $THREADS --un ${TRNA_UNALIGNED}.fq $TRNA_INDEX ${RRNA_UNALIGNED}.fq 2 > ${LOG}.log |\
samtools view -bS --threads $THREADS --reference ${TRNA_FASTA}.fa -o ${TRNA_BAM}.bam -
```

#### 4.3.3 Filter cDNA

``` bash
# Search which hairpin miRNA are present in cDNA data
# Get hairpin miRNA only for human
seqkit grep -r -p "hsa" ${HAIRPIN_MIRNA}.fa > ${HSA_HAIRPIN}.fa  

# blat <database> <query>
blat ${CDNA}.fa ${HSA_HAIRPIN}.fa ${OUTPUT}.psl

# Extract significant BLAT hits
```
