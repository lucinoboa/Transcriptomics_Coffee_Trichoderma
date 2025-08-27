# Transcriptomics_Coffee_Trichoderma
Coffee plants were inoculated with Trichoderma rifaii CT5, and RNA was extracted from leaf tissue at three time points: prior to inoculation (Day 0), one day post-inoculation (Day 1), and two days post-inoculation (Day 3). Samples were collected from both inoculated plants (“Tricho”) and non-inoculated controls (“Control”).
This repository documents the complete workflow for differential expression analysis, including quality control, read trimming, sequence alignment, and quantification of gene expression.

## Step 1: Quality Control of Raw Data 
### Navigate to the folder that contains the raw ".fq.gz" files.  
```bash
cd /data2/lnoboa/tricho_transcriptomics/raw_data_tricho
```
### Activate the conda environment. 
```bash
conda activate transcriptomics_env
```
### Run FastQC to assess the quality of the sequences.  
```bash
fastqc *.fq.gz -o /data2/lnoboa/tricho_transcriptomics/fastqc_tricho/rawdata_fastqc
```
### Run MultiQC in the same folder. 
```bash
multiqc .
```
[Check the report](Results/multiqc_report_tricho-coffee-rawdata.html)

## Step 2: Trimming 
###Run the Trimmomatic tool to remove low quality reads and adapters. 
```bash
TRIMMOMATIC=/home/jupyter-alumno7/.conda/envs/rnseq/bin/trimmomatic
ADAPTERS=/home/jupyter-alumno7/.conda/envs/rnseq/share/trimmomatic/adapters/TruSeq3-PE.fa
THREADS=4
INPUT_DIR=/data2/lnoboa/tricho_transcriptomics/raw_data_tricho
OUTPUT_DIR=/data2/lnoboa/tricho_transcriptomics/trimmed_tricho

for f1 in $INPUT_DIR/*_1.fq.gz; do
    base=$(basename "$f1" _1.fq.gz)
    f2="$INPUT_DIR/${base}_2.fq.gz"
     
    $TRIMMOMATIC PE -threads $THREADS \
        "$f1" "$f2" \
        "$OUTPUT_DIR/${base}_1.trim.fq.gz" "$OUTPUT_DIR/${base}_1un.trim.fq.gz" \
        "$OUTPUT_DIR/${base}_2.trim.fq.gz" "$OUTPUT_DIR/${base}_2un.trim.fq.gz" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
  done
```

## Step 3: Quality Control of the Trimmed Sequences 
### Navigate to the folder that contains the trimmed ".fq.gz" files.  
```bash
cd /data2/lnoboa/tricho_transcriptomics/trimmed_tricho
```
### Activate the conda environment. 
```bash
conda activate transcriptomics_env
```
### Run FastQC to assess the quality of the sequences.  
```bash
fastqc *.fq.gz -o /data2/lnoboa/tricho_transcriptomics/fastqc_tricho/trimmed_fastqc
```
### Run MultiQC in the same folder. 
```bash
multiqc .
```
[Check the report](Results/multiqc_trichotranscriptomics_trimmed.html)

## Step 4: Alignment to the Reference Genome 
### Create a folder for results. 
```bash
cd /data2/lnoboa/tricho_transcriptomics/
mkdir mapping_results
```
### Upload the genome index. 
[Reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036785885.1/)
```bash
cd /data2/lnoboa/ref_genome_coffea
hisat2-build -p 4 GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna coffea_index
```
### Mapping of sequences. 
#### Alignment of the trimmed sequences to the reference genome. 
```bash
cd /data2/lnoboa/tricho_transcriptomics/trimmed_tricho
for sample in *_1.trim.fq.gz; do
    base=$(basename $sample _1.trim.fq.gz)
    hisat2 -p 4 \
        -x /data2/lnoboa/ref_genome_coffea/coffea_index \
        -1 ${base}_1.trim.fq.gz \
        -2 ${base}_2.trim.fq.gz \
        -S /data2/lnoboa/mapping_results/${base}.sam
done
```
### Convert SAM to BAM. 
```bash
cd /data2/lnoboa/tricho_transcriptomics/mapping_results
for f in *.sam; do
    base=$(basename "$f" .sam)
    samtools sort "$f" -o "${base}_sorted.bam"
```
## Step 5: Quantification of Mapped Reads. 
```bash
/data2/lnoboa/tricho_transcriptomics/mapping_results
featureCounts -p -t exon -g gene_id   -a /data2/lnoboa/ref_genome_coffea/genomic.gtf   -o counts_matrix.txt /data2/lnoboa/tricho_transcriptomics/mapping_results/*_sorted.bam
```
