# Running HISAT2 Alignment for Trimmed Files
High-quality reads are aligned to the reference genome using HISAT2. This step produces SAM files representing mapped read positions. The reference genome used for this project was Human Genome GRCh 38.


## 1. Download the Reference Files

Download the reference .fasta files using the wget link provided below and index them in order to use then for alignment. Use the wget link for the downlaod which can be downloaded from ensemble.

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

Index the reference files using the below command:

```
hisat2-build GRCh38.fasta GRCh38_index
```


## 2. Test if HISAT2 is Working

Run the following command to check if HISAT2 is functioning correctly on a single file:

```
hisat2 -p 4 -x hg38_index -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -S sample_output.sam
```

## OR

```
hisat2 -x grch38_index \
-1 ERR10851963_1_paired.fq.gz \
-2 ERR10851963_2_paired.fq.gz \
--threads 4 \
-S sample.bam
```

The threads correspond to the number of CPU cores used for the process which can be adjusted based on the computational power available.



## SORTING .SAM FILES AND CONVERSION TO .BAM

SAM files are converted to sorted and indexed BAM files using SAMtools. BAM files are more efficient for storage and further processing.

## 1. One-Liner Command: Convert, Sort, and Index in One Go

To perform all three steps (convert to BAM, sort, and index) in a single command:

```
for file in *.sam; do
    samtools view -bS "$file" | samtools sort -o "${file%.sam}.bam"
    samtools index "${file%.sam}.bam"
done
```

## 2. Running Commands Separately

If you prefer to do each step individually, use the following:

### a. Convert SAM to BAM:

```
for file in *.sam; do
    samtools view -bS "$file" -o "${file%.sam}.bam"
done
```

### b. Sort BAM Files:

```
for file in *.bam; do
    samtools sort "$file" -o "${file%.bam}_sorted.bam"
done
```

### c. Index Sorted BAM Files:

```
for file in *_sorted.bam; do
    samtools index "$file"
done
```
