
# Running Trim-Galore for FastQC Results
Adapters and low-quality bases are removed using Trim Galore. This step improves read quality and enhances the accuracy of alignment.

## 1. Running Trim Galore on Single-End Files

Basic usage for trimming a single-end FASTQ file:

```
trim_galore --fastqc sample.fastq.gz
```


## 2. Running Trim Galore on Paired-End Files

```
trim_galore --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

The specific quality and length to be trimmed can also be specified by adding additional variables:
```
trim_galore --paired --quality 20 --length 30 --fastqc sample_R1.fastq.gz sample_R2.fastq.gz

```

Batch-processing all `.fastq.gz` files in a directory (used in this project):

```
for file in *_R1.fastq.gz
do
    base=$(basename "$file" "_R1.fastq.gz")

    file1="${base}_R1.fastq.gz"
    file2="${base}_R2.fastq.gz"

    echo "Trimming $file1 and $file2..."

    trim_galore --paired --quality 20 --length 30 --fastqc "$file1" "$file2"
done
```



