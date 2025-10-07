# Running FastQC for All Files and MultiQC
Quality control is performed to assess the read quality, GC content, Per base sequence content, adapter content, etc. To perform quality control using FastQC and summarize the results using MultiQC.


Step 1: Navigate to the directory containing the sequencing files
```
cd /path/to/your/sequencing/files
```
Step 2: Run FastQC on all files in the directory
```
fastqc *
```
The above command generates the quality control files containing attributes such as GC content, Per base sequence, Adapters, etc.
Store all the iutput files into a new folder using the command - 
```
mv *.html *.zip results_files/
```

NOTE: both the above can be donw in a single step using -
```
fastqc *.fastq.gz -o fastqc_results/
```

Step 3: Run MultiQC to aggregate FastQC results interactively
```
multiqc * --interactive
```
The multiqc command combines all the result files of fastqc and outputs a single file which eases the analysis of the files.
