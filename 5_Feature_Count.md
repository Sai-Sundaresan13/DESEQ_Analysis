Generation of a count matrix containing the number of raw reads using feature counts tool.


## 1. Download a reference .gtf file
The reference .gtf file conataining the gene structure is required to run feature count tool. this can be downloaded using a wget link or from ensemble.

## 2. Running Featurecounts Tool
The below command can be run to perform feature counts on all the .bam files (irrespective of the distinguishing criteria) and outputs a .txt file containing all the counts.

```
featureCounts -T 4 -p -a annotation.gtf -o counts.txt *.bam
```
