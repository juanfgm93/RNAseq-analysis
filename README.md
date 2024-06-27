# RNAseq-analysis
---
### Objective

The purpose of this repository is to estimate differentially expressed genes (DEGs) from a high-throughput sequencing experiment. DEGs were identified by analyzing read counts per gene across all experimental conditions using the DESeq2 Bioconductor package (v.1.20.0), which model count data with negative normal distributions.
---
### Context

Total RNA extracts were initially obtained from SW480 cells treated with siRNAs targeting either LUC, DIS3L2 or DIS3L2+TUTs. RNA libraries were prepared from three independent biological replicates for each condition and enriched for long RNAs (> 200 nucleotides) by poly(A) selection. Next, mRNA was fragmented and converted to first strand cDNA using reverse transcriptase and random primers. Then, libraries were run on an Illumina HiSeq 2500 sequencing platform generating around 20 million paired end reads with an average read length of 282 base pairs.
---
