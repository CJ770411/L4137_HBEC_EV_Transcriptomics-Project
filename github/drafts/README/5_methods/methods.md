## Scripts

This pipeline processes raw small-RNA sequencing data, performs miRNA differential expression and finally pathway enrichment analysis.

N.B. Pathway enrichment analysis was only performed on primary human bronchial epithelial cells (HBECs).

The following scripts are used for primary HBEC and Calu-3 analysis:

1. 01_raw_read_trimming.sh
2. 02_reads_qc.sh
3. 03_download_mirna_references.sh
4. 04_collapse_reads.sh
5. 04_collapse_reads.sh
6. 06_create_count_matrix.sh
7. 07_DE_analysis.sh

The following scripts are used for **only** primary HBEC analysis:

8. 08_FE_indirect.sh
9. 09_FE_direct.sh

## Methods

### 1. Sequencing read trimming

Raw sequencing data (FASTQ format files) is trimmed using Cutadapt to remove adapter contamination and isolate reads that are the approximate length of an miRNA molecule.

### 2. Quality control / RNA composition analysis

Quality of raw reads and trimmed reads (FASTQ format files) is determined using FastQC and the composition of RNA within the library is established using miRTrace. The outputs are compiled into a single HTML report using MultiQC.

### 3. Download miRNA reference files

Known precursor (hairpin) and mature *Homo sapiens* miRNA sequences are downloaded from miRBase for subsequent read mapping.

### 4. Read mapping and quantification

Trimmed reads FASTQ files are collapsed and then reads are mapped to known mature miRNA sequences and quantified. miRNA read counts are compiled into a TSV file to enable differential expression analysis using edgeR. 

### 5. Differential expression analysis

miRNA counts are used to establish differentially expressed (DE) miRNAs between untreated and treated conditions. This outputs a CSV file containing signifacntly DE miRNAs.

### 5. Pathway enrichment analysis (primary HBECs only)

Significant DE miRNAs are subjected to target-gene level and miRNA-level enrichment.

- Target-gene level - All data processing is performed within RStudio.
- miRNA-level - Lists of DE miRNAs (and background miRNAs) are generated in RStudio and outputted as a TXT file. The text was then copied and pasted into the relevant box on the miEAA 2.0 website to perform over-representation analysis and miRNA-set enrichment analysis.