#!/usr/bin/env bash

###==== Cutadapt ====###

# Run cutadapt to trim raw read data for each sample
for sample in "${RAW_READS[@]}"; do

    SAMPLE_NAME=$(basename "$sample" .fastq.gz)

    cutadapt \
        -j "$SLURM_CPUS_PER_TASK" \
        -a Realseq3P=TGGAATTCTC \
        -u1 \
        --discard-untrimmed \
        -m "$MIN_LENGTH" \
        -q 30 \
        -o "${OUTDIR_CUTADAPT}/${SAMPLE_NAME}_trimmed.fastq.gz" \
        "$sample"

done



###==== FastQC ====###

# Run FastQC to generate QC reports
fastqc \
"${TRIMMED_READS[@]}" \
-o "$OUTDIR_FASTQC" \
-t $SLURM_CPUS_PER_TASK




###==== miRTrace ====###

# Run miRTrace to generate trace QC reports
mirtrace \
trace \
-o "$OUTDIR_MIRTRACE_TRACE" \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" 

# Run miRTrace to generate QC reports
mirtrace \
qc \
-s hsa \
-o "$OUTDIR_MIRTRACE_QC" \
-f \
-t $SLURM_CPUS_PER_TASK \
"${TRIMMED_READS[@]}" 




###==== MultiQC ====###

# Run MultiQC to merge all individual sample QC reports into one report
multiqc \
-o "$OUTDIR_MULTIQC" \
"$INDIR_SAMPLE_REPORTS"




###==== Download references ====###

wget https://www.mirbase.org/download/${mirna_form}.fa




###==== miRDeep2 (map reads) ====###


mapper.pl "$TEMP_TRIMMED_FASTQ" -e -h -j -l "$MIN_LENGTH" -o $SLURM_CPUS_PER_TASK -m -s "$COLLAPSED_READS" 






###==== miRDeep2 (quantification) ====###

# Quantify read counts to produce count matrix
quantifier.pl \
-T $SLURM_CPUS_PER_TASK \
-p "$INFILE_HAIR_MIRNA" \
-m "$INFILE_MAT_MIRNA" \
-r "$INFILE_COLLAPSED_READS" \
-t hsa

#==============================================================================#
# End of script
#==============================================================================#
