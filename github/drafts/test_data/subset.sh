for i in $(seq 1 27); do
    seqtk sample -s$i SRR17084821.fastq.gz 1000 \
        > SRR17084821_subset_${i}.fastq
done

