#!/usr/bin/env bash

rule adapterTrimming:
    input:  "%s/fastq_raw/{filename}.fastq" % outPath
    output: temp("%s/fastq_trimmed/{filename}.fastq" % outPath)
    log:    "%s/logs/fastq_trimmed/{filename}.log" % outPath
    conda:  "envs/cutadapt.yml"
    threads: config['threads']['high']
    shell:
        "cutadapt {config[adapter]} --minimum-length {config[readMinLength]} --quality-cutoff {config[qualityCutoff]} --discard-untrimmed --cores={threads} -o '{output}' '{input}' > {log} 2>&1"

                "cutadapt -a Realseq3P=TGGAATTCTC -u1 --minimum-length 17 --quality-cutoff 30 --discard-untrimmed --cores=6 -o '{output}' '{input}' > {log} 2>&1"

rule mappingBowtieGenome:
    input:  "%s/mapping_bowtie_spikeins/{filename}.unmapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_genome/{filename}.map" % outPath,
            mapped = temp("%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_genome/{filename}.unmapped.fasta" % outPath),
            mappedSeqSizeDist = "%s/mapping_bowtie_genome/{filename}.mapped.ssd" % outPath,
            unmappedSeqSizeDist = "%s/mapping_bowtie_genome/{filename}.unmapped.ssd" % outPath
    log:    "%s/logs/mapping_bowtie_genome/{filename}.log" % outPath
    conda:  "envs/bowtie1.yml"
    threads: config['threads']['medium']
    shell:
        """
            bowtie --threads {threads} -f -k1 -v2 --fullref --un '{output.unmapped}' --al '{output.mapped}' {config[repoPath]}/{config[genomeID]}/{config[genomeVersion]}/bowtiedb/genome '{input}' > '{output.map}' 2> {log}
            cat '{output.mapped}' | awk '{{if(NR%2==1) {{printf "%s\t",$0}} else {{printf "%i\\n",length($1)}} }}' | cut -d "x" -f2 | awk '{{i[$2]+=$1}} END{{for(x in i){{print i[x]" "x}}}}' | sort -k2 -n > {output.mappedSeqSizeDist} 2> {log}
            cat '{output.unmapped}' | awk '{{if(NR%2==1) {{printf "%s\t",$0}} else {{printf "%i\\n",length($1)}} }}' | cut -d "x" -f2 | awk '{{i[$2]+=$1}} END{{for(x in i){{print i[x]" "x}}}}' | sort -k2 -n > {output.unmappedSeqSizeDist} 2> {log}
        """


rule mappingBowtieMirna:
    input:  "%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath
    output: map = "%s/mapping_bowtie_mirna/{filename}.map" % outPath,
            # mapped = temp("%s/mapping_bowtie_mirna/{filename}.mapped.fasta" % outPath),
            unmapped = temp("%s/mapping_bowtie_mirna/{filename}.unmapped.fasta" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_bowtie_mirna/{filename}.log" % outPath
    conda:
        "envs/bowtie1.yml"
    shell:
        "bowtie --threads {threads} -f -k1 --fullref --best -v1 --un '{output.unmapped}' {config[repoPath]}/mirbase/{config[miRBaseVersion]}/hairpin/bowtiedb/hairpin-{config[speciesCode]} '{input}' > '{output.map}' 2> {log}"

rule miRDeepPrep:
    input: "%s/fastq_trimmed/{filename}.fastq" % outPath
    output:
        collapsedReads = temp("%s/mapping_mirdeep2_miRNA/{filename}.fasta" % outPath)
#        readsVsGenome = temp("%s/mapping_mirdeep2_miRNA/{filename}.arf" % outPath)
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_mirdeep2_miRNA/{filename}.log" % outPath
    conda:
        "envs/mirdeep2.yml"
    shell:
        """
            mapper.pl '{input}' -e -h -j -l %s -o {threads} -m -s '{output.collapsedReads}' > '{log}' 2>&1
            rm -fr mapper_logs
            rm -fr bowtie.log
        """ % config['readMinLength']


rule mappingMiRDeep2MiRNA:
    input:
        collapsedReads = "%s/mapping_bowtie_genome/{filename}.mapped.fasta" % outPath
        # readsVsGenome = "%s/miRDeep2/genomemapped/{filename}.arf"
    output:
        map = "%s/mapping_mirdeep2_miRNA/{filename}/miRNAs_expressed_all_samples_default.csv" % outPath
    params:
        outputDir = "%s/mapping_mirdeep2_miRNA/{filename}" % outPath
    threads: config['threads']['medium']
    log:    "%s/logs/mapping_mirdeep2_miRNA/{filename}.log" % outPath
    conda:
        "envs/mirdeep2.yml"
    shell:
        """
        subDirs=$(awk -F"/" '{{print NF-1}}' <<< "{input.collapsedReads}")
        subDirPath=$(seq ${{subDirs}} | awk '{{printf "../"}}')
        mkdir -p '{params.outputDir}'
        cd '{params.outputDir}'
        quantifier.pl -p '{workflow.basedir}/{config[repoPath]}/mirbase/{config[miRBaseVersion]}/hairpin/uncompressed/hairpin-{config[speciesCode]}.dna.fa' \
            -y 'default' \
            -T {threads} \
            -d \
            -m '{workflow.basedir}/{config[repoPath]}/mirbase/{config[miRBaseVersion]}/mature/uncompressed/mature-{config[speciesCode]}.dna.fa' \
            -r "../${{subDirPath}}{input.collapsedReads}" > "../${{subDirPath}}{log}" 2>&1
        rm -fr expression_analyses/expression_analyses_default/*.ebwt
        """


rule filterCollapsedReads:
    input: "%s/fastq_collapsed/{filename}.fastq" % outPath
    output:temp("%s/fastq_collapsed_filtered/{filename}.fastq" % outPath)
    threads: 1
    run:
        shell("awk 'BEGIN {{FS = \"_x\" ; OFS = \"\\n\"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if ($2 >= {config[minReadCount]}) {{print header, seq, qheader, qseq}}}}' < '{input}' > '{output}'")

