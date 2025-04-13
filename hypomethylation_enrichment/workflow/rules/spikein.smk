# Map to spikein. Please change the parameters.
rule mapping_spikein_PE:
    input:
        r1="cutadapt_PE/trimmed_{sample}_R1.fq.gz",
        r2="cutadapt_PE/trimmed_{sample}_R2.fq.gz",
    output:
        bam="mapping_spikein/{sample}.spikein.bam",
        summary="mapping_spikein/summary/{sample}.summary",
    params:
        ref_spikein=config["ref_spikein"],
        tmp="mapping_spikein/tmp",
    threads: 12
    resources:
        mem_mb=6000,
    shell:
        """
        hisat2 -p {threads}  \
            -x {params.ref_spikein} -1 {input.r1} -2 {input.r2}\
            -t --no-spliced-alignment --maxins 500 --minins 0 --no-unal \
            --summary-file {output.summary} --new-summary \
            | samtools view -@ {threads} -Shub - \
            | samtools sort -T {params.tmp}/ -@ {threads} \
            -o {output.bam}
        """

rule featureCounts_spikein_all:
    input:
        bam=expand(
            "mapping_spikein/{sample}.spikein.bam", sample=SAMPLE
        ),
    output:
        counts="feature_counts/spikein/counts_spikein_all.txt",
        summary="feature_counts/spikein/counts_spikein_all.txt.summary",
    params:
        ref_spikein_all=config["ref_spikein_all"],
        featureCounts=config["featureCounts"],
    threads: 12
    resources:
        mem_mb=4000,
    shell:
        """
        {params.featureCounts} -T {threads} -F SAF --countReadPairs -p -a {params.ref_spikein_all} -o {output.counts} {input}
        """

rule featureCounts_spikein_lambda_window:
    input:
        bam=expand(
            "mapping_spikein/{sample}.spikein.bam", sample=SAMPLE
        ),
    output:
        counts="feature_counts/spikein/counts_spikein_lambda_window.txt",
        summary="feature_counts/spikein/counts_spikein_lambda_window.txt.summary",
    params:
        ref_spikein_lambda_window=config["ref_spikein_lambda_window"],
        featureCounts=config["featureCounts"],
    threads: 12
    resources:
        mem_mb=4000,
    shell:
        """
        {params.featureCounts} -T {threads} -F SAF --countReadPairs -p -a {params.ref_spikein_lambda_window} -o {output.counts} {input}
        """


# DO NOT DEDUP SPIKE-IN since the umi do not support a large coverage.
rule flag_sort_index_depth_spikein:
    input:
        "mapping_spikein/{sample}.spikein.bam",
    output:
        flagstat="mapping_spikein/flagstat/{sample}.spikein.flagstat",
        index="mapping_spikein/{sample}.spikein.bam.bai",
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.flagstat}
        samtools index {input}
        """


rule report_mapping_spikein:
    input:
        expand("mapping_spikein/summary/{sample}.summary", sample=SAMPLE),
        expand("mapping_spikein/flagstat/{sample}.spikein.flagstat", sample=SAMPLE),
    output:
        "report_qc/report_mapping_spikein.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -n {output} {input}"
