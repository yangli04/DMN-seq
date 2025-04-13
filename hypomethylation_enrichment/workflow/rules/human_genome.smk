rule mapping_human_genome_PE:
    input:
        r1="cutadapt_PE/trimmed_{sample}_R1.fq.gz",
        r2="cutadapt_PE/trimmed_{sample}_R2.fq.gz",
    output:
        bam="mapping_hg_genome/{sample}.hg_genome.bam",
        summary="mapping_hg_genome/summary/{sample}.summary",
    params:
        ref_hg_genome=config["ref_hg_genome"],
        tmp="mapping_hg_genome/tmp",
    threads: 12
    resources:
        mem_mb=6000,
    shell:
        """
        hisat2 -p {threads}  \
            -x {params.ref_hg_genome} -1 {input.r1} -2 {input.r2}\
            -t --no-spliced-alignment --maxins 5000 --minins 0 --no-unal \
            --summary-file {output.summary} --new-summary \
            | samtools view -@ {threads} -Shub - \
            | samtools sort -T {params.tmp}/ -@ {threads} \
            -o {output.bam}
        """


rule flag_sort_index_depth_hg_genome:
    input:
        bam="mapping_hg_genome/{sample}.hg_genome.bam",
    output:
        flagstat="mapping_hg_genome/flagstat/{sample}.hg_genome.flagstat",
        index="mapping_hg_genome/{sample}.hg_genome.bam.bai",
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.flagstat}
        samtools index {input}
        """


rule report_mapping_hg_genome:
    input:
        expand("mapping_hg_genome/summary/{sample}.summary", sample=SAMPLE),
        expand("mapping_hg_genome/flagstat/{sample}.hg_genome.flagstat", sample=SAMPLE),
    output:
        "report_qc/report_mapping_hg_genome.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -n {output} {input}"


rule dedup_hg_genome:
    input:
        bam="mapping_hg_genome/{sample}.hg_genome.bam",
        index="mapping_hg_genome/{sample}.hg_genome.bam.bai",
    output:
        dedup_bam="mapping_hg_genome/{sample}.hg_genome.dedup.bam",
        dedup_metrics="mapping_hg_genome/{sample}.hg_genome.dedup.metrics",
    params:
        picard=config["picard"],
        tmp="mapping_hg_genome/tmp",
    threads: 1
    resources:
        mem_mb=96000,
    shell:
        """
        java -Xmx96g -jar {params.picard} MarkDuplicates \
            I={input.bam} O={output.dedup_bam} M={output.dedup_metrics} \
            REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
            TMP_DIR={params.tmp}
        """
rule index_deduped_genome:
    input:
        dedup_bam="mapping_hg_genome/{sample}.hg_genome.dedup.bam",
    output:
        dedup_index="mapping_hg_genome/{sample}.hg_genome.dedup.bam.bai",
    shell:
        "samtools index {input.dedup_bam}"


rule featureCounts_hg_tss_2000bp:
    input:
        dedup_bam=expand(
            "mapping_hg_genome/{sample}.hg_genome.dedup.bam", sample=SAMPLE
        ),
    output:
        counts="feature_counts/hg38_tss_2000bp/counts_hg38_tss_2000bp.txt",
        summary="feature_counts/hg38_tss_2000bp/counts_hg38_tss_2000bp.txt.summary",
    params:
        ref_hg38_tss_2000bp=config["ref_hg38_tss_2000bp"],
        featureCounts=config["featureCounts"],
    threads: 12
    resources:
        mem_mb=4000,
    shell:
        """
        {params.featureCounts} -T {threads} -F SAF --countReadPairs -p -a {params.ref_hg38_tss_2000bp} -o {output.counts} {input}
        """


rule report_dedup_hg_genome:
    input:
        expand("mapping_hg_genome/{sample}.hg_genome.dedup.metrics", sample=SAMPLE),
    output:
        "report_qc/report_dedup_hg_genome.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -n {output} {input}"


rule report_featureCounts_hg_genome:
    input:
        expand("feature_counts/hg_genome/counts_hg_genome.txt.summary", sample=SAMPLE),
    output:
        "report_qc/report_featureCounts_hg_genome.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -n {output} {input}"

rule bam_coverage:
    input: bam="mapping_hg_genome/{sample}.hg_genome.dedup.bam",
        index="mapping_hg_genome/{sample}.hg_genome.dedup.bam.bai",
    output:
        coverage="deeptools/hg_genome/bam_coverage/{sample}.hg_genome.dedup.bw",
    shell:
        """
        bamCoverage -b {input.bam} -o {output.coverage} --ignoreDuplicates --normalizeUsing CPM -bs 20 --extendReads --effectiveGenomeSize 2913022398
        """

rule compute_matrix_scale_region_hg_genome:
    input:
        "deeptools/hg_genome/bam_coverage/{sample}.hg_genome.dedup.bw",
    output:
        matrix="deeptools/hg_genome/compute_matrix_scale_regions/{sample}.matrix_hg_genome.mat.gz",
    params:
        ref_hg_genome_gtf=config["ref_hg_genome_gtf"],
    shell:
        """
        computeMatrix scale-regions --missingDataAsZero --skipZeros --regionsFileName {params.ref_hg_genome_gtf} --scoreFileName {input} -o {output.matrix} --metagene --upstream 2000 --downstream 2000 --regionBodyLength 4000
        """

rule compute_matrix_reference_point_hg_genome:
    input:
        "deeptools/hg_genome/bam_coverage/{sample}.hg_genome.dedup.bw",
    output:
        matrix="deeptools/hg_genome/compute_matrix_reference_point/{sample}.matrix_hg_genome.mat.gz",
    params:
        ref_hg_genome_gtf=config["ref_hg_genome_gtf"],
    shell:
        """
        computeMatrix reference-point --regionsFileName {params.ref_hg_genome_gtf} --scoreFileName {input} -o {output.matrix} --referencePoint TSS --upstream 2000 --downstream 2000 --missingDataAsZero --skipZeros
        """


rule plot_heatmap_reference_point_hg_genome:
    input:
        matrix="deeptools/hg_genome/compute_matrix_reference_point/{sample}.matrix_hg_genome.mat.gz",
    output:
        heatmap="deeptools/hg_genome/plot_heatmap/{sample}.reference_point_heatmap_hg_genome.pdf",
    shell:
        """
        plotHeatmap -m {input.matrix} -o {output.heatmap}
        """

rule plot_heatmap_scale_region_hg_genome:
    input:
        matrix="deeptools/hg_genome/compute_matrix_scale_regions/{sample}.matrix_hg_genome.mat.gz",
    output:
        heatmap="deeptools/hg_genome/plot_heatmap/{sample}.scale_region_heatmap_hg_genome.pdf",
    shell:
        """
        plotHeatmap -m {input.matrix} -o {output.heatmap}
        """

