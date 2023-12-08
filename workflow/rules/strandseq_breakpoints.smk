

rule count_raw_sseq_alignments:
    input:
        bam=lambda wildcards: sample2mem(wildcards, 'bam'),
        bai=lambda wildcards: sample2mem(wildcards, 'bam.bai'),
        fastmap=sample2fastmap,
    output:
        mem = "sseq_alignment_counts/{sample}_sseq_mem_raw.csv",
        fastmap = "sseq_alignment_counts/{sample}_sseq_fastmap_raw.csv"
    params:
        script = get_script_path('R', 'count_alignments.snakemake.R')
    singularity: singularity_r_env
    conda: '../' + conda_r_env
    threads: 1
    resources:
        mem_mb=calc_mem(32),
        walltime=calc_walltime(1, 2)
    log: "log/count_raw_sseq_alignments_{sample}.log"
    benchmark: "benchmark/count_raw_sseq_alignments_{sample}.benchmark"
    shell:
        '''
        (Rscript --vanilla {params.script} \\
        --mem-alignment-bams {input.bam} \\
        --fastmap-alignments {input.fastmap} \\
        --threads {threads} \\
        --aggregate-alignments FALSE \\
        --output-mem {output.mem} \\
        --output-fastmap {output.fastmap}) > {log} 2>&1
        '''

#######################################
############ Breakpointing ############
#######################################

rule locate_sseq_breakpoints:
    input:
        lib_weights = "library_weights/{sample}_library_weights.csv",
        mem_alignments="sseq_alignment_counts/{sample}_sseq_mem_raw.csv",
        hap_counts='haplotype_marker_counts/{sample}_fudged_haplotype_marker_counts.csv'
    output:
       "breakpoints/{sample}_breakpoint_windows.csv"
    params:
        script = get_script_path('R','locate_sseq_breakpoints.snakemake.R')
    conda: '../envs/env_Renv.yaml'
    threads: 1
    resources:
        mem_mb=calc_mem(24),
        walltime=calc_walltime()
    log: "log/locate_sseq_breakpoints_{sample}.log"
    benchmark: "benchmark/locate_sseq_breakpoints_{sample}.benchmark"
    shell:
        '''
        (Rscript --vanilla {params.script} \\
        --lib-weights {input.lib_weights} \\
        --haplotype-marker-counts {input.hap_counts} \\
        --sseq-alignments {input.mem_alignments} \\
        --output {output}) > {log} 2>&1
        '''
