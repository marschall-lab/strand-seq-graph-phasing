

rule count_raw_sseq_alignments:
    input:
        sam=lambda wildcards: expand('temp/sams/{{sample}}/{lib}.mdup.filt.sam', lib=MAP_SAMPLE_TO_INPUT[wildcards.sample]['strandseq_libs']),
        fastmap=sample2fastmap
    output:
        mem = "sseq_alignment_counts/{sample}_sseq_mem_raw.csv",
        fastmap = "sseq_alignment_counts/{sample}_sseq_fastmap_raw.csv"
    params:
        script = get_script_path('R', 'count_alignments.snakemake.R')
    conda: '../envs/env_Renv.yaml'
    threads: 1
    resources:
        mem_mb=calc_mem(32),
        walltime=calc_walltime(1, 2)
    log: "log/count_raw_sseq_alignments_{sample}.log"
    benchmark: "benchmark/count_raw_sseq_alignments_{sample}.benchmark"
    shell:
        '''
        (Rscript --vanilla {params.script} \\
        --mem-alignment-bams {input.sam} \\
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
