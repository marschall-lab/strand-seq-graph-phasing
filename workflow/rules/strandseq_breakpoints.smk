

rule count_sseq_alignments_mem_unaggregated:
    input:
        aln=lambda wildcards: expand('temp/sams/{{sample}}/{lib}.mdup.filt.sam', lib=MAP_SAMPLE_TO_INPUT[wildcards.sample]['strandseq_libs']),
        unitig_lengths = 'temp/unitig_lengths/unitig_lengths_{sample}.tsv'
    output: "sseq_alignment_counts/{sample}_sseq_mem_unagg.csv",
    params:
        script = get_script_path('R', 'count_alignments_mem.snakemake.R')
    conda: '../envs/env_Renv.yaml'
    threads: 1
    resources:
        mem_mb=calc_mem(32),
        walltime=calc_walltime(1, 2)
    log: "log/count_sseq_alignments_mem_unaggregated_{sample}.log"
    benchmark: "benchmark/count_sseq_alignments_mem_unaggregated_{sample}.benchmark"
    shell:
        'Rscript --vanilla {params.script} '
        '--input {input.aln} '
        '--lengths {input.unitig_lengths} '
        '--threads {threads} '
        '--output {output} > {log} 2>&1'


#######################################
############ Breakpointing ############
#######################################

rule locate_sseq_breakpoints:
    input:
        lib_weights = "library_weights/{sample}_library_weights.csv",
        mem_alignments=rules.count_sseq_alignments_mem_unaggregated.output,
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
