################################
############ Config ############
################################
# TODO move this to main snakefile?
ref_name = pathlib.Path(reference).stem

# ##########################################################################################################################
# #########     Using reference genome only for evaluation- mapping/haplotagging overlap graph unitigs ##############
# ##########################################################################################################################

ref_out = 'reference/'+ref_name+'.homopolymer-compressed.fasta'
# print(ref_out)

rule homopolymer_compress_ref:
    input: reference
    output: ref_out
    params:
        script=get_script_path('python','homopolymer_compress_fasta.py')
    conda:'../envs/env_cl.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * 16 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    benchmark: "benchmark/homopolymer_compress_ref.benchmark"
    shell:
        '''
        seqtk hpc {input} > {output}
        '''

# TODO explore -H homopolymer compression option. Do both reads and ref need to be non-compressed? Does it not matter?
rule map_unitigs_to_ref:
    input:
        ref=lambda wildcards: ref_out if MAP_SAMPLE_TO_INPUT[wildcards.sample]['hpc'] else reference,
        unitigs="fasta/{sample}/{sample}_unitigs-hpc.fasta"
    output: "reference_alignments/{ref_name}/{sample}_{ref_name}_ref-aln.paf"
    conda: '../envs/env_cl.yaml'
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * 96 * attempt,
        walltime=lambda wildcards, attempt: f'{8 + attempt * attempt:02}:59:00'
    log: "log/map_unitigs_to_ref_{sample}_{ref_name}.log"
    threads: 12
    shell:
        "(time minimap2 -t {threads} --eqx -x asm20 {input.ref} {input.unitigs} > {output}) > {log} 2>&1"


rule bandage_annotations_with_reference:
    input:
        hmc="haplotype_marker_counts/{sample}_fudged_haplotype_marker_counts.csv",
        ref="reference_alignments/{ref_name}/{sample}_{ref_name}_ref-aln.paf",
    output:
        'rmc/{sample}_{ref_name}_reference_marker_counts.csv'
    conda: '../envs/env_Renv.yaml'
    resources:
        mem_mb=calc_mem(16),
        walltime=calc_walltime(1,0)
    params:
        script = get_script_path('R','make_bandage_annotations.snakemake.R'),
    log: 'log/bandage_annotations_{sample}_{ref_name}.log'
    shell:
        '''
        (Rscript --vanilla {params.script} \\
        --haplotype-marker-counts {input.hmc} \\
        --ref-alignment {input.ref} \\
        --output {output}) > {log} 2>&1
        '''
