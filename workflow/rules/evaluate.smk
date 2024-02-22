
if reference:
    ref_alns = expand("reference_alignments/{ref_name}/{sample}_{ref_name}_ref-aln.paf", sample=SAMPLES, ref_name=ref_name)
    out = 'plots/evaluation_plots_{ref_name}' + plot_output_suffix + '.pdf'
    log = 'log/evaluation_plots_{ref_name}' + plot_output_suffix + '.log'
else:
    ref_alns = [] # '.none.'
    out = 'plots/evaluation_plots_noref' + plot_output_suffix + '.pdf'
    log = 'log/evaluation_plots_noref' + plot_output_suffix + '.log'

rule evaluation_plots:
    input:
        hmc = expand("haplotype_marker_counts/{sample}_fudged_haplotype_marker_counts.csv", sample=SAMPLES),
        paths   = expand("rukki/{sample}/{sample}_rukki_paths.tsv", sample=SAMPLES),
        ref_alns = ref_alns
    output: out
    conda: '../envs/env_Renv.yaml'
    resources:
        mem_mb=calc_mem(16),
        walltime=calc_walltime(1, 0)
    params:
        script = get_script_path('R', 'evaluation_plots.snakemake.R'),
        dip_regex = dip_chroms_regex,
        hap_regex = hap_chroms_regex,
        hem_regex = hem_chroms_regex,
        hta_regex = hta_chroms_regex,
        plot_width = plot_width,
        plot_height = plot_height,
        slt = segment_length_threshold,
        samples = SAMPLES,
    log: log
    shell:
        '''
        (Rscript --vanilla {params.script} \\
        --included-diploid-chroms-regex "{params.dip_regex}" \\
        --included-haploid-chroms-regex "{params.hap_regex}" \\
        --included-hemiploid-chroms-regex "{params.hem_regex}" \\
        --hta-chroms-regex "{params.hta_regex}" \\
        --plot-width {params.plot_width} \\
        --plot-height {params.plot_height} \\
        --segment-length-threshold {params.slt} \\
        --samples {params.samples} \\
        --haplotype-marker-counts {input.hmc} \\
        --rukki-paths {input.paths} \\
        --ref-alignments {input.ref_alns} \\
        --output {output}) > {log} 2>&1
        '''
