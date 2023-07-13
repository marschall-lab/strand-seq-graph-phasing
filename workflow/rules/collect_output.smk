# From sample_table.smk
# MAP_SAMPLE_TO_INPUT
# SAMPLES

ALL_OUTPUT=None
def collect_output():

    # Mandatory Output: Haplotype marker counts
    all_output = [expand("haplotype_marker_counts/{sample}_haplotype_marker_counts.csv",sample=SAMPLES)]
    rukki_samples = set()
    for sample in SAMPLES:

        # Optional Output: Rukki paths
        coverage = MAP_SAMPLE_TO_INPUT[sample].get('coverage', None)
        if coverage is not None:
            rukki_samples.add(sample)

    all_output.append(expand("rukki/{sample}/{sample}_rukki_paths.gaf", sample=rukki_samples))

    # HiFiasm samples ~ yak databases
    hifi_samples = set()
    for sample in SAMPLES:
        assembler = MAP_SAMPLE_TO_INPUT[sample].get('assembler', None)
        if assembler == 'hifiasm':
            hifi_samples.add(sample)

    all_output.append(expand('yak/{sample}/{sample}_{hap}.yak',sample=hifi_samples, hap=['hap1', 'hap2']))

    global ALL_OUTPUT
    ALL_OUTPUT = all_output

collect_output()