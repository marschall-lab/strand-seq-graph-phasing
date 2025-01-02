[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13356328.svg)](https://doi.org/10.5281/zenodo.13356328)

# Graphasing: Phasing Diploid Genome Assembly Graphs with Single-Cell Strand Sequencing

## Getting Started

1. Clone the project directory.

```bash
git clone https://github.com/marschall-lab/strand-seq-graph-phasing.git
```
2. Create the working directory. This is where the intermediate and final outputs will be generated.

```bash
mkdir wd/
```

3. Create the conda environment that will be used to run snakemake. Keep this environment in the project directory. 

```bash
conda env create -p snakemake_runner/ -f strand-seq-graph-phasing/workflow/envs/env_snakemake.yaml
```

4. Adjust config file and sample sheet. Template config and sample sheets are `config/config.yaml` and `config/samples.tsv` respectively

5. Before running Graphasing, activate the environment created in step 3. Once activated, you should be ready to go!

```bash
conda activate snakemake_runner/
cd strand-seq-graph-phasing/
```

### Config and Sample Sheet

#### Sample Sheet

The sample sheet is a .tsv file with the following columns:

`sample: str` Sample ID. This will serve as the identifier for the sample outputs.	

`strandseq_dir: str` Path to directory containing Strand-seq .fasta.gz/fastq.gz files.	

`gfa: str` Path to unphased input assembly graph `.gfa` file

`coverage: str, NA` Path to coverage file used for rukki path calculations. The path coverage file is a tsv file with three columns 'node', 'coverage', and 'length'. For Verkko assemblies, this is the `assembly.hifi-coverage.csv` file. For hifiasm use "NA", and coverage will be extracted from the input `p_utg.gfa`.

`hpc: [TRUE, FALSE]` Is the unphased input assembly graph `gfa` homopolymer compressed?

`assembler: [verkko, hifiasm]` Which genome assembly program was used to create the input assembly graph. 

#### Config settings

`samples: str` Path to sample sheet. 

`scripts_dir: str` Path to folder containing Graphasing R scripts. This is the folder `scripts` in the project folder.

`reference: str (optional, default:None)` Path to reference genome. If provided, the input assembly will be aligned to the reference using `minimap2`, and used to supplement the output summaries and evaluation plots.

`segmentLengthThreshold: int` Filtration parameter. Unitigs in the input assembly graph `.gfa` less than `segmentLengthThreshold` in basepairs will be filtered out.

`tangleSegmentLengthThreshold: int (optional, default: 250000)` Filtration parameter. Unitigs in the input assembly graph `.gfa` less than `tangleSegmentLengthThreshold` in basepairs will be considered during large tangle removal.

`per_thread_memory: bool (optional, default:False)` If `True`, memory specified for each rule is divided by the number of threads.

##### Experimental Features: Settings

##### Input Switch Error Detection
`calc_breakpoints bool (optional, default:False)` If `True`, Graphasing will attempt to locate haplotype switches in the input assembly. See below for more detail.

## Example command

```bash
snakemake \\
-d ../wd/ \\
 --configfiles config/config.yaml \\
 --config samples=/path/to/samples.tsv  \\
 --restart-times 3
```
#### Explanation

`-d` Path to working directory

`--configfiles` Path to config .yaml

`--config samples=` Path to sample sheet. Overrides the config file.

`--restart-times` While the resource allocation for each rule generally works on the first attempt, some rules occasionally need more memory. 3 restarts generally works to provide enough memory for the workflow to finish.


## Generating Phased Assemblies.

### Verkko
Note that verkko, being implemented in snakemake, is very sensitive to timestamps, an issue which can lead to the phased paths input being ignored. See [this issue] (https://github.com/marbl/verkko/issues/302#issuecomment-2495672213) for guidance on how to maintain proper timestamps while generating phased assemblies.

To generate phased assemblies with [Verkko](https://github.com/marbl/verkko), the output .gaf file, found at `rukki/{sample}/{sample}_rukki_paths.gaf` can be passed to the `--paths` argument. At the same time, the folder containing the unphased assembly, as well as paths to the hifi and ONT reads used to constuct the unphased assembly, need to be input as well, eg:

```bash
verkko
# Other Arguments
...
# Same as for unphased assembly
--assembly {unphased_assembly_folder}
--hifi {hifi_reads}
--nano {nano_reads}
# Pass the paths .gaf here:
--paths {Graphasing_paths_gaf}
```

### hifiasm
To generate phased assemblies with [hifiasm](https://github.com/chhylp123/hifiasm), the output Yak databases files, found at `yak/{sample}_hap[12].yak` can be passed to the `-1` and `-2` arguments of hifiasm trio mode assembly. Read overlaps from the unphased assembly can be reused by passing the same output prefix as the unphased assembly, eg:

```bash
hifiasm -o {unphased_output_prefix} ... -1 {Graphasing_yak_database_1} -2 {Graphasing_yak_database_2}
```

See the [trio-phasing documentation](https://hifiasm.readthedocs.io/en/latest/trio-assembly.html) for more.

## Visualization with Bandage

Graphasing outputs [bandage-friendly](https://rrwick.github.io/Bandage/) haplotype marker summaries, which can be found in `haplotype_marker_counts/`. Each chromosome cluster is assigned two opposing colors, one for each haplotype, and each unitig is shaded according to the haplotype markers that align to it, with grey representing a 50/50 balance of markers. 

## Warnings

### Haploid Chromosomes may be assigned to the same haplotype.

Strand-seq alignments to haploid chromosomes have been observed to occasionally behave strangely, especially when they are already full phased in the input assembly .gfa file. Normally unitigs from haploid chromosomes display a haploid alignment signal and cluster together. However, it has been observed that haploid unitigs can display a diploid alignment signal. When this happens, it is possible for all haploid chromosomes, eg X and Y, to be assigned to the same haplotype.


## Experimental Features

### Input Switch Error Detection

When specified in the configuration file, Graphasing will attempt to identify haplotype switches within input unitigs. Roughly speaking, the switch error detection process convolves a step filter with the Strand-seq signal from Watson-Crick libraries, looking for signs of a sudden change, or "breakpoint peaks", in haplotype signal on a unitig. Breakpoints can only be calculated for unitigs that were assigned to a cluster during the phasing process, as knowledge of the inherited Strand-seq library state is necessary for the breakpoint calculation. Output can be found in the `breakpoints/` folder within the working directory. The output currently favors higher recall. More specific breakpoint detection can be achieved by filtering the output for homozygous state switches and big peaks. Alternatively, the breakpoint plots can be inspected to filter out spurious calls.

### Manual Adjustments
(more documentation to come)
It is possible to manually intercede and adjust results as multiple points within the phasing process. Though not strictly necesssary, better understanding of this process can be achieved by reading the Graphasing [methods section](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03409-1#Sec17)

After running Graphasing, within the `intermediate_output/{sample}/` directory there will be a collection of output files reflecting the state of the phasing process at a particular point. The files are listed here in order of when the files are computed:

* included_libraries.tsv
  * What Strand-seq libraries are used for the phasing calculations. This is after QC is run to filter out loq-quality libraries. Modifying this list will require use of the files in the `library_ids/` directory.
* prior_clusters.tsv
  * Unitig clustering by chromosome before the clustering algorithm is run. It is expected for all cluster assignments to be `NA`
* refined_clusters.tsv
  * Unitig clustering by chromosome before final cluster refinement
* final_clusters.tsv
  * Final unitig clustering by chromosome
* unitig_orientation.tsv
  * The detetcted orientation of each unitig. Note that two unitigs appearing to be of differing orientation may reflect the flexiblity with which the same connection between two nodes can be specified, rather than a genuine misoriented unitig.
* counting_methods.tsv
  * What method is used to calculate the phase vector for each chromosome cluster of unitigs.

To pass one of these files as input to Graphasing, create a directory named `phasing_input/` in the working directory and copy `intermediate_output/{sample}/` to `phasing_input/{sample}/`. Any files in this folder will be used in place of executing the corresponding calculations of the graphasing workflow during a subsequent run. One can manually edit these files to affect downstream calcuation. For example, manually editing  `final_clusters.tsv`, which is then used as input for unitig orientation detection, counting method decision, and ultimately the haplotype marker count output. Note that including downstream files may nullify the effect of upstream manual edits. For example, if both `refined_clusters.tsv` and `final_clusters.tsv` are included in `phasing_input/{sample}/`, then graphasing will use `final_clusters.tsv` as provided, and skip performing any calculations with  refined_clusters.tsv. 

Some plots are also included in `intermediate_output/{sample}/` that can help to evaluate the quality of the results (more documentation to come). The most useful (possibly) for the purpose of evaluation and making manual edits are the cosine similarity heatmaps. For more fragmented assembiles, the similarity barcode and heatmap plots will not output, as the plots become too large.

#### Manual Adjustments Summary

* run graphasing workflow
* copy `intermediate_output/{sample}/` to `phasing_input/{sample}/`
* Using plots for assistance, make edits to corresponding files. Delete files downstream of the edits to ensure edited files used.
* rerun graphasing
* double check current files in `intermediate_output/{sample}/` to ensure edits were used as expected. Repeat editing process as necessary