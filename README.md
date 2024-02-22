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
conda env create -p strand-seq-graph-phasing/env/snakemake_runner/ -f strand-seq-graph-phasing/workflow/envs/env_snakemake.yaml
```

4. Adjust config file and sample sheet. Template config and sample sheets are `config/config.yaml` and `config/samples.tsv` respectively

5. Before running Graphasing, activate the environment created in step 3. Once activated, you should be ready to go!

```bash
cd strand-seq-graph-phasing/
conda activate env/snakemake_runner/
```

### Difficulties installing R Bioconductor Packages

If there are difficulties in installing R Bioconductor packages in your execution environment, a [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) container may be used. You can use the [definition file](https://docs.sylabs.io/guides/latest/user-guide/definition_files.html) provided at `workflow/envs/Renv2_singularity_definition.def` to [build the required Singularity container](https://docs.sylabs.io/guides/latest/user-guide/build_a_container.html). Once built, copy the `.sif` contained into the working directory created in step 2.

### Config and Sample Sheet

#### Sample Sheet

The sample sheet is a .tsv file with the following columns:

`sample: str` Sample ID. This will serve as the identifier for the sample outputs.	

`strandseq_dir: str` Path to directory containing Strand-seq .fasta.gz/fastq.gz files.	

`gfa: str` Path to unphased input assembly graph `.gfa` file

`coverage: str, NA` Path to coverage file used for rukki path calculations. The path coverage file is a tsv file with three columns 'node', 'coverage', and 'length'. For Verkko assemblies, this is the `assembly.hifi-coverage.csv` file. For hifiasm use "NA", and coverage will be extracted from the input `p_utg.gfa`.

`hpc: [TRUE, FALSE]` Is the unphased input assembly graph `gfa` homopolymer compressed?

`assembler: [verkko, hifiasm]` Which genome assembly program was used to create the input assembly graph. 

`cluster_PAR_with_haploid: [TRUE, FALSE]` If TRUE, Graphasing will attempt to merge the PAR with any detected haploid chromosomes on the same connected component. Is safest if it can be expected that the only diploid component connected to the X and Y chromosomes will be the PAR. 

#### Config settings

`samples: str` Path to sample sheet. 

`scripts_dir: str` Path to folder containing Graphasing R and python scripts. This is the folder `scripts` in the project folder.

`reference: str (optional, default:None)` Path to reference genome. If provided, the input assembly will be aligned to the reference using `minimap2`, and used to supplement the output summaries and evaluation plots.

`segmentLengthThreshold: int` Filtration parameter. Unitigs in the input assembly graph `.gfa` less than `segmentLengthThreshold` in basepairs will be filtered out.

`per_thread_memory: bool (optional, default:False)` If `True`, memory specified for each rule is divided by the number of threads. 

These two parameters only apply to the rules which use R Bioconductor packages:

`deploy_offline: bool` If True, use the singularity to run the rules using R Bioconductor packages

`singularity_Renv:` Path to the R Singularty environment. Only needed if `deploy_offline: True`

##### Experimental Features

##### Breakpoints
`calc_breakpoints bool (optional, default:False)` If `True`, Graphasing will attempt to locate haplotype switches in the input assembly. Detected peaks, as well as plots, can be found in `breakpoints/`

##### Evaluation Plots
`plot_evaluation bool (optional, default:True)` If `True`, evaluation plots will be made for all input samples. If a reference is provided, evaluation plots will include reference-derived evaluations. Without a reference, some plots will be left blank. Evaluation plots can be found in `plots/`. 

`plot_width, plot_height int (optional, default:20,15)` 
Width and height of evaluation plots .pdf output file.

`plot_output_suffix str (optional, default:'')`
Suffix to output file name. May be useful if using multiple configs within the same working directory

The following options are only needed if a reference is provided for the evaluation plots:

`dip_chroms_regex: str (optional, default: '(chr|CHR)[0-9]+')` 
\
`hap_chroms_regex: str (optional, default: '(chr|CHR)[yY]'` 
\
`hem_chroms_regex: str (optional, default: '(chr|CHR)[xX]'` 
\
`hta_chroms_regex: str (optional, default: '(chr|CHR)(13|14|15|21|22)'`

These regexs control what is plotted, as well as additional plot annotations. Annotations are made for chromosomes that are diploid, haploid, hemiploid (sometimes diploid or haploid eg, chr X in humans), and hard-to-assemble. Only chromosomes matching the diploid, haploid, and hemiploid regexs are plotted. Defaults are set for human genomes.

## Example command

```bash
snakemake \\
-d ../wd/ \\
 --configfiles config/config.yaml \\
 --config samples=/path/to/samples.tsv  \\
 --use-singularity --singularity-args "-B /path/to/scripts/folder/:/path/to/scripts/folder/" \\
 --restart-times 3
```
#### Explanation

`-d` Path to working directory

`--configfiles` Path to config .yaml

`--config samples=` Path to sample sheet. Overrides the config file.

`--use-singularity --singularity-args "-B ..."` Used when using Singularity for the R Bioconductor environment. Unfortunately, the Singularity container can have trouble locating paths outside of the working directory, and thus may fail to locate the scripts folder. If that occurs, the path to the scripts folder needs to be explicitly bound in singularity: `-B /path/to/scripts/folder/:/path/to/scripts/folder/`.

`--restart-times` While the resource allocation for each rule generally works on the first attempt, some rules occasionally need more memory. 3 restarts generally works to provide enough memory for the workflow to finish.


## Generating Phased Assemblies.

### Verkko 
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

Graphasing outputs [bandage-friendly](https://rrwick.github.io/Bandage/) haplotype marker summaries, which can be found in `haplotype_marker_counts/`. Each chromosome cluster is assigned two opposing colors, one for each haplotype, and each unitig is shaded according to the haplotype markers that align to it, with grey representing a 50/50 balance of markers. If a reference genome was provided, then summaries, containing additional reference-alignment-derived annotations can be found in `rmc/`.

## Warnings

### Results are stochastic

An important step in the workflow is to cluster unitigs from the input assembly `.gfa` into groups corresponding to the same chromosome. Part of this step is achieved using a stochastic ensemble clustering algorithm from the `contiBAIT` package [(See figure 7.3)](https://open.library.ubc.ca/media/stream/pdf/24/1.0135595/1). Though `contiBAIT` parameters are set to attempt to limit variation in clustering, it is possible that the clustering of unitigs, as well as all downstream haplotype inference, will vary even if the pipeline is run multiple times with the same configuration.

### Haploid Chromosomes may be assigned to the same haplotype.

Strand-seq alignments to haploid chromosomes have been observed to occasionally behave strangely, especially when they are already full phased in the input assembly .gfa; while normally unitigs from haploid chromosomes would each display a haploid alignment signal and cluster together, it has been observed that some or all haploid can display a diploid alignment signal and cluster separately. When this happens, it is possible for all haploid chromosomes, eg X and Y, to be assigned to the same haplotype.