# strand-seq-graph-phasing

## Getting the workflow working on HHU Hilbert cluster:

Due to Hilbert's difficulties with R bioconductor packages (which attempt to connect to the internet during installation), the cluster cannot create the conda environment needed for `R` rules. Instead, a `.sif` [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) file, corresponding to the environment in `workflow/envs/env_Renv2.yaml`, has to be used.

### Step-by-step how Mir gets the workflow set up on HHU Hilbert

1. Clone project directory

2. Create the working directory. Mir has been working with the working directory at the same level as the project folder. Copy the Hilbert cluster profile into the working directory.

3. Copy the `R` singularity environment, `/gpfs/project/projects/medbioinf/projects/mihen108/env_Renv2.sif`, into the working directory.

4. Create the conda environment that will be used to run snakemake. The environment file Mir has been using can be found at `workflow/envs/env_snakemake.yaml`. Keep this environment somewhere in the project directory, EG: 

```bash
conda env create -p env/snakemake_runner/ -f workflow/envs/env_snakemake.yaml
```

Before running the pipeline, remember to activate this conda environment.

5. Adjust config file and sample sheet. Reference config yaml is `config/config.yaml`. Reference sample sheet is `config/samples-hgsvc-verkko14.tsv`

Hopefully the pipeline should work after this.

## Running the pipeline on HHU Hilbert

### Config and Sample Sheet
sample	strandseq_dir	gfa	coverage	hpc	assembler	expect_XY_separate

#### Sample sheet settings

`sample: str` Sample ID	

`strandseq_dir: str` Path to Strand-seq files	

`gfa: str` Path to input `.gfa`

`coverage: str, NA` Path to coverage file that will be injected into `gfa` for rukki path calculations. If no coverage file exists, use "NA"

`hpc: [TRUE, FALSE]` Is the `gfa` homopolymer compressed?

`assembler: [verkko, hifiasm]` which genome assembly program was used to create the gfa. 

`expect_XY_separate: [TRUE, FALSE]` Is the assembly graph such that it can be expected that the only diploid component connected to the X and Y chromosomes will be the PAR? From Mir's experience, this has been the case for verkko graphs, while for hifiasm graphs, most chromosomes have been joined together in a single component.

#### Config settings

`samples: str` Path to sample sheet. This is an optional parameter, as the path to the sample sheet can instead be entered during the snakemake command, which would allow the same config file to be used for multiple sample sheets.

`segmentLengthThreshold: int` Filtration parameter. Unitigs in the input `.gfa` less than `segmentLengthThreshold` in basepairs will be filtered out.

`scripts_dir: str` Path to folder containing R and python scripts. This is the folder `scripts` in the project folder.

`reference: str, None` Path to reference alignment. If provided, the assembly will be aligned to the reference using `minimap2`

The last two parameters only apply to the rules which run R scripts:

`deploy_offline: bool` If True, use the singularity to run the R script rules

`singularity_Renv:` Path to the R singularty environment. Only needed if `deploy_offline: True`

#### Running the pipeline on Hilbert: 

First, activate the snakemake running environment:

```bash
conda activate env/snakemake_runner/
```

Next, if the pipeline is being run with the R singularity container, load the `Singularity` module: 
```bash
module load Singularity
```

Example command:

```bash
snakemake \\
-d ../test_wd/ \\
 --configfiles config/config.yaml \\
 --config samples=/gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing/config/samples-hgsvc-verkko14.tsv  \\
 --profile ../test_wd/prf_SSGP_Mir_MT/ \\
 --use-singularity --singularity-args "-B /gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing/scripts/:/gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing/scripts/" \\
 --restart-times 3
```
#### Explanation

`-d` Path to working directory

`--configfiles` Path to config .yaml

`--config samples=` Path to sample sheet. Overrides the config file.

`--profile` Path to Hilbert cluster profile

`--use-singularity --singularity-args "-B ..."` Used when using Singularity for the R environments. Unfortunately, the singularity container can have trouble locating paths outside of the working directory, and thus may fail to locate the scripts folder. If that occurs, the path to the scripts folder needs to be explicitly bound in singularity: `-B /path/to/scripts/folder/:/path/to/scripts/folder/`.

`--restart-times` While the resource allocation for each rule generally works on the first attempt, some rules occasionally need more memory. 3 restarts generally works to provide enough memory for the pipeline to finish.

## Warnings

### Results are stochastic

One important step in the pipeline is to cluster unitigs from the `.gfa` into groups corresponding to the same chromosome. This step is achieved using a stochastic ensemble clustering algorithm from the `contiBAIT` package [(See figure 7.3)](https://open.library.ubc.ca/media/stream/pdf/24/1.0135595/1). Though `contiBAIT` parameters are set to attempt to limit variation in clustering, it is possible that the clustering of unitigs with less Strand-seq signal will vary even if the pipeline is run multiple times with the same configuration.

This can (rarely) result in **catastrophic misclustering**. One example Mir has encountered is the X/Y chromosome for sample NA18989. Over many tens of reruns, the output of the pipeline for NA18989 remains unchanged. However, on two different occasions, the X/Y chromosome ended up totally misclustered, with much of the X and Y chromosomes assigned to neither haplotype, or to the same haplotype. Another example is HG00171, is when large (several Mbp) centromeric unitigs from chromosome 4 were clustered with chromosome 12.

### X and Y may be assigned to the same haplotype.

Strand-seq alignments to the X and Y chromosomes can behave strangely when they are already full phased in the GFA; while normally unitigs from the X/Y chromosomes would display a "haploid" alignment signal and cluster together, it is not uncommon for fully phased X and Y chromosomes to display a "diploid" alignment signal and cluster separately. Mir has observed both Y ~ haploid and X ~ diploid separate clustering and X and Y both ~ diploid separate clustering. When this happens, it is possible for both the  X and Y chromosomes to be assigned to the same haplotype.