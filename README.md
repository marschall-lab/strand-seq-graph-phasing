# strand-seq-graph-phasing

## Getting the workflow working on HHU Hilbert cluster:

There are two difficulties with setting up this workflow on Hilbert: The `R` and `BubbleGun` conda environments. 

Due to Hilbert's difficulties with R bioconductor packages (which attempt to connect to the internet during installation), the cluster cannot create the conda environment needed for `R` rules. Instead, a `.sif` [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) file, corresponding to the environment in `workflow/envs/env_Renv2.yaml`, has to be used. Peter has already created this file, but it still needs to be accounted for during setup.  

`BubbleGun` is a package that is only available from `pip`, and not from conda. Accordingly, the `BubbleGun` environment has to be set-up as a special case described below.


### Step-by-step how Mir gets the workflow set up on HHU Hilbert

1. Clone project directory

2. Create the working directory. Mir has been working with the working directory at the same level as the project folder. Copy the Hilbert cluster profile into the working directory.

3. Copy the `R` singularity environment, `/gpfs/project/projects/medbioinf/projects/mihen108/env_Renv2.sif`, into the working directory.

4. Create the conda environment that will be used to run snakemake. The environment file Mir has been using can be found at `workflow/envs/env_snakemake.yaml`. Keep this environment somewhere in the project directory, EG: 

```bash
conda env create -p env/snakemake_runner/ -f workflow/envs/env_snakemake.yaml
```

Before running the pipeline, remember to activate this conda environment

5. Setup the `BubbleGun` conda environment:

   1. Run snakemake with `--conda-create-envs-only`
   2. `cd` to `wd/.snakemake/conda/`
   3. To figure out which conda environment corresponds to the 
      `BubbleGun` environment, investigate the .yaml files and look for which one contains the `BubbleGun` package. `conda activate` that environment.
   4. Use `pip` to install BubbleGun, according to [this link](https://wiki.hhu.de/display/HPC/Python),
   using the command:
```
PIP_CONFIG_FILE=/software/python/pip.conf pip install --user BubbleGun==1.1.3
```

6. Adjust config file and sample sheet. Reference config yaml is `config/config.yaml`. Reference sample sheet is `config/samples.tsv`

Hopefully pipeline should work after this.

## Running the pipeline on HHU Hilbert

__NOTE__: Peter has mentioned that, in general, `snakemake` appears to have trouble with relative paths. Therefore, all paths should be entered as absolute paths, which unfortunately can make commands quite long :confused:


### Config and Sample Sheet

#### Sample sheet settings

`sample: str` ID	

`gfa: str` Path to input `.gfa`

`strandseq_dir: str` Path to Strand-seq files	

`strandseq_suffix: str` Strand-seq file suffix. All Strand-seq file names are expected to be of the form `{lib}_{1,2}_{strandseq_suffix}`

#### Config settings

`samples: str` Path to sample sheet. This is an optional parameter, as the path to the sample sheet can instead be entered during the snakemake command, which would allow the same config file to be used for multiple sample sheets.

`segmentLengthThreshold: int` Filtration parameter. Unitigs in the input `.gfa` less than segmentLengthThreshold in basepairs will be filtered out.

`scripts_dir: str` Path to folder containing R and python scripts. This is the folder `scripts` in the project folder.

`reference: str, None` Path to reference alignment  

The last two parameters only apply to the rules which run R scripts:

`deploy_offline: bool` If True, use the singularity to run the R script rules

`singularity_Renv:` Path to the R singularty environment. Only needed if `deploy_offline: True`

#### Example Command: 

First, activate the snakemake running environment:

```bash
conda activate env/snakemake_runner/
```

Next, if the pipeline is being run with the R singularity container, load the `Singularity` module: 
```bash
module load Singularity
```

```bash
snakemake \\
-d ../wd/ \\ 
--configfiles config/config_HGSVC-renewal.yaml \\
--config samples=/gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing_testing/config/sample_sheet_NA24385_test.txt \\
--profile ../wd/prf_SSGP_Mir_MT/ \\
--use-singularity --singularity-args "-B /gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing_testing/scripts/:/gpfs/project/projects/medbioinf/projects/mihen108/strand-seq-graph-phasing_testing/scripts/" \\
--restart-times 2 
```
#### Explanation

`-d` Path to working directory

`--configfiles` Path to yaml

`--config samples=` Path to sample sheet. Overrides the config file.

`--profile` Path to Hilbert cluster profile

`--use-singularity --singularity-args "-B ..."` Used when using Singularity for the R environments. Unfortunately, when the singularity container can have trouble locating paths outside of the working directory, and thus may fail to locate the scripts folder. If that occurs, the path to the scripts folder needs to be explicitly bound in singularity: `-B /path/to/scripts/folder/:/path/to/scripts/folder/`.

`--restart-times` While the resource allocation for each rule generally works on the first attempt, some rules occasionally need more memory for. 2 restarts generally works to provide enough memory for the pipeline to finish.

## Warnings

### Results are stochastic!

One important step in the pipeline is to cluster unitigs from the `.gfa` into groups corresponding to the same chromosome. This step is achieved using a stochastic ensemble clustering algorithm from the `contiBAIT` package [(See figure 7.3)](https://open.library.ubc.ca/media/stream/pdf/24/1.0135595/1). The clustering is run with settings that will likely be consistent for large unitigs with high-quality Strand-seq signal. However, it is not unlikely that the clustering of smaller unitigs with less Strand-seq signal will vary if the pipeline is run multiple times with the same configuration.
