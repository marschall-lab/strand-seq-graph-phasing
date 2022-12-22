# strand-seq-graph-phasing

## Getting the workflow working on HHU Hilbert cluster:

There are two difficulties with setting up this workflow on Hilbert: The `R` and `BubbleGun` conda environments. 

Due to Hilbert's difficulties with R bioconductor packages (which attempt to connect to the internet during installation), the cluster cannot create the conda environment needed for `R` rules. Instead, a `.sif` [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) file, corresponding to the environment in `workflow/envs/env_Renv2.yaml`, has to be used. Peter has already created this file, but it still needs to be accounted for during setup.  

`BubbleGun` is a package that is only available from `pip`, and not from conda. Accordingly, the `BubbleGun` environment has to be set-up as a special case described below.


### Step-by-step how Mir gets the workflow working on HHU Hilbert

1. Clone project directory

2. Move working directory, `wd`, to wherever you want it to be located. Mir has been working with `wd` at the same level as the project directory. The `scripts` folder has to be in `wd`.

3. Copy Hilbert cluster profile into `wd`.

4. Copy `R` singularity environment, `/gpfs/project/projects/medbioinf/projects/mihen108/env_Renv2.sif`, into `wd`.

5. Create the conda environment that will be used to run snakemake. Keep this environment somewhere in the project directory, EG: 

```python
conda create -p strand-seq-graph-phasing/env/snakemake_runner/ snakemake
```

An environment file that contains the Python and Snakemake versions that Mir has been using to test can be found at `workflow/envs/env_snakemake.yaml`

6. load the `Singularity` module: `module load Singularity`

7. Setup the `BubbleGun` conda environment:

   1. Run snakemake with `--conda-create-envs-only`
   2. `cd` to `wd/.snakemake/conda/`
   3. To figure out which conda environment corresponds to the 
      `BubbleGun` environment, investigate the .yaml files and look for which one contains the `BubbleGun` package. `conda activate` that environment.
   4. Use `pip` to install BubbleGun, according to [this link](https://wiki.hhu.de/display/HPC/Python),
   using the command:
```
PIP_CONFIG_FILE=/software/python/pip.conf pip install --user BubbleGun==1.1.3
```

8. Adjust config file. Reference config is `HG002tester.yaml`

Hopefully pipeline should work after this.