# strand-seq-graph-phasing

Getting the workflow working on HHU hilbert:

1. clone project

2. Move working directory, `wd`, with scripts to wherever you want it to be located.
Folder `scripts` has to be in the working directory

3. Copy hilbert profile into `wd`.

4. Copy `/gpfs/project/projects/medbioinf/projects/mihen108/haploclust_Renv2.sif` into `wd`.

5. Create env which runs snakemake.

6. load `Singularity`

7. Setup the `BubbleGun` conda environment:

   1. Run snakemake with `--conda-create-envs-only`
   2. `cd` to `wd/.snakemake/conda/`
   3. Figure out which conda environemnt corresponds to the 
      `BubbleGun` environement and `conda activate` that environment
   4. Use `pip` to install BubbleGun, according to [this link](https://wiki.hhu.de/display/HPC/Python),
   using the command:
```
PIP_CONFIG_FILE=/software/python/pip.conf pip install --user BubbleGun==1.1.3
```

8. Adjust config file. Reference config is `HG002tester`

Hopefully pipeline should work after this.