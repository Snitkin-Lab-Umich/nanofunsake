# nanofunsake
A quality control workflow implemented in Snakemake to perform hybrid assembly for Candida auris sequencing data. It is a modified version of nanosake (https://github.com/Snitkin-Lab-Umich/Nanosake).

## Summary

This pipeline is very similar to [funQCD](https://github.com/Snitkin-Lab-Umich/funQCD) and [nanofunQC](https://github.com/Snitkin-Lab-Umich/nanofunQC) in terms of installation, setup, and output. The main difference is that [BWA](https://github.com/lh3/bwa) and [Polypolish](https://github.com/rrwick/Polypolish) are used to polish the long-read assembly produced by Medaka. This pipeline assumes that you've already run short-read and long-read assemblies on your data, and that both versions passed their respective QC checks.


## Installation

> Clone the github directory onto your system.

```

git clone https://github.com/Snitkin-Lab-Umich/nanofunsake

```

> Load Bioinformatics, snakemake, and singularity modules.

```

module load Bioinformatics snakemake singularity

```

## Setup

Similar to funQCD, you'll first need to run setup.py with the path to your raw reads:

```

python setup.py [path_to_reads]

```

Next, you'll need to update config/config.yaml with a new prefix, a new long_reads path, and a new short_reads path. The _reads fields should be the paths to the raw data used for the original assemblies, and the prefix field should be a unique name for the current run. Further options can be adjusted in the profile/config.yaml file.

## Running nanofunsake

> First, perform a dry run of nanofunsake by running the following command. This will show the steps and commands that the pipeline will execute in the real run, without actually executing anything. Remove the --quiet flag for a more detailed view.

```

snakemake -s workflow/nfqc.smk -p --configfile config/config.yaml --profile ./profile/ -n --quiet

```

> The snakemake options present in profile/config.yaml should be visible in the detailed dry run (such as memory and runtime for each rule). By default, --slurm is enabled in these options, and snakemake will submit jobs to the cluster using the account in your profile. If everything looks correct, start the run using a job script with minimal CPUs, moderate memory, and a long runtime. An example job script is provided in `run_nfs.job`.

```
#!/bin/bash

#SBATCH --mail-user=[your_email]
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=[your_account]
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=20:00:00

module load Bioinformatics snakemake singularity
snakemake -s workflow/nanofunsake.smk -p --configfile config/config.yaml --profile ./profile/

```



