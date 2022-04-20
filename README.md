# Standard RNA-seq workflow

Nextflow workflow to process single-end or paired-end raw RNA-seq data. It takes compressed fastq files as input and generates normalized signal tracks and a table of read counts. Also, it saves log files that can be parsed to calculate the number of reads filtered out in each step.

## Authors

Yuvia Alhelí PÉREZ-RICO

Léna CLERQUIN

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that we have used:

- nextflow (20.04.1)
- FastQC (v0.11.8)
- Trim Galore (0.6.3)
- STAR (2.7.2b)
- samtools (1.9)
- Picard Tools (2.20.8)
- featureCounts (v2.0.1)
- R (3.2.2)
- deepTools (3.1.3)

Note: I suggest to install all programs in a conda environment.

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script, the configuration file and the additional script within the bin folder, you will need to prepare the following files and change the paths in the configuration file accordingly:

- A STAR index for mapping.
- Gene annotations in GTF format.
- A comma-separated file indicating the sample name and paths to the fastq files (examples for single-end and paired-end libraries are available in the 'docs' folder). The name of the fastq files is relevant for parsing in the mapping process.

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/RNA-seq`

`nextflow run standard_RNA-seq.nf`

## Acknowledgements

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771.

## License

Licenced under the GNU general public license.

