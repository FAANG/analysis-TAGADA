# GENE-SWitCH project RNA-Seq analysis pipeline

## Requirements

The following dependencies are required:

- Python with [Pandas](https://pandas.pydata.org/docs/getting_started/install.html#installing-using-your-linux-distribution-s-package-manager)
- [Nextflow 19.10.0](https://www.nextflow.io/docs/latest/getstarted.html)
- [FastQC 0.11.8](https://github.com/s-andrews/FastQC)
- [Cutadapt 2.7](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [Trim Galore 0.6.5](https://github.com/FelixKrueger/TrimGalore)
- [STAR 2.7.3a](https://github.com/alexdobin/STAR)
- [StringTie 2.0.5](https://github.com/gpertea/stringtie)

## Usage

Execute this nextflow pipeline with:

```bash
./run rnaseq.nf --output <output directory> \
                --paired-reads <paired reads 1.fq> <paired reads 2.fq> ... \
                --single-reads <single reads.fq> ... \
                --mapped-reads <mapped reads.bam> ... \
                --reference <genome reference.fa> \
                --index <genome index input directory> \
                --annotation <genome annotation.gff> \
                --threads <number of threads to use> \
                [nextflow arguments]
```

The pipeline has the following requirements:
- Both `--output` and `--annotation` must be provided.
- At least one of `--paired-reads`, `--single-reads`, `--mapped-reads` must be provided.
- One of `--index` or `--reference` must be provided when `--paired-reads` or `--single-reads` are provided.

The pipeline creates the following directories in the output directory:
- `logs` contains a history file listing user commands, and log files from dependencies.
- `genome` contains the computed indexed genome, if a genome reference and unmapped reads were provided.
- `reads` contains computed trimmed reads and mapped reads, and quality control files.
- `counts` contains all quantified genes and transcripts in tab separated matrix files, and annotation files.
- `temp` contains temporary files that can be safely deleted.

For the pipeline to correctly generate output file names, use the following input file names:
- `name[_R{1,2}].{fq,fastq}[.gz]` for unmapped reads.
- `name.{sam,bam}` for mapped reads.

Using the `./run` wrapper script instead of the `nextflow run` command grants several benefits:
- All results and temporary files are written to the output directory, keeping the execution directory clean.
- All temporary files are deleted once the pipeline has successfully completed.
- The pipeline can be resumed from any directory with nextflow's `-resume` option.
- The pipeline's arguments follow GNU conventions.

## About

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
