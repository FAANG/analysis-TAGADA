# GENE-SWitCH project RNA-Seq analysis pipeline

## Requirements

The following dependencies are required:

- [Nextflow 19.10.0](https://www.nextflow.io/docs/latest/getstarted.html)
- [FastQC 0.11.8](https://github.com/s-andrews/FastQC)
- [Cutadapt 2.7](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [Trim Galore 0.6.5](https://github.com/FelixKrueger/TrimGalore)
- [STAR 2.7.3a](https://github.com/alexdobin/STAR)
- [StringTie 2.0.5](https://github.com/gpertea/stringtie)

## Usage

Execute this nextflow pipeline with:

```bash
./run rnaseq.nf --outdir <output directory> \
                --reads <reads1.fq> <reads2.fq> \
                --genome <genome.fa> \
                --annotation <annotation.gff> \
                [nextflow arguments]
```

Using the `./run` wrapper script instead of the `nextflow run` command provides several benefits:
- All results and temporary files are written to the output directory, keeping the execution directory clean.
- All temporary files are deleted once the pipeline has successfully completed.
- The pipeline can be resumed from any directory with nextflow's `-resume` option.
- The pipeline's arguments are parsed in accordance with the GNU conventions.

## About

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
