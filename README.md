# GENE-SWitCH project RNA-Seq analysis pipeline

## Requirements

The following dependencies are required:

- [Nextflow 20.01.0](https://www.nextflow.io/docs/latest/getstarted.html)
- [FastQC 0.11.9](https://github.com/s-andrews/FastQC)
- [Cutadapt 2.9](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [Trim Galore 0.6.5](https://github.com/FelixKrueger/TrimGalore)
- [STAR 2.7.3a](https://github.com/alexdobin/STAR)
- [Samtools 1.10](https://github.com/samtools/samtools)
- [StringTie 2.1.1](https://github.com/gpertea/stringtie)
- Python with [Pandas](https://pandas.pydata.org/docs/getting_started/install.html#installing-using-your-linux-distribution-s-package-manager)


## Usage

Execute this nextflow pipeline with:

    ./run rnaseq.nf [arguments]

The pipeline accepts the following arguments:

    --output <directory>            Output directory                  Required
    --reads <'*.{fq,bam}'>          Input reads glob pattern          Required
    --genome <genome.fa>            Input genome sequence file        Required if raw --reads and no --index
    --index <directory>             Input genome index directory      Required if raw --reads and no --genome
    --annotation <annotation.gff>   Input reference annotation file   Required
    --metadata <metadata.tsv>       Input metadata file               Required if --merge
    --merge <factor1,factor2>       Factor(s) to merge mapped reads   Optional
    --direction <rf|fr>             Direction of reads                Optional
    --threads <number>              Number of threads                 Optional
    -resume                         Resume pipeline                   Optional

For the pipeline to correctly generate output file names, use the following input file names:
- `prefix[_R{1,2}].{fq,fastq}[.gz]` for `fastq` raw reads.
- `prefix.bam` for `bam` mapped reads.

The pipeline executes the following processes:
1. __Control__ raw reads quality with [FastQC](https://github.com/s-andrews/FastQC).  
   Outputs quality reports to `output/quality/raw`.
2. __Trim__ adaptators from reads with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).  
   Outputs quality reports to `output/quality/trimmed`.
3. __Index__ genome sequence wih [STAR](https://github.com/alexdobin/STAR).  
   Outputs indexed genome to `output/index`.
4. __Map__ reads to indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   Outputs mapped reads to `output/maps`.
5. __Merge__ mapped reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merging factors](#merging-factors) section for details.
6. __Assemble__ transcripts and __combine__ them into a new assembly annotation with [StringTie](https://github.com/gpertea/stringtie).  
   Outputs the new assembly annotation to `output/annotation`.
7. __Count__ genes and transcripts with [StringTie](https://github.com/gpertea/stringtie), and __format__ them into tabulated files.  
   Outputs TPM counts and average per-base read coverage to `output/counts`.  
   Counts are given for the reference and assembly annotations separately.


## Merging factors

Use the `--merge` and `--metadata` options together to merge multiple mapped reads.  
This results in genes and transcripts being counted by __factors__ rather than by __inputs__.

#### Examples

Given the following tabulated metadata file:

    prefix    diet     tissue
    A         corn     liver
    B         corn     liver
    C         wheat    liver
    D         wheat    muscle

With the following options:

    --reads '{A.fq,B.fq,C.fq,D.bam}' --metadata metadata.tsv --merge diet

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet.
- __C__ and __D__ mapped reads will be merged, resulting in gene and transcript counts for the __wheat__ diet.


With the following options:

    --reads '{A.fq,B.fq,C.fq,D.bam}' --metadata metadata.tsv --merge diet,tissue

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet and __liver__ tissue pair.
- __C__ mapped reads will be left untouched, resulting in gene and transcript counts for the __wheat__ diet and __liver__ tissue pair.
- __D__ mapped reads will be left untouched, resulting in gene and transcript counts for the __wheat__ diet and __muscle__ tissue pair.


## Nextflow wrapper script

Using the `./run` wrapper script instead of the `nextflow run` command grants several benefits:
- All results and temporary files are written to the output directory, keeping the execution directory clean.
- All temporary files are deleted once the pipeline has successfully completed.
- The pipeline can be resumed from any directory with nextflow's `-resume` option.


## About this project

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
