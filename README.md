# TAGADA: Transcripts And Genes Assembly, Deconvolution, Analysis

TAGADA is a Nextflow pipeline that processes RNA-Seq data. It parallelizes multiple tasks to control reads quality, map reads to a reference genome, assemble new transcripts to create a novel annotation, and quantify genes and transcripts.

## Table of Contents

- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Nextflow options](#nextflow-options)
  - [Pipeline options](#pipeline-options)
  - [Merge factors](#merge-factors)
- [Workflow](#workflow)
- [About](#about)

## Dependencies

To use this pipeline you will need:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 21.04.1
- [Docker](https://docs.docker.com/engine/install/) >= 19.03.2 or [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) >= 3.7.3

## Usage

A small dataset is provided to test this pipeline. To try it out, use this command:

    nextflow run FAANG/analysis-TAGADA -profile test,docker -revision 1.1.0 --output directory

### Nextflow Options

The pipeline is written in Nextflow, which provides the following default options.  
All Nextflow options are [documented here](https://www.nextflow.io/docs/latest/cli.html#run).

| Option | Parameters | Description | Requirement |
|--------|--------------|-------------|-------------|
| __`-profile`__ | `profile1,profile2,[...]` | Profile(s) to use when<br>running the pipeline.<br>Specify the profiles that<br>fit your infrastructure<br>among `singularity`,<br>`docker`, `kubernetes`,<br>`slurm`. | Required |
| __`-revision`__ | `version` | Version of the pipeline<br>to launch. | Optional |
| __`-resume`__ | | Resume the pipeline<br>from the last completed<br>process. | Optional |
| __`-work-dir`__ | | Work directory where<br>all temporary files are<br>written. | Optional |

### Pipeline options

| Option | Parameters | Description | Requirement |
|--------|--------------|-------------|-------------|
| __`--output`__ | `directory` | Output directory where<br>all results are written. | Required |
| __`--reads`__ | `'path/to/reads/*'` | Input `fastq` file(s)<br>and/or `bam` file(s).<br><br>For single-end reads,<br>your files must end with:<br>`.fq[.gz]`<br><br>For paired-end reads,<br>your files must end with:<br>`_[R]{1,2}.fq[.gz]`<br><br>For mapped reads,<br>your files must end with:<br>`.bam`<br><br>If the files are numerous,<br>you may provide a `.txt`<br>sheet with one path or url<br>per line. | Required |
| __`--annotation`__ | `annotation.gtf[.gz]` | Input reference<br>annotation file or url. | Required |
| __`--genome`__ | `genome.fa[.gz]` | Input genome<br>sequence file or url. | Required |
| __`--index`__ | `directory[.tar.gz]` | Input genome index<br>directory or url. | Optional, to skip<br>genome indexing. |
| __`--metadata`__ | `metadata.tsv` | Input tabulated<br>metadata file or url. | Required if<br>`--assemble-by` or<br>`--quantify-by`<br>are provided. |
| __`--assemble-by`__ | `factor1,factor2,[...]` | Factor(s) defining groups<br>in which transcripts are<br>assembled. Mapped<br>reads of identical factors<br>are merged and each<br>resulting merge group is<br>processed individually.<br>See the [merge factors](#merge-factors)<br>section for details. | Optional |
| __`--quantify-by`__ | `factor1,factor2,[...]` | Factor(s) defining groups<br>in which transcripts are<br>quantified. Mapped<br>reads of identical factors<br>are merged and each<br>resulting merge group is<br>processed individually.<br>See the [merge factors](#merge-factors)<br>section for details. | Optional |
| __`--skip-assembly`__ | | Skip transcripts<br>assembly with StringTie. | Optional |
| __`--skip-feelnc`__ | | Skip detection of long<br>non-coding RNAs<br>with FEELnc. | Optional |
| __`--feelnc-args`__ | `'--mode shuffle [...]'` | Custom arguments to<br>pass to FEELnc's<br>[coding potential](https://github.com/tderrien/FEELnc#2--feelnc_codpotpl) script<br>when detecting long<br>non-coding RNAs. | Optional |
| __`--max-cpus`__ | `16` | Maximum number of<br>CPU cores that can be<br>used for each process.<br>This is a limit, not the<br>actual number of<br>requested CPU cores. | Optional |
| __`--max-memory`__ | `64GB` | Maximum memory that<br>can be used for each<br>process. This is a limit,<br>not the actual amount<br>of alloted memory. | Optional |
| __`--max-time`__ | `12h` | Maximum time that can<br>be spent on each<br>process. This is a limit<br>and has no effect on the<br>duration of each process.| Optional |

### Merge factors

Transcripts assembly and quantification can be done by __factors__ instead of by __input__. When using the `--assemble-by` and `--quantify-by` options, mapped reads of identical factors are merged and each resulting merge group is processed individually.

Factors are specified in the metadata file. This file consists of tab-separated values describing your inputs. The first column must contain file names without extensions. There is no restriction on column names or number of columns.

#### Examples

Given the following tabulated metadata file:

    input    tissue     stage
    A        liver      30 days
    B        liver      30 days
    C        liver      60 days
    D        muscle     60 days

And the following inputs:

    A_R1.fq
    A_R2.fq
    B.fq.gz
    C.bam
    D.fq

With the following arguments:

    --assemble-by tissue --quantify-by stage

- __A__ and __B__ and __C__ mapped reads will be merged for transcripts assembly in the __liver__ tissue.
- __D__ mapped reads will be left alone for transcripts assembly in the __muscle__ tissue.
- __A__ and __B__ mapped reads will be merged for quantification in the __30 days__ stage.
- __C__ and __D__ mapped reads will be merged for quantification in the __60 days__ stage.

With the following arguments:

    --quantify-by tissue,stage

- all mapped reads will be left alone for transcripts assembly in each individual input.
- __A__ and __B__ mapped reads will be merged for quantification in the __liver__ at __30 days__.
- __C__ mapped reads will be left alone for quantification in the __liver__ at __60 days__.
- __D__ mapped reads will be left alone for quantification in the __muscle__ at __60 days__.


## Workflow

The pipeline executes the following processes:
1. Control reads quality with [FastQC](https://github.com/s-andrews/FastQC).
2. Trim adapters with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).
3. Index the genome with [STAR](https://github.com/alexdobin/STAR).  
   The indexed genome is saved to `output/index`.
4. Map reads to the indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   Mapped reads are saved to `output/maps` in `.bam` files.
5. Compute genome coverage with [Bedtools](https://github.com/arq5x/bedtools2).  
   Coverage information is saved to `output/coverage` in `.bed` files.
6. Merge mapped reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merge factors](#merge-factors) section for details.
7. Assemble transcripts and create a novel annotation with [StringTie](https://github.com/gpertea/stringtie).  
   The novel annotation is saved to `output/annotation` in a `.gtf` file.
8. Detect long non-coding RNAs with [FEELnc](https://github.com/tderrien/FEELnc).  
   The annotation saved to `output/annotation` is updated with the results.
9. Quantify genes and transcripts with [StringTie](https://github.com/gpertea/stringtie).  
   Counts and TPM matrices are saved to `output/quantification` in `.tsv` files.
10. Aggregate quality controls into a report with [MultiQC](https://github.com/ewels/MultiQC).  
    The report is saved to `output/control` in a `.html` file.


## About

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
