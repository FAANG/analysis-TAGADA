# GENE-SWitCH project RNA-Seq analysis pipeline


## Installation

To use this pipeline, simply clone or download this repository, and install the dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 20.01.0
- [Docker](https://docs.docker.com/engine/install/) >= 19.03.2 or [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) >= 3.4


## Usage

Execute this nextflow pipeline with:

    ./run rnaseq.nf [arguments]

The `./run` launcher script replaces the `nextflow run` command and grants these benefits:
- Options can receive multiple space-separated parameters.
- Long options are preceded by double dashes, following GNU conventions.
- Temporary files and logs are written to the output directory, keeping the execution directory clean.
- Temporary files are deleted after the pipeline has successfully completed.
- The pipeline can be resumed from any directory with the `--resume` option.

### Arguments

| Option | Parameter(s) | Description | Requirement |
|--------|--------------|-------------|-------------|
| __`--profile`__ | `<profile1>` `<profile2>` `...` | Profile(s) to use when running the<br>pipeline. Specify the profiles that fit<br>your infrastructure among `slurm`,<br>`singularity`, `docker`. | Required |
| __`--output`__ | `<directory>` | Output directory where all temporary<br>files, logs, and results are written. | Required |
| __`--reads`__ | `<reads.fq>` `<*.bam>` `...` | Input `fastq` file(s) and/or `bam` file(s).<br><br>For single-end reads, name your files:<br>`prefix.{fq,fastq}[.gz]`<br><br>For paired-end reads, name your files:<br>`prefix_R{1,2}.{fq,fastq}[.gz]`<br><br>For mapped reads, name your files:<br>`prefix.bam` | Required |
| __`--annotation`__ | `<annotation.gff>` | Input reference annotation file. | Required |
| __`--genome`__ | `<genome.fa>` | Input genome sequence file. | Required if `fastq`<br>files are provided<br>and `--index` is<br>absent. |
| __`--index`__ | `<directory>` | Input genome index directory.<br>Overrides `--genome`. | Required if `fastq`<br>files are provided<br>and `--genome` is<br>absent. |
| __`--metadata`__ | `<metadata.tsv>` | Input tabulated metadata file. | Required if `--merge`<br>is provided. |
| __`--merge`__ | `<factor1>` `<factor2>` `...` | Factor(s) to merge reads files. See<br>the [merge factors](https://github.com/FAANG/proj-gs-rna-seq#merge-factors) section for details. | Optional |
| __`--direction`__ | `<rf\|fr>` | Direction of reads. Either `rf` or `fr`. | Optional |
| __`--max-cpus`__ | `<16>` | Maximum number of CPU cores that<br>can be used for each process. This<br>is a limit, not the actual number of<br>requested CPU cores. | Optional |
| __`--max-memory`__ | `<64GB>` | Maximum memory that can be used<br>for each process. This is a limit, not<br>the actual amount of alloted memory. | Optional |
| __`--max-time`__ | `<12h>` | Maximum time that can be spent<br>on each process. This is a limit and<br>has no effect on the duration of each<br>process.| Optional |
| __`--resume`__ | | Resume the pipeline after interruption.<br>Previously completed processes will<br>be skipped. | Optional |

### Merge factors

Use the `--merge` and `--metadata` options together to merge reads files after trimming and mapping. This results in genes and transcripts being counted by __factor__ rather than by __input file__.

The metadata file consists of tab-separated values describing your input files. The first column must contain input file prefixes without extensions. There is no restriction on column names or number of columns.

#### Examples

Given the following tabulated metadata file:

    input    diet      tissue
    A        corn      liver
    B        corn      liver
    C        wheat     liver
    D        wheat     muscle

With the following arguments:

    --reads A.fq B.fq C.fq D.bam --metadata metadata.tsv --merge diet

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet.
- __C__ and __D__ mapped reads will be merged, resulting in gene and transcript counts for the __wheat__ diet.

With the following arguments:

    --reads A.fq B.fq C.fq D.bam --metadata metadata.tsv --merge diet tissue

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet and __liver__ tissue pair.
- __C__ mapped reads will be left alone, resulting in gene and transcript counts for the __wheat__ diet and __liver__ tissue pair.
- __D__ mapped reads will be left alone, resulting in gene and transcript counts for the __wheat__ diet and __muscle__ tissue pair.


## Workflow

The pipeline executes the following processes:
1. __Control__ reads quality with [FastQC](https://github.com/s-andrews/FastQC).  
   Outputs quality reports to `output/quality/raw`.
2. __Trim__ adaptators from reads with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).  
   Outputs quality reports to `output/quality/trimmed`.
3. __Index__ genome sequence wih [STAR](https://github.com/alexdobin/STAR).  
   Outputs indexed genome to `output/index`.
4. __Map__ reads to indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   Outputs mapped reads to `output/maps`.
5. __Merge__ mapped reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merge factors](#merge-factors) section for details.
6. __Assemble__ transcripts and __combine__ them into a new assembly annotation with [StringTie](https://github.com/gpertea/stringtie).  
   Outputs the new assembly annotation to `output/annotation`.
7. __Count__ genes and transcripts with [StringTie](https://github.com/gpertea/stringtie), and __format__ them into tabulated files.  
   Outputs TPM counts and average per-base read coverage to `output/counts`.  
   Counts are given for the reference and assembly annotations separately.


## About this project

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
