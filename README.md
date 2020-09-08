# GENE-SWitCH project RNA-Seq analysis pipeline


## Installation

To use this pipeline, simply clone or download this repository, and install the dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 20.01.0
- [Docker](https://docs.docker.com/engine/install/) >= 19.03.2 or [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) >= 3.4


## Usage

Execute this nextflow pipeline with:

    ./nextflow-run proj-gs-rna-seq [arguments]

The `./nextflow-run` launcher script replaces the `nextflow run` command and grants these benefits:
- Options can receive multiple space-separated parameters and unquoted globs.
- Long options are preceded by double dashes, following GNU conventions.
- Temporary files and logs are written to the output directory, keeping the execution directory clean.
- Temporary files are deleted after the pipeline has successfully completed.
- The pipeline can be resumed from any directory with the `--resume` option.

### Arguments

| Option | Parameter(s) | Description | Requirement |
|--------|--------------|-------------|-------------|
| __`--profile`__ | `<profile1>` `<profile2>` `...` | Profile(s) to use when<br>running the pipeline.<br>Specify the profiles that<br>fit your infrastructure<br>among `singularity`,<br>`docker`, `kubernetes`,<br>`slurm`. | Required |
| __`--output`__ | `<directory>` | Output directory where<br>all temporary files, logs,<br>and results are written. | Required |
| __`--reads`__ | `<reads.fq>` `<*.bam>` `...` | Input `fastq` file(s)<br>and/or `bam` file(s).<br><br>For single-end reads,<br>name your files:<br>`name.fq[.gz]`<br><br>For paired-end reads,<br>name your files:<br>`name_R{1,2}.fq[.gz]`<br><br>For mapped reads,<br>name your files:<br>`name.bam` | Required |
| __`--annotation`__ | `<annotation.gff>` | Input reference<br>annotation file. | Required |
| __`--genome`__ | `<genome.fa>` | Input genome<br>sequence file. | Required if `fastq`<br>files are provided<br>and `--index` is<br>absent. |
| __`--index`__ | `<directory>` | Input genome index<br>directory. Overrides<br>`--genome`. | Required if `fastq`<br>files are provided<br>and `--genome` is<br>absent. |
| __`--metadata`__ | `<metadata.tsv>` | Input tabulated<br>metadata file. | Required if `--merge`<br>is provided. |
| __`--merge`__ | `<factor1>` `<factor2>` `...` | Factor(s) to merge<br>mapped reads. See<br>the [merge factors](https://github.com/FAANG/proj-gs-rna-seq#merge-factors)<br>section for details. | Optional |
| __`--max-cpus`__ | `<16>` | Maximum number of<br>CPU cores that can be<br>used for each process.<br>This is a limit, not the<br>actual number of<br>requested CPU cores. | Optional |
| __`--max-memory`__ | `<64GB>` | Maximum memory that<br>can be used for each<br>process. This is a limit,<br>not the actual amount<br>of alloted memory. | Optional |
| __`--max-time`__ | `<12h>` | Maximum time that can<br>be spent on each<br>process. This is a limit<br>and has no effect on the<br>duration of each process.| Optional |
| __`--resume`__ | | Resume the pipeline after<br>an interruption. Previously<br>completed processes will<br>be skipped. | Optional |
| __`--keep-temp`__ | | Do not delete temporary<br>files upon completion. | Optional |

### Merge factors

Use the `--merge` and `--metadata` options together to merge mapped reads. This results in genes and transcripts being quantified by __factors__ rather than by __inputs__.

The metadata file consists of tab-separated values describing your inputs. The first column must contain file names without extensions. There is no restriction on column names or number of columns.

#### Examples

Given the following tabulated metadata file:

    input    diet      tissue
    A        corn      liver
    B        corn      liver
    C        wheat     liver
    D        wheat     muscle

With the following arguments:

    --reads A_R1.fq A_R2.fq B.fq C.fq D.bam --metadata metadata.tsv --merge diet

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet.
- __C__ and __D__ mapped reads will be merged, resulting in gene and transcript counts for the __wheat__ diet.

With the following arguments:

    --reads A_R1.fq A_R2.fq B.fq C.fq D.bam --metadata metadata.tsv --merge diet tissue

- __A__ and __B__ mapped reads will be merged, resulting in gene and transcript counts for the __corn__ diet and __liver__ tissue pair.
- __C__ mapped reads will be left alone, resulting in gene and transcript counts for the __wheat__ diet and __liver__ tissue pair.
- __D__ mapped reads will be left alone, resulting in gene and transcript counts for the __wheat__ diet and __muscle__ tissue pair.


## Workflow

The pipeline executes the following processes:
1. Control reads __quality__ with [FastQC](https://github.com/s-andrews/FastQC).  
   Outputs quality reports to `output/quality/raw`.
2. __Trim__ adaptators from reads with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).  
   Outputs quality reports to `output/quality/trimmed`.
3. Estimate __overhang__ length of splice junctions, and __index__ the genome sequence with [STAR](https://github.com/alexdobin/STAR).  
   Outputs indexed genome to `output/index`.
4. __Map__ reads to indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   Outputs mapped reads to `output/maps`.
5. Estimate __direction__ and __length__ of mapped reads, and compute genome __coverage__ with [Bedtools](https://github.com/arq5x/bedtools2).  
   Outputs bedGraph files to `output/coverage`.
6. __Merge__ mapped reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merge factors](#merge-factors) section for details.
7. __Assemble__ transcripts and __combine__ them into a new assembly annotation with [StringTie](https://github.com/gpertea/stringtie).  
   Outputs the new assembly annotation to `output/assembly`.
8. __Quantify__ genes and transcripts with [StringTie](https://github.com/gpertea/stringtie), and __format__ them into tabulated files.  
   Outputs TPM values and read counts to `output/quantification` for the reference and assembly annotations.


## About this project

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
