# TAGADA: Transcripts And Genes Assembly, Deconvolution, Analysis

## Installation

To use this pipeline, simply clone or download this repository, and install the dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 21.04.1
- [Docker](https://docs.docker.com/engine/install/) >= 19.03.2 or [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) >= 3.7.3

## Usage

Try out this nextflow pipeline with:

    ./nextflow-run FAANG/analysis-TAGADA --revision 1.0.1 --output directory --profile test docker

The `./nextflow-run` launcher script replaces the `nextflow run` command and grants these benefits:
- Options can receive multiple space-separated parameters and unquoted globs.
- Long options are preceded by double dashes, following GNU conventions.
- Temporary files and logs are written to the output directory, keeping the execution directory clean.
- Temporary files can be automatically deleted after the pipeline has successfully completed.
- The pipeline can be resumed from any directory with the `--resume` option.

### Arguments

| Option | Parameter(s) | Description | Requirement |
|--------|--------------|-------------|-------------|
| __`--profile`__ | `profile1` `profile2` `...` | Profile(s) to use when<br>running the pipeline.<br>Specify the profiles that<br>fit your infrastructure<br>among `singularity`,<br>`docker`, `kubernetes`,<br>`slurm`. | Required |
| __`--output`__ | `directory` | Output directory where<br>all temporary files, logs,<br>and results are written. | Required |
| __`--reads`__ | `reads.fq` `*.bam` `...` | Input `fastq` file(s)<br>and/or `bam` file(s).<br><br>For single-end reads,<br>your files must end with:<br>`.fq[.gz]`<br><br>For paired-end reads,<br>your files must end with:<br>`_[R]{1,2}.fq[.gz]`<br><br>For mapped reads,<br>your files must end with:<br>`.bam`<br><br>You may also provide urls<br>of files to be downloaded.<br><br>If the files are numerous,<br>you may provide a `.txt`<br>sheet with one path or url<br>per line. | Required |
| __`--annotation`__ | `annotation.gtf[.gz]` | Input reference<br>annotation file or url. | Required |
| __`--genome`__ | `genome.fa[.gz]` | Input genome<br>sequence file or url. | Required |
| __`--index`__ | `directory[.tar.gz]` | Input genome index<br>directory or url. | Optional to skip<br>genome indexing. |
| __`--metadata`__ | `metadata.tsv` | Input tabulated<br>metadata file or url. | Required if `--merge`<br>is provided. |
| __`--merge`__ | `factor1` `factor2` `...` | Factor(s) to merge<br>mapped reads. See<br>the [merge factors](https://github.com/FAANG/analysis-TAGADA#merge-factors)<br>section for details. | Optional |
| __`--max-cpus`__ | `16` | Maximum number of<br>CPU cores that can be<br>used for each process.<br>This is a limit, not the<br>actual number of<br>requested CPU cores. | Optional |
| __`--max-memory`__ | `64GB` | Maximum memory that<br>can be used for each<br>process. This is a limit,<br>not the actual amount<br>of alloted memory. | Optional |
| __`--max-time`__ | `12h` | Maximum time that can<br>be spent on each<br>process. This is a limit<br>and has no effect on the<br>duration of each process.| Optional |
| __`--resume`__ | | Preserve temporary files<br>and resume the pipeline<br>from the last completed<br>process. If this option is<br>absent, temporary files<br>will be deleted upon<br>completion, and the<br>pipeline will not be<br>resumable. | Optional |
| __`--feelnc-args`__ | `'--mode shuffle ...'` | Custom arguments to<br>pass to FEELnc's<br>[coding potential](https://github.com/tderrien/FEELnc#2--feelnc_codpotpl) script<br>when detecting long<br>non-coding RNAs. | Optional |
| __`--skip-feelnc`__ | | Skip the detection of long<br>non-coding RNAs with FEELnc. | Optional |

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

- __A__ and __B__ mapped reads will be merged, resulting in counts for the __corn__ diet and __liver__ tissue pair.
- __C__ mapped reads will be left alone, resulting in counts for the __wheat__ diet and __liver__ tissue pair.
- __D__ mapped reads will be left alone, resulting in counts for the __wheat__ diet and __muscle__ tissue pair.


## Workflow

The pipeline executes the following main processes:
1. Control reads __quality__ with [FastQC](https://github.com/s-andrews/FastQC).
2. __Trim__ adaptators from reads with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).
3. Estimate __overhang__ length of splice junctions, and __index__ the genome sequence with [STAR](https://github.com/alexdobin/STAR).  
   The indexed genome is saved to `output/index`.
4. __Map__ reads to the indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   The mapped reads are saved to `output/maps`.
5. Estimate __direction__ and __length__ of mapped reads, and compute genome __coverage__ with [Bedtools](https://github.com/arq5x/bedtools2).  
   The coverage files are saved to `output/coverage`.
6. __Merge__ mapped reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merge factors](#merge-factors) section for details.
7. __Assemble__ transcripts and __combine__ them into a novel annotation with [StringTie](https://github.com/gpertea/stringtie).  
   The novel annotation is saved to `output/annotation`.
8. __Detect__ long non-coding RNAs with [FEELnc](https://github.com/tderrien/FEELnc).  
   The annotations of long non-coding RNAs are saved to `output/annotation`.
9. __Quantify__ genes and transcripts with [StringTie](https://github.com/gpertea/stringtie), and __format__ them into tabulated files.  
   The TPM values and read counts for each annotation are saved to `output/quantification`.
10. Aggregates quality controls into a __report__ with [MultiQC](https://github.com/ewels/MultiQC).  
    The report is saved to `output/control`.


## About this project

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.
