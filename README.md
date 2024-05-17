# TAGADA: Transcript And Gene Assembly, Deconvolution, Analysis

TAGADA is a Nextflow pipeline that processes RNA-Seq data. It parallelizes multiple tasks to control reads quality, align reads to a reference genome, assemble new transcripts to create a novel annotation, and quantify genes and transcripts.


## Table of contents

- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Nextflow options](#nextflow-options)
  - [Input and output options](#input-and-output-options)
  - [Merge options](#merge-options)
  - [Assembly options](#assembly-options)
  - [Skip options](#skip-options)
  - [Resources options](#resources-options)
- [Custom resources](#custom-resources)
  - [Example configuration](#example-configuration)
- [Metadata](#metadata)
  - [Example metadata](#example-metadata)
- [Merging inputs](#merging-inputs)
  - [Merging inputs by a single factor](#merging-inputs-by-a-single-factor)
  - [Merging inputs by an intersection of factors](#merging-inputs-by-an-intersection-of-factors)
- [Workflow and results](#workflow-and-results)
- [Novel annotation](#novel-annotation)
- [Funding](#funding)
- [Citing](#citing)


## Dependencies

To use this pipeline you will need:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 21.04.1
- [Docker](https://docs.docker.com/engine/install/) >= 19.03.2 or [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) >= 3.7.3


## Usage

A small dataset is provided to test this pipeline. To try it out, use this command:

    nextflow run FAANG/analysis-TAGADA -profile test,docker -revision 2.1.3 --output directory

### Nextflow options

The pipeline is written in Nextflow, which provides the following default options:

<table>
  <thead>
    <tr>
      <th width=222px>Option</th>
      <th width=220px>Example</th>
      <th width=215px>Description</th>
      <th width=153px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>-profile</code></strong></td>
      <td nowrap><code>profile1,profile2,etc.</code></td>
      <td>Profile(s) to use when running the pipeline. Specify the profiles that fit your infrastructure among <code>singularity</code>, <code>docker</code>, <code>kubernetes</code>, <code>slurm</code>.</td>
      <td align=center>Required</td>
    </tr>
    <tr>
      <td nowrap><strong><code>-config</code></strong></td>
      <td nowrap><code>custom.config</code></td>
      <td>
        Configuration file tailored to your infrastructure and dataset.<br><br>
        To find a configuration file for your infrastructure, browse <a href="https://github.com/nf-core/configs/tree/master/conf">nf-core configs</a>.<br><br>
        Some large datasets require more computing resources than the pipeline defaults. To specify custom resources for specific processes, see the <a href="#custom-resources">custom resources</a> section.
      </td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>-revision</code></strong></td>
      <td nowrap><code>version</code></td>
      <td>Version of the pipeline to launch.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>-work-dir</code></strong></td>
      <td nowrap><code>directory</code></td>
      <td>Work directory where all temporary files are written.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>-resume</code></strong></td>
      <td nowrap></td>
      <td>Resume the pipeline from the last completed process.</td>
      <td align=center>Optional</td>
    </tr>
  </tbody>
</table>

For more Nextflow options, see [Nextflow's documentation](https://www.nextflow.io/docs/latest/cli.html#run).

### Input and output options

<table>
  <thead>
    <tr>
      <th width=222px>Option</th>
      <th width=220px>Example</th>
      <th width=215px>Description</th>
      <th width=153px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>--output</code></strong></td>
      <td nowrap><code>directory</code></td>
      <td>Output directory where all results are written.</td>
      <td align=center>Required</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--reads</code></strong></td>
      <td nowrap><code>'path/to/reads/*'</code></td>
      <td>
        Input <code>fastq</code> file(s) and/or <code>bam</code> file(s).<br><br>
        For single-end reads, your files must end with:<br><code>.fq[.gz]</code><br><br>
        For paired-end reads, your files must end with:<br><code>_[R]{1,2}.fq[.gz]</code><br><br>
        <code>TAGADA</code>code> will automatically infer read size and strandedness of libraries, but all libraries must have the same<br><br>
        For aligned reads, your files must end with:<br><code>.bam</code><br><br>
        If the provided path includes a wildcard character like <code>*</code>, you must enclose it with quotes to prevent Bash glob expansion, as per <a href="https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters">Nextflow's requirements</a>.<br><br>
        If the files are numerous, you may provide a <code>.txt</code> sheet with one path or url per line.</td>
      <td align=center>Required</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--annotation</code></strong></td>
      <td nowrap><code>annotation.gtf</code></td>
      <td>Input reference<br>annotation file or url.
        Be careful this file should contain both exon and transcript rows and should include gene_id and transcript_id in the 9th field.</td>
      <td align=center>Required</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--genome</code></strong></td>
      <td nowrap><code>genome.fa</code></td>
      <td>Input genome<br>sequence file or url.</td>
      <td align=center>Required</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--index</code></strong></td>
      <td nowrap><code>directory</code></td>
      <td>Input genome index<br>directory or url.</td>
      <td align=center>Optional, to<br>skip genome indexing</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--metadata</code></strong></td>
      <td nowrap><code>metadata.tsv</code></td>
      <td>Input tabulated<br>metadata file or url.</td>
      <td align=center>Required if<br><code>--assemble-by</code><br>or<br><code>--quantify-by</code><br>are provided</td>
    </tr>
  </tbody>
</table>

### Merge options

<table>
  <thead>
    <tr>
      <th width=222px>Option</th>
      <th width=220px>Example</th>
      <th width=215px>Description</th>
      <th width=153px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>--assemble-by</code></strong></td>
      <td nowrap><code>factor1,factor2,etc.</code></td>
      <td>Factor(s) defining groups in which transcripts are assembled. Aligned reads of identical factors are merged and each resulting merge group is processed individually. See the <a href="#merging-inputs">merging inputs</a> section for details.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--quantify-by</code></strong></td>
      <td nowrap><code>factor1,factor2,etc.</code></td>
      <td>Factor(s) defining groups in which transcripts are quantified. Aligned reads of identical factors are merged and each resulting merge group is processed individually. See the <a href="#merging-inputs">merging inputs</a> section for details.</td>
      <td align=center>Optional</td>
    </tr>
  </tbody>
</table>

### Assembly options

<table>
  <thead>
    <tr>
      <th width=262px>Option</th>
      <th width=180px>Example</th>
      <th width=268px>Description</th>
      <th width=100px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>--min-transcript-occurrence</code></strong></td>
      <td nowrap><code>2</code></td>
      <td>After transcripts assembly, rare novel transcripts that appear in few assembly groups are removed from the final novel annotation. By default, if a transcript occurs in less than <code>2</code> assembly groups, it is removed. If there is only one assembly group, this option defaults to <code>1</code>.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--min-monoexonic-occurrence</code></strong></td>
      <td nowrap><code>2</code></td>
      <td>If specified, rare novel monoexonic transcripts are filtered according to the provided threshold. Otherwise, this option takes the value of<br><code>--min-transcript-occurrence</code>.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--min-transcript-tpm</code></strong></td>
      <td nowrap><code>0.1</code></td>
      <td>After transcripts assembly, novel transcripts with low TPM values in every assembly group are removed from the final novel annotation. By default, if a transcript's TPM value is lower than <code>0.1</code> in every assembly group, it is removed.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--min-monoexonic-tpm</code></strong></td>
      <td nowrap><code>1</code></td>
      <td>If specified, novel monoexonic transcripts with low TPM values are filtered according to the provided threshold. Otherwise, this option takes the value of<br><code>--min-transcript-tpm * 10</code>.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--coalesce-transcripts-with</code></strong></td>
      <td nowrap><code>tmerge</code></td>
      <td>Tool used to coalesce transcripts assemblies into a non-redundant set of transcripts for the novel annotation. Can be <code>tmerge</code> or <code>stringtie</code>. Defaults to <code>tmerge</code>.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--tmerge-args</code></strong></td>
      <td nowrap><code>'--endFuzz 10000'</code></td>
      <td>Custom arguments to pass to <a href="https://github.com/julienlag/tmerge#options">tmerge</a> when coalescing transcripts.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--feelnc-filter-args</code></strong></td>
      <td nowrap><code>'--size 200'</code></td>
      <td>Custom arguments to pass to FEELnc's <a href="https://github.com/tderrien/FEELnc#1--feelnc_filterpl">filter</a> script when detecting long non-coding transcripts.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--feelnc-codpot-args</code></strong></td>
      <td nowrap><code>'--mode shuffle'</code></td>
      <td>Custom arguments to pass to FEELnc's <a href="https://github.com/tderrien/FEELnc#2--feelnc_codpotpl">coding potential</a> script when detecting long non-coding transcripts.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--feelnc-classifier-args</code></strong></td>
      <td nowrap><code>'--window 10000'</code></td>
      <td>Custom arguments to pass to FEELnc's <a href="https://github.com/tderrien/FEELnc#3--feelnc_classifierpl">classifier</a> script when detecting long non-coding transcripts.</td>
      <td align=center>Optional</td>
    </tr>
  </tbody>
</table>

### Skip options

<table>
  <thead>
    <tr>
      <th width=222px>Option</th>
      <th width=220px>Example</th>
      <th width=215px>Description</th>
      <th width=153px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>--skip-assembly</code></strong></td>
      <td nowrap></td>
      <td>Skip transcripts assembly with StringTie and skip all subsequent processes working with a novel annotation.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--skip-lnc-detection</code></strong></td>
      <td nowrap></td>
      <td>Skip detection of long non-coding transcripts in the novel annotation with FEELnc.</td>
      <td align=center>Optional</td>
    </tr>
  </tbody>
</table>

### Resources options

<table>
  <thead>
    <tr>
      <th width=222px>Option</th>
      <th width=220px>Example</th>
      <th width=215px>Description</th>
      <th width=153px>Required</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td nowrap><strong><code>--max-cpus</code></strong></td>
      <td nowrap><code>16</code></td>
      <td>Maximum number of CPU cores that can be used for each process. This is a limit, not the actual number of requested CPU cores.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--max-memory</code></strong></td>
      <td nowrap><code>64GB</code></td>
      <td>Maximum memory that can be used for each process. This is a limit, not the actual amount of allotted memory.</td>
      <td align=center>Optional</td>
    </tr>
    <tr>
      <td nowrap><strong><code>--max-time</code></strong></td>
      <td nowrap><code>24h</code></td>
      <td>Maximum time that can be spent on each process. This is a limit and has no effect on the duration of each process.</td>
      <td align=center>Optional</td>
    </tr>
  </tbody>
</table>

## Custom resources

With large datasets, some [workflow processes](#workflow-and-results) may require more computing resources than the pipeline defaults. To customize the amount of resources allotted to specific processes, add a [process scope](https://www.nextflow.io/docs/edge/config.html#scope-process) to your configuration file. Resources provided in the configuration file override the [resources options](#resources-options).

### Example configuration

    -config custom.config

`custom.config`

    process {

      withName: TRIMGALORE_trim_adapters {
        cpus = 8
        memory = 18.GB
        time = 36.h
      }

      withName: STAR_align_reads {
        cpus = 16
        memory = 64.GB
        time = 2.d
      }

    }


## Metadata

Using `--metadata`, you may provide a file describing your inputs with tab-separated factors. The first column must contain file names without file type extensions or paired-end suffixes. There are no constraints on column names or number of columns.

### Example metadata

    --reads reads.txt --metadata metadata.tsv

`reads.txt`

    path/to/A_R1.fq
    path/to/A_R2.fq
    path/to/B.fq.gz
    path/to/C.bam
    path/to/D.fq

`metadata.tsv`

    input    tissue     stage
    A        liver      30 days
    B        liver      30 days
    C        liver      60 days
    D        muscle     60 days


## Merging inputs

When using `--assemble-by` and/or `--quantify-by`, your inputs are merged into experiment groups that share common factors. With `--assemble-by`, transcripts assembly is done individually for each assembly group, and consensus transcripts are kept to generate a novel annotation. With `--quantify-by`, quantification values are given individually for each quantification group.

### Merging inputs by a single factor

    --assemble-by tissue --quantify-by stage

<table>
  <thead align=center>
    <tr>
      <th colspan=3 width=250px>Metadata</th>
      <th rowspan=2 width=190px>Transcripts assembly<br>by tissue</th>
      <th rowspan=2 width=190px>Annotation</th>
      <th rowspan=2 width=190px>Quantification<br>by stage</th>
    </tr>
    <tr>
      <th>input</th>
      <th>tissue</th>
      <th>stage</th>
    </tr>
  </thead>
  <tbody align=center>
    <tr>
      <td>A</td>
      <td>liver</td>
      <td>30 days</td>
      <td rowspan=3>A, B, C<br>↓<br>liver</td>
      <td rowspan=4>liver, muscle<br>↓<br>novel annotation</td>
      <td rowspan=2>A, B<br>↓<br>30 days</td>
    </tr>
    <tr>
      <td>B</td>
      <td>liver</td>
      <td>30 days</td>
    </tr>
    <tr>
      <td>C</td>
      <td>liver</td>
      <td>60 days</td>
      <td rowspan=2>C, D<br>↓<br>60 days</td>
    </tr>
    <tr>
      <td>D</td>
      <td>muscle</td>
      <td>60 days</td>
      <td>D<br>↓<br>muscle</td>
    </tr>
  </tbody>
</table>

### Merging inputs by an intersection of factors

    --assemble-by tissue,stage

<table>
  <thead align=center>
    <tr>
      <th colspan=3 width=250px>Metadata</th>
      <th rowspan=2 width=190px>Transcripts assembly<br>by tissue and stage</th>
      <th rowspan=2 width=190px>Annotation</th>
      <th rowspan=2 width=190px>Quantification<br>by input</th>
    </tr>
    <tr>
      <th>input</th>
      <th>tissue</th>
      <th>stage</th>
    </tr>
  </thead>
  <tbody align=center>
    <tr>
      <td>A</td>
      <td>liver</td>
      <td>30 days</td>
      <td rowspan=2>A, B<br>↓<br>liver at 30 days</td>
      <td rowspan=4>liver at 30 days,<br>liver at 60 days,<br>muscle at 60 days<br>↓<br>novel annotation</td>
      <td>A</td>
    </tr>
    <tr>
      <td>B</td>
      <td>liver</td>
      <td>30 days</td>
      <td>B</td>
    </tr>
    <tr>
      <td>C</td>
      <td>liver</td>
      <td>60 days</td>
      <td>C<br>↓<br>liver at 60 days</td>
      <td>C</td>
    </tr>
    <tr>
      <td>D</td>
      <td>muscle</td>
      <td>60 days</td>
      <td>D<br>↓<br>muscle at 60 days</td>
      <td>D</td>
    </tr>
  </tbody>
</table>


## Workflow and results

The pipeline executes the following processes:

1. `FASTQC_control_reads`  
   Control reads quality with [FastQC](https://github.com/s-andrews/FastQC).

2. `TRIMGALORE_trim_adapters`  
   Trim adapters with [Trim Galore](https://github.com/FelixKrueger/TrimGalore).

3. `STAR_index_genome`  
   Index genome with [STAR](https://github.com/alexdobin/STAR).  
   The indexed genome is saved to `output/index`.

4. `STAR_align_reads`  
   Align reads to the indexed genome with [STAR](https://github.com/alexdobin/STAR).  
   Aligned reads are saved to `output/alignment` in `.bam` files.

5. `BEDTOOLS_compute_coverage`  
   Compute genome coverage with [Bedtools](https://github.com/arq5x/bedtools2).  
   Coverage information is saved to `output/coverage` in `.bed` files.

6. `SAMTOOLS_merge_reads`  
   Merge aligned reads by factors with [Samtools](https://github.com/samtools/samtools).  
   See the [merging inputs](#merging-inputs) section for details.

7. `STRINGTIE_assemble_transcripts`  
   Assemble transcripts in each individual assembly group with [StringTie](https://github.com/gpertea/stringtie).

8. `TAGADA_filter_transcripts`  
   Filter rare transcripts that appear in few assembly groups and poorly-expressed transcripts with low TPM values.

9. `STRINGTIE_coalesce_transcripts` or `TMERGE_coalesce_transcripts`  
   Create a novel annotation with [StringTie](https://github.com/gpertea/stringtie) or [Tmerge](https://github.com/julienlag/tmerge).  
   The novel annotation is saved to `output/annotation` in a `.gtf` file.

10. `FEELNC_classify_transcripts`  
   Detect long non-coding transcripts with [FEELnc](https://github.com/tderrien/FEELnc).  
   The annotation saved to `output/annotation` is updated with the results.

11. `STRINGTIE_quantify_expression`  
    Quantify genes and transcripts with [StringTie](https://github.com/gpertea/stringtie).  
    Counts and TPM matrices are saved to `output/quantification` in `.tsv` files.

12. `MULTIQC_generate_report`  
    Aggregate quality controls into a report with [MultiQC](https://github.com/ewels/MultiQC).  
    The report is saved to `output/control` in a `.html` file.


## Novel annotation

The novel annotation contains information from [StringTie](https://github.com/gpertea/stringtie), [Tmerge](https://github.com/julienlag/tmerge), and [FEELnc](https://github.com/tderrien/FEELnc). It is provided in gtf format with exon, transcript and gene rows. Row attributes vary depending on which tool was used to coalesce transcripts.

<br>

    --coalesce-transcripts-with tmerge

- `gene_id`  
  All rows. The Tmerge `gene_id` starting with LOC.

- `ref_gene_id`  
  All rows. A comma-separated list of reference annotation `gene_id` when a Tmerge transcript is made of at least one reference transcript, otherwise a dot.

- `transcript_id`  
  Exon and transcript rows. The Tmerge `transcript_id` starting with TM, unless the transcript is exactly identical to a reference transcript, in which case the reference annotation `transcript_id` is provided.

- `tmerge_tr_id`  
  Exon and transcript rows. Optional. A comma-separated list of Tmerge `transcript_id` if the current `transcript_id` is from the reference annotation, to list which initial Tmerge transcripts it is made of.

- `transcript_biotype`  
  Exon and transcript rows. Optional. The reference annotation `transcript_biotype` of the `transcript_id`.

- `feelnc_biotype`  
  Exon and transcript rows. Optional. The transcript biotype determined by FEELnc (lncRNA, mRNA, noORF, or TUCp) if the transcript has been classified.

- `contains`, `contains_count`, `3p_dists_to_3p`, `5p_dists_to_5p`, `flrpm`, `longest`, `longest_FL_supporters`, `longest_FL_supporters_count`, `mature_RNA_length`, `meta_3p_dists_to_5p`, `meta_5p_dists_to_5p`, `rpm`, `spliced`  
  Transcript rows. Attributes provided by Tmerge.

<br>

    --coalesce-transcripts-with stringtie

- `gene_id`  
  All rows. The StringTie `gene_id` starting with MSTRG.

- `ref_gene_id`  
  All rows. Optional. The reference annotation `gene_id`.

- `ref_gene_name`  
  All rows. Optional. The reference annotation `gene_name`.

- `transcript_id`  
  Exon and transcript rows. The StringTie `transcript_id` starting with MSTRG, unless the transcript is exactly identical to a reference transcript, in which case the reference annotation `transcript_id` is provided.

- `transcript_biotype`  
  Exon and transcript rows. Optional. The reference annotation `transcript_biotype` of the `transcript_id`.

- `feelnc_biotype`  
  Exon and transcript rows. Optional. The transcript biotype determined by FEELnc (lncRNA, mRNA, noORF, or TUCp) if the transcript has been classified.

- `exon_number`  
  Exon rows. The StringTie `exon_number` starting from 1 within a given transcript.


## Funding

The GENE-SWitCH project has received funding from the European Union’s [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) research and innovation program under Grant Agreement No 817998.

This repository reflects only the listed contributors views. Neither the European Commission nor its Agency REA are responsible for any use that may be made of the information it contains.

## Citing

If you use TAGADA in a publication, please cite [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10578202/pdf/lqad089.pdf):

Kurylo C, Guyomar C, Foissac S, Djebali S. TAGADA: a scalable pipeline to improve genome annotations with RNA-seq data. NAR Genomics and Bioinformatics. 2023 Dec 1;5(4):lqad089.
