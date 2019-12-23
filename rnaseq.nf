// Parse parameters
if (!params.outdir || !params.reads || !params.genome || !params.annotation) {
  exit 1, """Required options:
  --outdir <directory>                   Output directory
  --reads <reads1.fq> <reads2.fq>        Reads files
  --genome <genome.fa>                   Reference genome file
  --annotation <annotation.gff>          Genome annotation file
  """
}

Channel.fromPath(
  params.reads.tokenize(','),
  checkIfExists: true
).into {
  reads_to_control;
  reads_to_trim;
}

Channel.fromPath(
  params.genome,
  checkIfExists: true
).set {
  genome_to_index;
}

Channel.fromPath(
  params.annotation,
  checkIfExists: true
).into {
  annotation_to_index;
  annotation_to_assemble;
}

// Control reads quality
process control {
  tag 'FastQC'

  publishDir "${params.outdir}/quality/raw", mode: 'copy'

  input:
    path reads from reads_to_control.collect()

  output:
    file "*_fastqc.html"

  script:
    """
    fastqc $reads
    """
}

// Remove adaptators from reads
process trim {
  tag 'Trim Galore'

  publishDir "${params.outdir}", mode: 'copy', saveAs: {
    filename ->
    if (filename.contains("_trimming_report.txt")) "logs/trim_galore/$filename"
    else if (filename.contains("_fastqc.html")) "quality/trimmed/$filename"
    else if (filename.contains(".fq.gz")) "reads/trimmed/$filename"
  }

  input:
    path reads from reads_to_trim.collect()

  output:
    file "*_trimming_report.txt"
    file "*_fastqc.html"
    file "*_val_?.fq.gz" into reads_to_map

  script:
    """
    trim_galore $reads --paired --fastqc --gzip
    """
}

// Index reference genome
process index {
  tag 'STAR'

  input:
    path genome from genome_to_index.collect()
    path annotation from annotation_to_index.collect()

  output:
    file "index" into index_to_map

  script:
    """
    mkdir index
    STAR --runThreadN 10 \\
         --runMode genomeGenerate \\
         --genomeDir index \\
         --sjdbGTFfile $annotation \\
         --genomeFastaFiles $genome
    """
}

// Map reads to reference genome
process map {
  tag 'STAR'

  publishDir "${params.outdir}", mode: 'copy', saveAs: {
    filename ->
    if (filename.contains("Log.")) "logs/star/$filename"
    else if (filename.contains(".out.")) "mapping/$filename"
  }

  input:
    path reads from reads_to_map.collect()
    path index from index_to_map.collect()

  output:
    file "Log.*"
    file "*.out.tab"
    file "*.out.bam" into mapping_to_assemble

  script:
    """
    STAR --runThreadN 10 \\
         --readFilesCommand zcat \\
         --outSAMtype BAM SortedByCoordinate \\
         --genomeDir $index \\
         --readFilesIn $reads
    """
}

// Assemble and quantify transcripts
process assemble {
  tag 'StringTie'

  publishDir "${params.outdir}/assembly", mode: 'copy'

  input:
    path annotation from annotation_to_assemble.collect()
    path mapping from mapping_to_assemble.collect()

  output:
    file "assembled_transcripts.gff"

  script:
    """
    stringtie $mapping -G $annotation -o assembled_transcripts.gff
    """
}
