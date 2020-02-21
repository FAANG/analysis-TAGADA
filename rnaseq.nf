// ---------- PARAMETERS

output = params.containsKey('output') ? params.'output' : ''
paired_reads = params.containsKey('paired-reads') ? params.'paired-reads' : ''
single_reads = params.containsKey('single-reads') ? params.'single-reads' : ''
mapped_reads = params.containsKey('mapped-reads') ? params.'mapped-reads' : ''
reference = params.containsKey('reference') ? params.'reference' : ''
index = params.containsKey('index') ? params.'index' : ''
annotation = params.containsKey('annotation') ? params.'annotation' : ''
THREADS = params.containsKey('threads') ? params.'threads' : 1

error = ''

if (!output) {
  error += 'No output given.\n'
}

if (!annotation) {
  error += 'No annotation given.\n'
}

if (!paired_reads && !single_reads && !mapped_reads) {
  error += 'No reads or mapped reads given.\n'
}

if ((paired_reads || single_reads) && (!reference && !index)) {
  error += 'No genome reference or index given.\n'
}

if (error) {
  exit 1, error + '\n' + '''Options:
  --output <directory>                          Output directory
  --paired-reads <reads1.fq> <reads2.fq> ...    Paired reads file(s)
  --single-reads <reads.fq> ...                 Single reads file(s)
  --mapped-reads <map.bam> ...                  Mapped reads file(s)
  --reference <reference.fa>                    Genome reference file
  --index <directory>                           Genome index input directory
  --annotation <annotation.gff>                 Genome annotation file
  --threads <number>                            Number of threads
  '''
}

if (paired_reads) {
  Channel.fromPath(
    paired_reads.tokenize(','),
    checkIfExists: true
  ).into {
    paired_reads_to_control
    paired_reads_to_trim
  }
}

if (single_reads) {
  Channel.fromPath(
    single_reads.tokenize(','),
    checkIfExists: true
  ).into {
    single_reads_to_control
    single_reads_to_trim
  }
}

if (paired_reads || single_reads) {
  Channel.empty().concat(
    paired_reads ? paired_reads_to_control : Channel.empty(),
    single_reads ? single_reads_to_control : Channel.empty()
  ).set {
    reads_to_control
  }

  Channel.empty().concat(
    paired_reads ?
      paired_reads_to_trim.buffer(size: 2).combine(Channel.of('--paired')) :
      Channel.empty(),
    single_reads ?
      single_reads_to_trim.combine(Channel.of('')).combine(Channel.of('')) :
      Channel.empty()
  ).set {
    reads_to_trim
  }
}

if (mapped_reads) {
  Channel.fromPath(
    mapped_reads.tokenize(','),
    checkIfExists: true
  ).into {
    mapped_reads_to_assemble
    mapped_reads_to_quantify
  }
}

if (index) {
  Channel.fromPath(
    index,
    checkIfExists: true
  ).set {
    index_to_map
  }
}

if (reference) {
  Channel.fromPath(
    reference,
    checkIfExists: true
  ).set {
    reference_to_index
  }
}

Channel.fromPath(
  annotation,
  checkIfExists: true
).into {
  annotation_to_index
  annotation_to_assemble
  annotation_to_merge
}

// ---------- WORKFLOW

if (paired_reads || single_reads) {

  // Control reads quality
  process control {
    tag 'FastQC'

    publishDir "$output/reads/raw", mode: 'copy'

    input:
      path reads from reads_to_control.collect()

    output:
      file '*_fastqc.html'

    script:
      """
      fastqc $reads
      """
  }

  // Trim adaptators from reads
  process trim {
    tag 'Trim Galore'

    publishDir "$output", mode: 'copy', saveAs: {
      filename ->
      if (filename.endsWith('_trimming_report.txt')) "logs/trim_galore/$filename"
      else if (filename.endsWith('_fastqc.html')) "reads/trimmed/$filename"
      else if (filename.endsWith('.fq.gz')) {
        filename = filename.minus(~/_paired$/)
        "reads/trimmed/$filename"
      }
    }

    input:
      tuple path(reads1), path(reads2), val(paired) from reads_to_trim

    output:
      file '*_trimming_report.txt'
      file '*_fastqc.html'
      file '*_trimmed.fq.gz_paired' optional true into paired_reads_to_map
      file '*_trimmed.fq.gz' optional true into single_reads_to_map

    script:
      """
      trim_galore $reads1 $reads2 $paired --cores $THREADS --fastqc --gzip

      for f in *_val_?.fq.gz; do
        [ -f "\$f" ] || break
        outfile="\${f%_val_?.fq.gz}"
        mv "\$f" "\$outfile"_trimmed.fq.gz_paired
      done

      for f in *_val_?_fastqc.html; do
        [ -f "\$f" ] || break
        outfile="\${f%_val_?_fastqc.html}"
        mv "\$f" "\$outfile"_trimmed_fastqc.html
      done

      for f in *_trimming_report.txt; do
        [ -f "\$f" ] || break
        outfile="\${f%_trimming_report.txt}"
        outfile="\${outfile%.gz}"
        outfile="\${outfile%.fq}"
        outfile="\${outfile%.fastq}"
        mv "\$f" "\$outfile"_trimming_report.txt
      done
      """
  }

  Channel.empty().concat(
    paired_reads ?
      paired_reads_to_map.flatten().buffer(size: 2) :
      Channel.empty(),
    single_reads ?
      single_reads_to_map.flatten().combine(Channel.of('')) :
      Channel.empty()
  ).set {
    reads_to_map
  }

  if (!index) {

    // Index reference genome
    process index {
      tag 'STAR'

      publishDir "$output/genome", mode: 'copy'

      input:
        path reference from reference_to_index
        path annotation from annotation_to_index

      output:
        file 'index' into index_to_map

      script:
        """
        mkdir index
        STAR --runThreadN $THREADS \\
             --runMode genomeGenerate \\
             --genomeDir index \\
             --sjdbGTFfile $annotation \\
             --genomeFastaFiles $reference
        """
    }
  }

  // Map reads to reference genome
  process map {
    tag 'STAR'

    publishDir "$output", mode: 'copy', saveAs: {
      filename ->
      if (filename.endsWith('.out')) "logs/star/$filename"
      else if (filename.endsWith('.out.tab')) "logs/star/$filename"
      else if (filename.endsWith('.bam')) "reads/mapped/$filename"
    }

    input:
      tuple path(reads1), path(reads2), path(index) from reads_to_map.combine(index_to_map)

    output:
      file '*.out'
      file '*.out.tab'
      file '*.bam' into maps_to_assemble
      file '*.bam' into maps_to_quantify

    script:
      """
      outfile=$reads1
      outfile="\${outfile%_paired}"
      outfile="\${outfile%_trimmed.fq.gz}"
      outfile="\${outfile%_R?}"

      STAR --runThreadN $THREADS \\
           --readFilesCommand zcat \\
           --outSAMtype BAM SortedByCoordinate \\
           --genomeDir $index \\
           --readFilesIn $reads1 $reads2 \\
           --outFileNamePrefix "\$outfile".

      mv *.bam "\$outfile".bam
      """
  }
}

Channel.empty().concat(
  (single_reads || paired_reads) ? maps_to_assemble : Channel.empty(),
  mapped_reads ? mapped_reads_to_assemble : Channel.empty()
).set {
  maps_to_assemble
}

Channel.empty().concat(
  (single_reads || paired_reads) ? maps_to_quantify : Channel.empty(),
  mapped_reads ? mapped_reads_to_quantify : Channel.empty()
).set {
  maps_to_quantify
}

// Assemble transcripts
process assemble {
  tag 'StringTie'

  input:
    path annotation from annotation_to_assemble
    each path(map) from maps_to_assemble

  output:
    file '*.gff' into assemblies_to_merge

  script:
    """
    stringtie $map -p $THREADS -G $annotation -o "\$RANDOM".gff
    """
}

// Merge assemblies
process merge {
  tag 'StringTie'

  input:
    path annotation from annotation_to_merge
    path assemblies from assemblies_to_merge.collect()

  output:
    file '*.gff' into annotation_to_quantify

  script:
    """
    stringtie --merge $assemblies -G $annotation -o merged.gff
    """
}

// Quantify genes
process quantify {
  tag 'StringTie'

  publishDir "$output/genes", mode: 'copy', saveAs: {
    filename ->
    if (filename.endsWith('.gff')) "$filename"
  }

  input:
    path annotation from annotation_to_quantify
    each path(map) from maps_to_quantify

  output:
    file '*.gff'
    file '*.tsv' into genes_to_format

  script:
    """
    outfile="$map"
    outfile="\${outfile%.*}"
    stringtie $map -e \
              -p $THREADS \
              -G $annotation \
              -o "\$outfile".gff \
              -A "\$outfile".tsv
    """
}

// Format results
process format {
  tag 'Bash'

  publishDir "$output/genes", mode: 'copy'

  input:
    path genes from genes_to_format.collect()

  output:
    file '*.tsv'

  script:
    """
    files=( $genes )

    echo "" > all_genes.tsv

    # Rename columns, sort rows, join file
    for (( f = \$(( \${#files[@]} -1 )); f >= 0; f-- )); do
      awk '
        BEGIN {
          FS = "\\t"
          OFS = FS
        }
        NR == 1 {
          gsub(/.*\\//, "", FILENAME)
          gsub(/\\.[^\\.]*\$/, "", FILENAME)
          for (i = 1; i <= NF; i++) {
            if (\$i == "Coverage") \$i = "Coverage " FILENAME
            else if (\$i == "FPKM") \$i = "FPKM " FILENAME
            else if (\$i == "TPM") \$i = "TPM " FILENAME
          }
          print
          next
        }
        {
          print | "sort -k 1,1"
        }
      ' "\${files[f]}" \\
      | join --header -a 1 -a 2 -t \$'\\t' -o auto - all_genes.tsv > temp.tsv
      mv temp.tsv all_genes.tsv
    done

    # Merge redundant columns
    awk '
      BEGIN {
        FS = "\\t"
        OFS = FS
      }
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if (!merge[i]) {
            merge[i] = i
            if (i == 1) printf "%s",\$i
            else printf "%s%s",OFS,\$i

            for (j = i+1; j <= NF; j++) {
              if (\$j == \$i) merge[j] = i
            }
          }
        }
        print ""
        next
      }
      {
        for (i = 1; i <= NF; i++) {
          if (merge[i] == i) {
            for (j = i; j <= NF; j++) {
              if (merge[j] == i && \$j != "") {
                if (i == 1) printf "%s",\$j
                else printf "%s%s",OFS,\$j
                break
              }
              if (j == NF && i != 1) printf "%s",OFS
            }
          }
        }
        print ""
      }
    ' all_genes.tsv > temp.tsv
    mv temp.tsv all_genes.tsv
    """
}
