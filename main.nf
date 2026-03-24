#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import nextflow.util.MemoryUnit
import nextflow.util.Duration

// CHECK PARAMETERS ------------------------------------------------------------

// Validation
// ---------------------------
def validateConfig(cfg) {
    def error = ''

    // Required inputs
    if (!cfg.output)      error += '\nNo --output provided'
    if (!cfg.reads)       error += '\nNo --reads provided'
    if (!cfg.genome)      error += '\nNo --genome provided'
    if (!cfg.annotation)  error += '\nNo --annotation provided'

    // Deprecated parameters
    if (params.containsKey('merge'))       error += '\nReplace deprecated --merge with --assemble_by and --quantify_by'
    if (params.containsKey('feelnc_args')) error += '\nReplace deprecated --feelnc_args with --feelnc_codpot_args'

    // Logical conflicts
    if (cfg.assemble_by && cfg.skip_assembly) {
        error += '\nCannot use --assemble_by and --skip_assembly together'
    }

    if ((cfg.assemble_by || cfg.quantify_by) && !cfg.metadata) {
        error += '\nNo --metadata provided for --assemble_by or --quantify_by'
    }

    // Allowed values
    if (!['tmerge', 'stringtie'].contains(cfg.coalesce_transcripts_with)) {
        error += '\nInvalid --coalesce_transcripts_with must be "tmerge" or "stringtie"'
    }

    // Warn about ignored parameters
    if (cfg.skip_filtering &&
        (cfg.min_transcript_tpm || cfg.min_monoexonic_occurrence ||
        cfg.min_transcript_occurrence || cfg.min_monoexonic_tpm)) {
        log.warn "Ignoring filtering parameters because --skip_filtering has been set"
    }

    if (error) exit(1, error)
}

// INCLUDE WORKFLOWS -----------------------------------------------------------

include { CREATE_CHANNELS } from './workflows/create_channels.nf'
include { MERGE_READS } from './workflows/merge_reads.nf'

// INCLUDE MODULES -------------------------------------------------------------

include {
  SAMTOOLS_sort_reads
  SAMTOOLS_control_flags
  SAMTOOLS_control_alignments
  SAMTOOLS_control_contigs
} from './modules/samtools.nf'

include {
  TRIMGALORE_trim_adapters
} from './modules/trimgalore.nf'

include {
  FASTQC_control_reads
} from './modules/fastqc.nf'

include {
  TAGADA_estimate_reads
  TAGADA_filter_transcripts
  TAGADA_format_quantifications
  TAGADA_cluster_expression
  TAGADA_control_expression
  TAGADA_control_annotation
} from './modules/tagada.nf'

include {
  STAR_index_genome
  STAR_align_reads
} from './modules/star.nf'

include {
  BEDTOOLS_compute_coverage
} from './modules/bedtools.nf'

include {
  STRINGTIE_assemble_transcripts
  STRINGTIE_coalesce_transcripts
  STRINGTIE_quantify_expression
} from './modules/stringtie.nf'

include {
  TMERGE_coalesce_transcripts
} from './modules/tmerge.nf'

include {
  FEELNC_classify_transcripts
} from './modules/feelnc.nf'

include {
  FEATURECOUNTS_control_exons
} from './modules/featurecounts.nf'

include {
  MULTIQC_generate_report
} from './modules/multiqc.nf'

// WORKFLOW --------------------------------------------------------------------

workflow {

  def cfg = [
    output                     : params.get('output', ''),
    reads                      : params.get('reads', ''),
    genome                     : params.get('genome', ''),
    annotation                 : params.get('annotation', ''),
    index                      : params.get('index', ''),
    metadata                   : params.get('metadata', ''),

    assemble_by                : (params.get('assemble_by', '')).split(',').collect{ s -> s.trim() }.findAll{ s -> s },
    quantify_by                : (params.get('quantify_by', '')).split(',').collect{ s -> s.trim() }.findAll{ s -> s },

    min_transcript_occurrence  : params.get('min_transcript_occurrence', ''),
    min_monoexonic_occurrence  : params.get('min_monoexonic_occurrence', ''),
    min_transcript_tpm         : params.get('min_transcript_tpm', ''),
    min_monoexonic_tpm         : params.get('min_monoexonic_tpm', ''),

    coalesce_transcripts_with  : params.get('coalesce_transcripts_with', 'tmerge'),
    tmerge_args                : params.get('tmerge_args', ''),
    feelnc_filter_args         : params.get('feelnc_filter_args', ''),
    feelnc_codpot_args         : params.get('feelnc_codpot_args', ''),
    feelnc_classifier_args     : params.get('feelnc_classifier_args', ''),

    skip_assembly              : params.get('skip_assembly', false),
    skip_filtering             : params.get('skip_filtering', false),
    skip_lnc_detection         : params.get('skip_lnc_detection', false),

    transcript_biotype_field   : params.get('transcript_biotype_field', 'transcript_biotype'),

    max_cpus                   : params.get('max_cpus', 16) as int,
    max_memory                 : params.get('max_memory', 64.GB) as MemoryUnit,
    max_time                   : params.get('max_time', 24.h) as Duration
  ]

  validateConfig(cfg)
    
  CREATE_CHANNELS(cfg)

  // each [prefix, [fastqs]]
  channel_raw_reads =
    CREATE_CHANNELS.out.channel_raw_reads

  // each [prefix, bam]
  channel_aligned_reads =
    CREATE_CHANNELS.out.channel_aligned_reads

  // one index
  channel_index =
    CREATE_CHANNELS.out.channel_index

  // one genome
  channel_genome =
    CREATE_CHANNELS.out.channel_genome

  // one annotation
  channel_reference_annotation =
    CREATE_CHANNELS.out.channel_reference_annotation

  // each [prefix, id, [column: value, column: value]]
  channel_assembly_metadata =
    CREATE_CHANNELS.out.channel_assembly_metadata

  // each [prefix, id, [column: value, column: value]]
  channel_quantification_metadata =
    CREATE_CHANNELS.out.channel_quantification_metadata

  // each [reports]
  channel_reports =
    Channel.empty()
  // SORT ALIGNED READS --------------------------------------------------------

  // each [prefix, bam] => [prefix, bam]
  SAMTOOLS_sort_reads(channel_aligned_reads)

  // each [prefix, bam]
  channel_aligned_reads =
    SAMTOOLS_sort_reads.out.map({ output ->
      ['prefix': output[0], 'bam': output[1]]
    })

  // TRIM ADAPTERS FROM RAW READS ----------------------------------------------

  // each [prefix, [fastqs]] => [prefix, [fastqs]] & [reports]
  TRIMGALORE_trim_adapters(channel_raw_reads)

  // each [prefix, [fastqs]]
  channel_trimmed_reads =
    TRIMGALORE_trim_adapters.out.results.map({ output ->
      [
        'prefix': output[0],
        'fastqs': output[1] instanceof List ? output[1].sort({ fastq ->
          (fastq =~ /[\._ ][Rr]?([12]) \(.+?\)\.fq\.gz$/).with({ match ->
            match ? match[0][1] : fastq
          })
        }) : [output[1]]
      ]
    })

  // each [reports]
  channel_reports =
    channel_reports.mix(TRIMGALORE_trim_adapters.out.reports)

  // CONTROL RAW AND TRIMMED READS QUALITY -------------------------------------

  // each [prefix, fastq, suffix] => report
  FASTQC_control_reads(
    Channel.empty().mix(
      channel_raw_reads.map({ raw -> raw + ['suffix': 'raw'] }),
      channel_trimmed_reads.map({ trimmed -> trimmed + ['suffix': 'trimmed'] })
    ).map({ reads ->
      (0..reads['fastqs'].size()-1).inject([], { result, n ->
        result + [
          [
            'prefix': reads['prefix'],
            'fastq': reads['fastqs'][n],
            'suffix': reads['suffix']
          ]
        ]
      })
    }).flatten()
  )

  // each [reports]
  channel_reports =
    channel_reports.mix(FASTQC_control_reads.out.reports)

  // INDEX GENOME --------------------------------------------------------------

  if (!params.index) {

    // one genome & annotation & [fastqs] => index
    STAR_index_genome(
      channel_genome,
      channel_reference_annotation,
      channel_trimmed_reads.map({ trimmed ->
        trimmed['fastqs']
      }).flatten().collect()
    )

    // one index
    channel_index =
      STAR_index_genome.out
  }

  // ALIGN TRIMMED READS -------------------------------------------------------

  // each [prefix, [fastqs], index] => [prefix, bam] & report
  STAR_align_reads(
    channel_trimmed_reads.map({ trimmed ->
      [trimmed['prefix'], trimmed['fastqs']]
    }).combine(channel_index)
  )

  // each [prefix, bam]
  channel_aligned_reads =
    channel_aligned_reads.mix(
      STAR_align_reads.out.results.map({ output ->
        ['prefix': output[0], 'bam': output[1]]
      })
    )

  // each [reports]
  channel_reports =
    channel_reports.mix(STAR_align_reads.out.reports)

  // CONTROL ALIGNMENT QUALITY -------------------------------------------------

  // each [prefix, bam] => report
  SAMTOOLS_control_flags(channel_aligned_reads)
  SAMTOOLS_control_alignments(channel_aligned_reads)
  SAMTOOLS_control_contigs(channel_aligned_reads)

  // each [reports]
  channel_reports =
    channel_reports.mix(
      SAMTOOLS_control_flags.out,
      SAMTOOLS_control_alignments.out,
      SAMTOOLS_control_contigs.out
    )

  // GET READ LENGTHS AND DIRECTIONS -------------------------------------------

  // each [prefix, bam, annotation] => [prefix, bam, length, direction]
  TAGADA_estimate_reads(
    channel_aligned_reads.map({ aligned ->
      [aligned['prefix'], aligned['bam']]
    }).combine(channel_reference_annotation)
  )

  // each [prefix, bam, length, direction]
  channel_aligned_reads =
    TAGADA_estimate_reads.out.map({ output ->
      [
        'prefix': output[0],
        'bam': output[1],
        'length': output[2].toInteger(),
        'direction': output[3]
      ]
    })

  channel_aligned_reads.toSortedList({ a, b ->
    a['prefix'] <=> b['prefix']
  }).view({ aligned_reads ->
    if (aligned_reads.size() == 0) return ''
    s = aligned_reads.size() > 1 ? 's' : ''
    info = '\nEstimated read length' + s + ' and direction' + s + ':\n  '
    info += aligned_reads.collect({ aligned ->
      aligned['prefix'] +
      ' (' + aligned['length'] + ')' +
      ' (' + aligned['direction'] + ')'
    }).join('\n  ')
    return info
  })

  // COMPUTE COVERAGE ----------------------------------------------------------

  // each [prefix, bam, direction]
  BEDTOOLS_compute_coverage(
    channel_aligned_reads.map({ aligned ->
      [aligned['prefix'], aligned['bam'], aligned['direction']]
    })
  )

  // MERGE ALIGNED READS -------------------------------------------------------

  MERGE_READS(
    channel_aligned_reads,
    channel_assembly_metadata,
    channel_quantification_metadata
  )

  // each [id, bam, length, direction]
  channel_assembly_reads =
    MERGE_READS.out.channel_assembly_reads

  // each [id, bam, length, direction]
  channel_quantification_reads =
    MERGE_READS.out.channel_quantification_reads

  // ASSEMBLE TRANSCRIPTS ------------------------------------------------------

  if (!cfg.skip_assembly) {

    // each [id, bam, direction, annotation] => assembly
    STRINGTIE_assemble_transcripts(
      channel_assembly_reads.map({ group ->
        [group['id'], group['bam'], group['direction']]
      }).combine(channel_reference_annotation)
    )

    STRINGTIE_assemble_transcripts.out.view()
    ch_assembly = STRINGTIE_assemble_transcripts.out.collect()
    ch_assembly.view()
  

    if (!cfg.skip_filtering){
      TAGADA_filter_transcripts(ch_assembly, cfg)
      ch_filtered_assembly = TAGADA_filter_transcripts.out.results
    } else {
      ch_filtered_assembly = ch_assembly
    }

    if (cfg.coalesce_transcripts_with == 'stringtie') {

      // one [assemblies] & annotation => annotation
      STRINGTIE_coalesce_transcripts(
        ch_filtered_assembly,
        channel_reference_annotation,
        cfg
      )

      // one annotation
      channel_novel_annotation =
        STRINGTIE_coalesce_transcripts.out

    } else {

      // one [assemblies] & annotation => annotation
      TMERGE_coalesce_transcripts(
        ch_filtered_assembly,
        channel_reference_annotation,
        cfg
      )

      // one annotation
      channel_novel_annotation =
        TMERGE_coalesce_transcripts.out
    }

  } else {

    // none annotation
    channel_novel_annotation =
      Channel.empty()
  }

  // DETECT LONG NON-CODING TRANSCRIPTS ----------------------------------------

  if (!cfg.skip_lnc_detection && !cfg.skip_assembly) {

    // one genome & reference_annotation & novel_annotation => [reports]
    FEELNC_classify_transcripts(
      channel_genome,
      channel_reference_annotation,
      channel_novel_annotation,
      cfg
    )

    // each [reports]
    channel_reports =
      channel_reports.mix(FEELNC_classify_transcripts.out.reports)
  }

  // QUANTIFY GENES AND TRANSCRIPTS --------------------------------------------

  // each [id, bam, length, direction, annotation, type] => [quantifications]
  STRINGTIE_quantify_expression(
    channel_quantification_reads.map({ group ->
      [group['id'], group['bam'], group['length'], group['direction']]
    }).combine(
      Channel.empty().mix(
        channel_reference_annotation.combine(Channel.of('reference')),
        channel_novel_annotation.combine(Channel.of('novel'))
      )
    )
  )

  // each [quantifications]
  channel_quantifications =
    STRINGTIE_quantify_expression.out.flatten().map({ quantification ->
      [
        (quantification.getName() =~ /\.([^\.]+)\.tsv$/)[0][1],
        quantification
      ]
    }).groupTuple().map({ quantifications -> quantifications[1] })

  // each [quantifications] => quantification
  TAGADA_format_quantifications(channel_quantifications)

  // each quantification
  channel_quantifications =
    TAGADA_format_quantifications.out

  // CONTROL EXPRESSION --------------------------------------------------------

  // one genes_tpm_quantification & quantification_metadata
  TAGADA_cluster_expression(
    channel_quantifications.filter({ quantification ->
      quantification.getName() == 'reference_genes_TPM.tsv'
    }),
    channel_quantification_metadata.collect()
  )

  // one genes_tpm_quantification & genes_counts_quantification => [reports]
  TAGADA_control_expression(
    channel_quantifications.filter({ quantification ->
      quantification.getName() == 'reference_genes_TPM.tsv'
    }),
    channel_quantifications.filter({ quantification ->
      quantification.getName() == 'reference_genes_counts.tsv'
    })
  )

  // each [reports]
  channel_reports =
    channel_reports.mix(TAGADA_control_expression.out.reports)

  // CONTROL ANNOTATION --------------------------------------------------------

  if (!cfg.skip_assembly) {

    // one reference_annotation & novel_annotation &
    // transcripts_tpm_quantification & genes_tpm_quantification => [reports]
    TAGADA_control_annotation(
      channel_reference_annotation,
      channel_novel_annotation,
      channel_quantifications.filter({ quantification ->
        quantification.getName() == 'novel_transcripts_TPM.tsv'
      }),
      channel_quantifications.filter({ quantification ->
        quantification.getName() == 'novel_genes_TPM.tsv'
      })
    )

    // each [reports]
    channel_reports =
      channel_reports.mix(TAGADA_control_annotation.out.reports)
  }

  // CONTROL EXONS -------------------------------------------------------------

  // each [bams] & [annotation, type] => report
  FEATURECOUNTS_control_exons(
    channel_aligned_reads.map({ aligned -> aligned['bam'] }).collect(),
    Channel.empty().mix(
      channel_reference_annotation.combine(Channel.of('reference')),
      channel_novel_annotation.combine(Channel.of('novel'))
    )
  )

  // each [reports]
  channel_reports =
    channel_reports.mix(FEATURECOUNTS_control_exons.out.reports)

  // GENERATE HTML REPORT ------------------------------------------------------

  // one [reports] & config & custom_images & metadata
  MULTIQC_generate_report(
    channel_reports.flatten().collect(),
    Channel.fromPath(
      projectDir + '/assets/multiqc/config.yaml',
      checkIfExists: true
    ),
    Channel.fromPath(
      projectDir + '/assets/multiqc/custom_images.tsv',
      checkIfExists: true
    ),
    cfg.metadata ? Channel.fromPath(cfg.metadata) : []
  )
}
