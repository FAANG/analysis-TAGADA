// INCLUDE MODULES -------------------------------------------------------------

include {
  GZIP_decompress as GZIP_decompress_index
  GZIP_decompress as GZIP_decompress_genome
  GZIP_decompress as GZIP_decompress_annotation
} from '../modules/gzip.nf'

// DEFINE FUNCTIONS ------------------------------------------------------------

def set_metadata_groups(channel_metadata, factors) {
  return channel_metadata.map({ row ->
    [
      'prefix': row['prefix'],
      'columns': row['columns'].subMap(factors).findAll({ it -> it.value })
    ]
  }).map({ row ->
    [
      'prefix': row['prefix'],
      'id': factors.size() > 0 ? (
        row['columns'].size() == factors.size() ?
        row['columns'].values().join('_') : '?'
      ) : row['prefix'],
      'columns': row['columns']
    ]
  })
}

def log_metadata_groups(channel_metadata, task) {
  channel_metadata.map({ row ->
    [row['id'], row['prefix']]
  }).groupTuple().map({ row ->
    ['id': row[0], 'prefixes': row[1]]
  }).toSortedList({ a, b -> a['id'] <=> b['id'] }).view({ rows ->
    if (rows.size() == 0) return ''
    s = rows.size() > 1 ? 's' : ''
    info = '\nThe following ' + task + ' group' + s + ' will be created:\n  '
    info += rows.collect({ row ->
      row['id'] + ':  ' + row['prefixes'].join('  ')
    }).join('\n  ')
    return info
  })
}

// WORKFLOW --------------------------------------------------------------------

workflow CREATE_CHANNELS {

  emit:

    channel_raw_reads
    channel_aligned_reads
    channel_index
    channel_genome
    channel_reference_annotation
    channel_assembly_metadata
    channel_quantification_metadata

  main:

    // SPLIT READS INTO R1 R2 SINGLE ALIGNED -----------------------------------

    channel_reads =
      Channel.fromPath(
        params.reads,
        checkIfExists: true
      )

    if (params.reads.endsWith('.txt')) {
      channel_reads = channel_reads.splitText().map({ line ->
        file(line.strip(), checkIfExists: true)
      })
    }

    channel_reads =
      channel_reads.map({ path ->
        filename = path.getName()
        return (
          filename =~ /^(.+?)(?:[\._ ][Rr]?([12]))?(?:\.(fastq|fq|gz|bam))+$/
        ).with({ it ->
          matches() ? [
            'prefix': it[0][1],
            'R': it[0][2] ? it[0][2].toInteger() : null,
            'bam': it[0][3] == 'bam',
            'path': path,
            'filename': filename
          ] : filename
        })
      }).branch({ it ->
        invalid: it instanceof String
        aligned: it['bam']
        r1: it['R'] == 1
        r2: it['R'] == 2
        single: true
      })

    // CHECK FILE NAMES --------------------------------------------------------

    channel_reads.invalid.toList().subscribe({ invalid ->
      if (invalid.size() > 0) {
        s = invalid.size() > 1 ? 's' : ''
        error = '\nWrong format for file name' + s + ':\n  '
        error += invalid.join('\n  ')
        exit(1, error)
      }
    })

    // CHECK UNPAIRED AND DUPLICATED READS -------------------------------------

    channel_incoherent_reads =
      Channel.empty().concat(
        channel_reads.r1.map({ r1 ->
          [r1['prefix'], r1['filename'], 'r1']
        }),
        channel_reads.r2.map({ r2 ->
          [r2['prefix'], r2['filename'], 'r2']
        }),
        channel_reads.single.map({ single ->
          [single['prefix'], single['filename'], 'single']
        }),
        channel_reads.aligned.map({ aligned ->
          [aligned['prefix'], aligned['filename'], 'aligned']
        })
      ).groupTuple().filter({ pairing ->
        pairing[2] != ['r1', 'r2'] &&
        pairing[2] != ['single'] &&
        pairing[2] != ['aligned']
      }).map({ pairing -> pairing[1] })

    channel_incoherent_reads.filter({ filenames ->
      filenames.size() == 1
    }).map({ filenames ->
      filenames[0]
    }).toList().subscribe({ unpaired ->
      if (unpaired.size() > 0) {
        s = unpaired.size() > 1 ? 's' : ''
        error = '\nNo pair' + s + ' for file' + s + ':\n  '
        error += unpaired.join('\n  ')
        exit(1, error)
      }
    })

    channel_incoherent_reads.filter({ filenames ->
      filenames.size() > 1
    }).toList().subscribe({ duplicated ->
      if (duplicated.size() > 0) {
        error = '\nDuplicated inputs:\n  '
        error += duplicated.collect({ filenames ->
          filenames.join('  ')
        }).join('\n  ')
        exit(1, error)
      }
    })

    // CREATE READS CHANNELS ---------------------------------------------------

    channel_raw_reads =
      channel_reads.r1.map({ r1 ->
        [r1['prefix'], r1['path']]
      }).join(
        channel_reads.r2.map({ r2 ->
          [r2['prefix'], r2['path']]
        })
      ).map({ paired ->
        ['prefix': paired[0], 'fastqs': [paired[1], paired[2]]]
      }).mix(
        channel_reads.single.map({ single ->
          ['prefix': single['prefix'], 'fastqs': [single['path']]]
        })
      )

    channel_aligned_reads =
      channel_reads.aligned.map({ aligned ->
        ['prefix': aligned['prefix'], 'bam': aligned['path']]
      })

    // CREATE METADATA CHANNEL -------------------------------------------------

    if (!params.metadata) {

      channel_metadata =
        Channel.empty().mix(
          channel_raw_reads,
          channel_aligned_reads
        ).map({ row ->
          ['prefix': row['prefix'], 'columns': [:]]
        })

    } else {

      channel_metadata =
        Channel.fromPath(
          params.metadata,
          checkIfExists: true
        ).splitCsv(header: true, sep: '\t').map({ row ->
          [
            'prefix': row.values()[0],
            'columns': row.subMap(row.keySet().collect()[1..-1])
          ]
        })

      // CHECK DUPLICATED METADATA ROWS ----------------------------------------

      channel_metadata.map({ row ->
        [row['prefix'], null]
      }).groupTuple().filter({ row ->
        row[1].size() > 1
      }).map({ row -> row[0] }).toList().subscribe({ duplicated_rows ->
        if (duplicated_rows.size() > 0) {
          s = duplicated_rows.size() > 1 ? 's' : ''
          error = '\nDuplicated metadata row' + s + ':\n  '
          error += duplicated_rows.join('\n  ')
          exit(1, error)
        }
      })

      // KEEP METADATA ROWS MATCHING INPUT READS -------------------------------

      channel_metadata =
        Channel.empty().concat(
          channel_metadata.map({ row ->
            [row['prefix'], false, row['columns']]
          }),
          channel_raw_reads.mix(channel_aligned_reads).map({ reads ->
            [reads['prefix'], true]
          })
        ).groupTuple().filter({ row ->
          row[1].findAll().size() > 0
        }).map({ row ->
          ['prefix': row[0], 'columns': row.size() == 3 ? row[2][0] : null]
        })

      // CHECK MISSING METADATA ROWS -------------------------------------------

      channel_metadata.filter({ row ->
        !row['columns']
      }).map({ row ->
        row['prefix']
      }).flatten().toList().subscribe({ missing_rows ->
        if (missing_rows.size() > 0) {
          s = missing_rows.size() > 1 ? 's' : ''
          error = '\nMissing metadata row' + s + ':\n  '
          error += missing_rows.join('\n  ')
          exit(1, error)
        }
      })

      channel_metadata =
        channel_metadata.filter({ row -> row['columns'] })

      // CHECK MISSING METADATA COLUMNS ----------------------------------------

      channel_metadata.first().map({ row ->
        row['columns'].keySet().collect()
      }).map({ column_names ->
        (
          params.assemble_by + params.quantify_by
        ).sort().unique().findAll({ factor ->
          !(factor in column_names)
        })
      }).subscribe({ missing_columns ->
        if (missing_columns.size() > 0) {
          s = missing_columns.size() > 1 ? 's' : ''
          error = '\nMissing metadata column' + s + ':\n  '
          error += missing_columns.join('\n  ')
          exit(1, error)
        }
      })

      // CHECK MISSING METADATA VALUES -----------------------------------------

      channel_metadata.map({ row ->
        [
          'prefix': row['prefix'],
          'empty': (
            params.assemble_by + params.quantify_by
          ).sort().unique().findAll({ factor ->
            factor in row['columns'].keySet() && !row['columns'][factor]
          })
        ]
      }).filter({ row ->
        row['empty'].size() > 0
      }).toList().subscribe({ invalid_rows ->
        if (invalid_rows.size() > 0) {
          s = invalid_rows.collect({ row ->
            row['empty']
          }).flatten().toSet().size() > 1 ? 's' : ''
          s = invalid_rows.size() > 1 || s ? 's' : ''
          error = '\nMissing metadata value' + s + ':\n  '
          error += invalid_rows.collect({ row ->
            row['prefix'] + '  ' + row['empty'].join('  ')
          }).join('\n  ')
          exit(1, error)
        }
      })
    }

    // DETERMINE MERGE GROUPS --------------------------------------------------

    channel_assembly_metadata =
      set_metadata_groups(channel_metadata, params.assemble_by)

    if (params.assemble_by) {
      log_metadata_groups(channel_assembly_metadata, 'assembly')
    }

    channel_quantification_metadata =
      set_metadata_groups(channel_metadata, params.quantify_by)

    if (params.quantify_by) {
      log_metadata_groups(channel_quantification_metadata, 'quantification')
    }

    // DECOMPRESS INDEX --------------------------------------------------------

    if (params.index) {

      channel_index =
        Channel.fromPath(
          params.index,
          type: 'dir',
          checkIfExists: true
        )

      if (params.index.endsWith('.gz') || params.index.endsWith('.tar')) {

        // index => index
        GZIP_decompress_index(channel_index)

        channel_index =
          GZIP_decompress_index.out
      }

    } else {

      channel_index =
        Channel.empty()
    }

    // DECOMPRESS GENOME -------------------------------------------------------

    channel_genome =
      Channel.fromPath(
        params.genome,
        checkIfExists: true
      )

    if (params.genome.endsWith('.gz')) {

      // genome => genome
      GZIP_decompress_genome(channel_genome)

      channel_genome =
        GZIP_decompress_genome.out
    }

    // DECOMPRESS ANNOTATION ---------------------------------------------------

    channel_reference_annotation =
      Channel.fromPath(
        params.annotation,
        checkIfExists: true
      )

    if (params.annotation.endsWith('.gz')) {

      // annotation => annotation
      GZIP_decompress_annotation(channel_reference_annotation)

      channel_reference_annotation =
        GZIP_decompress_annotation.out
    }
}
