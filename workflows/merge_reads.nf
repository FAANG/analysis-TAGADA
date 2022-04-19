// INCLUDE MODULES -------------------------------------------------------------

include {
  SAMTOOLS_merge_reads as SAMTOOLS_merge_reads_for_assembly
  SAMTOOLS_merge_reads as SAMTOOLS_merge_reads_for_quantification
} from '../modules/samtools.nf'

// DEFINE FUNCTIONS ------------------------------------------------------------

def group_reads(channel_aligned_reads, channel_metadata) {
  return channel_aligned_reads.map({ aligned ->
    [
      aligned['prefix'],
      aligned['bam'],
      aligned['length'],
      aligned['direction']
    ]
  }).join(
    channel_metadata.map({ row ->
      [row['prefix'], row['id'], row['columns']]
    })
  ).map({ pairing ->
    [pairing[4], pairing[1], pairing[2], pairing[3]]
  }).groupTuple().map({ group ->
    [
      'id': group[0],
      'bams': group[1],
      'lengths': group[2],
      'directions': group[3]
    ]
  })
}

def check_grouped_lengths(channel_grouped_reads, task) {
  channel_grouped_reads.filter({ group ->
    group['lengths'].subsequences().findAll({ lengths ->
      lengths.size() == 2
    }).collect({ lengths ->
      (lengths[0] - lengths[1]).abs()
    }).findAll({ difference ->
      difference > 10
    }).size() > 0
  }).toList().subscribe({ differing_lengths ->
    if (differing_lengths.size() > 0) {
      error = '\nIdentical read lengths required for ' + task + ':\n  '
      error += differing_lengths.collect({ group ->
        group['id'] + ':  ' + [
          group['bams'],
          group['lengths']
        ].transpose().collect({ it ->
          it[0].getName() + ' (' + it[1] + ')'
        }).join('  ')
      }).join('\n  ')
      exit(1, error)
    }
  })
}

def check_grouped_directions(channel_grouped_reads, task) {
  channel_grouped_reads.filter({ group ->
    group['directions'].toSet().size() > 1
  }).toList().subscribe({ differing_directions ->
    if (differing_directions.size() > 0) {
      error = '\nIdentical read directions required for ' + task + ':\n  '
      error += differing_directions.collect({ group ->
        group['id'] + ':  ' + [
          group['bams'],
          group['directions']
        ].transpose().collect({ it ->
          it[0].getName() + ' (' + it[1] + ')'
        }).join('  ')
      }).join('\n  ')
      exit(1, error)
    }
  })
}

// WORKFLOW --------------------------------------------------------------------

workflow MERGE_READS {

  take:

    channel_aligned_reads
    channel_assembly_metadata
    channel_quantification_metadata

  emit:

    channel_assembly_reads
    channel_quantification_reads

  main:

    // each [id, [bams], [lengths], [directions]]
    channel_assembly_reads =
      group_reads(
        channel_aligned_reads,
        channel_assembly_metadata
      )

    check_grouped_directions(
      channel_assembly_reads,
      'assembly'
    )

    // each [id, [bams], length, direction]
    channel_assembly_reads =
      channel_assembly_reads.map({ group ->
        [
          group['id'],
          group['bams'].size() == 1 &&
          (group['bams'][0].getName() =~ /^(.+?)\.bam$/)[0][1] == group['id'] ?
          group['bams'][0] : group['bams'],
          (group['lengths'].sum()/group['lengths'].size()).toInteger(),
          group['directions'][0]
        ]
      }).branch({ group ->
        grouped: group[1] instanceof List
        single: true
      })

    // each [id, [bams], length, direction] => [id, bam, length, direction]
    SAMTOOLS_merge_reads_for_assembly(
      channel_assembly_reads.grouped
    )

    // each [id, bam, length, direction]
    channel_assembly_reads =
      SAMTOOLS_merge_reads_for_assembly.out.mix(
        channel_assembly_reads.single
      ).map({ output ->
        [
          'id': output[0],
          'bam': output[1],
          'length': output[2],
          'direction': output[3]
        ]
      })

    // each [id, [bams], [lengths], [directions]]
    channel_quantification_reads =
      group_reads(
        channel_aligned_reads,
        channel_quantification_metadata
      )

    check_grouped_lengths(
      channel_quantification_reads,
      'quantification'
    )

    check_grouped_directions(
      channel_quantification_reads,
      'quantification'
    )

    // each [id, [bams], length, direction]
    channel_quantification_reads =
      channel_quantification_reads.map({ group ->
        [
          group['id'],
          group['bams'].size() == 1 &&
          (group['bams'][0].getName() =~ /^(.+?)\.bam$/)[0][1] == group['id'] ?
          group['bams'][0] : group['bams'],
          (group['lengths'].sum()/group['lengths'].size()).toInteger(),
          group['directions'][0]
        ]
      }).branch({ group ->
        grouped: group[1] instanceof List
        single: true
      })

    // each [id, [bams], length, direction] => [id, bam, length, direction]
    SAMTOOLS_merge_reads_for_quantification(
      channel_quantification_reads.grouped
    )

    // each [id, bam, length, direction]
    channel_quantification_reads =
      SAMTOOLS_merge_reads_for_quantification.out.mix(
        channel_quantification_reads.single
      ).map({ output ->
        [
          'id': output[0],
          'bam': output[1],
          'length': output[2],
          'direction': output[3]
        ]
      })
}
