process GZIP_decompress {

  input:
    path(compressed)

  output:
    path(decompressed)

  shell:
    match = (compressed.getName() =~ /^(.+?)((\.tar)?(\.gz)?)$/)
    decompressed = match[0][1]
    compression = match[0][2]
    if (compression.startsWith('.tar'))
      '''
      tar -xf !{compressed} -C '!{decompressed}'
      '''
    else if (compression.startsWith('.gz'))
      '''
      gzip -c -d !{compressed} > '!{decompressed}'
      '''
    else
      '''
      '''
}
