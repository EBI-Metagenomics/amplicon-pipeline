
def dada2_input_preparation_function( concat_input, reads_qc, cutadapt_channel ) {

    // Monster combination of operators that does multiple things:
    // 1. Groups tuples by sample, splitting meta into variable regions and variable regions size with the expected number of elements
    // 1. being the number of variable regions identified for that sample (var_regions_size) so it doesn't wait until all items are emitted to continue
    // 1. This link explains how it works: https://training.nextflow.io/advanced/grouping/#grouping-using-submap

    // 2. Adds in QCd reads from both fastp and cutadapt, which get grouped using groupkey with the expected number of elements of 2 so that
    // 2. Nextflow doesn't wait until all items are emitted to continue

    // 3. Does a long map that reorganises the channel into:
    // 3.A Cleans up the meta dictionary as some duplications happen during the different groupings
    // 3.B Containing either the fastq from fastp if no primer trimming was done, or the fastq from cutadapt if primer trimming was done
    // 3.B This is done by checking if the cutadapt fastq is empty or not, if it's empty it's because no primers needed to be trimmed
    // 3.B How this is done also depends on if the sample is single or paired-end which added complexity

    dada2_input = concat_input
                  .map{ meta, std_primer, auto_primer ->                                                                    // 1.
                    key = groupKey(meta.subMap('id', 'single_end'), meta['var_regions_size'])                               // 1.
                    [ key, meta['var_region'], meta['var_regions_size'] ]                                                   // 1.
                   }                                                                                                        // 1.
                  .groupTuple(by: 0)                                                                                        // 1.
                  .join(reads_qc, by: 0)                                                                                    // 2.
                  .mix(cutadapt_channel)                                                                                    // 2.
                  .map { meta, var_region, var_regions_size, reads ->                                                       // 2.
                    [ groupKey(meta.subMap('id', 'single_end'), 2), var_region, var_regions_size, reads ]                   // 2.
                  }                                                                                                         // 2.
                  .groupTuple(by:0)                                                                                         // 2.
                  .map{ meta, var_region, var_regions_size, reads ->                                                        // 3.
                    final_var_region = var_region.unique()[0]                                                               // 3.A
                    final_var_region_size = var_regions_size[0][0]                                                          // 3.A
                    final_meta = meta + ['var_region': final_var_region, 'var_regions_size': final_var_region_size]         // 3.A

                    fastp_reads = reads[0]                                                                                  // 3.B
                    cutadapt_reads = reads[1]                                                                               // 3.B
                                                                                                                            // 3.B
                    cutadapt_read_size = 0                                                                                  // 3.B
                    if (meta.single_end == true){                                                                           // 3.B
                        cutadapt_read_size = cutadapt_reads.size()                                                          // 3.B
                    }                                                                                                       // 3.B
                    else{                                                                                                   // 3.B
                        cutadapt_read_size = cutadapt_reads[0].size()                                                       // 3.B
                    }                                                                                                       // 3.B
                                                                                                                            // 3.B
                    final_reads = ''                                                                                        // 3.B
                    if (meta['single_end'] ){                                                                               // 3.B
                        if (cutadapt_read_size > 0){                                                                        // 3.B
                            final_reads = cutadapt_reads                                                                    // 3.B
                        }                                                                                                   // 3.B
                        else{                                                                                               // 3.B
                            final_reads = fastp_reads                                                                       // 3.B
                        }                                                                                                   // 3.B
                    }                                                                                                       // 3.B
                    else{                                                                                                   // 3.B
                        if (cutadapt_read_size > 0){                                                                        // 3.B
                            final_reads = cutadapt_reads                                                                    // 3.B
                        }                                                                                                   // 3.B
                        else{                                                                                               // 3.B
                            final_reads = fastp_reads                                                                       // 3.B
                        }                                                                                                   // 3.B
                    }                                                                                                       // 3.B
                    [ final_meta, final_reads ]                                                                             // FINAL OUTPUT
                }

    return dada2_input
}