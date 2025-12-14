process MERGE_BAMS {
    tag "$out_basename"

    input:
    path in_bams
    val sample_name
    val out_basename

    output:
    path "${out_basename}.bam", emit: out_bam
    path "${out_basename}.bam.flagstat.txt", emit: flagstat

    script:
    """
    mem_in_mb=\$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    # If list has >1 file, merge. Else just copy (or symlink logic handled by python script)
    # Nextflow 'in_bams' will be a list of files if multiple are passed.
    
    # We need to construct the command arguments correctly. 
    # Nextflow automatically handles input file staging.
    
    read_utils.py merge_bams ${in_bams} merged.bam --JVMmemory="\$mem_in_mb"m --loglevel DEBUG

    # Remap SM values
    if [ -n "${sample_name}" ]; then
        samtools view -H merged.bam | perl -n -e'/SM:(\\S+)/ && print "SM\\t\$1\\t${sample_name}\\n"' | sort | uniq >> reheader_table.txt
    fi

    if [ -s reheader_table.txt ]; then
        read_utils.py reheader_bam merged.bam reheader_table.txt "${out_basename}.bam" --loglevel DEBUG
    else
        mv merged.bam "${out_basename}.bam"
    fi

    samtools flagstat "${out_basename}.bam" > "${out_basename}.bam.flagstat.txt"
    """
}
