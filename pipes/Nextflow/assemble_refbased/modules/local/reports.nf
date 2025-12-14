process ALIGNMENT_METRICS {
    tag "$out_basename"

    input:
    tuple val(out_basename), path(aligned_bam)
    path ref_fasta
    path primers_bed
    val min_coverage

    output:
    path "${out_basename}.raw_wgs_metrics.txt", emit: wgs_metrics
    path "${out_basename}.alignment_metrics.txt", emit: alignment_metrics
    path "${out_basename}.insert_size_metrics.txt", emit: insert_size_metrics
    path "${out_basename}.ampliconstats.txt", emit: amplicon_stats

    script:
    """
    # Picard needs dictionary
    python3 <<CODE
    import shutil
    import util.file
    with util.file.fastas_with_sanitized_ids("${ref_fasta}", use_tmp=True) as sanitized_fastas:
        shutil.copyfile(sanitized_fastas[0], 'reference.fasta')
    CODE
    
    # Picard create dict
    picard CreateSequenceDictionary -R reference.fasta

    # Collect WGS Metrics
    picard CollectRawWgsMetrics \
        -R reference.fasta \
        -I "${aligned_bam}" \
        -O picard_raw.raw_wgs_metrics.txt
    
    grep -v '#' picard_raw.raw_wgs_metrics.txt | grep . | head -2 > "${out_basename}.raw_wgs_metrics.txt"

    # Collect Alignment Metrics
    picard CollectAlignmentSummaryMetrics \
        -R reference.fasta \
        -I "${aligned_bam}" \
        -O picard_raw.alignment_metrics.txt
    
    grep -v '#' picard_raw.alignment_metrics.txt | grep . | head -4 > "${out_basename}.alignment_metrics.txt"

    # Collect Insert Size Metrics
    picard CollectInsertSizeMetrics \
        -I "${aligned_bam}" \
        -O picard_raw.insert_size_metrics.txt \
        -H picard_raw.insert_size_metrics.pdf \
        --INCLUDE_DUPLICATES true
    
    grep -v '#' picard_raw.insert_size_metrics.txt | grep . | head -2 > "${out_basename}.insert_size_metrics.txt"

    # Amplicon stats if bed provided
    touch "${out_basename}.ampliconstats.txt"
    if [ -f "${primers_bed}" ]; then
        cat "${primers_bed}" | sort -k 1,1 -k 4,4 -t \$'\\t' > primers-sorted.bed
        samtools ampliconstats -s -@ \$(nproc) \
            -d ${min_coverage} \
            -o "${out_basename}.ampliconstats.txt" \
            primers-sorted.bed "${aligned_bam}"
    fi
    """
}

process PLOT_COVERAGE {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(aligned_reads_bam)

    output:
    path "${sample_name}.coverage_plot.pdf", emit: coverage_plot
    path "${sample_name}.coverage_plot.txt", emit: coverage_tsv
    path "mean_coverage", emit: mean_coverage

    script:
    """
    samtools view -c ${aligned_reads_bam} > reads_aligned
    if [ "\$(cat reads_aligned)" != "0" ]; then
        samtools index -@ \$(nproc) "${aligned_reads_bam}"
        
        reports.py plot_coverage \
            "${aligned_reads_bam}" \
            "${sample_name}.coverage_plot.pdf" \
            --outSummary "${sample_name}.coverage_plot.txt" \
            --plotFormat pdf \
            --plotTitle "${sample_name} coverage plot" \
            --loglevel=DEBUG
    else
        touch "${sample_name}.coverage_plot.pdf" "${sample_name}.coverage_plot.txt"
    fi

    # Mean coverage calculation
    samtools view -H ${aligned_reads_bam} | perl -n -e'/^@SQ.*LN:(\\d+)/ && print "\$1\\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" > assembly_length
    samtools view ${aligned_reads_bam} | cut -f10 | tr -d '\\n' | wc -c > bases_aligned
    python -c "print(float(\$(cat bases_aligned))/\$(cat assembly_length)) if \$(cat assembly_length)>0 else print(0)" > mean_coverage
    """
}
