process LOFREQ {
    tag "$out_basename"

    input:
    tuple val(out_basename), path(aligned_bam)
    path reference_fasta

    output:
    path "${out_basename}.vcf", emit: report_vcf

    script:
    """
    # Local copies because sometimes CWD is writeable but input is not (Nextflow handles staging, but scripts might assume write access to dir)
    cp ${aligned_bam} aligned.bam
    
    # Sanitize fasta IDs if needed (copying logic from WDL)
    python3 <<CODE
    import shutil
    import util.file
    with util.file.fastas_with_sanitized_ids("${reference_fasta}", use_tmp=True) as sanitized_fastas:
        shutil.copyfile(sanitized_fastas[0], 'reference.fasta')
    CODE

    if [ \$(grep -v '^>' reference.fasta | tr -d '\\nNn' | wc -c) == "0" ]; then
        touch "${out_basename}.vcf"
    else
        samtools faidx reference.fasta
        samtools index aligned.bam
        
        lofreq call \
            -f reference.fasta \
            -o "${out_basename}.vcf" \
            aligned.bam
    fi
    """
}
