process ALIGN_READS {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(reads_bam)
    path reference_fasta
    val aligner
    path novocraft_license
    val skip_mark_dupes
    val aligner_options

    output:
    tuple val(sample_name), path("${sample_name}.mapped.bam"), emit: aligned_bam
    path "${sample_name}.all.bam", emit: all_bam
    path "${sample_name}.all.bam.flagstat.txt", emit: flagstat
    path "${sample_name}.mapped_fastqc.html", emit: fastqc_html
    path "${sample_name}.mapped_fastqc.zip", emit: fastqc_zip
    
    script:
    def skip_mark_dupes_arg = skip_mark_dupes ? '--skipMarkDupes' : ''
    def license_arg = novocraft_license ? "--NOVOALIGN_LICENSE_PATH=${novocraft_license}" : ''
    def options_arg = aligner_options ? "--aligner_options \"${aligner_options}\"" : ''
    """
    cp ${reference_fasta} assembly.fasta
    read_utils.py index_fasta_picard assembly.fasta
    read_utils.py index_fasta_samtools assembly.fasta

    # Calculate memory for JVM (approx 90% of available)
    mem_in_mb=\$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    read_utils.py align_and_fix \
        "${reads_bam}" \
        assembly.fasta \
        --outBamAll "${sample_name}.all.bam" \
        --outBamFiltered "${sample_name}.mapped.bam" \
        --aligner ${aligner} \
        ${options_arg} \
        ${skip_mark_dupes_arg} \
        --JVMmemory "\$mem_in_mb"m \
        ${license_arg} \
        --loglevel=DEBUG

    samtools view -h -F 260 "${sample_name}.all.bam" | samtools flagstat - > "${sample_name}.all.bam.flagstat.txt"
    reports.py fastqc "${sample_name}.mapped.bam" "${sample_name}.mapped_fastqc.html" --out_zip "${sample_name}.mapped_fastqc.zip"
    """
}

process IVAR_TRIM {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(aligned_bam)
    path trim_coords_bed

    output:
    tuple val(sample_name), path("${sample_name}.trimmed.bam"), emit: trimmed_bam
    path "IVAR_TRIM_PCT", emit: trim_pct
    path "IVAR_TRIM_COUNT", emit: trim_count

    script:
    """
    if [ -f "${trim_coords_bed}" ]; then
        ivar trim -e \
            -b ${trim_coords_bed} \
            -m 30 -s 4 -q 20 \
            -i ${aligned_bam} -p trim | tee IVAR_OUT
        samtools sort -@ \$(nproc) -m 1000M -o "${sample_name}.trimmed.bam" trim.bam
    else
        cp ${aligned_bam} "${sample_name}.trimmed.bam"
        echo "Trimmed primers from 0% (0) of reads." > IVAR_OUT
    fi
    
    PCT=\$(grep "Trimmed primers from" IVAR_OUT | perl -lape 's/Trimmed primers from (\\S+)%.*/\$1/')
    if [[ \$PCT == -* ]]; then echo 0; else echo \$PCT; fi > IVAR_TRIM_PCT
    grep "Trimmed primers from" IVAR_OUT | perl -lape 's/Trimmed primers from \\S+% \\((\\d+)\\).*/\$1/' > IVAR_TRIM_COUNT
    """
}

process REFINE_ASSEMBLY {
    tag "$out_basename"

    input:
    tuple val(out_basename), path(reads_aligned_bam)
    path reference_fasta
    val min_coverage
    val major_cutoff

    output:
    path "${out_basename}.fasta", emit: assembly_fasta
    path "${out_basename}.sites.vcf.gz", emit: sites_vcf
    path "assembly_length", emit: len
    path "assembly_length_unambiguous", emit: len_unambig

    script:
    """
    # Calculate memory for JVM
    mem_in_mb=\$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    ln -s ${reference_fasta} assembly.fasta
    ln -s ${reads_aligned_bam} temp_markdup.bam
    samtools index -@ \$(nproc) temp_markdup.bam 

    assembly.py refine_assembly \
        assembly.fasta \
        temp_markdup.bam \
        refined.fasta \
        --already_realigned_bam=temp_markdup.bam \
        --outVcf "${out_basename}.sites.vcf.gz" \
        --min_coverage ${min_coverage} \
        --major_cutoff ${major_cutoff} \
        --JVMmemory "\$mem_in_mb"m \
        --loglevel=DEBUG

    if [ -s refined.fasta ]; then
        fasta-trim-terminal-ambigs.pl --3rules \
        --ten 4 \
        --fifty 15 \
        --maxfrac 0.9 \
        refined.fasta > trimmed.fasta
    else
        touch trimmed.fasta
    fi

    file_utils.py rename_fasta_sequences \
        trimmed.fasta "${out_basename}.fasta" "${out_basename}"

    grep -v '^>' trimmed.fasta | tr -d '\\n' | wc -c | tee assembly_length
    grep -v '^>' trimmed.fasta | tr -d '\\nNn' | wc -c | tee assembly_length_unambiguous
    """
}

process RUN_DISCORDANCE {
    tag "$out_basename"

    input:
    tuple val(out_basename), path(reads_aligned_bam)
    path reference_fasta
    val min_coverage

    output:
    path "${out_basename}.discordant.vcf", emit: discordant_vcf
    path "num_concordant", emit: num_concordant
    path "num_discordant_snps", emit: num_discordant_snps

    script:
    """
    read_utils.py --version | tee VERSION

    python3 <<CODE
    import tools.samtools
    header = tools.samtools.SamtoolsTool().getHeader("${reads_aligned_bam}")
    rgids = [[x[3:] for x in h if x.startswith('ID:')][0] for h in header if h[0]=='@RG']
    with open('readgroups.txt', 'wt') as outf:
        for rg in rgids:
            outf.write(rg+'\\t'+rg+'\\n')
    CODE

    if [ \$(grep -v '^>' "${reference_fasta}" | tr -d '\\nNn' | wc -c) != "0" ]; then
        samtools faidx "${reference_fasta}"
    fi

    bcftools mpileup \
        -G readgroups.txt -d 10000 -a "FORMAT/DP,FORMAT/AD" \
        -q 1 -m 2 -Ou \
        -f "${reference_fasta}" "${reads_aligned_bam}" \
        | bcftools call \
        -P 0 -m --ploidy 1 \
        --threads \$(nproc) \
        -Ov -o everything.vcf

    cat everything.vcf | bcftools filter -e "FMT/DP<${min_coverage}" -S . > filtered.vcf
    cat filtered.vcf | bcftools filter -i "MAC>0" > "${out_basename}.discordant.vcf"

    bcftools filter -i 'MAC=0' filtered.vcf | bcftools query -f '%POS\\n' | wc -l | tee num_concordant
    bcftools filter -i 'TYPE="snp"'  "${out_basename}.discordant.vcf" | bcftools query -f '%POS\\n' | wc -l | tee num_discordant_snps
    """
}
