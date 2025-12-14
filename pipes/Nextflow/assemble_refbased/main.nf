nextflow.enable.dsl=2

include { ALIGN_READS; IVAR_TRIM; REFINE_ASSEMBLY; RUN_DISCORDANCE } from './modules/local/assembly.nf'
include { MERGE_BAMS } from './modules/local/read_utils.nf'
include { LOFREQ } from './modules/local/intrahost.nf'
include { ALIGNMENT_METRICS; PLOT_COVERAGE } from './modules/local/reports.nf'

workflow {
    // Input channels
    ch_reads_bam = Channel.fromPath(params.reads_unmapped_bams)
                          .map { file -> tuple(file.simpleName, file) }
    
    ch_ref_fasta = Channel.value(file(params.reference_fasta))
    ch_trim_bed = params.trim_coords_bed ? file(params.trim_coords_bed) : file('NO_TRIM_BED')

    // Optional input handling
    ch_novocraft_license = params.novocraft_license ? file(params.novocraft_license) : []
    
    // Step 1: Align reads to reference (Scatter)
    ALIGN_READS(
        ch_reads_bam,
        ch_ref_fasta,
        params.aligner,
        ch_novocraft_license,
        params.skip_mark_dupes,
        "" // aligner_options
    )

    // Step 2: Trim primers (Scatter)
    IVAR_TRIM(
        ALIGN_READS.out.aligned_bam,
        ch_trim_bed
    )

    // Step 3: Merge BAMs (Gather)
    IVAR_TRIM.out.trimmed_bam
        .map { it[1] }
        .collect()
        .set { ch_all_trimmed_bams }

    MERGE_BAMS(
        ch_all_trimmed_bams,
        params.sample_name,
        params.sample_name + ".align_to_ref.trimmed"
    )

    // Prepare channel for downstream steps acting on the merged BAM
    MERGE_BAMS.out.out_bam
        .map { bam -> tuple(params.sample_name, bam) }
        .set { ch_merged_bam_tuple }

    // Step 4: Intrahost Variants (Reference)
    LOFREQ(
        ch_merged_bam_tuple,
        ch_ref_fasta
    )

    // Step 5: Metrics & Reports
    ALIGNMENT_METRICS(
        ch_merged_bam_tuple,
        ch_ref_fasta,
        ch_trim_bed,
        params.min_coverage
    )

    RUN_DISCORDANCE(
        ch_merged_bam_tuple,
        ch_ref_fasta,
        params.min_coverage + 1
    )

    PLOT_COVERAGE(
        ch_merged_bam_tuple
    )

    // Step 6: Refine Assembly (Call Consensus)
    REFINE_ASSEMBLY(
        ch_merged_bam_tuple,
        ch_ref_fasta,
        params.min_coverage,
        params.major_cutoff
    )

    // Step 7: Align to Self (New Assembly)
    // Combine original READS with NEW reference
    // We need to spread the new assembly to all original read BAMs
    // ch_reads_bam is [sample_chunk, bam] coverage
    // REFINE_ASSEMBLY.out.assembly_fasta is [path] (since it's not a tuple output, just path)
    
    // We want to align EACH chunk to the NEW reference.
    
    ch_reads_bam.combine(REFINE_ASSEMBLY.out.assembly_fasta)
        .set { ch_reads_and_new_ref }
    
    ALIGN_TO_SELF(
        ch_reads_and_new_ref,
        params.aligner,
        ch_novocraft_license,
        params.skip_mark_dupes,
        ""
    )

    // Step 8: Merge Self-Aligned BAMs
    ALIGN_TO_SELF.out.aligned_bam
        .map { it[1] }
        .collect()
        .set { ch_all_self_bams }

    MERGE_BAMS_SELF(
        ch_all_self_bams,
        params.sample_name,
        params.sample_name + ".merge_align_to_self"
    )
    
    // Prepare channel for downstream self-related steps
    MERGE_BAMS_SELF.out.out_bam
        .map { bam -> tuple(params.sample_name, bam) }
        .set { ch_merged_self_bam_tuple }

    // Step 9: Intrahost Variants (Self)
    LOFREQ_SELF(
        ch_merged_self_bam_tuple,
        REFINE_ASSEMBLY.out.assembly_fasta
    )
    
    // Step 10: Plot Coverage (Self)
    PLOT_COVERAGE_SELF(
        ch_merged_self_bam_tuple
    )
}

// We define a separate process call for ALIGN_TO_SELF because it uses different input structures
// OR we can alias ALIGN_READS. In DSL2 we can alias.
process ALIGN_TO_SELF {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(reads_bam), path(new_ref_fasta)
    val aligner
    path novocraft_license
    val skip_mark_dupes
    val aligner_options

    output:
    tuple val(sample_name), path("${sample_name}.mapped.bam"), emit: aligned_bam

    script:
    def skip_mark_dupes_arg = skip_mark_dupes ? '--skipMarkDupes' : ''
    def license_arg = novocraft_license ? "--NOVOALIGN_LICENSE_PATH=${novocraft_license}" : ''
    def options_arg = aligner_options ? "--aligner_options \"${aligner_options}\"" : ''
    """
    cp ${new_ref_fasta} assembly.fasta
    read_utils.py index_fasta_picard assembly.fasta
    read_utils.py index_fasta_samtools assembly.fasta
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
    """
}

// Aliasing MERGE_BAMS for clarity
process MERGE_BAMS_SELF {
    tag "$out_basename"
    input:
    path in_bams
    val sample_name
    val out_basename
    output:
    path "${out_basename}.bam", emit: out_bam
    script:
    """
    mem_in_mb=\$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
    read_utils.py merge_bams ${in_bams} merged.bam --JVMmemory="\$mem_in_mb"m --loglevel DEBUG
    mv merged.bam "${out_basename}.bam"
    """
}

// Aliasing LOFREQ for clarity
process LOFREQ_SELF {
    tag "$out_basename"
    input:
    tuple val(out_basename), path(aligned_bam)
    path reference_fasta
    output:
    path "${out_basename}.vcf", emit: report_vcf
    script:
    """
    cp ${aligned_bam} aligned.bam
    cp ${reference_fasta} reference.fasta
    samtools faidx reference.fasta
    samtools index aligned.bam
    lofreq call -f reference.fasta -o "${out_basename}.vcf" aligned.bam
    """
}

// Aliasing PLOT_COVERAGE for clarity
process PLOT_COVERAGE_SELF {
    tag "$sample_name"
    input:
    tuple val(sample_name), path(aligned_reads_bam)
    output:
    path "${sample_name}.coverage_plot.pdf", emit: coverage_plot
    script:
    """
    samtools view -c ${aligned_reads_bam} > reads_aligned
    if [ "\$(cat reads_aligned)" != "0" ]; then
        samtools index -@ \$(nproc) "${aligned_reads_bam}"
        reports.py plot_coverage "${aligned_reads_bam}" "${sample_name}.coverage_plot.pdf" --outSummary "${sample_name}.coverage_plot.txt" --plotFormat pdf --plotTitle "${sample_name} coverage plot"
    else
        touch "${sample_name}.coverage_plot.pdf"
    fi
    """
}
