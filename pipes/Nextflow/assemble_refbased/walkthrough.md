# WDL to Nextflow Conversion Walkthrough: `assemble_refbased`

I have successfully converted the `assemble_refbased` WDL workflow to Nextflow.

## Directory Structure
The new Nextflow pipeline is located at `pipes/Nextflow/assemble_refbased`:

```
pipes/Nextflow/assemble_refbased/
├── nextflow.config       # Configuration (Docker profiles, resources)
├── main.nf               # Main workflow definition
└── modules/
    └── local/
        ├── assembly.nf   # Alignment, Trimming, Refinement, Discordance
        ├── intrahost.nf  # Lofreq variant calling
        ├── reports.nf    # Metrics and Plots
        └── read_utils.nf # BAM Merging
```

## Prerequisites

1.  **Docker Desktop**: Must be installed and running.
    -   Verify with `docker ps`.
2.  **Java 11+**: Nextflow requires a recent Java version.
    -   Recommended: OpenJDK 17 or 19.
    -   If your default Java is old (e.g., Java 8), set `JAVA_HOME` before running.

## Running the Pipeline

```bash
cd pipes/Nextflow/assemble_refbased

# Set JAVA_HOME (Example for Temurin 19 on macOS)
export JAVA_HOME="/Library/Java/JavaVirtualMachines/temurin-19.jdk/Contents/Home"

# Run the pipeline
# Note: The first run will take time to download Docker images (~7GB).
nextflow run main.nf \
    -profile docker \
    --reads_unmapped_bams "/Users/xg203/workspace/antigravity/viral-pipelines/test/input/G5012.3.testreads.bam" \
    --reference_fasta "/Users/xg203/workspace/antigravity/viral-pipelines/test/input/ebov-makona.fasta" \
    --sample_name "test_sample" \
    --outdir "./results" \
    -resume
```

## Verification Results
- **Logic**: All WDL tasks mapped to Nextflow processes. Channel logic handles scatter/gather patterns.
- **Execution**: Verified that Nextflow successfully connects to Docker and pulls the required images (e.g., `quay.io/broadinstitute/viral-core`).
- **Output**: Results will appear in the `./results` directory.

## Next Steps
- Expand coverage to other workflows in `pipes/WDL/workflows/`.
