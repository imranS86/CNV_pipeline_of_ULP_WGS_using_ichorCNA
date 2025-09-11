# Copy Number Variant (CNV) Workflow for Ultra Low-Pass WGS Using ichorCNA
This is my working bash script in command line for processing ultra low-pass whole genome sequencing (WGS) data to detect copy number alterations (CNA's) using the `ichorCNA` tool. I have added details so that anyone can use it starting from processing the raw WGS fastq files to get the genome-wide plots for CNAs. ichorCNA predicts large scale copy number alterations and estimate the fraction of tumor in cell-free DNA of an ultra-low-pass whole genome sequencing sample (ULP-WGS, 0.1x coverage). It uses a probabilistic model which segment the genome (based on a selected bin size) and predict large scale CNA's.
## Prerequisites
To run this pipeline, you will need to install the following tools:
* [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): A quality control tool for raw sequence data.
* [`fastp`](https://github.com/OpenGene/fastp) : An all-in-one FASTQ preprocessor for quality control, filtering, and trimming.
* [`MultiQC`](https://github.com/MultiQC/MultiQC/) : A tool that aggregates results from multiple bioinformatics analyses into a single report.
* [`bwa`](http://bio-bwa.sourceforge.net/) : A tool for aligning sequencing reads to a reference genome.
* [`samtools`](http://www.htslib.org/) : A suite of tools for manipulating and analyzing alignment files.
* [`sambamba`](https://lomereiter.github.io/sambamba/) : A tool for manipulating alinment files.
* [`ichorCNA`](https://github.com/broadinstitute/ichorCNA) : An R package for estimating tumor fraction and predicting large-scale CNVs.
# Step-by-Step Pipeline
##  Pipeline setup and directory structure
To make this pipeline reproducible, we will use a clear directory structure and variable paths for each major step to make your code clean and avoid any confusion.
```bash
# Project directories (environment variables for general paths)
PROJECT_DIR="/path/to/your/project_name"
RAW_DATA_DIR="$PROJECT_DIR/raw_data"
QC_DIR="$PROJECT_DIR/qc_reports"
TRIMMED_DATA_DIR="$PROJECT_DIR/trimmed_data"
ALIGNMENT_DIR="$PROJECT_DIR/alignment"
WIG_DIR="$PROJECT_DIR/wig_files"
RESULTS_DIR="$PROJECT_DIR/results"
REF_DIR="$PROJECT_DIR/ichorCNA_references"

# Reference files (place your files in the $REF_DIR)
REFERENCE_GENOME="$REF_DIR/hg38.fa.gz" # make sure this folder contain all the indexed files of reference genome by BWA-mem
#Reference files for ichorCNA
GC_WIG="$REF_DIR/gc_hg38_1000kb.wig"
MAP_WIG="$REF_DIR/map_hg38_1000kb.wig"
CENTROMERE_FILE="$REF_DIR/GRCh38.GCA_000001405.2_centromere_acen.txt"
NORMAL_PANEL="$REF_DIR/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"

# ichorCNA script files (adjust these paths based on your system)
ICHORCNA_SCRIPTS="/path/to/ichorCNA/scripts"
HMMCOPY_BIN="/path/to/hmmcopy_utils-master/bin"

# Create all necessary directories (folders)
mkdir -p "$RAW_DATA_DIR"
mkdir -p "$QC_DIR/raw"
mkdir -p "$TRIMMED_DATA_DIR"
mkdir -p "$QC_DIR/trimmed"
mkdir -p "$ALIGNMENT_DIR/raw_bams"
mkdir -p "$ALIGNMENT_DIR/sorted_bams"
mkdir -p "$ALIGNMENT_DIR/filtered_bams"
mkdir -p "$ALIGNMENT_DIR/duplicate_removed_bams"
mkdir -p "$WIG_DIR"
mkdir -p "$RESULTS_DIR"
mkdir -p "$REF_DIR"
```
Move all your raw fastq files for all your samples into the `raw_data` directory. For example, KT3_S29_R1_001.fastq.gz and KT3_S29_R2_001.fastq.gz (forward and reverese fastq files) of a sample respectively. 
## Quality Control on Raw FASTQ files
Before we start the analysis it is important to assess the quality of our raw fastq files(forward and reverse end) for each sample before trimming. We use FASTQC tool to generate quality reports. 
```bash
# Navigate to the raw_data directory
cd $RAW_DATA_DIR

# Run FastQC on all zipped FASTQ files and output reports to the 'raw' QC directory
echo "Running FastQC on raw FASTQ files..."
fastqc *.fastq.gz -o $QC_DIR/raw/
```
 `fastqc *.fastq.gz`: This command runs FastQC on all files in the current directory that end with ".fastq.gz".
 
`-o $QC_DIR/raw/`: This option specifies the output directory for the FastQC reports for each file.
## Trimming adaptor and low quality sequences
Raw reads contain low-quality base calls and adapter sequences attached during sequencing, which can interfere with downstream analysis. The `fastp` tool efficiently remove the adaptor sequences and low quality reads.
```bash
# Navigate to the raw_data directory
cd $RAW_DATA_DIR

# Use a for loop to process multiple paired-end samples at once. This makes it easier to apply the script to multiple files intead of running files seperately.
echo "Starting fastp trimming..."
for R1_FILE in *_R1_001.fastq.gz; do
    R2_FILE=${R1_FILE/_R1_001/_R2_001}
    BASE_NAME=$(basename ${R1_FILE} _R1_001.fastq.gz)
    
    fastp -i ${R1_FILE} -I ${R2_FILE} \
          -o $TRIMMED_DATA_DIR/${BASE_NAME}_R1_trimmed.fastq.gz \
          -O $TRIMMED_DATA_DIR/${BASE_NAME}_R2_trimmed.fastq.gz \
done
```
Here `-i` and `-I` demonstrate the input files for Read 1 and Read 2, respectively. And `-o` and `-O`options show the output files for Read 1 and Read 2 after trimming.
## Quality Control after trimming with FastQC + MultiQC
Once the fastq files have been trimmed, its important to run FastQC again to see the trimming has improved the overall quality of the reads and adaptor sequences have been removed. `multiqc` combines all the FastQC html reports into a single user friendly html report to scan the overall quality of all the files.
```bash
# Navigate to the'trimmed_data directory
cd $TRIMMED_DATA_DIR

# Run FastQC on all trimmed gzipped FASTQ files
echo "Running FastQC on trimmed FASTQ files..."
fastqc *.fastq.gz -o $QC_DIR/trimmed/

# merge all the fastqc html reports with MultiQC. multiqc will scan all the html files in the current folder
echo "Generating MultiQC report..."
multiqc $QC_DIR/raw/ $QC_DIR/trimmed/ -o $QC_DIR/multiqc_report
```
Here `$QC_DIR/raw/` and `$QC_DIR/trimmed/` are the input directories containing fastqc reports to be summarized, and `-o $QC_DIR/multiqc_report` mentioned to download the output directory for the MultiQC report.
## Alignment with BWA
This step aligns the trimmed fastq files with the reference genome (hg38) using `bwa-mem` aligner. `bwa mem` is a fast and accurate alignment algorithm for longer reads (>70bp).The reference genome must be indexed before alignment. If you have not indexed your reference genome you can do this by `bwa index $REFERENCE_GENOME`. Before indexing, download the hg38 reference genome fasta file to REFERENCE_GENOME folder.
```bash
# Navigate to the alignment directory
cd $ALIGNMENT_DIR/raw_bams

# Use a for loop to process the alignment for multiple samples in a single run.
echo "Starting bwa-mem alignment..."
for R1_FILE in $TRIMMED_DATA_DIR/*_R1_trimmed.fastq.gz; do
    R2_FILE=${R1_FILE/_R1_trimmed/_R2_trimmed}
    BASE_NAME=$(basename ${R1_FILE} _R1_trimmed.fastq.gz)
    
    # Run bwa mem, pipe output to samtools view for conversion to BAM file
    bwa mem -t 12 -M $REFERENCE_GENOME $R1_FILE $R2_FILE | \
    samtools view -b -o ${BASE_NAME}.bam -

    echo "Alignment complete for sample ${BASE_NAME}."
done
```
Since i use my personal computer with 16 GB RAM so i use thread option `-t 12` which refers to number of threads use to speed up the process.
`-M` marks shorter split hits as secondary alignments, which is recommended by bwa-mem (for Picard and GATK compatibility). `|` is the pipe symbol which redirects the standard output (`stdout`) of `bwa-mem` to the standard input (`stdin`) of `samtools view`. This pipe operator avoids creating a large intermediate SAM file. `samtools view -b`converts the SAM format output from bwa-mem into the compressed BAM format.`-o ${BASE_NAME}.bam` specifies the output BAM file name.`-` tells `samtools view` to read from stdin.
## Sorting and indexing BAM files
After alignment, the BAM files are sorted and indexed to be useful for downstream analysis. Sorting improves access and downstream compatibility and indexing is important for read retrieval by region. 
```bash
# Navigate to the raw BAMs directory
cd $ALIGNMENT_DIR/raw_bams

echo "Sorting and indexing BAM files..."
for BAM_FILE in *.bam; do
    BASE_NAME=$(basename ${BAM_FILE} .bam)
    
    # Sort the BAM file by coordinate
    samtools sort -@ 16 -o $ALIGNMENT_DIR/sorted_bams/${BASE_NAME}_sorted.bam ${BAM_FILE}
    
    # Index the sorted BAM file for fast lookup
    samtools index $ALIGNMENT_DIR/sorted_bams/${BASE_NAME}_sorted.bam
   
    echo "Processing complete for BAM ${BAM_FILE}."
done
```
`samtools sort` sorts the BAM file by genomic coordinates. `-@ 16` uses 16 threads for sorting, which speeds up the process. `-o`specifies the output file name for the sorted BAM.`samtools index` creates an index BAM file(`.bai`) extension, which allows for quick retrieval of alignments in specific genomic regions. 

At this step you can also View the mapping results to explore the alignment summary (e.g., total reads, mapped reads, properly paired reads, un paired reads, singleton, duplicates) using `samtools flagstat` command on sorted BAM files.
## Filter Properly Paired, Mapped Reads
We used `sambamba view` to filter out reads that are not properly paired, unmapped, or are secondary alignments and keep only properly paired and mapped reads. These artifacts can interfere with accurate CNV detection.
```bash
# Navigate to the sorted BAMs directory
cd $ALIGNMENT_DIR/sorted_bams

echo "Filtering reads with sambamba view..."
for f in *_sorted.bam; do
    base=$(basename "$f" _sorted.bam)
    sambamba view -t 16 -f bam -F "proper_pair and not (unmapped or secondary_alignment)" -o "$ALIGNMENT_DIR/filtered_bams/${base}_filtered.bam" "$f"
done
```
`-f bam` specifies BAM output format. `-F "proper_pair and not (unmapped or secondary_alignment)"` is the filter expression. `proper_pair` keeps only paired-end reads that aligned correctly, with the expected orientation and insert size.`not (unmapped or secondary_alignment)` excludes reads that are either unmapped or represent secondary alignments, ensuring we only use the best possible alignment for each read.
## Remove PCR Duplicates
We use `sambamba markdup` to identify and remove PCR duplicates. PCR duplicates are reads that start at the same genomic location and are likely artifacts of the sequencing process rather than true biological signal. Dupliates can bias copy number estimation, also suggested by ichorCNA documentation to remove PCR duplicates.
```bash
# Navigate to the filtered BAMs directory
cd $ALIGNMENT_DIR/filtered_bams

echo "Removing duplicates with sambamba markdup..."
for f in *_filtered.bam; do
    base=$(basename "$f" _filtered.bam)
    sambamba markdup -t 16 -r "$f" "$ALIGNMENT_DIR/duplicate_removed_bams/${base}_duprem.bam"
done
```
`-r` removes duplicate reads instead of just marking them. This is often preferred for CNV analysis and suggested in ichorCHA guidelines. `"$f"` is the input BAM file and `"$ALIGNMENT_DIR/duplicate_removed_bams/${base}_duprem.bam"` demonstrate the output BAM file with duplicates removed.
# ichorCNA workflow
Now we will run the ichorCNA workflow which will use duplicate removed BAM files in the previous step to call CNV variants. 
For the ichorCNA steps, you will need the following files:
* HMMcopy_utils: This package contains the `readCounter` utility, which is required for preparing the WIG files for `ichorCNA`. WIG files are read count files generated from duplicate_removed_bam file as input.
* ichorCNA scripts: Download the core scripts from the ichorCNA github repository [scripts](https://github.com/broadinstitute/ichorCNA/tree/master/scripts/) including `runIchorCNA.R` and `createPanelOfNormals.R`, are needed.
* Reference files: You will also need a GC content WIG file, a mappability WIG file, a centromere file, and a panel of normals (PoN) file. These files should correspond to your reference genome (e.g., hg38) and the window size (1Mb or 500 kb) you plan to use. You can download all these files from ichorCNA github repository [extdata](https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata/).
* Duplicate-removed BAM files: The input for this section of the pipeline is the final, duplicate-removed BAM file generated in the previous steps.

Note: ichorCNA also provide its own PON created from normal samples for 1 Mb (e.g HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds) and 500 kb and also give a script to create your own PON from normal samples. You can use this default PON and can also explore PON from your samples.
## Count reads in fixed windows
The first step before running ichorCNA pipeline is creating WIG files which will be used as an input in the ichorCNA. This step involves counting the reads within fixed-size genomic windows or bin size across the reference genome. This is done using the `readCounter` function from the `HMMcopy_utils` package that generate WIG file. You can select the bin size either 1 Mb or 500 kb based on your requirements. A kep point to mention is that for making WIG files you will need input BAM file (i-e duplicate_removed BAM file in previous step), along with its indexed BAM file in the same folder for `readCounter` function. A handy thing is that even if you don't have an index file in duplicate_removed BAM folder, `readCounter` by default first generate an index file. This index file will be updated (as you have done filtering and duplication on your sorted BAM files so the index file used in this step should be made after these steps). You can either make an index file after removing duplicate step or use the index file generated in `readCounter` step. 
```bash
# Navigate to the directory containing the HMMcopy_utils bin file
cd "$HMMCOPY_BIN"

# Loop through all duplicate-removed BAM files to create WIG files
echo "Counting reads in fixed windows..."
for BAM_FILE in "$ALIGNMENT_DIR/duplicate_removed_bams"/*_duprem.bam; do
    BASE_NAME=$(basename "${BAM_FILE}"_duprem.bam)
    ./readCounter --window 1000000 --quality 20 \
        --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
        "$BAM_FILE" \
        > "$WIG_DIR/${BASE_NAME}.wig"
done
```
`--window 1000000` is window size of 1,000,000 base pairs (1Mb), commonly used filter for low-pass WGS, ichorCNA also suggested to use 1Mb size initially. `--quality 20` filters a minimum mapping quality score of 20. Reads with a score below this will be ignored, helping to remove poor-quality alignments. `--chromosome` lists the chromosomes to be included in the analysis. This ensures consistency and avoids any issues with different chromosome notations, "chr1" notation is used to align with hg38 notation. `>` redirects the standard output of `readCounter` to a new `.wig` file.
## Run ichorCNA for CNV analysis
This is the core ichorCNA analysis step. ichorCNA uses the WIG files, along with reference files and a panel of normals, to estimate tumor fraction, ploidy, and copy number alteration plots. This set of parameters is a good starting point if you expect a moderate to high tumor fraction (e.g., > 50%). It explores a range of normal contamination values from 0.5 to 0.9 in `--normal` filter below. You can optimize this filter if there is very low tumor fraction for example `--normal "c(0.95, 0.99, 0.995, 0.999)`. You can find further details regarding fine tuning for very low tumor fraction in the ichorCNA vignette. This scrip iterate through the generated WIG files, runs the `runIchorCNA.R` script file (R file downloaded from ichorCNA github repository), and the file paths required for reference files including GC content file (gc_hg38_1000kb.wig), mappability file (map_hg38_1000kb.wig),centromeres (GRCh38.GCA_000001405.2_centromere_acen.txt), and a PON (HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds) from [ichorCNA github repository](https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata/) for 1Mb window size for hg38 assembly. Results for each sample are saved in a dedicated output directory containing genome-wide plots for CNVs and other files (details for each result file is explained in the ichorCNA's vignette.
```bash
# Navigate to the ichorCNA scripts directory
cd "$ICHORCNA_SCRIPTS"

echo "Running ichorCNA with lenient parameters..."
for WIG_FILE in "$WIG_DIR"/*.wig; do
    BASE_NAME=$(basename "$WIG_FILE" .wig)
    Rscript runIchorCNA.R \
        --id "${BASE_NAME}" \
        --WIG "${WIG_FILE}" \
        --ploidy "c(2)" \
        --normal "c(0.5,0.6,0.7,0.8,0.9)" \
        --maxCN 3 \
        --gcWig "${GC_WIG}" \
        --mapWig "${MAP_WIG}" \
        --centromere "${CENTROMERE_FILE}" \
        --normalPanel "${NORMAL_PANEL}" \
        --includeHOMD False \
        --chrs "c(1:22)" \
        --chrTrain "c(1:22)" \
        --estimateNormal True \
        --estimatePloidy True \
        --estimateScPrevalence False \
        --scStates "c()" \
        --txnE 0.9999 \
        --txnStrength 10000 \
        --outDir "$RESULTS_DIR/${BASE_NAME}"
done
```
`--id`: Sample ID.
`--WIG`: Path to the WIG file from the `readCounter` step.
`--ploidy "c(2)"`: The expected ploidy. `ichorCNA` can estimate this if `estimatePloidy` is set to `True`. `c(2)` is the diploid state.
`--normal "c(0.5,0.6,0.7,0.8,0.9)"`: A vector of normal contamination fractions to test.
`--maxCN 3`: Maximum copy number of copy number variation to consider for e.g Deletion, duplication etc. 
`--gcWig`, `--mapWig`, `--centromere`: Path to the reference files.
`--normalPanel`: Path to the panel of normals (`.rds` file).
`--includeHOMD False`: Don't include homozygous deletions in the model.
`--chrs "c(1:22)"`: Specifies the chromosomes to analyze.
`--chrTrain "c(1:22)"`: Specifies the chromosomes to use for model training.
`--estimateNormal True`: Estimates the normal fraction.
`--estimatePloidy True`: Estimates the ploidy.
`--txnE` and `--txnStrength`: These parameters control the HMM's transition probabilities, balancing sensitivity and specificity of CNV calls.
`--outDir`: Output folder that will contain results based on wig file names.
    
