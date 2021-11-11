import re
import snakemake.io

# Genome URLs: ensembl
genomeURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa.gz"
gtfURL = "http://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf"
cdnaURL = "http://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa"
cdsURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz"


# Directories
fastqDir = "fastq/"
genomeDir = "genome/"
starIndexDir = "starindex/"
salmonIndexDir = "salmonindex/"

# Genome files
genomeFile = genomeDir + "Mus_musculus.GRCm38.dna_rm.primary_assembly.fa"
gtfFile = genomeDir + "Mus_musculus.GRCm38.93.gtf"
cdsFile = genomeDir + "Mus_musculus.GRCm38.cds.all.fa"
transcriptomeFile = genomeDir + "Mus_musculus.GRCm38.cdna.all.fa"
genomeIndex = genomeFile + ".fai"
genomeSizeFile = genomeFile + ".txt"


# Files to test if an index was created
starIndexTestFile = starIndexDir + "chrName.txt"
salmonIndexTestFile = salmonIndexDir + "indexing.log"

DIR = "/home/omar/Downloads/Data/ObOb/rawdata/"
(SAMPLES,LANES, READS, ) = snakemake.io.glob_wildcards(DIR+ "fastq/{sample}_{lane}_{read}_001.fastq.gz")

READS = ["R1", "R2"]

# Standardize filenames
SAMPLESsum  = list(set([re.sub("(ob_ob|WT)_(ND|HFD|HDF)_\d_(\D\D).*", "\\1_\\2_\\3",x) for x in SAMPLES]))
SAMPLESsum = list(set([re.sub("HDF", "HFD" ,x ) for x in SAMPLESsum]))
#print("FASTQ Number %d" % len(SAMPLES))
#print("Conditions Number %d" % len(SAMPLESsum))

# Lists for rule all
QCS = expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=READS)
MULTIQC = ["multiqc/report.html"]
QUANTS = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
COUNTS = "counts/TPMaverage.csv"
INIT= expand("models/context/model_{sample}.mat", sample = SAMPLESsum)

rule all: # requires all outputs as input
     input:
         QCS,
         MULTIQC,
         QUANTS,
         COUNTS,
         INIT

####################################################################
# Merge lanes
rule combine_lanes:
    # Add file one by one to the output target
    output:
        r1 = "fastq/{sample}_R1.fastq.gz",
        r2 = "fastq/{sample}_R2.fastq.gz"
    shell:
        """
        cat {DIR}/fastq/{wildcards.sample}_L*_R1_001.fastq.gz > {output.r1}
        cat {DIR}/fastq/{wildcards.sample}_L*_R2_001.fastq.gz > {output.r2}
        """
####################################################################
# Quality control
rule fastqc:
    input:
        R1 = "fastq/{sample}_R1.fastq.gz",
        R2 = "fastq/{sample}_R2.fastq.gz"
    output:
        "qc/{sample}_R1_fastqc.html",
        "qc/{sample}_R2_fastqc.html"
    threads: 2
    shell:
        "fastqc -o qc --threads {threads} -f fastq {input.R1} {input.R2}"

# https://github.com/steveped/ngsReports/tree/ef0da938c6dcd677c24895dc0f25cf27eed489d5
# To be added after each step
# rule ngsReports:
#     input: expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=PAIRS)
#     output: "qc/ngsReports_Fastqc.html"
#     shell:
#         "Rscript Rscript/ngs_reports.R fastq"

rule multiqc:
    input: expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=READS)
    output: "multiqc/report.html"
    shell:
        """
        mkdir -p multiqc
        multiqc -f --filename report --outdir multiqc qc
        """

####################################################################
# Load genome files

rule load_db_files:
    output: genomeFile, gtfFile, transcriptomeFile, cdsFile
    shell:
        """
        wget {genomeURL} -O - | gunzip -c > {genomeFile}
        wget {gtfURL} -O - | gunzip -c > {gtfFile}
        wget {cdnaURL} -O - | gunzip -c > {transcriptomeFile}
        wget {cdsURL} -O - | gunzip -c > {cdsFile}
        """

####################################################################
# Salmon

rule salmon_index:
    input: transcriptomeFile
    output: salmonIndexTestFile
    log: "logs/salmon_index.log"
    shell:
        "salmon index -t {input} -i {salmonIndexDir} &> {log}"


rule salmon_quant:
    input:
        R1 = fastqDir + "{sample}_R1.fastq.gz",
        R2 = fastqDir + "{sample}_R2.fastq.gz",
        testfile = salmonIndexTestFile
    output: "salmon/{sample}/quant.sf"
    params:
        prefix = "salmon/{sample}"
    threads: 12
    log: "logs/salmon_{sample}_quant.log"
    shell:
        """
        salmon quant \
        --index {salmonIndexDir} \
        --libType A \
        --numBootstraps 100 \
        --threads {threads} \
        -1 {input.R1} -2 {input.R2} \
        --output {params.prefix} &> {log}
        """


rule extract_TPM:
    input: QUANTS
    output: COUNTS
    shell:
        """
        Rscript tx2median.R --genome {gtfFile} --directory salmon --out {COUNTS}
        """

rule context_models:
    input: COUNTS
    output: INIT
    params:
        samples = SAMPLESsum
    run:
        for i in range(0, len(params.samples)):
            AMP = params.samples[i]
            shell("""/home/omar/NoApps/matlab/bin/matlab  -nosplash -nodesktop -r "tInit('%s');quit;" """ % (AMP))
