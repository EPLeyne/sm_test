# To run:
# module load miniconda3/4.3.24
# snakemake -j                (-j flag runs script on multiple cores)

# To dryrun:
# module load miniconda3/4.3.24
# snakemake -n

# To create a DAG file in dryrun:
# module load miniconda3/4.3.24
# snakemake -n --dag | dot -Tsvg > dag.svg 

temp_loc = "/scratch1/ley015"
IDS = [1,2]
PUs = ['P','U']

sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]#,
#            "CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002"]

#samples = {f[:74] for f in os.listdir(".") if f.endswith("fastq.gz")}

MAX_THREADS = 32
trim_path = '/apps/trimmomatic/0.38/trimmomatic-0.38.jar'

rule all:
    input:
        expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.zip", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
        expand("{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.html", temp_loc = temp_loc, sample = sample, id = IDS, pu = PUs),
        # expand("{temp_loc}/reports/trimmed_reads/{sample}_1P.fq.gz", temp_loc = temp_loc, sample = sample),
        # expand("{temp_loc}/reports/trimmed_reads/{sample}_1U.fq.gz", temp_loc = temp_loc, sample = sample),
        # expand("{temp_loc}/reports/trimmed_reads/{sample}_2P.fq.gz", temp_loc = temp_loc, sample = sample),
        # expand("{temp_loc}/reports/trimmed_reads/{sample}_2U.fq.gz", temp_loc = temp_loc, sample = sample),
        expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.zip", sample = sample, temp_loc = temp_loc, id = IDS),
        expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.html", sample = sample, temp_loc = temp_loc, id = IDS),

rule fastqc_raw:
    input:
        fastq = expand("test_data/{sample}_R{id}.fastq.gz", sample = sample, id = IDS)
    output:
        zip1 = expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.zip", sample = sample, temp_loc = temp_loc, id = IDS),
        html1 = expand("{temp_loc}/reports/raw_reads/qc/{sample}_R{id}_fastqc.html", sample = sample, temp_loc = temp_loc, id = IDS),
    threads:
        MAX_THREADS
    shell:
        """
        module load fastqc/0.11.8
        mkdir -p {temp_loc}/reports/raw_reads/
        fastqc -t {threads} {input.fastq} -o {temp_loc}/reports/raw_reads/qc/
        """
 
rule trim:
    input:
#        fastq = expand("test_data/{sample}.fastq.gz", sample = sample)
        forward = "test_data/{sample}_R1.fastq.gz",
        reverse = "test_data/{sample}_R2.fastq.gz",
    output:
        forward_paired = "{temp_loc}/reports/trimmed_reads/{sample}_1P.fq.gz",
        forward_unpaired = "{temp_loc}/reports/trimmed_reads/{sample}_1U.fq.gz",
        reverse_paired = "{temp_loc}/reports/trimmed_reads/{sample}_2P.fq.gz",
        reverse_unpaired = "{temp_loc}/reports/trimmed_reads/{sample}_2U.fq.gz",        
    threads:
        MAX_THREADS
    shell:
        """
        module load trimmomatic/0.38               # Had to remove -trimlog as it was not creating the correct 
        java -jar {trim_path} PE -phred33 \
        {input.forward} {input.reverse} \
        {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
        LEADING:3 \
        TRAILING:3 \
        ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10
        module unload trimmomatic/0.38
        """

rule fastqc_trimmed:
    input:
        fq1 = expand("{temp_loc}/reports/trimmed_reads/{sample}_{id}{pu}.fq.gz", id = IDS, temp_loc = temp_loc, sample = sample, pu = PUs),
    output:
        zip2 = "{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.zip",
        html2 = "{temp_loc}/reports/trimmed_reads/qc/{sample}_{id}{pu}_fastqc.html",
    threads:
        32
    shell:
        """
        module load fastqc/0.11.8
        mkdir -p {temp_loc}/reports/trimmed_reads/qc/
        fastqc -t 32 {input.fq1} -o {temp_loc}/reports/trimmed_reads/qc/
        module unload fastqc/0.11.8
        """

