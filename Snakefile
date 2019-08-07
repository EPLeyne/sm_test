temp_loc = "/scratch1/ley015"

sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002_R1",
            "CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002_R2"]

MAX_THREADS = 32

rule all:
  input:
      expand("{temp_loc}/reports/raw_reads/{sample}_fastqc.zip", sample = sample, temp_loc = temp_loc),
      expand("{temp_loc}/reports/raw_reads/{sample}_fastqc.html", sample = sample, temp_loc = temp_loc)

rule fastqc:
    input:
        fastq = expand("test_data/{sample}.fastq.gz", sample = sample)
    output:
        zip1 = expand("{temp_loc}/reports/raw_reads/{sample}_fastqc.zip", sample = sample, temp_loc = temp_loc),
        html = expand("{temp_loc}/reports/raw_reads/{sample}_fastqc.html", sample = sample, temp_loc = temp_loc),
    threads:
        MAX_THREADS
    shell:
        """
        module load fastqc/0.11.8
        fastqc -t {threads} {input.fastq} -o {temp_loc}/reports/raw_reads/
        module unload fastqc/0.11.8
        """
 

# Below works!

# sample = ["CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002_R1",
#             "CA73YANXX_8_161220_BPO--000_Other_TAAGGCGA-CTCTCTAT_R_161128_SHADIL_LIB2500_M002_R2"]

# MAX_THREADS = 32

# rule all:
#   input:
#       expand("reports/raw_reads/{sample}_fastqc.zip", sample = sample),
#       expand("reports/raw_reads/{sample}_fastqc.html", sample = sample)

# rule fastqc:
#     input:
#         fastq = expand("test_data/{sample}.fastq.gz", sample = sample)
#     output:
#         zip1 = expand("reports/raw_reads/{sample}_fastqc.zip", sample = sample),
#         html = expand("reports/raw_reads/{sample}_fastqc.html", sample = sample),
#     threads:
#         MAX_THREADS
#     shell:
#         """
#         module load fastqc/0.11.8
#         fastqc -t {threads} {input.fastq} -o {temp_loc}/reports/raw_reads/
#         module unload fastqc/0.11.8
#         """
 