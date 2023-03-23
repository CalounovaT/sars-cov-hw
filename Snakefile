
rule all:
    input:
        "pipeline.finished"

rule index_reference:
  """
  Indexes the reference sequence.
  """
  input:
    reference = "data/reference.fasta"
  output:
    "data/reference.fasta.bwt"
  log:
    "logs/index_reference.log"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    bwa index {input.reference} &> {log}
    """

rule align_reads:
    """
    Aligns the reads to the reference sequence.
    """
    input:
        reference = "data/reference.fasta",
        reference_index = "data/reference.fasta.bwt",
        r1_reads = "data/{sample}.R1.fastq.gz",
        r2_reads = "data/{sample}.R2.fastq.gz"
    output:
        "data/{sample}_align.bam"
    log:
        "logs/align_reads_{sample}.log"
    params:
        memory="4"
    threads:
        6
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.r1_reads} {input.r2_reads} | samtools sort -@ {threads} -o {output} &> {log}
        """

rule index_align:
    """
    Indexes the alignment.
    """
    input:
        "data/{sample}_align.bam"
    output:
        "data/{sample}_align.bam.bai"
    log:
        "logs/index_align_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        samtools index {input} &> {log}
        """

rule call_variants:
    """
    Calls variants from the alignment.
    """
    input:
        reference = "data/reference.fasta",
        alignment = "data/{sample}_align.bam",
        alignment_index = "data/{sample}_align.bam.bai"
    output:
        "data/{sample}_variants.vcf.gz"
    log:
        "logs/call_variants_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        bcftools mpileup -d 10000 -f {input.reference} {input.alignment} | bcftools call --ploidy 1 -mv -Ob -o {output} &> {log}
        """

rule index_variants:
    """
    Indexes the variants.
    """
    input:
        "data/{sample}_variants.vcf.gz"
    output:
        "data/{sample}_variants.vcf.gz.csi"
    log:
        "logs/index_variants_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        bcftools index {input} &> {log}
        """

rule create_consensus:
    """
    Creates a consensus sequence from the variants.
    """
    input:
        reference = "data/reference.fasta",
        variants = "data/{sample}_variants.vcf.gz",
        variants_index = "data/{sample}_variants.vcf.gz.csi"
    output:
        "data/{sample}_consensus.fasta"
    log:
        "logs/create_consensus_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        cat {input.reference} | bcftools consensus -o {output} {input.variants} &> {log}
        """

def get_samples(wildcards):
    indeces = glob_wildcards("data/{sample}.R1.fastq.gz").sample
    completed = expand("data/{sample}_consensus.fasta",sample=indeces)
    return completed

rule finish:
    input:
        get_samples
    output:
        "pipeline.finished"
    log:
        "logs/finish.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        touch {output}
        """