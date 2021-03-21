
from snakemake.remote import FTP
import os
import pandas as pd

FTP = FTP.RemoteProvider()


configfile: "config.yaml"

samples = pd.read_csv(config['samplesheet']).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")





rule all:
    input:
        expand(["ref/annotation.chr{chrom}.gtf",
                "ref/genome.chr{chrom}.fa"], chrom=config["chrom"]),
        expand("reads/{sample}.chr{chrom}_R{group}.fastq.gz", 
               group=[1, 2], sample=samples["sample"], chrom=config["chrom"])


rule annotation:
    input:
        FTP.remote("ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz", static=True, keep_local=True)
    output:
        "ref/annotation.chr{chrom}.gtf"
    shell:
        "zgrep -P ^{wildcards.chrom} {input} > {output}"


rule genome:
    input:
        FTP.remote("ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.{chrom}.fa.gz", static=True, keep_local=True)
    output:
        "ref/genome.chr{chrom}.fa"
    shell:
        "gzip -d -c {input} > {output}"

rule transcriptome:
    output:
        "ref/transcriptome.chr{chrom}.fa"
    shell:
        """wget -O {output} 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Filter name = "chromosome_name" value = "21"/><Attribute name = "ensembl_transcript_id" /><Attribute name = "cdna" /></Dataset></Query>'"""


rule reads:
    input:
        bam = lambda wildcards: os.path.join(config['bamdir'], wildcards.sample + ".sorted.bam"),
        bai = lambda wildcards: os.path.join(config['bamdir'], wildcards.sample + ".sorted.bam.bai")
    output:
        "reads/{sample}.chr{chrom}_R1.fastq.gz",
        "reads/{sample}.chr{chrom}_R2.fastq.gz"
    params:
        # url=config["bam"],
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools bam2fq -1 {output[0]} -2 {output[1]} "
        # "<(samtools view -b -s{params.seed}.2 {input.bam} chr{wildcards.chrom})"
        "<(samtools view -b {input.bam} chr{wildcards.chrom} | samtools sort -n)"



    
rule samtools_index:
    input:
        "{filename}.bam"
    output:
        "{filename}.bam.bai"
    shell:
        "samtools index {input}"


rule sam_to_bam_sort:
    input:
        "{filename}.bam"
    output:
        "{filename}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"


