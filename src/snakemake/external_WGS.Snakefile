configfile: "config/config.yml"

DATASETS = list(config["external_WGS"]["dataset_to_library_sra"])

genome = config["external_WGS"]["genome"]

#
#  Default targets
#

rule all:
  input:
    lambda wildcards:
      rules.all_multiqc.input,
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.vcf.gz.tbi",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.vcf.gz.tbi",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.stats",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.filter.vcf.gz.tbi",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.filter.stats",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.filter.tsv.gz",
      "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.indels.filter.vcf.gz.tbi",
      "data/external_WGS/fasta/" + genome + "_DGRP-57_DGRP-439_maternal.fa.gz",
      "data/external_WGS/fasta/" + genome + "_DGRP-57_DGRP-439_paternal.fa.gz"

rule all_multiqc:
  input:
    "data/external_WGS/qc/multiqc_report.html"

#
#  SRA data download
#

rule download_sra:
  output:
    "data/external_WGS/fastq/{sample}_1.fastq.gz",
    "data/external_WGS/fastq/{sample}_2.fastq.gz"
  conda:
    "../../env/sra-tools.yaml"
  shell:
    "fastq-dump --split-files --origfmt --gzip -O data/external_WGS/fastq {wildcards.sample}"

#
#  FASTQ quality control
#

rule fastqc:
  input:
    "data/external_WGS/fastq/{file}.fastq.gz"
  output:
    "data/external_WGS/qc/{file}_fastqc.zip"
  conda:
    "../../env/fastqc.yaml"
  shell:
    """
    fastqc -o data/external_WGS/qc {input}
    """

def library_fastq_qcfiles(wildcards):
  return ["data/external_WGS/qc/" + library_sra + "_" + readindex + "_fastqc.zip"
    for libraries in config["external_WGS"]["dataset_to_library_sra"].values()
      for library_sra in libraries
        for readindex in ["1", "2"]]

rule multiqc:
  input:
    library_fastq_qcfiles
  output:
    "data/external_WGS/qc/multiqc_report.html"
  conda:
    "../../env/multiqc.yaml"
  shell:
    """
    cd data/external_WGS/qc
    rm -rf multiqc_data/
    multiqc --interactive .
    """

#
#  WGS read mapping, merging sequencing libraries for the same replicate
#

rule bwa_mem:
  input:
    fastq1 = "data/external_WGS/fastq/{file}_1.fastq.gz",
    fastq2 = "data/external_WGS/fastq/{file}_2.fastq.gz"
  output:
    temp("data/external_WGS/bam/all_reads/" + genome + "_{file}.bam")
  params:
    fasta_genome = config["genomes"][genome]["fasta"]
  threads:
    16
  conda:
    "../../env/bwa_samtools.yaml"
  shell:
    """
    bwa mem -t {threads} {params.fasta_genome} {input.fastq1} {input.fastq2} | samtools view -bT {params.fasta_genome} - > {output}
    """

def dataset_to_library_sra_bamfiles(wildcards):
  if wildcards.dataset in config["external_WGS"]["dataset_to_library_sra"]:
    library_sras = config["external_WGS"]["dataset_to_library_sra"][wildcards.dataset]
    return ["data/external_WGS/bam/all_reads/" + genome + "_" + sra + ".bam" for sra in library_sras]
  else:
    return "/missing_fastq"

rule samtools_cat:
  input:
    lambda wildcards: dataset_to_library_sra_bamfiles(wildcards)
  output:
    temp("data/external_WGS/bam/merged_reads/" + genome + "_{dataset}.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools cat {input} -o {output}
    """

#
#  Sort reads and mark duplicates using Picard, index the resulting .BAM files
#

rule picard_sort:
  input:
    "data/external_WGS/bam/merged_reads/{dataset}.bam"
  output:
    temp("data/external_WGS/bam/merged_reads/{dataset}.sort.bam")
  params:
    tmpdir = config["paths"]["tmpdir"]
  conda:
    "../../env/picard-slim.yaml"
  shell:
    """
    picard -Xmx32g SortSam -INPUT {input} -OUTPUT {output} -SORT_ORDER coordinate -TMP_DIR {params.tmpdir}
    """

rule picard_markdup:
  input:
    "data/external_WGS/bam/merged_reads/{dataset}.sort.bam"
  output:
    bam = "data/external_WGS/bam/merged_reads/{dataset}.sort.markdup.bam",
    metrics = "data/external_WGS/bam/merged_reads/{dataset}.sort.markdup.metrics.txt"
  params:
    tmpdir = config["paths"]["tmpdir"]
  conda:
    "../../env/picard-slim.yaml"
  shell:
    """
    picard MarkDuplicates -INPUT {input} -OUTPUT {output.bam} -METRICS_FILE {output.metrics} -TMP_DIR {params.tmpdir}
    """

rule samtools_index:
  input:
    "data/external_WGS/bam/merged_reads/{dataset}.bam"
  output:
    "data/external_WGS/bam/merged_reads/{dataset}.bam.bai"
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools index {input}
    """

#
#  Call variants using freebayes, filter them and calculate statistics
#

rule freebayes_DGRP57_DGRP439:
  input:
    bam_DGRP57 = "data/external_WGS/bam/merged_reads/" + genome + "_DGRP-57.sort.markdup.bam",
    bai_DGRP57 = "data/external_WGS/bam/merged_reads/" + genome + "_DGRP-57.sort.markdup.bam.bai",
    bam_DGRP439 = "data/external_WGS/bam/merged_reads/" + genome + "_DGRP-439.sort.markdup.bam",
    bai_DGRP439 = "data/external_WGS/bam/merged_reads/" + genome + "_DGRP-439.sort.markdup.bam.bai"
  output:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.vcf.gz"
  params:
    fasta_genome = config["genomes"][genome]["fasta"]
  conda:
    "../../env/freebayes.yaml"
  shell:
    """
    bamaddrg \
      -b {input.bam_DGRP57} -s DGRP-57 \
      -b {input.bam_DGRP439} -s DGRP-439 \
      | freebayes -f {params.fasta_genome} --stdin --no-population-priors \
      | vcffilter -f 'QUAL > 20' \
      | bgzip > {output}
    """

rule tabix:
  input:
    "data/external_WGS/vcf/{dataset}.vcf.gz"
  output:
    "data/external_WGS/vcf/{dataset}.vcf.gz.tbi"
  conda:
    "../../env/freebayes.yaml"
  shell:
    """
    tabix {input}
    """

rule vcfsnps:
  input:
    "data/external_WGS/vcf/{dataset}.vcf.gz"
  output:
    "data/external_WGS/vcf/{dataset}.SNPs.vcf.gz"
  conda:
    "../../env/freebayes.yaml"
  shell:
    """
    bgzip -cd {input} \
      | vcfallelicprimitives -kg \
      | tr '|' '/' \
      | vcfsnps \
      | vcfuniq \
      | vcfbiallelic \
      | bgzip > {output}
    """

rule vcfindels:
  input:
    "data/external_WGS/vcf/{dataset}.vcf.gz"
  output:
    "data/external_WGS/vcf/{dataset}.indels.vcf.gz"
  conda:
    "../../env/freebayes.yaml"
  shell:
    """
    bgzip -cd {input} \
      | vcfallelicprimitives -kg \
      | tr '|' '/' \
      | vcfindels \
      | vcfuniq \
      | bgzip > {output}
    """

rule vcf_statistics:
  input:
    "data/external_WGS/vcf/{dataset}.vcf.gz"
  output:
    "data/external_WGS/vcf/{dataset}.stats"
  conda:
    "../../env/htslib.yaml" # provides bgzip
  shell:
    """
    bgzip -cd {input} \
      | grep -v '^##' \
      | awk '(NR == 1)' \
      | cut -f 10- \
      | sed -e 's/^/count\t/' > {output}
    bgzip -cd {input} \
      | grep -v '^#' \
      | cut -f 10- \
      | sed -e 's/:[^\t]*//g' \
      | sort \
      | uniq -c >> {output}
    """

rule vcf_filter_DGRP57_DGRP439:
  input:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.{type}.vcf.gz"
  output:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.{type}.filter.vcf.gz"
  conda:
    "../../env/htslib.yaml" # provides bgzip
  shell:
    """
    # take only the variants that are homozygous in DGRP-57 and homozygous but different in DGRP-439
    bgzip -cd {input} \
      | awk -F '\t' '$0 ~ /^#/ || ($10 ~ /^0\/0:/ && $11 ~ /^1\/1:/) || ($10 ~ /^1\/1:/ && $11 ~ /^0\/0:/)' \
      | bgzip > {output}
    """

#
#  Extract SNPs in a custom format for haplotype separation
#

rule vcf_extract_DGRP57_DGRP439:
  input:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.filter.vcf.gz"
  output:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.filter.tsv.gz"
  conda:
    "../../env/htslib.yaml" # provides bgzip
  shell:
    """
    bgzip -cd {input} \
      | grep -v '^##' \
      | awk 'BEGIN {{ FS="\t"; OFS="\t" }} \
        {{
          if ($0 ~ /^#/) print "chrom", "pos", $10, $11;
          else if ($10 ~ /^0\/0:/ && $11 ~ /^1\/1:/) print $1, $2, $4, $5;
          else if ($10 ~ /^1\/1:/ && $11 ~ /^0\/0:/) print $1, $2, $5, $4;
        }}' \
      | bgzip > {output}
    """

#
#  Extract haplotype-specific genomes ('alternate references') in FASTA format
#  (taking only the SNPs to keep genomic coordinates invariant)
#

rule vcf_consensus_DGRP57_DGRP439:
  input:
    "data/external_WGS/vcf/" + genome + "_DGRP-57_DGRP-439.SNPs.vcf.gz"
  output:
    "data/external_WGS/fasta/" + genome + "_DGRP-57_DGRP-439_{haplotype}.fa.gz"
  wildcard_constraints:
    haplotype = "maternal|paternal"
  params:
    fasta_genome = config["genomes"][genome]["fasta"],
    requested_haplotype = lambda wildcards: "DGRP-57" if wildcards.haplotype == "maternal" else "DGRP-439"
  conda:
    "../../env/vcftools.yaml"
  shell:
    """
    cat {params.fasta_genome} | vcf-consensus -s {params.requested_haplotype} {input} | gzip > {output}
    """
