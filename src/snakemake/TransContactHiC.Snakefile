configfile: "config/config.yml"

DATASETS = list(config["haplotype_HiC"]["datasets"]) + list(config["haplotype_HiC"]["combined_dataset_to_dataset"])

wildcard_constraints:
  resolution = "\d+|rs",
  phasing_mode = "phased|fully_phased"

#
#  Default targets
#

rule all:
  input:
    lambda wildcards:
      rules.all_multiqc.input + \
      rules.all_trans_pairs.input + \
      rules.all_reads_in_selected_regions.input + \
      rules.all_hicexplorer.input + \
      rules.all_mcool.input

rule all_multiqc:
  input:
    "data/external_HiC/qc/multiqc_report.html",
    "data/hicexplorer/qc/multiqc_report.html"

rule all_trans_pairs:
  input:
    expand("data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_WE_2-4h_{replicate}.pairs.{ext}.tsv.gz",
      replicate = ['Rep1', 'Rep2'], ext = ['phased', 'fully_phased']),
    expand("data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_PnM_{replicate}.pairs.{ext}.tsv.gz",
      replicate = ['Rep1', 'Rep2'], ext = ['phased', 'fully_phased'])

rule all_reads_in_selected_regions:
  input:
    expand("data/haplotype_HiC/bam/reads_in_selected_regions/HiC_D_mel_DGRP-57_DGRP-439_WE_2-4h_{replicate}.pairs.2R_22592501-22595000_{hapl1}.2R_22650001-22652500_{hapl2}.sort.bam.bai",
      replicate = ['Rep1', 'Rep2'], hapl1 = ['DGRP-57', 'DGRP-439'], hapl2 = ['DGRP-57', 'DGRP-439']),
    expand("data/haplotype_HiC/bam/reads_in_selected_regions/HiC_D_mel_DGRP-57_DGRP-439_WE_2-4h_{replicate}.pairs.2R_22592501-22595000_{hapl1}.2R_22650001-22652500_{hapl2}.sort.bam.bai",
      replicate = ['Rep1', 'Rep2'], hapl1 = ['DGRP-57', 'DGRP-439'], hapl2 = ['DGRP-57', 'DGRP-439'])

rule all_hicexplorer:
  input:
    list([
      [
        "data/hicexplorer/h5/" + dataset + "_" + resolution + ".h5",
        "data/hicexplorer/h5/" + dataset + "_" + resolution + "_corrected.h5",
      ] for resolution in ["5000", "rs"]
    ] for dataset in DATASETS)

rule all_mcool:
  input:
    list([
      "data/hicexplorer/cool/" + dataset + "_5000.mcool",
    ] for dataset in DATASETS)

#
#  SRA data download
#

rule download_sra:
  output:
    "data/external_HiC/fastq/SRR{sample}_1.fastq.gz",
    "data/external_HiC/fastq/SRR{sample}_2.fastq.gz"
  conda:
    "../../env/sra-tools.yaml"
  shell:
    """
    fastq-dump --split-files --origfmt --gzip -O data/external_HiC/fastq SRR{wildcards.sample}
    """

#
#  FASTQ quality control
#

rule fastqc:
  input:
    "data/external_HiC/fastq/{file}.fastq.gz"
  output:
    "data/external_HiC/qc/{file}_fastqc.zip"
  conda:
    "../../env/fastqc.yaml"
  shell:
    """
    fastqc -o data/external_HiC/qc {input}
    """

def library_fastq_qcfiles(wildcards):
  zipfiles_fastq = ["data/external_HiC/fastq/" + library_fastq
    for libraries in config["haplotype_HiC"]["dataset_to_library_fastq"].values()
      for library_fastq_pair in libraries
        for library_fastq in library_fastq_pair.values()]
  zipfiles_sra = ["data/external_HiC/qc/" + library_sra + "_" + readindex + "_fastqc.zip"
    for libraries in config["haplotype_HiC"]["dataset_to_library_sra"].values()
      for library_sra in libraries
        for readindex in ["1", "2"]]
  return zipfiles_fastq + zipfiles_sra

rule multiqc_fastqc:
  input:
    library_fastq_qcfiles
  output:
    "data/external_HiC/qc/multiqc_report.html"
  conda:
    "../../env/multiqc.yaml"
  shell:
    """
    cd data/external_HiC/qc
    rm -rf multiqc_data/
    multiqc --interactive .
    """

#
#  Map Hi-C reads, merging sequencing libraries for the same replicate
#

def dataset_to_genome(dataset):
  dataset_wo_replicate = dataset.split('_Rep')[0]
  return config["haplotype_HiC"]["dataset_to_genome"][dataset_wo_replicate]

def dataset_to_fasta_genome(dataset):
  genome = dataset_to_genome(dataset)
  return config["genomes"][genome]["fasta"]

def dataset_to_fastq(wildcards):
  if wildcards.dataset in config["haplotype_HiC"]["dataset_to_library_fastq"]:
    library_fastqs = [pair["read" + wildcards.read]
      for pair in config["haplotype_HiC"]["dataset_to_library_fastq"][wildcards.dataset]]
    return ["data/external_HiC/fastq/" + library_fastq for library_fastq in library_fastqs]
  elif wildcards.dataset in config["haplotype_HiC"]["dataset_to_library_sra"]:
    library_sras = config["haplotype_HiC"]["dataset_to_library_sra"][wildcards.dataset]
    return ["data/external_HiC/fastq/" + sra + "_" + wildcards.read + ".fastq.gz" for sra in library_sras]
  else:
    return "/missing_fastq"

rule bwa_mem:
  input:
    dataset_to_fastq
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}_{read}.bam")
  wildcard_constraints:
    read = "1|2",
    dataset = "HiC_D_mel_DGRP-57_DGRP-439_WE_2-4h_Rep1|HiC_D_mel_DGRP-57_DGRP-439_WE_2-4h_Rep2|HiC_D_mel_DGRP-57_DGRP-439_PnM_Rep1|HiC_D_mel_DGRP-57_DGRP-439_PnM_Rep2"
  params:
    fasta_genome = lambda wildcards: dataset_to_fasta_genome(wildcards.dataset)
  threads:
    16
  conda:
    "../../env/bwa_samtools.yaml"
  shell:
    """
    zcat {input} | bwa mem -E50 -L0 -5 -t {threads} {params.fasta_genome} - | samtools view -bT {params.fasta_genome} - > {output}
    """

#
#  Merge Hi-C reads from separate files for read1 and read2 to single files
#

rule samtools_namesort:
  input:
    "data/external_HiC/bam/all_reads/{dataset}.bam"
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}.ns.bam")
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    rm -f {output}.tmp.*
    samtools sort -n -@ {threads} -m 4G -O bam {input} -o {output}
    """

rule set_bam_flag_read1:
  input:
    "data/external_HiC/bam/all_reads/{dataset}_1.ns.bam"
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}_1.ns.addFlag.bam")
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/py/set_bam_flag.py -i {input} 65 -o {output} # read paired (0x1), first in pair (0x40)
    """

rule set_bam_flag_read2:
  input:
    "data/external_HiC/bam/all_reads/{dataset}_2.ns.bam"
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}_2.ns.addFlag.bam")
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/py/set_bam_flag.py -i {input} 129 -o {output} # read paired (0x1), second in pair (0x80)
    """

rule samtools_merge:
  input:
    "data/external_HiC/bam/all_reads/{dataset}_1.ns.addFlag.bam",
    "data/external_HiC/bam/all_reads/{dataset}_2.ns.addFlag.bam"
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}.merge.bam")
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools merge -n -@ {threads} {output} {input}
    """

rule samtools_fixmate:
  input:
    "data/external_HiC/bam/all_reads/{dataset}.merge.bam"
  output:
    temp("data/external_HiC/bam/all_reads/{dataset}.merge.fixmate.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools fixmate -p {input} {output}
    """

#
#  Remove PCR and optical duplicates, annotate and separate haplotypes
#

rule picard_rmdup_DGRP57_DGRP439:
  input:
    "data/external_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.bam"
  output:
    bam = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.bam",
    metrics = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.metrics.txt"
  log:
    "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.log"
  params:
    tmpdir = config["paths"]["tmpdir"]
  conda:
    "../../env/picard-slim.yaml"
  shell:
    """
    picard -Xmx32g MarkDuplicates -INPUT {input} -OUTPUT {output.bam} -REMOVE_DUPLICATES true -READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' -ASSUME_SORT_ORDER queryname -VALIDATION_STRINGENCY SILENT -METRICS_FILE {output.metrics} --MAX_RECORDS_IN_RAM 5000000 -TMP_DIR {params.tmpdir} &> {log}
    """

#
#  Process Hi-C reads using TransContactHiC
#

rule pairtools_parse:
  input:
    "data/haplotype_HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.bam"
  output:
    "data/haplotype_HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.pairs.gz"
  params:
    genome = lambda wildcards: dataset_to_genome(wildcards.dataset),
    fasta_genome = lambda wildcards: config["genomes"][dataset_to_genome(wildcards.dataset)]["fasta"]
  conda:
    "../../env/pairtools.yaml"
  shell:
    """
    pairtools parse {input} \
      --assembly {params.genome} \
      --chroms-path {params.fasta_genome}.fai \
      --walks-policy all \
      --add-columns pos5,pos3 \
      --drop-sam \
      --no-flip \
      --output {output}
    """

rule TransContactHiC_DGRP57_DGRP439:
  input:
    bam = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.bam",
    pairs = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.pairs.gz"
  output:
    pairs_tsv = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.pairs.{phasing_mode}.tsv.gz",
    pairs_stats = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.pairs.{phasing_mode}.stats"
  log:
    "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.pairs.{phasing_mode}.log"
  params:
    SNPs = config["haplotype_HiC"]["DGRP57_DGRP439_SNPs"],
    mode_params = lambda wildcards: "--fully-phased" if wildcards.phasing_mode == "fully_phased" else ""
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/py/TransContactHiC.py \
      --variants {params.SNPs} \
      --input-bam {input.bam} \
      --input-pairs {input.pairs} \
      --output-tsv {output.pairs_tsv} \
      --output-stats {output.pairs_stats} \
      {params.mode_params} &> {log}
    """

#
#  Extract the reads corresponding to contacts between specified haplotypes in selected regions
#  (also sort and index the resulting .BAM files)
#

rule reads_in_selected_regions_DGRP57_DGRP439:
  input:
    bam = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.bam",
    pairs = "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.pairs.gz"
  output:
    pairs_tsv = "data/haplotype_HiC/bam/reads_in_selected_regions/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.pairs.{chrom1}_{start1}-{end1}_{hapl1}.{chrom2}_{start2}-{end2}_{hapl2}.tsv.gz",
    bam = "data/haplotype_HiC/bam/reads_in_selected_regions/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.pairs.{chrom1}_{start1}-{end1}_{hapl1}.{chrom2}_{start2}-{end2}_{hapl2}.bam"
  wildcard_constraints:
    start1 = "\d+",
    end1 = "\d+",
    hapl1 = "[^.]+",
    start2 = "\d+",
    end2 = "\d+",
    hapl2 = "[^.]+"
  params:
    SNPs = config["haplotype_HiC"]["DGRP57_DGRP439_SNPs"],
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/py/TransContactHiC.py \
      --variants {params.SNPs} \
      --input-bam {input.bam} \
      --input-pairs {input.pairs} \
      --output-tsv {output.pairs_tsv} \
      --output-bam {output.bam} \
      --region {wildcards.chrom1}:{wildcards.start1}-{wildcards.end1}:{wildcards.hapl1} \
      --region {wildcards.chrom2}:{wildcards.start2}-{wildcards.end2}:{wildcards.hapl2}
    """

rule samtools_sort:
  input:
    "data/haplotype_HiC/bam/reads_in_selected_regions/{dataset}.bam"
  output:
    "data/haplotype_HiC/bam/reads_in_selected_regions/{dataset}.sort.bam"
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    rm -f {output}.tmp.*
    samtools sort -@ {threads} -m 4G -O bam {input} -o {output}
    """

rule samtools_index:
  input:
    "data/haplotype_HiC/bam/reads_in_selected_regions/{dataset}.bam"
  output:
    "data/haplotype_HiC/bam/reads_in_selected_regions/{dataset}.bam.bai"
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools index {input}
    """

#
#  Extract same-haplotype Hi-C reads for obtaining haplotype-specific contact maps
#

rule extract_haplotype_DGRP57_DGRP439:
  input:
    "data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{dataset}_{replicate}.merge.fixmate.rmdup.bam"
  output:
    temp("data/haplotype_HiC/bam/all_reads/HiC_D_mel_DGRP-57_DGRP-439_{haplotype}_{dataset}_{replicate}.merge.fixmate.rmdup.bam")
  wildcard_constraints:
    haplotype = "maternal|paternal"
  params:
    SNPs = config["haplotype_HiC"]["DGRP57_DGRP439_SNPs"],
    requested_haplotype = lambda wildcards: "DGRP-57" if wildcards.haplotype == "maternal" else "DGRP-439"
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/py/extract_haplotype.py --variants {params.SNPs} --input {input} --output {output} --requested-haplotype {params.requested_haplotype}
    """

ruleorder:
  extract_haplotype_DGRP57_DGRP439 > bwa_mem

#
#  Split Hi-C reads from single files to separate files for read1 and read2
#

rule samtools_extract_read1:
  input:
    "data/haplotype_HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.bam"
  output:
    temp("data/haplotype_HiC/bam/all_reads/{dataset}_R1.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools view -b -f 64 -o {output} {input}
    """

rule samtools_extract_read2:
  input:
    "data/haplotype_HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.bam"
  output:
    temp("data/haplotype_HiC/bam/all_reads/{dataset}_R2.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools view -b -f 128 -o {output} {input}
    """

#
#  Process Hi-C reads using HiCExplorer
#

def dataset_to_restrictionEnzyme(dataset):
  dataset_wo_replicate = dataset.split('_Rep')[0]
  return config["haplotype_HiC"]["dataset_to_restrictionEnzyme"][dataset_wo_replicate]

def restrictionEnzyme_to_restrictionSequence(restrictionEnzyme):
  return config["haplotype_HiC"]["restrictionEnzyme_to_restrictionSequence"][restrictionEnzyme]

def dataset_to_restrictionSequence(dataset):
  return config["haplotype_HiC"]["restrictionEnzyme_to_restrictionSequence"][dataset_to_restrictionEnzyme(dataset)]

def restrictionEnzyme_to_danglingSequence(restrictionEnzyme):
  return config["haplotype_HiC"]["restrictionEnzyme_to_danglingSequence"][restrictionEnzyme]

def dataset_to_danglingSequence(dataset):
  return config["haplotype_HiC"]["restrictionEnzyme_to_danglingSequence"][dataset_to_restrictionEnzyme(dataset)]

rule hicFindRestSite:
  output:
    "results/rest_site_positions_{restrictionEnzyme}_{genome}.bed"
  wildcard_constraints:
    restrictionEnzyme = "[^_]+"
  params:
    fasta_genome = lambda wildcards: config["genomes"][wildcards.genome]["fasta"],
    restrictionSequence = lambda wildcards: '"' + restrictionEnzyme_to_restrictionSequence(wildcards.restrictionEnzyme).replace(" ", "|") + '"'
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicFindRestSite --fasta {params.fasta_genome} --searchPattern {params.restrictionSequence} -o {output}
    """

rule hicBuildMatrix:
  input:
    bam1 = "data/haplotype_HiC/bam/all_reads/{dataset}_R1.bam",
    bam2 = "data/haplotype_HiC/bam/all_reads/{dataset}_R2.bam",
    bed = lambda wildcards: "results/rest_site_positions_" + dataset_to_restrictionEnzyme(wildcards.dataset) + "_" + dataset_to_genome(wildcards.dataset) + ".bed"
  output:
    bam = "data/hicexplorer/bam/filtered_reads/{dataset}_{resolution}.bam",
    h5 = "data/hicexplorer/h5/{dataset}_{resolution}.h5",
    qc = "data/hicexplorer/qc/{dataset}_{resolution}/hicQC.html"
  params:
    restrictionSequence = lambda wildcards: dataset_to_restrictionSequence(wildcards.dataset),
    danglingSequence = lambda wildcards: dataset_to_danglingSequence(wildcards.dataset),
    resolution_spec = lambda wildcards:
      "" if wildcards.resolution == "rs" else "--binSize " + wildcards.resolution
  threads:
    8
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicBuildMatrix --samFiles {input.bam1} {input.bam2} --outBam {output.bam} --outFileName {output.h5} \
      --restrictionSequence {params.restrictionSequence} \
      --danglingSequence {params.danglingSequence} \
      --restrictionCutFile {input.bed} \
      --skipDuplicationCheck \
      {params.resolution_spec} \
      --QCfolder data/hicexplorer/qc/{wildcards.dataset}_{wildcards.resolution} \
      --threads {threads}
    """

def combined_dataset_to_dataset_h5(wildcards):
  if wildcards.dataset in config["haplotype_HiC"]["combined_dataset_to_dataset"]:
    datasets = config["haplotype_HiC"]["combined_dataset_to_dataset"][wildcards.dataset]
    return ["data/hicexplorer/h5/" + dataset + "_" + wildcards.resolution + ".h5" for dataset in datasets]
  else:
    return "/missing_dataset"

rule hicSumMatrices:
  input:
    combined_dataset_to_dataset_h5
  output:
    "data/hicexplorer/h5/{dataset}_{resolution}.h5"
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicSumMatrices -m {input} -o {output}
    """

ruleorder:
  hicSumMatrices > hicBuildMatrix

def dataset_to_chromosomes(dataset):
  genome = dataset_to_genome(dataset)
  return config["genomes"][genome]["chromosomes"]

rule hicCorrectMatrix:
  input:
    "data/hicexplorer/h5/{dataset}_{resolution}.h5"
  output:
    png = "data/hicexplorer/qc/{dataset}_{resolution}_corrected.png",
    h5 = "data/hicexplorer/h5/{dataset}_{resolution}_corrected.h5"
  params:
    chromosomes = lambda wildcards: dataset_to_chromosomes(wildcards.dataset),
    filterThreshold = config["haplotype_HiC"]["filterThreshold"]
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    # note that the first invocation of hicCorrectMatrix modifies the input file(!)
    hicCorrectMatrix diagnostic_plot --chromosomes {params.chromosomes} -m {input} -o {output.png}
    hicCorrectMatrix correct --chromosomes {params.chromosomes} -m {input} --filterThreshold {params.filterThreshold} -o {output.h5}
    """

# when combining replicates: first sum raw count matrices, then normalize the resulting matrix
ruleorder:
  hicCorrectMatrix > hicSumMatrices

rule hicConvertFormat_to_cool:
  input:
    "data/hicexplorer/h5/{dataset}.h5"
  output:
    temp("data/hicexplorer/cool/{dataset}.cool")
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicConvertFormat --matrices {input} --outFileName {output} --inputFormat h5 --outputFormat cool
    """

rule multiqc_hicexplorer:
  input:
    ["data/hicexplorer/qc/" + dataset + "_rs/hicQC.html" for dataset in config["haplotype_HiC"]["datasets"]]
  output:
    "data/hicexplorer/qc/multiqc_report.html"
  params:
    directories = [dataset + "_rs" for dataset in config["haplotype_HiC"]["datasets"]]
  conda:
    "../../env/multiqc.yaml"
  shell:
    """
    cd data/hicexplorer/qc
    rm -rf multiqc_haplotype_HiC_data/
    multiqc --interactive {params.directories}
    """

#
#  Convert HiCExplorer output to HiGlass .mcool files
#

rule cooler_zoomify:
  input:
    "data/hicexplorer/cool/{dataset}_{resolution}.cool"
  output:
    "data/hicexplorer/cool/{dataset}_{resolution}.mcool"
  wildcard_constraints:
    resolution = "1000|5000"
  params:
    resolution_spec = lambda wildcards: "1000N" if wildcards.resolution == "1000" else "5000,10000N"
  conda:
    "../../env/cooler.yaml"
  shell:
    """
    cooler zoomify --resolutions {params.resolution_spec} --balance {input}
    """
