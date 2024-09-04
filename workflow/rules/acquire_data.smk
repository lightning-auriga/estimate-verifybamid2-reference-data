rule download_reference_data:
    """
    Conditionally acquire reference data files from either remote
    source (e.g. S3) or somewhere else on local filesystem.

    All the interesting stuff happens in the input mapping function,
    which uses the path contents between 'reference_data' and the actual
    filename to determine which configured reference file to pull.

    In some instances, the remote source will be gzipped but the local version
    won't be; in that instance, decompress midflight.

    I'm giving up on remotes, because the FTP one is unusable with conda env creation.
    """
    output:
        step1=temp("results/reference_genome/reference.fasta.staging"),
        step2="results/reference_genome/reference.fasta",
    benchmark:
        "results/performance_benchmarks/download_reference_data/results.tsv"
    params:
        genome=config["reference-fasta"],
    conda:
        "../envs/awscli.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        'if [[ "{params}" == "s3://"* ]] ; then aws s3 cp {params} {output.step1} ; '
        'elif [[ "{params}" == "http://"* ]] || [[ "{params}" == "https://"* ]] || [[ "{params}" == "ftp://"* ]] ; then wget -O {output.step1} {params} ; '
        "else cp {params} {output.step1} ; fi ; "
        'if [[ "{params}" = *".gz" ]] ; then cat {output.step1} | gunzip -c > {output.step2} ; '
        "else cp {output.step1} {output.step2} ; fi"


rule index_fasta:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.fai",
    benchmark:
        "results/performance_benchmarks/index_fasta/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "samtools faidx {input}"


rule download_input_vcf:
    output:
        temp("results/input_data/{filename}.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/download_input_vcf/{filename}.tsv"
    params:
        target=lambda wildcards: manifest.loc[wildcards.filename, "vcf"],
    conda:
        "../envs/awscli.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        'if [[ "{params}" == "s3://"* ]] ; then aws s3 cp {params} {output} ; '
        'elif [[ "{params}" == "http://"* ]] || [[ "{params}" == "https://"* ]] || [[ "{params}" == "ftp://"* ]] ; then wget -O {output} {params} ; '
        "else cp {params} {output} ; fi"


rule index_vcf:
    """
    Use tabix to generate tbi file for vcf.gz input
    """
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    benchmark:
        "results/performance_benchmarks/index_vcf/{prefix}.tsv"
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "tabix -p vcf {input}"
