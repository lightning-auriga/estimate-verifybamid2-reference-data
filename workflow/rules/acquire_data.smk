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
    params:
        genome=config["reference-fasta"],
    conda:
        "../envs/awscli.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
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
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
    shell:
        "samtools faidx {input}"


rule download_input_vcf:
    output:
        temp("results/input_data/{filename}.vcf.gz"),
    params:
        target=lambda wildcards: manifest.loc[wildcards.filename, "vcf"],
    conda:
        "../envs/awscli.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
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
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
    shell:
        "tabix -p vcf {input}"
