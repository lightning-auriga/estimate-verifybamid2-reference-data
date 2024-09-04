rule annotate_and_shrink:
    input:
        vcf="results/input_data/{filename}.vcf.gz",
        tbi="results/input_data/{filename}.vcf.gz.tbi",
    output:
        vcf="results/annotate_and_shrink/{filename}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["bcftools"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bcftools annotate --threads {threads} --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -x ^FORMAT/GT -Ou {input.vcf} | "
        "bcftools view --threads {threads} -M2 -q 0.05 -v snps -f 'PASS,' -Oz -o {output}"


rule report_markers:
    input:
        vcf="results/annotate_and_shrink/{filename}.vcf.gz",
        tbi="results/annotate_and_shrink/{filename}.vcf.gz.tbi",
    output:
        tsv="results/report_markers/{filename}.tsv",
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bcftools query -f '%ID\\n' {input.vcf} > {output.tsv}"


rule select_downsampled_markers:
    input:
        expand("results/report_markers/{filename}.tsv", filename=manifest.index),
    output:
        "results/select_downsampled_markers/target_variants.tsv",
    params:
        variant_count=config["downsampled-variant-count"],
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "cat {input} | shuf -n {params.variant_count} -o {output}"


rule subset_vcf:
    input:
        vcf="results/annotate_and_shrink/{filename}.vcf.gz",
        tbi="results/annotate_and_shrink/{filename}.vcf.gz.tbi",
        targets="results/select_downsampled_markers/target_variants.tsv",
    output:
        vcf="results/subset_vcf/{filename}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["default"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bcftools view -i 'ID=@{input.targets}' -Oz -o {output.vcf} {input.vcf}"


rule combine_vcfs:
    input:
        expand("results/subset_vcf/{filename}.vcf.gz", filename=manifest.index),
    output:
        "results/combine_vcfs/analysis-ready.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["bcftools"]["partition"], config_resources["partitions"]
        ),
    shell:
        "bcftools concat --threads {threads} --naive -Oz -o {output} {input}"
