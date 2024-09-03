rule annotate_and_shrink:
    input:
        vcf="results/input_data/{filename}.vcf.gz",
        tbi="results/input_data/{filename}.vcf.gz.tbi",
    output:
        vcf="results/annotate_and_shrink/{filename}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
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
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
    shell:
        "bcftools query -f '%ID\\n' {input.vcf} > {output.tsv}"


rule select_downsampled_markers:
    input:
        expand("results/report_markers/{filename}.tsv", filename=manifest.index),
    output:
        "results/select_downsampled_markers/target_variants.tsv",
    params:
        variant_count=config["downsampled-variant-count"],
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
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
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
    shell:
        "bcftools view -i 'ID=@{input.targets}' -Oz -o {output.vcf} {input.vcf}"


rule combine_vcfs:
    input:
        expand("results/subset_vcf/{filename}.vcf.gz", filename=manifest.index),
    output:
        "results/combine_vcfs/analysis-ready.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=3800,
        slurm_partition="spotshort",
    shell:
        "bcftools concat --threads {threads} --naive -Oz -o {output} {input}"
