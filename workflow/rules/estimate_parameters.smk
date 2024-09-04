rule estimate_verify_parameters:
    input:
        vcf="results/combine_vcfs/analysis-ready.vcf.gz",
        tbi="results/combine_vcfs/analysis-ready.vcf.gz.tbi",
        fasta="results/reference_genome/reference.fasta",
        fai="results/reference_genome/reference.fasta.fai",
    output:
        ud=temp("results/combine_vcfs/analysis-ready.vcf.gz.UD"),
        mu=temp("results/combine_vcfs/analysis-ready.vcf.gz.mu"),
        bed=temp("results/combine_vcfs/analysis-ready.vcf.gz.bed"),
        v=temp("results/combine_vcfs/analysis-ready.vcf.gz.V"),
    conda:
        "../envs/verifybamid2.yaml"
    threads: config_resources["verifybamid2"]["threads"]
    resources:
        mem_mb=config_resources["verifybamid2"]["memory"],
        slurm_partition=rc.select_partition(
            config_resources["verifybamid2"]["partition"], config_resources["partitions"]
        ),
    shell:
        "verifybamid2 --RefVCF {input.vcf} --Reference {input.fasta}"


localrules:
    copy_verify_files,


rule copy_verify_files:
    input:
        ud="results/combine_vcfs/analysis-ready.vcf.gz.UD",
        mu="results/combine_vcfs/analysis-ready.vcf.gz.mu",
        bed="results/combine_vcfs/analysis-ready.vcf.gz.bed",
        v="results/combine_vcfs/analysis-ready.vcf.gz.V",
    output:
        ud="results/verify_output/" + config["dataset-name"] + ".UD",
        mu="results/verify_output/" + config["dataset-name"] + ".mu",
        bed="results/verify_output/" + config["dataset-name"] + ".bed",
        v="results/verify_output/" + config["dataset-name"] + ".V",
    threads: 1
    shell:
        "cp {input.ud} {output.ud} && "
        "cp {input.mu} {output.mu} && "
        "cp {input.bed} {output.bed} && "
        "cp {input.v} {output.v}"
