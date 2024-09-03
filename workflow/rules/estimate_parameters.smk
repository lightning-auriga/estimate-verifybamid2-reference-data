rule estimate_verify_parameters:
    input:
        vcf="results/combine_vcfs/analysis-ready.vcf.gz",
        tbi="results/combine_vcfs/analysis-ready.vcf.gz.tbi",
        fasta="results/reference_genome/reference.fasta",
        fai="results/reference_genome/reference.fasta.fai",
    output:
        "results/combine_vcfs/analysis-ready.vcf.gz.UD",
        "results/combine_vcfs/analysis-ready.vcf.gz.mu",
        "results/combine_vcfs/analysis-ready.vcf.gz.bed",
        "results/combine_vcfs/analysis-ready.vcf.gz.V",
    conda:
        "../envs/verifybamid2.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        slurm_partition="comp",
    shell:
        "verifybamid2 --RefVCF {input.vcf} --Reference {input.fasta}"
