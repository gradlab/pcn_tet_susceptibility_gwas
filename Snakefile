localrules:
    all,
    filter_significant,


rule all:
    input:
        "data/figures/manhattan_plot.pdf"


rule gubbins:
    input:
        "data/gwas/{abx}/alignments/gc_{abx}_pseudogenomes.fasta",
    output:
        "data/gwas/{abx}/gubbins/gc_{abx}.final_tree.tre",
    params:
        directory="data/gwas/{abx}/gubbins/"
    resources:
        cpus=12,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
        time=lambda wildcards, attempt: attempt * 8000,
    conda:
        "conda_envs/gubbins.yml"
    shell:
        """
        mkdir -p {params.directory}
        run_gubbins.py --threads {resources.cpus} --prefix data/gwas/{abx}/gubbins/gc_{abx} {input}
        """


rule unitigs:
    input:
        strain_list="data/gwas/{abx}/unitigs/strain_list.txt",
    output:
        "data/gwas/{abx}/unitigs/unitigs/unitigs.txt",
    params:
        directory="data/gwas/{abx}/unitigs/unitigs/",
    resources:
        cpus=4,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 480,
    conda:
        "conda_envs/unitig-counter.yml"
    shell:
        """
        unitig-counter -strains {input.strain_list} -output unitigs/ -nb-cores {resources.cpus}
        mv unitigs/* {params.directory}
        rmdir unitigs/
        """


rule similarity_matrix:
    input:
        tree="data/gwas/{abx}/gubbins/gc_{abx}.final_tree.tre",
    output:
        matrix="data/gwas/{abx}/gubbins/similarity_matrix.txt",
    resources:
        cpus=1,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 30,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        "python software/pyseer/scripts/phylogeny_distance.py --lmm {input.tree} > {output.matrix}"


rule lmm_gwas:
    input:
        pheno="data/gwas/{abx}/{abx}_mics.tsv",
        unitigs="data/gwas/{abx}/unitigs/unitigs/unitigs.txt",
        similarity="data/gwas/{abx}/gubbins/similarity_matrix.txt",
        metadata="data/gwas/{abx}/{abx}_metadata.tsv"
    output:
        patterns="data/gwas/{abx}/unitig_patterns.txt",
        significance="data/gwas/{abx}/{abx}_unitig_significance.txt",
    params:
        pheno_column="{abx}_log2"
    resources:
        cpus=1,
        mem_mb=lambda wildcards, attempt: attempt * 10000,
        time=lambda wildcards, attempt: attempt * 360,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        pyseer --lmm --phenotypes {input.pheno} --similarity {input.similarity} --uncompressed --kmers {input.unitigs} --phenotype-column {params.pheno_column} --covariates {input.metadata} --use-covariates 2 3 4 --output-patterns {output.patterns} --cpu 1 > {output.significance}
        """


rule pyseer_count_patterns:
    input:
        "data/gwas/{abx}/unitig_patterns.txt",
    output:
        "data/gwas/{abx}/significance_limits.txt",
    resources:
        cpus=1,
        mem_mb=1000,
        time=10,
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        python software/pyseer/scripts/count_patterns.py {input} > {output}
        """


rule filter_significant:
    input:
        script="scripts/filter_significant_unitigs.R",
        limit="data/gwas/{abx}/significance_limits.txt",
        unitig_significance="data/gwas/{abx}/{abx}_unitig_significance.txt",
    output:
        "data/gwas/{abx}/{abx}_unitig_significance_filtered.txt",
    shell:
        """
        load_R
        load_R_packages
        Rscript scripts/filter_significant_unitigs.R {input.limit} {input.unitig_significance} {output}
        """

rule annotate_unitigs:
    input:
        unitig_significance="data/gwas/{abx}/{abx}_unitig_significance.txt",
        reference="data/gwas/references.txt"
    output:
        "data/gwas/{abx}/{abx}_unitig_annotated_WHO_N.txt",
    resources:
        cpus=1,
        mem_mb=1000,
        time=10
    conda:
        "conda_envs/pyseer.yml"
    shell:
        """
        annotate_hits_pyseer {input.unitig_significance} {input.reference} {output}
        """

rule manhattan_plot:
    input:
        expand("data/gwas/{abx}/{abx}_unitig_annotated_WHO_N.txt", abx=["pcn","tet"]),
        expand("data/gwas/{abx}/significance_limits.txt", abx=["pcn", "tet"])
    output:
        "data/figures/manhattan_plot.pdf"
    shell:
        """
        load_R
        load_R_packages
        Rscript manhattan_plot.R
        """