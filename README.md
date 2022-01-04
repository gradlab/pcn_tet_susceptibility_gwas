# pcn_tet_susceptibility_gwas
Data and code for "Genome wide association studies define minimal set of loci for prediction of penicillin and tetracycline susceptibility in Neisseria gonorrhoeae"

Accessions for the global dataset can be found in `data/prediction/pcn_tet_genotype_phenotype.tsv`, and accessions for the validation dataset can be found in `data/validation/2021-06-24_reimche_gisp_alleles.tsv`.

Files > 100MB unavailable through GitHub repository (e.g. assemblies, alignments, unitig and kmer presence/absence). These are available upon request from mortimer@hsph.harvard.edu

## conda_envs/
Conda environments for GWAS software

## data/

### gwas/

#### references.txt
File describing reference for Manhattan plot

#### pcn/
PCN conditional GWAS input and output

* gubbins/ : contains newick formatted tree file from Gubbins
* kmer_significance_limits.txt: Bonferroni corrected significance limits for kmer based GWAS
* pcn_kmer_filtered.txt: pyseer output with unitig significance, due to file size, kmers have been filtered to those passing significance threshold
* pcn_metadata.tsv: metadata used for covariates in GWAS
* pcn_mics.tsv: penicillin MICs used as phenotype in GWAS
* pcn_unitig_annotated_WHO_N.txt: unitigs annotated with position in WHO-N reference genome
* pcn_unitig_signfiicance.txt: pyseer output with unitig significance
* significance_limits.txt: Bonferroni corrected significance limits for unitig based GWAS
* unitig_patterns.txt: Unitig presence absence patterns used to calculate p-value thershold

#### tet/
TET conditional GWAS input and output

* gubbins/ : contains newick formatted tree file from Gubbins
* kmer_significance_limits.txt: Bonferroni corrected significance limits for kmer based GWAS
* tet_kmer_filtered.txt: pyseer output with unitig significance, due to file size, kmers have been filtered to those passing significance threshold
* tet_metadata.tsv: metadata used for covariates in GWAS
* tet_mics.tsv: tetracycline MICs used as phenotype in GWAS
* tet_unitig_annotated_WHO_N.txt: unitigs annotated with position in WHO-N reference genome
* tet_unitig_signfiicance.txt: pyseer output with unitig significance
* significance_limits.txt: Bonferroni corrected significance limits for unitig based GWAS
* unitig_patterns.txt: Unitig presence absence patterns used to calculate p-value thershold


### prediction/
Input files required for predicting PCN and TET susceptibility based on identified loci. Used by `susceptibility_prediction.R`

* pcn_tet_genotype_phenotype.tsv
* qc.txt

### validation/
Input files required for predicting PCN and TET susceptibility based on identified loci in validationd dataset. Used by `validation.R`

* 2021-06-24_reimche_gisp_alleles.tsv
* reimche_metadata.tsv

## scripts/
Scripts for processing GWAS output and predicting phenotype from genotype.
* filter_significant_unitigs.R
* manhattan_plot.R
* plot_qc.R
* susceptibility_prediction.R
* validation.R

## Snakefile
Snakefile describing GWAS pipeline
