# pcn_tet_susceptibility_gwas
Data and code for "Genome wide association studies define minimal set of loci for prediction of penicillin and tetracycline susceptibility in Neisseria gonorrhoeae"

Files > 100MB unavailable through GitHub repository (e.g. assemblies, alignments, unitigs). These are available upon request from mortimer@hsph.harvard.edu

## conda_envs/
Conda environments for GWAS software

## data/

### gwas/

#### references.txt
File describing reference for Manhattan plot

#### pcn/
PCN conditional GWAS input and output

#### tet/
TET conditional GWAS input and output

### prediction/
Input files required for predicting PCN and TET susceptibility based on identified loci. Used by `susceptibility_prediction.R`

* pcn_tet_genotype_phenotype.tsv

### validation/
Input files required for predicting PCN and TET susceptibility based on identified loci in validationd dataset. Used by `validation.R`

* 2021-06-24_reimche_gisp_alleles.tsv
* reimche_metadata.tsv

## scripts/
Scripts for processing GWAS output and predicting phenotype from genotype.
* filter_significant_unitigs.R
* manhattan_plot.R
* susceptibility_prediction.R
* validation.R

## Snakefile
Snakefile describing GWAS pipeline
