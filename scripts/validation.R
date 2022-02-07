library(tidyverse)
library(cowplot)
library(infer)

# read in data
validation_metadata <- read_tsv("data/validation/reimche_metadata.tsv")
validation_resistance <- read_tsv("data/validation/2021-06-24_reimche_gisp_alleles.tsv")

# combine data

all_data <- validation_resistance %>%
    left_join(validation_metadata)

# add genotype interpretation
all_data <- all_data %>%
    mutate(PCN_genotype = if_else(penA_01 == 1 & blaTEM == 0, "S", "R")) %>%
    mutate(TET_genotype = if_else(RpsJ_57 == "V" & tetM == 0, "S", "R"))

# filter to valid TET MICs and categorize based on breakpoints
all_data <- all_data %>% 
    mutate(TET = str_replace_all(tetracycline, "[><= ]", "")) %>%
    mutate(TET = as.numeric(TET)) %>%
    mutate(TET_interpretation = case_when(TET <= 0.25 ~ "S", TET > 0.25 & TET < 2 ~ "I", TET >= 2 ~ "R")) %>%
    mutate(PCN = str_replace_all(penicillin, "[><= ]", "")) %>%
    mutate(PCN = as.numeric(PCN)) %>%
    mutate(PCN_interpretation = case_when(PCN <= 0.064 ~ "S", PCN > 0.064 & PCN < 2 ~ "I", PCN >= 2 ~ "R"))

# get genotype counts for all data
all_data %>% count(PCN_genotype) %>%
    pivot_wider(names_from = PCN_genotype, values_from = n) %>%
    mutate(PCN_total = S + R, S_percent = S*100/PCN_total, R_percent = R*100/PCN_total) %>%
    print(n = 30)

all_data %>% count(TET_genotype) %>%
    pivot_wider(names_from = TET_genotype, values_from = n) %>%
    mutate(TET_total = S + R, S_percent = S*100/TET_total, R_percent = R*100/TET_total) %>%
    print(n = 30)

all_data %>% filter(is.na(TET_genotype)) %>% select(wgs_id, tetracycline, TET_genotype, RpsJ_57, tetM)

# calculate sensitivity and specificity

# PCN

PCN_total_susceptible_phenotype <- all_data %>% filter(PCN_interpretation == "S") %>% nrow()
PCN_total_nonresistant_phenotype <- all_data %>% filter(PCN_interpretation != "R") %>% nrow()
PCN_total_susceptible_genotype <- all_data %>% filter(PCN_genotype == "S") %>% nrow()

# PCN sensitivity & specificity susceptible
PCN_true_positives <- all_data %>% filter(PCN_interpretation == "S" & PCN_genotype == "S") %>% nrow()
PCN_false_negatives <- all_data %>% filter(PCN_interpretation == "S" & PCN_genotype != "S") %>% nrow()
PCN_true_negatives <- all_data %>% filter(PCN_interpretation != "S" & PCN_genotype == "R") %>% nrow()
PCN_false_positives <- all_data %>% filter(PCN_interpretation != "S" & PCN_genotype == "S") %>% nrow()

PCN_sensitivity_susceptible <- (all_data %>% filter(PCN_interpretation == "S" & PCN_genotype == "S") %>% nrow())/PCN_total_susceptible_phenotype
PCN_specificity_susceptible <- (all_data %>% filter(PCN_interpretation != "S" & PCN_genotype == "R") %>% nrow())/(all_data %>% filter(PCN_interpretation != "S") %>% nrow())

PCN_sensitivity_susceptible + sqrt((PCN_sensitivity_susceptible*(1-PCN_sensitivity_susceptible))/PCN_true_positives)
PCN_sensitivity_susceptible - sqrt((PCN_sensitivity_susceptible*(1-PCN_sensitivity_susceptible))/PCN_true_positives)

PCN_specificity_susceptible + sqrt((PCN_specificity_susceptible*(1-PCN_specificity_susceptible))/PCN_true_positives)
PCN_specificity_susceptible - sqrt((PCN_specificity_susceptible*(1-PCN_specificity_susceptible))/PCN_true_positives)

# PCN sensitivity & specificity non-resistant

PCN_sensitivity_nonresistant <- (all_data %>% filter(PCN_interpretation != "R" & PCN_genotype == "S") %>% nrow())/PCN_total_nonresistant_phenotype
PCN_specificity_nonresistant <- (all_data %>% filter(PCN_interpretation == "R" & PCN_genotype == "R") %>% nrow())/(all_data %>% filter(PCN_interpretation == "R") %>% nrow())

PCN_true_positives <- all_data %>% filter(PCN_interpretation != "R" & PCN_genotype == "S") %>% nrow()
PCN_false_negatives <- all_data %>% filter(PCN_interpretation != "R" & PCN_genotype != "S") %>% nrow()
PCN_true_negatives <- all_data %>% filter(PCN_interpretation == "R" & PCN_genotype == "R") %>% nrow()
PCN_false_positives <- all_data %>% filter(PCN_interpretation == "R" & PCN_genotype == "S") %>% nrow()

PCN_sensitivity_nonresistant + sqrt((PCN_sensitivity_nonresistant*(1-PCN_sensitivity_nonresistant))/PCN_true_positives)
PCN_sensitivity_nonresistant - sqrt((PCN_sensitivity_nonresistant*(1-PCN_sensitivity_nonresistant))/PCN_true_positives)

PCN_specificity_nonresistant + sqrt((PCN_specificity_nonresistant*(1-PCN_specificity_nonresistant))/PCN_true_positives)
PCN_specificity_nonresistant - sqrt((PCN_specificity_nonresistant*(1-PCN_specificity_nonresistant))/PCN_true_positives)

PCN_total_susceptible_phenotype
PCN_total_nonresistant_phenotype
PCN_total_susceptible_genotype

PCN_sensitivity_susceptible
PCN_specificity_susceptible

PCN_sensitivity_nonresistant
PCN_specificity_nonresistant

# TET
TET_total_susceptible_phenotype <- all_data %>% filter(TET_interpretation == "S") %>% nrow()
TET_total_nonresistant_phenotype <- all_data %>% filter(TET_interpretation != "R") %>% nrow()
TET_total_susceptible_genotype <- all_data %>% filter(TET_genotype == "S") %>% nrow()

# TET sensitivity & specificity susceptible

TET_true_positives <- all_data %>% filter(TET_interpretation == "S" & TET_genotype == "S") %>% nrow()
TET_false_negatives <- all_data %>% filter(TET_interpretation == "S" & TET_genotype == "R") %>% nrow()
TET_true_negatives <- all_data %>% filter(TET_interpretation != "S" & TET_genotype == "R") %>% nrow()
TET_false_positives <- all_data %>% filter(TET_interpretation != "S" & TET_genotype == "S") %>% nrow()

TET_sensitivity_susceptible <- (all_data %>% filter(TET_interpretation == "S" & TET_genotype == "S") %>% nrow())/TET_total_susceptible_phenotype
TET_specificity_susceptible <- (all_data %>% filter(TET_interpretation != "S" & TET_genotype == "R") %>% nrow())/(all_data %>% filter(TET_interpretation != "S") %>% nrow())

TET_sensitivity_susceptible + sqrt((TET_sensitivity_susceptible*(1-TET_sensitivity_susceptible))/TET_true_positives)
TET_sensitivity_susceptible - sqrt((TET_sensitivity_susceptible*(1-TET_sensitivity_susceptible))/TET_true_positives)

TET_specificity_susceptible + sqrt((TET_specificity_susceptible*(1-TET_specificity_susceptible))/TET_true_positives)
TET_specificity_susceptible - sqrt((TET_specificity_susceptible*(1-TET_specificity_susceptible))/TET_true_positives)

# TET sensitivity & specificity non-resistant

TET_sensitivity_nonresistant <- (all_data %>% filter(TET_interpretation != "R" & TET_genotype == "S") %>% nrow())/TET_total_nonresistant_phenotype
TET_specificity_nonresistant <- (all_data %>% filter(TET_interpretation == "R" & TET_genotype == "R") %>% nrow())/(all_data %>% filter(TET_interpretation == "R") %>% nrow())

TET_true_positives <- all_data %>% filter(TET_interpretation != "R" & TET_genotype == "S") %>% nrow()
TET_false_negatives <- all_data %>% filter(TET_interpretation != "R" & TET_genotype == "R") %>% nrow()
TET_true_negatives <- all_data %>% filter(TET_interpretation == "R" & TET_genotype == "R") %>% nrow()
TET_false_positives <- all_data %>% filter(TET_interpretation == "R" & TET_genotype == "S") %>% nrow()

TET_sensitivity_nonresistant + sqrt((TET_sensitivity_nonresistant*(1-TET_sensitivity_nonresistant))/TET_true_positives)
TET_sensitivity_nonresistant - sqrt((TET_sensitivity_nonresistant*(1-TET_sensitivity_nonresistant))/TET_true_positives)

TET_specificity_nonresistant + sqrt((TET_specificity_nonresistant*(1-TET_specificity_nonresistant))/TET_true_positives)
TET_specificity_nonresistant - sqrt((TET_specificity_nonresistant*(1-TET_specificity_nonresistant))/TET_true_positives)

TET_total_susceptible_phenotype
TET_total_nonresistant_phenotype
TET_total_susceptible_genotype

TET_sensitivity_susceptible
TET_specificity_susceptible

TET_sensitivity_nonresistant
TET_specificity_nonresistant

# plasmid only prediction

# PCN
PCN_total_susceptible_phenotype <- all_data %>% filter(PCN_interpretation == "S") %>% nrow()
PCN_total_nonresistant_phenotype <- all_data %>% filter(PCN_interpretation != "R") %>% nrow()
PCN_total_susceptible_genotype <- all_data %>% filter(blaTEM == 0) %>% nrow()

PCN_true_positives <- all_data %>% filter(PCN_interpretation == "S" & blaTEM == 0) %>% nrow()
PCN_false_negatives <- all_data %>% filter(PCN_interpretation == "S" & blaTEM == 1) %>% nrow()
PCN_true_negatives <- all_data %>% filter(PCN_interpretation != "S" & blaTEM == 1) %>% nrow()
PCN_false_positives <- all_data %>% filter(PCN_interpretation != "S" & blaTEM == 0) %>% nrow()

PCN_sensitivity_susceptible <- (all_data %>% filter(PCN_interpretation == "S" & blaTEM == 0) %>% nrow())/PCN_total_susceptible_phenotype
PCN_specificity_susceptible <- (all_data %>% filter(PCN_interpretation != "S" & blaTEM == 1) %>% nrow())/(all_data %>% filter(PCN_interpretation != "S") %>% nrow())

PCN_sensitivity_nonresistant <- (all_data %>% filter(PCN_interpretation != "R" & blaTEM == 0) %>% nrow())/PCN_total_nonresistant_phenotype
PCN_specificity_nonresistant <- (all_data %>% filter(PCN_interpretation == "R" & blaTEM == 1) %>% nrow())/(all_data %>% filter(PCN_interpretation == "R") %>% nrow())

PCN_true_positives <- all_data %>% filter(PCN_interpretation != "R" & blaTEM == 0) %>% nrow()
PCN_false_negatives <- all_data %>% filter(PCN_interpretation != "R" & blaTEM == 1) %>% nrow()
PCN_true_negatives <- all_data %>% filter(PCN_interpretation == "R" & blaTEM == 1) %>% nrow()
PCN_false_positives <- all_data %>% filter(PCN_interpretation == "R" & blaTEM == 0) %>% nrow()

PCN_total_susceptible_phenotype
PCN_total_nonresistant_phenotype
PCN_total_susceptible_genotype

PCN_sensitivity_susceptible
PCN_specificity_susceptible

PCN_sensitivity_nonresistant
PCN_specificity_nonresistant

# TET
TET_total_susceptible_phenotype <- all_data %>% filter(TET_interpretation == "S") %>% nrow()
TET_total_nonresistant_phenotype <- all_data %>% filter(TET_interpretation != "R") %>% nrow()
TET_total_susceptible_genotype <- all_data %>% filter(tetM == 0) %>% nrow()

TET_true_positives <- all_data %>% filter(TET_interpretation == "S" & tetM == 0) %>% nrow()
TET_false_negatives <- all_data %>% filter(TET_interpretation == "S" & tetM == 1) %>% nrow()
TET_true_negatives <- all_data %>% filter(TET_interpretation != "S" & tetM == 1) %>% nrow()
TET_false_positives <- all_data %>% filter(TET_interpretation != "S" & tetM == 0) %>% nrow()

TET_sensitivity_susceptible <- (all_data %>% filter(TET_interpretation == "S" & tetM == 0) %>% nrow())/TET_total_susceptible_phenotype
TET_specificity_susceptible <- (all_data %>% filter(TET_interpretation != "S" & tetM == 1) %>% nrow())/(all_data %>% filter(TET_interpretation != "S") %>% nrow())

TET_sensitivity_nonresistant <- (all_data %>% filter(TET_interpretation != "R" & tetM == 0) %>% nrow())/TET_total_nonresistant_phenotype
TET_specificity_nonresistant <- (all_data %>% filter(TET_interpretation == "R" & tetM == 1) %>% nrow())/(all_data %>% filter(TET_interpretation == "R") %>% nrow())

TET_true_positives <- all_data %>% filter(TET_interpretation != "R" & tetM == 0) %>% nrow()
TET_false_negatives <- all_data %>% filter(TET_interpretation != "R" & tetM == 1) %>% nrow()
TET_true_negatives <- all_data %>% filter(TET_interpretation == "R" & tetM == 1) %>% nrow()
TET_false_positives <- all_data %>% filter(TET_interpretation == "R" & tetM == 0) %>% nrow()

TET_total_susceptible_phenotype
TET_total_nonresistant_phenotype
TET_total_susceptible_genotype

TET_sensitivity_susceptible
TET_specificity_susceptible

TET_sensitivity_nonresistant
TET_specificity_nonresistant


# utility in different populations

chisq_test(all_data, sexual_behavior ~ PCN_genotype)
chisq_test(all_data, sexual_behavior ~ TET_genotype)

msm <- all_data %>% filter(sexual_behavior == "MSM")
chisq_test(msm, race ~ PCN_genotype)
chisq_test(msm, race ~ TET_genotype)

msw <- all_data %>% filter(sexual_behavior == "MSW")
chisq_test(msw, race ~ PCN_genotype)
chisq_test(msw, race ~ TET_genotype)

all_data %>% count(sexual_behavior, race, PCN_genotype) %>%
    pivot_wider(names_from = PCN_genotype, values_from = n, values_fill = 0) %>%
    mutate(PCN_total = S + R, S_percent = S*100/PCN_total, R_percent = R*100/PCN_total) %>%
    view()

all_data %>% count(sexual_behavior, race, TET_genotype) %>%
    pivot_wider(names_from = TET_genotype, values_from = n, values_fill = 0) %>%
    mutate(TET_total = S + R, S_percent = S*100/TET_total, R_percent = R*100/TET_total) %>%
    view()
