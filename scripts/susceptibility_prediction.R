library(tidyverse)
library(cowplot)

# global data
all_data <- read_tsv("../data/prediction/pcn_tet_genotype_phenotype.tsv",
                            col_types = cols(penicillin = col_character()))

# validation data
validation_metadata <- read_tsv("../data/validation/reimche_metadata.tsv")
validation_resistance <- read_tsv("../data/validation/2021-06-24_reimche_gisp_alleles.tsv")

# combine data

validation_all_data <- validation_resistance %>%
    left_join(validation_metadata) %>%
    mutate(reference = "Reimche2021") %>%
    select(wgs_id, reference, tetracycline, penicillin, RpsJ_57, penA_01, tetM, blaTEM) %>%
    mutate(PCN_genotype = if_else(penA_01 == 1 & blaTEM == 0, "S", "R")) %>%
    mutate(TET_genotype = if_else(RpsJ_57 == "V" & tetM == 0, "S", "R")) %>%
    mutate(TET = str_replace_all(tetracycline, "[><= ]", "")) %>%
    mutate(TET = as.numeric(TET)) %>%
    mutate(PCN = str_replace_all(penicillin, "[><= ]", "")) %>%
    mutate(PCN = as.numeric(PCN)) %>%
    mutate(TET_interpretation = case_when(TET <= 0.25 ~ "S", TET > 0.25 & TET < 2 ~ "I", TET >= 2 ~ "R")) %>%
    mutate(PCN_interpretation = case_when(PCN <= 0.064 ~ "S", PCN > 0.064 & PCN < 2 ~ "I", PCN >= 2 ~ "R")) %>%
    mutate(dataset = "Validation")

# add genotype interpretation
all_data <- all_data %>%
    mutate(PCN_genotype = if_else(penA_01 == 1 & blaTEM == 0, "S", "R")) %>%
    mutate(TET_genotype = if_else(RpsJ_57 == "V" & tetM == 0, "S", "R"))

# filter to valid TET MICs and categorize based on breakpoints
metadata_tet <- all_data %>% 
    select(-penicillin) %>%
    filter(tetracycline != "<=4", tetracycline != "<=4.0", tetracycline != "<=8.0") %>%
    mutate(TET = str_replace_all(tetracycline, "[><= ]", "")) %>%
    mutate(TET = as.numeric(TET)) %>%
    mutate(TET_interpretation = case_when(TET <= 0.25 ~ "S", TET > 0.25 & TET < 2 ~ "I", TET >= 2 ~ "R")) %>%
    drop_na()

# filter to valid PCN MICs and categorized based on breakpoints (0.06 breakpoint interpretated to include 0.063/0.064)
metadata_pcn <- all_data %>%
    select(-tetracycline) %>%
    mutate(PCN = str_replace_all(penicillin, "[><= ]", "")) %>%
    mutate(PCN = as.numeric(PCN)) %>%
    mutate(PCN_interpretation = case_when(PCN <= 0.064 ~ "S", PCN > 0.064 & PCN < 2 ~ "I", PCN >= 2 ~ "R")) %>%
    drop_na()

# get genotype counts for all data by dataset
all_data %>% count(reference, PCN_genotype) %>%
    pivot_wider(names_from = PCN_genotype, values_from = n) %>%
    mutate(PCN_total = S + R, S_percent = S*100/PCN_total, R_percent = R*100/PCN_total) %>%
    print(n = 30)

all_data %>% count(reference, TET_genotype) %>%
    pivot_wider(names_from = TET_genotype, values_from = n) %>%
    mutate(TET_total = S + R, S_percent = S*100/TET_total, R_percent = R*100/TET_total) %>%
    print(n = 30)

# calculate sensitivity and specificity

# PCN
PCN_total_susceptible_phenotype <- metadata_pcn %>% filter(PCN_interpretation == "S") %>% nrow()
PCN_total_nonresistant_phenotype <- metadata_pcn %>% filter(PCN_interpretation != "R") %>% nrow()
PCN_total_susceptible_genotype <- all_data %>% filter(PCN_genotype == "S") %>% nrow()

PCN_sensitivity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation == "S" & PCN_genotype == "S") %>% nrow())/PCN_total_susceptible_phenotype
PCN_specificity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation != "S" & PCN_genotype == "R") %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation != "S") %>% nrow())

PCN_sensitivity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation != "R" & PCN_genotype == "S") %>% nrow())/PCN_total_nonresistant_phenotype
PCN_specificity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation == "R" & PCN_genotype == "R") %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation == "R") %>% nrow())

# confidence intervals susceptible breakpoint

PCN_true_positives <- metadata_pcn %>% filter(PCN_interpretation == "S" & PCN_genotype == "S") %>% nrow()
PCN_false_negatives <- metadata_pcn %>% filter(PCN_interpretation == "S" & PCN_genotype != "S") %>% nrow()
PCN_true_negatives <- metadata_pcn %>% filter(PCN_interpretation != "S" & PCN_genotype == "R") %>% nrow()
PCN_false_positives <- metadata_pcn %>% filter(PCN_interpretation != "S" & PCN_genotype == "S") %>% nrow()

PCN_sensitivity_susceptible + 1.96*sqrt((PCN_sensitivity_susceptible*(1-PCN_sensitivity_susceptible))/PCN_true_positives)
PCN_sensitivity_susceptible - 1.96*sqrt((PCN_sensitivity_susceptible*(1-PCN_sensitivity_susceptible))/PCN_true_positives)

PCN_specificity_susceptible + 1.96*sqrt((PCN_specificity_susceptible*(1-PCN_specificity_susceptible))/PCN_true_positives)
PCN_specificity_susceptible - 1.96*sqrt((PCN_specificity_susceptible*(1-PCN_specificity_susceptible))/PCN_true_positives)

# confidence intervals non-resistant breakpoint
PCN_true_positives <- metadata_pcn %>% filter(PCN_interpretation != "R" & PCN_genotype == "S") %>% nrow()
PCN_false_negatives <- metadata_pcn %>% filter(PCN_interpretation != "R" & PCN_genotype != "S") %>% nrow()
PCN_true_negatives <- metadata_pcn %>% filter(PCN_interpretation == "R" & PCN_genotype == "R") %>% nrow()
PCN_false_positives <- metadata_pcn %>% filter(PCN_interpretation == "R" & PCN_genotype == "S") %>% nrow()

PCN_sensitivity_nonresistant + 1.96*sqrt((PCN_sensitivity_nonresistant*(1-PCN_sensitivity_nonresistant))/PCN_true_positives)
PCN_sensitivity_nonresistant - 1.96*sqrt((PCN_sensitivity_nonresistant*(1-PCN_sensitivity_nonresistant))/PCN_true_positives)

PCN_specificity_nonresistant + 1.96*sqrt((PCN_specificity_nonresistant*(1-PCN_specificity_nonresistant))/PCN_true_positives)
PCN_specificity_nonresistant - 1.96*sqrt((PCN_specificity_nonresistant*(1-PCN_specificity_nonresistant))/PCN_true_positives)

PCN_total_susceptible_phenotype
PCN_total_nonresistant_phenotype
PCN_total_susceptible_genotype

PCN_sensitivity_susceptible
PCN_specificity_susceptible

PCN_sensitivity_nonresistant
PCN_specificity_nonresistant

# TET
TET_total_susceptible_phenotype <- metadata_tet %>% filter(TET_interpretation == "S") %>% nrow()
TET_total_nonresistant_phenotype <- metadata_tet %>% filter(TET_interpretation != "R") %>% nrow()
TET_total_susceptible_genotype <- all_data %>% filter(TET_genotype == "S") %>% nrow()

TET_sensitivity_susceptible <- (metadata_tet %>% filter(TET_interpretation == "S" & TET_genotype == "S") %>% nrow())/TET_total_susceptible_phenotype
TET_specificity_susceptible <- (metadata_tet %>% filter(TET_interpretation != "S" & TET_genotype == "R") %>% nrow())/(metadata_tet %>% filter(TET_interpretation != "S") %>% nrow())

TET_sensitivity_nonresistant <- (metadata_tet %>% filter(TET_interpretation != "R" & TET_genotype == "S") %>% nrow())/TET_total_nonresistant_phenotype
TET_specificity_nonresistant <- (metadata_tet %>% filter(TET_interpretation == "R" & TET_genotype == "R") %>% nrow())/(metadata_tet %>% filter(TET_interpretation == "R") %>% nrow())

# confidence intervals susceptible
TET_true_positives <- metadata_tet %>% filter(TET_interpretation == "S" & TET_genotype == "S") %>% nrow()
TET_false_negatives <- metadata_tet %>% filter(TET_interpretation == "S" & TET_genotype == "R") %>% nrow()
TET_true_negatives <- metadata_tet %>% filter(TET_interpretation != "S" & TET_genotype == "R") %>% nrow()
TET_false_positives <- metadata_tet %>% filter(TET_interpretation != "S" & TET_genotype == "S") %>% nrow()

TET_sensitivity_susceptible + 1.96*sqrt((TET_sensitivity_susceptible*(1-TET_sensitivity_susceptible))/TET_true_positives)
TET_sensitivity_susceptible - 1.96*sqrt((TET_sensitivity_susceptible*(1-TET_sensitivity_susceptible))/TET_true_positives)

TET_specificity_susceptible + 1.96*sqrt((TET_specificity_susceptible*(1-TET_specificity_susceptible))/TET_true_positives)
TET_specificity_susceptible - 1.96*sqrt((TET_specificity_susceptible*(1-TET_specificity_susceptible))/TET_true_positives)

# confidence intervals non-resistant

TET_true_positives <- metadata_tet %>% filter(TET_interpretation != "R" & TET_genotype == "S") %>% nrow()
TET_false_negatives <- metadata_tet %>% filter(TET_interpretation != "R" & TET_genotype == "R") %>% nrow()
TET_true_negatives <- metadata_tet %>% filter(TET_interpretation == "R" & TET_genotype == "R") %>% nrow()
TET_false_positives <- metadata_tet %>% filter(TET_interpretation == "R" & TET_genotype == "S") %>% nrow()

TET_sensitivity_nonresistant + 1.96*sqrt((TET_sensitivity_nonresistant*(1-TET_sensitivity_nonresistant))/TET_true_positives)
TET_sensitivity_nonresistant - 1.96*sqrt((TET_sensitivity_nonresistant*(1-TET_sensitivity_nonresistant))/TET_true_positives)

TET_specificity_nonresistant + 1.96*sqrt((TET_specificity_nonresistant*(1-TET_specificity_nonresistant))/TET_true_positives)
TET_specificity_nonresistant - 1.96*sqrt((TET_specificity_nonresistant*(1-TET_specificity_nonresistant))/TET_true_positives)

TET_total_susceptible_phenotype
TET_total_nonresistant_phenotype
TET_total_susceptible_genotype

TET_sensitivity_susceptible
TET_specificity_susceptible

TET_sensitivity_nonresistant
TET_specificity_nonresistant

# plasmid only prediction

# PCN
PCN_total_susceptible_phenotype <- metadata_pcn %>% filter(PCN_interpretation == "S") %>% nrow()
PCN_total_nonresistant_phenotype <- metadata_pcn %>% filter(PCN_interpretation != "R") %>% nrow()
PCN_total_susceptible_genotype <- metadata_pcn %>% filter(blaTEM == 0) %>% nrow()

PCN_true_positives <- metadata_pcn %>% filter(PCN_interpretation == "S" & blaTEM == 0) %>% nrow()
PCN_false_negatives <- metadata_pcn %>% filter(PCN_interpretation == "S" & blaTEM == 1) %>% nrow()
PCN_true_negatives <- metadata_pcn %>% filter(PCN_interpretation != "S" & blaTEM == 1) %>% nrow()
PCN_false_positives <- metadata_pcn %>% filter(PCN_interpretation != "S" & blaTEM == 0) %>% nrow()

PCN_sensitivity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation == "S" & blaTEM == 0) %>% nrow())/PCN_total_susceptible_phenotype
PCN_specificity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation != "S" & blaTEM == 1) %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation != "S") %>% nrow())

PCN_sensitivity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation != "R" & blaTEM == 0) %>% nrow())/PCN_total_nonresistant_phenotype
PCN_specificity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation == "R" & blaTEM == 1) %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation == "R") %>% nrow())

PCN_true_positives <- metadata_pcn %>% filter(PCN_interpretation != "R" & blaTEM == 0) %>% nrow()
PCN_false_negatives <- metadata_pcn %>% filter(PCN_interpretation != "R" & blaTEM == 1) %>% nrow()
PCN_true_negatives <- metadata_pcn %>% filter(PCN_interpretation == "R" & blaTEM == 1) %>% nrow()
PCN_false_positives <- metadata_pcn %>% filter(PCN_interpretation == "R" & blaTEM == 0) %>% nrow()

PCN_total_susceptible_phenotype
PCN_total_nonresistant_phenotype
PCN_total_susceptible_genotype

PCN_sensitivity_susceptible
PCN_specificity_susceptible

PCN_sensitivity_nonresistant
PCN_specificity_nonresistant

# TET
TET_total_susceptible_phenotype <- metadata_tet %>% filter(TET_interpretation == "S") %>% nrow()
TET_total_nonresistant_phenotype <- metadata_tet %>% filter(TET_interpretation != "R") %>% nrow()
TET_total_susceptible_genotype <- metadata_tet %>% filter(tetM == 0) %>% nrow()

TET_true_positives <- metadata_tet %>% filter(TET_interpretation == "S" & tetM == 0) %>% nrow()
TET_false_negatives <- metadata_tet %>% filter(TET_interpretation == "S" & tetM == 1) %>% nrow()
TET_true_negatives <- metadata_tet %>% filter(TET_interpretation != "S" & tetM == 1) %>% nrow()
TET_false_positives <- metadata_tet %>% filter(TET_interpretation != "S" & tetM == 0) %>% nrow()

TET_sensitivity_susceptible <- (metadata_tet %>% filter(TET_interpretation == "S" & tetM == 0) %>% nrow())/TET_total_susceptible_phenotype
TET_specificity_susceptible <- (metadata_tet %>% filter(TET_interpretation != "S" & tetM == 1) %>% nrow())/(metadata_tet %>% filter(TET_interpretation != "S") %>% nrow())

TET_sensitivity_nonresistant <- (metadata_tet %>% filter(TET_interpretation != "R" & tetM == 0) %>% nrow())/TET_total_nonresistant_phenotype
TET_specificity_nonresistant <- (metadata_tet %>% filter(TET_interpretation == "R" & tetM == 1) %>% nrow())/(metadata_tet %>% filter(TET_interpretation == "R") %>% nrow())

TET_true_positives <- metadata_tet %>% filter(TET_interpretation != "R" & tetM == 0) %>% nrow()
TET_false_negatives <- metadata_tet %>% filter(TET_interpretation != "R" & tetM == 1) %>% nrow()
TET_true_negatives <- metadata_tet %>% filter(TET_interpretation == "R" & tetM == 1) %>% nrow()
TET_false_positives <- metadata_tet %>% filter(TET_interpretation == "R" & tetM == 0) %>% nrow()

TET_total_susceptible_phenotype
TET_total_nonresistant_phenotype
TET_total_susceptible_genotype

TET_sensitivity_susceptible
TET_specificity_susceptible

TET_sensitivity_nonresistant
TET_specificity_nonresistant

# calculate sensitivity and specificity of only plasmid determinants

# PCN
PCN_total_susceptible_phenotype <- metadata_pcn %>% filter(PCN_interpretation == "S") %>% nrow()
PCN_total_nonresistant_phenotype <- metadata_pcn %>% filter(PCN_interpretation != "R") %>% nrow()
PCN_total_susceptible_genotype <- all_data %>% filter(blaTEM == 0) %>% nrow()

PCN_sensitivity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation == "S" & blaTEM == 0) %>% nrow())/PCN_total_susceptible_phenotype
PCN_specificity_susceptible <- (metadata_pcn %>% filter(PCN_interpretation != "S" & blaTEM == 1) %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation != "S") %>% nrow())

PCN_sensitivity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation != "R" & blaTEM == 0) %>% nrow())/PCN_total_nonresistant_phenotype
PCN_specificity_nonresistant <- (metadata_pcn %>% filter(PCN_interpretation == "R" & blaTEM == 1) %>% nrow())/(metadata_pcn %>% filter(PCN_interpretation == "R") %>% nrow())

PCN_total_susceptible_phenotype
PCN_total_nonresistant_phenotype
PCN_total_susceptible_genotype

PCN_sensitivity_susceptible
PCN_specificity_susceptible

PCN_sensitivity_nonresistant
PCN_specificity_nonresistant

# TET

TET_total_susceptible_phenotype <- metadata_tet %>% filter(TET_interpretation == "S") %>% nrow()
TET_total_nonresistant_phenotype <- metadata_tet %>% filter(TET_interpretation != "R") %>% nrow()
TET_total_susceptible_genotype <- all_data %>% filter(tetM == 0) %>% nrow()

TET_sensitivity_susceptible <- (metadata_tet %>% filter(TET_interpretation == "S" & tetM == 0) %>% nrow())/TET_total_susceptible_phenotype
TET_specificity_susceptible <- (metadata_tet %>% filter(TET_interpretation != "S" & tetM == 1) %>% nrow())/(metadata_tet %>% filter(TET_interpretation != "S") %>% nrow())

TET_sensitivity_nonresistant <- (metadata_tet %>% filter(TET_interpretation != "R" & tetM == 0) %>% nrow())/TET_total_nonresistant_phenotype
TET_specificity_nonresistant <- (metadata_tet %>% filter(TET_interpretation == "R" & tetM == 1) %>% nrow())/(metadata_tet %>% filter(TET_interpretation == "R") %>% nrow())

TET_total_susceptible_phenotype
TET_total_nonresistant_phenotype
TET_total_susceptible_genotype

TET_sensitivity_susceptible
TET_specificity_susceptible

TET_sensitivity_nonresistant
TET_specificity_nonresistant

# combine global and validation data

metadata_pcn <- metadata_pcn %>%
    mutate(dataset = "Global") %>%
    select(-penicillin) %>%
    rows_insert(validation_all_data %>% select(-tetracycline, -penicillin, -TET, -TET_interpretation))

metadata_tet <- metadata_tet %>%
    mutate(dataset = "Global") %>%
    select(-tetracycline) %>%
    rows_insert(validation_all_data %>% select(-tetracycline, -penicillin, -PCN, -PCN_interpretation)) %>%
    filter(!is.na(TET_genotype))

# plot

pcn_plot <- metadata_pcn %>% ggplot(aes(x = fct_relevel(PCN_genotype, "S"), y = PCN, fill = PCN_genotype)) +
    geom_violin(bw=0.8, draw_quantiles=c(0.5), trim=FALSE) +
    theme_bw() + 
    xlab("Genotype") +
    ylab(expression(paste("PCN MIC (", mu, "g/mL)", ))) +
    scale_fill_brewer() + 
    #theme(legend.position = "none") +
    #scale_x_discrete(labels = c("penA_01 present \nand blaTEM absent", "penA_01 absent \nor blaTEM present")) +
    scale_y_continuous(trans = 'log2', 
                       limits = c(0.002, 256),
                       breaks = c(2^-9, 2^-7, 2^-5, 2^-3, 2^-1, 2^1, 2^3, 2^5, 2^7), 
                       labels = c(0.002, 0.008, 0.03, 0.125, 0.5, 2, 8, 32, 128)) +
    geom_hline(yintercept = 0.06, linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    facet_wrap(~dataset)


tet_plot <- metadata_tet %>% ggplot(aes(x = fct_relevel(TET_genotype, "S"), y = TET, fill = TET_genotype)) + 
    geom_violin(bw=0.8, draw_quantiles=c(0.5), trim=FALSE) +
    theme_bw() + 
    xlab("Genotype") +
    ylab(expression(paste("TET MIC (", mu, "g/mL)", ))) +
    scale_fill_brewer(palette = "Purples") + 
    #theme(legend.position = "none") + 
    #scale_x_discrete(labels = c("rpsJ WT\nand tetM absent", "rpsJ V57M\n or tetM present")) +
    scale_y_continuous(trans = 'log2', 
                       limits = c(0.002, 256),
                       breaks = c(2^-9, 2^-7, 2^-5, 2^-3, 2^-1, 2^1, 2^3, 2^5, 2^7), 
                       labels = c(0.002, 0.008, 0.03, 0.125, 0.5, 2, 8, 32, 128)) +
    geom_hline(yintercept = 0.25, linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    facet_wrap(~dataset)

plot_grid(pcn_plot, tet_plot, labels = c("A","B"), nrow = 2) 

