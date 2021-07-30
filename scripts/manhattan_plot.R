library(tidyverse)
library(cowplot)

limits_tet <- read_tsv("data/gwas/pcn/significance_limits.txt", col_names = F)
p_threshold_tet <- limits_tet[[2,2]]

limits_pcn <- read_tsv("data/gwas/pcn/significance_limits.txt", col_names = F)
p_threshold_pcn <- limits_pcn[[2,2]]


unitigs <- read_tsv("data/gwas/tet/tet_unitig_annotated_WHO_N.txt", 
                    col_names = c("unitig", "af", "filter_p", "lrt_p", "beta", "beta_std_err", "h2", "notes", "annotation"))

unitigs <- unitigs %>% mutate(annotation = ifelse(is.na(annotation), notes, annotation))

unitigs <- unitigs %>% separate(annotation, c("annotation1", "annotation2", "annotation3", "annotation4", "annotation5", "annotation6"), sep = ',') %>%
  pivot_longer(annotation1:annotation6, names_to = "annotation_number", values_to = "annotation", values_drop_na = T) %>%
  separate(annotation, c("reference", "location"), sep = ":") %>%
  separate(location, c("coordinates", "gene"), sep = ";") %>%
  separate(coordinates, c("start", "stop"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>% 
  mutate(replicon = ifelse(reference == "WHO_N", "chromosome", "plasmid"))

p_chromosome <- ggplot(unitigs %>% mutate(significant = if_else(lrt_p < p_threshold_tet & beta < 0, T, F)) %>% filter(replicon == "chromosome"), 
            aes(x = start, y = -log10(lrt_p), color = significant)) +
  geom_point() + 
  scale_color_manual(values = c("#D3D3D3", "#756bb1")) +
  theme_minimal() +
  xlab("Genomic position") +
  ylab("Significance (-log10(p))") +
  geom_hline(yintercept = -log10(p_threshold_tet), linetype = "dashed") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(reference), scales = "free_x", space = "free_x")

p_plasmid <- ggplot(unitigs %>% mutate(significant = if_else(lrt_p < p_threshold_tet & beta < 0, T, F)) %>% filter(replicon == "plasmid"), 
                    aes(x = start, y = -log10(lrt_p), color = significant)) +
  geom_point() + 
  scale_color_manual(values = c("#D3D3D3", "#756bb1")) +
  theme_minimal() +
  xlab("Genomic position") +
  ylab("Significance (-log10(p))") +
  geom_hline(yintercept = -log10(p_threshold_tet), linetype = "dashed") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(reference), scales = "free_x")

tet <- plot_grid(p_chromosome, p_plasmid, ncol = 1)

unitigs_pcn <- read_tsv("data/gwas/pcn/pcn_unitig_annotated_WHO_N.txt", 
                    col_names = c("unitig", "af", "filter_p", "lrt_p", "beta", "beta_std_err", "h2", "notes", "annotation"))

unitigs_pcn <- unitigs_pcn %>% mutate(annotation = ifelse(is.na(annotation), notes, annotation))

unitigs_pcn <- unitigs_pcn %>% separate(annotation, c("annotation1", "annotation2", "annotation3", "annotation4", "annotation5", "annotation6"), sep = ',') %>%
  pivot_longer(annotation1:annotation6, names_to = "annotation_number", values_to = "annotation", values_drop_na = T) %>%
  separate(annotation, c("reference", "location"), sep = ":") %>%
  separate(location, c("coordinates", "gene"), sep = ";") %>%
  separate(coordinates, c("start", "stop"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>% 
  mutate(replicon = ifelse(reference == "WHO_N", "chromosome", "plasmid"))

p_chromosome_pcn <- ggplot(unitigs_pcn %>% mutate(significant = if_else(lrt_p < p_threshold_pcn & beta < 0, T, F)) %>% filter(replicon == "chromosome"), 
                       aes(x = start, y = -log10(lrt_p), color = significant)) +
  geom_point() + 
  scale_color_manual(values = c("#D3D3D3", "#3182bd")) +
  theme_minimal() +
  xlab("Genomic position") +
  ylab("Significance (-log10(p))") +
  geom_hline(yintercept = -log10(p_threshold_tet), linetype = "dashed") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(reference), scales = "free_x", space = "free_x")

p_plasmid_pcn <- ggplot(unitigs_pcn %>% mutate(significant = if_else(lrt_p < p_threshold_pcn & beta < 0, T, F)) %>% filter(replicon == "plasmid"), 
                    aes(x = start, y = -log10(lrt_p), color = significant)) +
  geom_point() + 
  scale_color_manual(values = c("#D3D3D3", "#3182bd")) +
  theme_minimal() +
  xlab("Genomic position") +
  ylab("Significance (-log10(p))") +
  geom_hline(yintercept = -log10(p_threshold_tet), linetype = "dashed") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(reference), scales = "free_x")

pcn <- plot_grid(p_chromosome_pcn, p_plasmid_pcn, ncol = 1)

manhattan_plot <- plot_grid(pcn, tet, nrow = 2, labels = c("A", "B"))
ggsave("data/figures/manhattan_plot.pdf", manhattan_plot, width = 6, height = 9, units = "in")
