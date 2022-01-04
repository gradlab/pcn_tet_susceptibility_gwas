library(tidyverse)
library(cowplot)

qc <- read_tsv("data/prediction/qc.txt")

p_coverage <- qc %>%
  ggplot(aes(x = reference_coverage)) +
  scale_x_continuous(trans='log10') +
  geom_histogram(bins = 100) +
  xlab("NCCP11945 coverage") +
  theme_minimal() 

p_length <- qc %>%
  ggplot(aes(x = `Total length`)) +
  geom_histogram(bins = 100) +
  xlab("Assembly length") +
  theme_minimal() 

p_n50 <- qc %>%
  ggplot(aes(x = N50)) +
  scale_x_continuous(trans='log10') +
  geom_histogram(bins = 100) +
  xlab("Assembly N50") +
  theme_minimal() 

plot_grid(p_coverage, p_length, p_n50, ncol = 1, labels = "AUTO")
