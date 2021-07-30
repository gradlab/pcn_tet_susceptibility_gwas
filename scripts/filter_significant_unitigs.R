library(tidyverse)

arguments <- commandArgs(trailingOnly = T)

limits <- read_tsv(arguments[1], col_names = F)
p_threshold <- limits[[2,2]]
unitigs <- read_tsv(arguments[2])
unitigs.filtered <- unitigs %>% filter(`lrt-pvalue` < p_threshold)
write_tsv(unitigs.filtered, arguments[3])
