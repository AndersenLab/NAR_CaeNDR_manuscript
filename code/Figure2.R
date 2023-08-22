library(tidyverse)

# set the working dir
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load the functions
source("code/functions.R")

# Make a manhatten plot from NemaScan Output
Figure2A <- manPlot(dataPath = "data/processed_trait_AGGREGATE_mapping_inbred.tsv",
                indepTestPath = "data/processed_trait_total_independent_tests.txt",
                xlab = "Genomic position (Mb)",
                ylab = expression(-log[10](italic(p))),
                title = NULL,
                sigThresh = NULL)

# make a pxg plot to fit
Figure2B <- pxgPlot(dataPath = "data/processed_trait_AGGREGATE_mapping_inbred.tsv") +
  labs(y = "Trait value")

# make a finemapping plot
fine.df <- data.table::fread("data/processed_trait_V_1757246-4295865_bcsq_genes_inbred.tsv")
Figure2C <- finePlotGene(dataPath = "data/processed_trait_V_1757246-4295865_bcsq_genes_inbred.tsv")
Figure2C_tweek <- Figure2C +
  xlim(2.5, 3) +
  ylim(9, max(Figure2C$data$VARIANT_LOG10p)+1) +
  theme(legend.position = c(0.15, 0.8))

# put them together.
figure2bc <- cowplot::plot_grid(Figure2B, Figure2C_tweek, rel_widths = c(1,2), align = "vh", labels = c("B", "C"), vjust = 1.25)
figure2bc_v2 <- cowplot::plot_grid(Figure2B, Figure2C_tweek, align = "vh", labels = c("B", "C"), vjust = 1.25)

fig2 <- cowplot::plot_grid(Figure2A, figure2bc, labels = c("A", ""), ncol = 1, vjust = 1.25)# align = "h", axis = "bl")
fig2_v2 <- cowplot::plot_grid(Figure2A, figure2bc_v2, labels = c("A", ""), ncol = 1, vjust = 1.25)# align = "h", axis = "bl")
#fig2

cowplot::ggsave2(fig2, filename = "plots/Figure2_v2.png", width = 8.5, height = 8.5)
cowplot::ggsave2(fig2_v2, filename = "plots/Figure2_v3.png", width = 8.5, height = 8.5)
