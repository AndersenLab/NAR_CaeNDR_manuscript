library(tidyverse)

# set the working dir
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# load the functions
source("code/functions.R")

# load Evans et al. 2021 abamectin traits
at <- data.table::fread("/projects/b1059/projects/Tim/NAR_CaeNDR_manuscript/data/Evans2021_S2_abamectin_traits.csv", drop = 1) %>%
  dplyr::select(-condition) %>%
  tidyr::pivot_wider(names_from = trait, values_from = phenotype) %>%
  dplyr::filter(!is.na(mean.TOF))

# export the mean.TOF trait to map with CaeNDR
mean.TOF <- at %>%
  dplyr::select(strain, mean.TOF) %>%
  rio::export(., file = "/projects/b1059/projects/Tim/NAR_CaeNDR_manuscript/data/mean.TOF.tsv")

# recreate the Manhattan plot from the NemaScan output
Figure2A <- manPlot(dataPath = "data/processed_mean_TOF_AGGREGATE_mapping_loco.tsv",
                    indepTestPath = "data/processed_mean_TOF_total_independent_tests.txt",
                    xlab = "Genomic position (Mb)",
                    ylab = expression(-log[10](italic(p))),
                    title = NULL,
                    sigThresh = NULL)

# restict the aggregate mapping to the chr V 16mb QTL for plotting
data.table::fread("data/processed_mean_TOF_AGGREGATE_mapping_loco.tsv") %>%
  dplyr::filter(peak_id == 3) %>%
  rio::export(., file = "data/processed_mean_TOF_AGGREGATE_mapping_loco_V_16175482.tsv")

# make a pxg plot to fit
Figure2B <- pxgPlot(dataPath = "data/processed_mean_TOF_AGGREGATE_mapping_loco_V_16175482.tsv") +
  labs(y = "Abamectin response")

# make a fine mapping plot
# read in data
fine.df <- data.table::fread("data/processed_mean_TOF_V_15840330-18896645_bcsq_genes_loco.tsv")

# add labels for glc-1
fine.df.labs <- fine.df %>%
  dplyr::mutate(label = case_when(GENE_NAME == "glc-1" ~ "glc-1",
                                  TRUE ~ NA_character_)) %>%
  dplyr::group_by(GENE_NAME) %>%
  dplyr::arrange(desc(VARIANT_LOG10p)) %>%
  dplyr::mutate(n.variants = n()/2,
                variant.id = rep(1:(n()/2), each = 2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label.2 = case_when(GENE_NAME == "glc-1" & variant.id == 1 ~ "italic('glc-1')",
                                    TRUE ~ NA_character_)) %>%
  dplyr::distinct(label.2, .keep_all = T) %>%
  dplyr::filter(!is.na(label.2)) %>%
  dplyr::mutate(POS2 = POS+2500,
                VARIANT_LOG10p2 = VARIANT_LOG10p+0.025)

# plot without labels
Figure2C <- finePlotGene(dataPath = "data/processed_mean_TOF_V_15840330-18896645_bcsq_genes_loco.tsv")

# tweek it with labels for glc-1 and zoom in to the action
Figure2C_tweek <- Figure2C +
  xlim(16, 16.5) +
  ylim(4, 12.5) + #max(Figure2C$data$VARIANT_LOG10p)+1
  theme(legend.position = c(0.15, 0.8)) +
  ggrepel::geom_text_repel(data = fine.df.labs, aes(x = POS2/1e6, y = VARIANT_LOG10p2, label = label.2),
                           min.segment.length = 0,
                           nudge_x = 0.03,
                           nudge_y = 0.8,
                           arrow = arrow(length = unit(0.01, "npc")),
                           parse = T)

# put them together.
figure2bc <- cowplot::plot_grid(Figure2B, Figure2C_tweek, rel_widths = c(1,2), align = "vh", labels = c("B", "C"), vjust = 1.25)

fig2 <- cowplot::plot_grid(Figure2A, figure2bc, labels = c("A", ""), ncol = 1, vjust = 1.25)# align = "h", axis = "bl")

cowplot::ggsave2(fig2, filename = "plots/Fig2.png", width = 8.5, height = 8.5)
cowplot::ggsave2(fig2, filename = "plots/Fig2.pdf", width = 8.5, height = 8.5)
