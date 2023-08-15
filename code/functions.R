library(tidyverse)

#==============================================================================#
# make a custom Manhattan plot
# datPath = the path to an aggregate .tsv mapping file from NemaScan
# indepTestPath = the path to the independent tests for EIGEN threshold from NeamScan: Genotype_Matrix/total_independent_tests.txt
# xlab = ggplot x label default is: "Genomic position (Mb)"
# ylab = ggplot y label default is: expression(-log[10](italic(p)))
# title = Whatever you want to title the plot: default is: NULL
# sigThresh = The significance threshold to use for plotting, defualt is: NULL which shows BF and EIGEN,
#   other options are "BF" and "EIGEN" to show just one or the other.
#==============================================================================#
manPlot <- function(dataPath,
                    indepTestPath,
                    xlab = "Genomic position (Mb)",
                    ylab = expression(-log[10](italic(p))),
                    title = NULL,
                    sigThresh = NULL) {
  
  # pull the nemascan data
  data <- data.table::fread(dataPath) %>%
    dplyr::filter(CHROM != "MtDNA") %>%
    dplyr::distinct(marker, .keep_all = T)
  
  # get the BF threshold
  BF <- data %>% 
    dplyr::group_by(trait) %>% 
    dplyr::filter(log10p != 0) %>% 
    dplyr::distinct(marker, log10p) %>%
    dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
    dplyr::ungroup() %>%
    dplyr::select(BF) %>%
    unique(.) %>%
    dplyr::slice(1) %>% # BF can be slightly different between loco and inbred... but just plot one (5.46 v 5.47...)
    as.numeric()
  
  # load independent tests result file to get EIGEN sig threshold
  total_independent_tests <- read.table(indepTestPath, quote="\"", comment.char="", stringsAsFactors=FALSE)
  EIGEN <- -log10(0.05/total_independent_tests[[1]])
  
  # set sig thresholds
  if(is.null(sigThresh)){
    sig.colors <- c("red","#EE4266", "black")
    names(sig.colors) <- c("BF","EIGEN", "NONSIG")
    # set the sig thresholds for plotting
    sig.df <- tibble::tibble(name = c("BF", "EIGEN"), sig = c(BF, EIGEN)) 
  }else{ 
    if(sigThresh == "BF"){
      sig.colors <- c("red", "black")
      names(sig.colors) <- c("BF", "NONSIG")
      sig.df <- tibble::tibble(name = "BF", sig = BF)
    }
    if(sigThresh == "EIGEN"){
      sig.colors <- c("red", "black")
      names(sig.colors) <- c("EIGEN", "NONSIG")
      sig.df <- tibble::tibble(name = "EIGEN", sig = EIGEN)
    }
    if(!is.null(sigThresh) & !(sigThresh %in%  c("BF", "EIGEN"))){
      message('Please specify a valid sigThresh: NULL, "BF", or "EIGEN"')
      stop()
    }
  }
  
  # add sig colors
  sig.data <- data %>%
    dplyr::mutate(sig = dplyr::case_when(log10p > BF ~ "BF",
                                         log10p > EIGEN & log10p <= BF ~ "EIGEN",
                                         TRUE ~ "NONSIG"))
  # Plot it
  man.plot <- ggplot() + 
    theme_bw() + 
    geom_point(data = sig.data, 
               mapping = aes(x = POS/1000000, 
                             y = log10p,
                             colour = sig,
                             alpha = sig)) +
    scale_alpha_manual(values = c("BF" = 1, "EIGEN" = 1, "user" = 1, "NONSIG" = 0.25)) +
    scale_colour_manual(values = sig.colors) + 
    scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
    geom_hline(data = sig.df, aes(yintercept = sig, linetype = name)) + 
    scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3, "user" = 2)) +
    labs(x = xlab,
         y = ylab) + #expression(-log[10](italic(p)))
    theme(legend.position = "none", 
          panel.grid = element_blank()) + 
    facet_grid(~CHROM, scales = "free") 
  
    # add a title if needed
    if(is.null(title)){
      return(man.plot)
    } else {
      man.plot2 <- man.plot +
        ggtitle(title)
      return(man.plot2)
    }
  }

#==============================================================================#
# dataPath = Nemascan Aggregate mapping .tsv file path. 
# Make a PxG plot like NemaScan does
#==============================================================================#
pxgPlot <- function(dataPath = "data/processed_trait_AGGREGATE_mapping_inbred.tsv") {
  # get the data
  processed_mapping <- read.delim(dataPath, stringsAsFactors=FALSE) %>%
    dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
    dplyr::select(-marker) %>%
    tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)
  
  nested.pxg.dat <- processed_mapping %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::select(CHROM, marker, trait, startPOS, peakPOS, endPOS, AF1, value, strain, allele, peak_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(startPOS = startPOS/1000000,
                  peakPOS = peakPOS/1000000,
                  endPOS = endPOS/1000000) %>%
    # dplyr::left_join(.,sweep.chrom.pivot) %>% # don't have this file yet
    # dplyr::group_by(trait, peak_id) %>%
    # dplyr::recode(allele, "-1" = "REF", "1" = "ALT") %>%
    dplyr::mutate(allele = dplyr::case_when(allele == "-1" ~ "REF",
                                            allele == "1" ~ "ALT",
                                            TRUE ~ "NA"),
                  allele = factor(allele, levels = c("REF", "ALT")))
  
  strains.of.interest <- c("PD1074", "N2", "CB4856", "RC301", "MY16", 
                           "ECA396", "ECA36", "XZ1516", "ECA248", "AB1", 
                           "CB4507", "CB4858", "CB4855", "CB4852", "MY1", 
                           "JU319", "JU345", "JU400", "PB306", "PX174", "PX179")
  
  plot <- nested.pxg.dat %>%
    dplyr::filter(allele != "NA" | !is.na(allele)) %>%
    dplyr::mutate(SOI = strain %in% strains.of.interest,
                  SOI.3 = dplyr::case_when(strain %in% c("N2", "PD1074") ~ "N2",
                                           strain == "CB4856" ~ "CB",
                                           strain %in% strains.of.interest ~ "special",
                                           TRUE ~ "other"),
                  SOI.2 = if_else(SOI == TRUE, true = strain, false = "")) %>%
    droplevels() %>%
    dplyr::arrange(SOI.2) %>%
    ggplot2::ggplot(mapping = aes(x = allele, y = value, text = SOI.2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_violin(aes(fill = allele), alpha = 0.5, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_fill_manual(values = c("REF" = "#A79F92", "ALT" = "mediumpurple4"), guide = FALSE) +
    ggnewscale::new_scale("fill") +
    ggplot2::geom_point(aes(fill = SOI.3, size = SOI), position = ggbeeswarm::position_beeswarm(), shape = 21) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "special" ="red", "other" ="grey50"), guide = FALSE) +
    ggplot2::scale_size_manual(values = c(1.5,2.5), guide = FALSE) +
    ggrepel::geom_text_repel(aes(label = SOI.2),
                             colour = "black", position = ggbeeswarm::position_beeswarm(),
                             max.overlaps = Inf,
                             min.segment.length = 0) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(y = "Trait Value",
                  x = "Genotype") +
    ggplot2::facet_grid(~marker)
  
  # return
  return(plot)
  
}

#==============================================================================#
# dataPath = a path to the NemaScan finemapping output bcsq_genes.tsv file 
# Make a fineMapping plot
#==============================================================================#
finePlotGene <- function(dataPath) {
  # load the data
  genes_in_region <- data.table::fread(dataPath)
 
  # get gene level data
  gene_df <- genes_in_region %>%
  dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
  dplyr::distinct(WBGeneID, GENE_NAME, PEAK_MARKER, CHROM, STRAND, TRANSCRIPTION_START_POS, 
                  TRANSCRIPTION_END_POS, START_POS, END_POS, VARIANT_LOG10p) %>%
  dplyr::group_by(WBGeneID) %>%
  dplyr::mutate(VARIANT_LOG10p = max(VARIANT_LOG10p, na.rm = T)) %>%
  dplyr::distinct()
  
# report gene counts
  #message(glue::glue("There are {gene_df")
# find the peak variant position
peak_variant <- as.numeric(strsplit(unique(gene_df$PEAK_MARKER), split = ":")[[1]][2])

# set the intergenic impacts
variant_df <- genes_in_region %>%
  dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
  dplyr::distinct(CHROM, POS, GENE_NAME, WBGeneID, VARIANT_LOG10p, VARIANT_IMPACT)

variant_df$VARIANT_IMPACT[is.na(variant_df$VARIANT_IMPACT)] <- "Intergenic"

# get the start, end, and chrom for finemap
xs <- unique(gene_df$START_POS)
xe <- unique(gene_df$END_POS)
cq <- unique(gene_df$CHROM)

# get an offset for the y axis to plot correctly based on max logP
max_logp <- unique(max(variant_df$VARIANT_LOG10p, na.rm = T))/150


gene_plot <- ggplot(gene_df) +
  aes(text = paste0(GENE_NAME, "\n", WBGeneID)) +
  geom_vline(aes(xintercept = peak_variant/1e6),
             linetype=3, color = "cyan")+
  geom_segment(aes(x = ifelse(STRAND == "+", TRANSCRIPTION_START_POS/1e6, TRANSCRIPTION_END_POS/1e6),
                   xend = ifelse(STRAND == "+", TRANSCRIPTION_END_POS/1e6, TRANSCRIPTION_START_POS/1e6),
                   y = VARIANT_LOG10p,
                   yend = VARIANT_LOG10p),
               arrow = arrow(length = unit(5, "points")), size = 1) +
  geom_segment(aes(x = POS/1e6,
                   xend = POS/1e6,
                   y = VARIANT_LOG10p+max_logp,
                   yend = VARIANT_LOG10p-max_logp,
                   color = VARIANT_IMPACT), data = variant_df) +
  scale_color_manual(values = c("MODIFIER" = "gray50",
                                "LOW" = "gray30",
                                "MODERATE" = "orange",
                                "HIGH" = "red",
                                "Linker" = "gray80", 
                                "Intergenic" = "gray80"),
                     breaks = c("HIGH", "MODERATE", "LOW", "MODIFIER", "Intergenic"),
                     name = "EFFECT")+
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  # y = expression(-log[10](italic(p))))+
  theme_bw() +
  xlim(c(xs/1e6, xe/1e6)) +
  theme(legend.position = "top",
        panel.grid = element_blank())
  
# return the plot
  return(gene_plot)
}

