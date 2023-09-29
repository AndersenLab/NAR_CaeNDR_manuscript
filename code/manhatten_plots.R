library(tidyverse)

# set dir relative to repo
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# read data
dat <- data.table::fread("data/processed_length_TRAIT_AGGREGATE_mapping_loco.tsv") %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::distinct(marker, .keep_all = T)

# Plot it
man_plot1 <- ggplot(dat,aes(x = POS/1000000, 
               y = log10p,
               fill = as.factor(aboveBF))) + 
  theme_classic() + 
  geom_rect(aes(xmin = startPOS/1000000, 
                xmax = endPOS/1000000, 
                ymin = 0, 
                ymax = Inf, 
                fill = "blue"),
            color = "blue",
            fill = "cyan",
            linetype = 2, 
            alpha=.3,
            linewidth = 0)+ 
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = c("black","red")) + 
  #scale_x_continuous(limits = c(1986,2014), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  #scale_fill_manual(values = c("black","yellow")) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,12.5)) +
  geom_hline(yintercept = unique(dat$BF), linetype = 2) +
  #labs(x = "", y = "Marker significance") +
       #y = expression(-log[10](italic(p)))) +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        strip.background = element_rect(color="black", fill='grey', linewidth=1, linetype="solid"),
        strip.text = element_text(color="black", face = "bold", size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        #axis.title.y = element_text(color="black", face = "bold", size = 18),
        axis.ticks.y = element_blank()) + 
  facet_grid(~CHROM, scales = "free_x", space = "free") 
man_plot1

# save it
ggsave(man_plot1, filename = "plots/manhatten_plot_2.png", width = 6.6, height = 4.75)
