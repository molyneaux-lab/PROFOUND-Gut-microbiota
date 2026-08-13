setwd("~/Documents/GitHub/PROFOUND-Baseline/")
renv::load()
renv::restore()
library(dplyr)
library(DT)
library(mia)
library(tidyverse)
library(ggsci)

nejm_colors <- pal_nejm()(8)
scales::show_col(nejm_colors)
Palette_2 = c("#BC3C29FF", "#0072B5FF")

Palette_10 = c("#003f5a", # deep navy (blue-teal)
               "#0072B5", # blue
               "#6F99AD", # slate blue
               "#6ea6a4", # muted turquoise
               "#20854E", # green-teal
               "#8ab184", # sage
               "#C4644A", # muted pink
               "#BC3C29", # red-brown
               "#de6600", # orange
               "#E18727", # gold
               "#FFDC91", # light gold
               "#fec682" # peach
)

Palette_diverging <- colorRampPalette(Palette_10)
Palette_diverging <- Palette_diverging(40)

set.seed(2026)

## PROFOUND DATA 
PROFOUND <- read_tsv("data/Genus-Absolute-Abundance-Counts-Transposed.txt")
PROFOUND <- column_to_rownames(PROFOUND, "...1")
lib_sizes <- rowSums(PROFOUND)

df_lib_sizes <- data.frame(
  PatientID = names(lib_sizes),
  LibrarySize = lib_sizes)

PROFOUND  =  PROFOUND/rowSums(PROFOUND)*100
PROFOUND = as.data.frame(t(PROFOUND))
PROFOUND = rownames_to_column(PROFOUND, "Taxa")
PROFOUND = PROFOUND %>%
  mutate(
    Taxa = str_replace_all(Taxa, "\\[", ""),
    Taxa = str_replace_all(Taxa, "\\]", "")) %>%
  dplyr::rename("PFND164" = "PFND164_amplified") %>% # remove sample that wasn't amplified to keep same
  dplyr::select(-MTN008) # failed sample
PROFOUND = column_to_rownames(PROFOUND, "Taxa")

pfnd_metadata = read_csv("data/16SSequenced_complete-metadata.csv")
pfnd_metadata = pfnd_metadata %>%
  dplyr::select(PatientID, Sex, Age_at_recruitment, Diagnosis, Smoking_history, Reflux_treatment) %>%
  filter(PatientID %in% colnames(PROFOUND)) #%>%

PROFOUND = PROFOUND %>%
  dplyr::select(all_of(pfnd_metadata$PatientID))

# Merged df
PROFOUND = rownames_to_column(PROFOUND, "Taxa")

PROFOUND <- column_to_rownames(PROFOUND, "Taxa")
abund_table <- as.data.frame(t(PROFOUND)) 

#### PCOA ####
library(vegan)
library(glue)
library(aplot)
library(ggplotify)

perm_diagnosis = adonis2(abund_table ~ Diagnosis, permutations = 9999, method = "robust.aitchison", data = pfnd_metadata)
permanova_p = perm_diagnosis$`Pr(>F)`

adonis2(abund_table ~ Reflux_treatment, permutations = 999,
        method = "robust.aitchison", data = pfnd_metadata)
adonis2(abund_table ~ Sex, permutations = 999,
        method = "robust.aitchison", data = pfnd_metadata)
adonis2(abund_table ~ Age_at_recruitment, permutations = 999,
        method = "robust.aitchison", data = pfnd_metadata)

abund_robust.aitchison = vegdist(abund_table, "robust.aitchison")
RA_pcoa = cmdscale(abund_robust.aitchison, k=2, eig=T)
RA_pcoa_eig <- RA_pcoa$eig
RA_total_variance <- sum(RA_pcoa_eig[RA_pcoa_eig > 0])

percentage_explained_pco1 <- (RA_pcoa_eig[1] / RA_total_variance) * 100
percentage_explained_pco2 <- (RA_pcoa_eig[2] / RA_total_variance) * 100
RA_pcoa_coord = RA_pcoa$points
colnames(RA_pcoa_coord) = c("PCoA1", "PCoA2")

plot.data <- cbind(pfnd_metadata, RA_pcoa_coord)

pacbio_pcoa = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Diagnosis, fill = Diagnosis),
             size = 3) + 
  stat_ellipse(aes(fill = Diagnosis, color = Diagnosis),
               alpha=0.2, level = 0.95, geom = "polygon") +
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_2) + 
  scale_colour_manual(values = Palette_2) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Healthy vs. IPF patients",
       subtitle = glue("Permanova: p={round(permanova_p, 3)}\n"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face="italic"),
    legend.position.inside = c(0.1,0.90), legend.title = element_text(size=18, face = "bold", colour = "firebrick"),
    legend.text = element_text(size=16, face = "italic"), legend.position = "inside",
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

pacbio_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA1,
                   fill = Diagnosis, alpha=0.2),
               show.legend = F) + coord_flip() +
  #scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(x = "none", y = "none") +
  scale_fill_manual(values = Palette_2) 

pacbio_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA2,
                   fill = Diagnosis, alpha=0.2),
               show.legend = F) +
  #scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  guides(x = "none", y = "none") +
  scale_fill_manual(values = Palette_2)

pacbio_pcoa_fig = pacbio_pcoa %>%
  aplot::insert_bottom(pacbio_pcoa1, height = 0.1) %>%
  aplot::insert_right(pacbio_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

pacbio_pcoa_fig

## Strep general
PROFOUND_species = read_tsv("data/Species-Absolute-Abundance-Counts.txt")
PROFOUND_species = column_to_rownames(PROFOUND_species, "Taxa")
PROFOUND_species = as.data.frame(t(PROFOUND_species))
PROFOUND_species  =  PROFOUND_species/rowSums(PROFOUND_species)*100
PROFOUND_species = as.data.frame(t(PROFOUND_species))
PROFOUND_species = rownames_to_column(PROFOUND_species, "Taxa")
PROFOUND_species = PROFOUND_species %>%
  mutate(
    Taxa = str_replace_all(Taxa, "\\[", ""),
    Taxa = str_replace_all(Taxa, "\\]", "")) %>%
  dplyr::rename("PFND164" = "PFND164_amplified") %>% # remove sample that wasn't amplified to keep same
  dplyr::select(-MTN008)
PROFOUND_species = column_to_rownames(PROFOUND_species, "Taxa")
pfnd_strep = PROFOUND_species %>%
  dplyr::filter(grepl("Streptococcus", rownames(.))) 
pfnd_strep = as.data.frame(t(pfnd_strep))

top = pfnd_strep[,order(colSums(pfnd_strep),decreasing=TRUE)]
N = 20
taxa_list = colnames(top)[1:N]
N = length(taxa_list)
top = data.frame(top[,colnames(top) %in% taxa_list])
top

other = pfnd_strep[,order(colSums(pfnd_strep),decreasing=TRUE)]
N = dim(pfnd_strep)[2]
taxa_list2 = colnames(other)[21:N]
N = length(taxa_list2)
Others = data.frame(other[,colnames(other) %in% taxa_list2])
Others = rowSums(Others)
Others = as.data.frame(Others)
Others

strep_df_other = cbind(top,Others)
total_abundance = colSums(strep_df_other)

## Stacked bar plot
strep_long_all = pfnd_strep %>%
  rownames_to_column("PatientID") %>%
  dplyr::left_join(pfnd_metadata[c("PatientID", "Reflux_treatment")]) %>%
  reshape2::melt(id.vars = c("PatientID", "Reflux_treatment"), 
       variable.name = "Species", 
       value.name = "value") %>%
  mutate(Diagnosis = ifelse(grepl("^PFND", PatientID), "IPF", "Healthy"))

strep_df_summary  = strep_long_all %>%
  group_by(Species, Diagnosis) %>%
  summarise(mean_value = mean(value),
            median_value = median(value),
            sd = sd(value),
            q1 = quantile(value, 0.25),
            q3 = quantile(value, 0.75))

total_abundance = colSums(pfnd_strep)
ordered_taxa = names(sort(total_abundance, decreasing = T))
strep_df_summary$Species = factor(strep_df_summary$Species, levels = ordered_taxa)

strep_comparison = strep_long_all %>%
  group_by(PatientID, Diagnosis) %>%
  summarise(TotalAbundance = sum(value, na.rm = TRUE)) %>% ungroup() 

strep_comparison %>%
  group_by(Diagnosis) %>%
  summarise(mean_value = mean(TotalAbundance),
            median_value = median(TotalAbundance),
            sd = sd(TotalAbundance),
            q1 = quantile(TotalAbundance, 0.25),
            q3 = quantile(TotalAbundance, 0.75))

wilcox.test(TotalAbundance ~ Diagnosis, data = strep_comparison)

strep_species = unique(strep_long_all$Species)

strep_df_other = dplyr::left_join(pfnd_metadata, strep_df_other %>% rownames_to_column("PatientID"))

wilcox_results_strep <- list()
for (species in strep_species) {
  clean_species <- gsub(" ", ".", species)
  if (clean_species %in% colnames(strep_df_other)) {
    formula_str <- as.formula(paste0("`", clean_species, "` ~ Diagnosis"))
    test_result <- wilcox.test(formula_str, data = strep_df_other)
    wilcox_results_strep[[species]] <- test_result$p.value
  } else {
    warning(paste("Column not found:", clean_species))
  }
}

results_df <- data.frame(
  Species = names(wilcox_results_strep),
  P_value = unlist(wilcox_results_strep)
)

Palette_many_smooth <- colorRampPalette(Palette_10)(10)

library(ggtext)
ggplot(strep_df_summary %>% filter(median_value>0), aes(x=Species, y=median_value)) + 
  geom_bar(aes(y = median_value, x = Species, fill = Species),
           stat="identity") + facet_wrap(Diagnosis~.)+
  geom_errorbar(aes(x=Species, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5) +
  theme_classic() + scale_fill_manual(values = Palette_10) + 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        strip.text = element_text(size=18),
        axis.text.y = element_text(face = "italic"),
        plot.caption = element_markdown(size = 12, color = "grey30", hjust = 0, margin = margin(t = 15)),
        plot.title = element_markdown(size=18),
        legend.position = "none") + coord_flip() + 
  labs(title = "*Streptococcus* species relative abundance in Healthy and IPF metagenomes",
       x = "", y = "Median relative abundance (%,IQR)",
       caption = "Only showing *Streptococcus* species with a median relative abundance of > 0%")

## Streptococcus heatmap ##
abund_raw <- read_tsv("data/Species-Absolute-Abundance-Counts.txt")

abund_mat <- abund_raw %>%
  column_to_rownames(var = colnames(abund_raw)[1]) %>%
  as.matrix()

strep_mat <- abund_mat[grepl("Streptococcus", rownames(abund_mat), ignore.case = TRUE), ]

# Ensure sample order matches metadata
metadata <- pfnd_metadata %>% dplyr::filter(PatientID %in% colnames(strep_mat))
strep_mat <- strep_mat[, metadata$PatientID]

strep_mat <- strep_mat[rowSums(strep_mat) > 0, ]

strep_log <- log10(strep_mat + 1e-5)  # pseudocount to handle zeros

strep_log[strep_mat == 0] <- NA

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(min(strep_log, na.rm = TRUE), 
             mean(c(min(strep_log, na.rm = TRUE), max(strep_log, na.rm = TRUE))),  # midpoint
             max(strep_log, na.rm = TRUE)),
  colors = c("white", "#6F99AD", "#BC3C29")
)
# Group annotation colours
group_colours <- c("Healthy" = "#BC3C29FF", "IPF" = "#0072B5FF")

col_annotation <- HeatmapAnnotation(
  Diagnosis = metadata$Diagnosis,
  col = list(Diagnosis = group_colours),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Diagnosis = list(title = "Diagnosis")
  )
)

ht <- Heatmap(
  strep_log,
  col = col_fun,
  na_col = "white",
  
  # Annotations
  top_annotation = col_annotation,
  cluster_columns = FALSE,      
  cluster_rows = F,             
  show_column_dend = FALSE,
  show_row_dend = TRUE,
  
  # Labels
  column_title = NULL,
  row_title = "Species",
  row_title_side = "left",
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9, fontface = "italic", fontfamily = "sans"),
  column_names_gp = gpar(fontsize = 8),
  
  # Legend
  heatmap_legend_param = list(
    title = "log10(rel. abund.)",
    at = c(min(strep_log, na.rm = TRUE), 
           mean(c(min(strep_log, na.rm = TRUE), max(strep_log, na.rm = TRUE))),
           max(strep_log, na.rm = TRUE)),
    labels = round(c(min(strep_log, na.rm = TRUE), 
                     mean(c(min(strep_log, na.rm = TRUE), max(strep_log, na.rm = TRUE))),
                     max(strep_log, na.rm = TRUE)), 3),
    legend_height = unit(3, "cm")),
  column_split = metadata$Diagnosis,
  column_gap = unit(3, "mm"),
  width = unit(0.6 * ncol(strep_log), "cm"),
  height = unit(0.5 * nrow(strep_log), "cm"))

ht
