setwd("~/Documents/GitHub/PROFOUND-Baseline/")

library(tidyverse)
library(reshape2)
library(Maaslin2)
library(svglite)
library(ggplotify)
library(patchwork)

# --- Load path abundances from Humann3 --- #
#### Unstratified ####
pathabund = read_tsv("data/all_pathabundance_relabund_unstratified.tsv")

pathabund = pathabund %>%
  dplyr::filter(!grepl("UNMAPPED", `# Pathway`)) %>%
  dplyr::filter(!grepl("UNINTEGRATED", `# Pathway`)) %>%
  column_to_rownames("# Pathway")

colnames(pathabund) = str_split_i(colnames(pathabund), "_", 1)
pathabund = pathabund %>% dplyr::rename("PFND164" = "PFND164-amplified")

# --- Load metadata --- #
pfnd_metadata = read_csv("data/16SSequenced_complete-metadata.csv")
pfnd_metadata = pfnd_metadata %>%
  dplyr::select(PatientID, Sex, Age_at_recruitment, Diagnosis, Smoking_history, Reflux_treatment,
                Baseline_ppFVC, Baseline_ppDLCO, FVC_change_1y,
                Survival_death, Death) %>%
  filter(PatientID %in% colnames(pathabund)) %>%
  column_to_rownames("PatientID")

# --- Maaslin2 --- #
library(pheatmap)
fit_data <- Maaslin2(
  input_data = pathabund,
  input_metadata = pfnd_metadata,
  output = "data/Maaslin2-humann3-unstratified",
  fixed_effects = c("Diagnosis"),
  reference=c("Diagnosis,Healthy"),
  standardize = FALSE, cores=1)

all_results <- fit_data$results
maaslin2_names = make.names(rownames(pathabund))
original_names = rownames(pathabund)
name_map <- tibble(
  maaslin2_name = maaslin2_names,
  original_name = original_names)

all_results_labelled = all_results %>%
  dplyr::left_join(name_map, by = c("feature" = "maaslin2_name"))

# individual samples
sig_pathways_results <- all_results_labelled %>%
  dplyr::filter(pval <= 0.05) %>%
  dplyr::filter(!grepl("\\|", original_name)) %>%
  dplyr::filter(metadata == "Diagnosis") %>% 
  #filter(value != "Healthy") %>%
  # select top 20
  mutate(abs_coef = abs(coef)) %>%
  arrange(desc(abs_coef)) %>% head(50) %>% # select top 50
  dplyr::select(original_name, coef, pval, qval) %>%
  dplyr::rename(Pathway = original_name,
                Coef_IPF_vs_Healthy = coef,
                Qvalue = qval)

sig_pathway_names <- sig_pathways_results$Pathway

# Manually renamed the pathways with MetaCyc database for better interpretation. 
# curated June 2026.
classify_pathway <- function(description) {
  case_when(
    # Amino acid biosynthesis
    grepl("L-methionine|L-cysteine|L-glutamine|L-phenylalanine|L-tyrosine|L-lysine|L-aspartate|L-asparagine|seleno-amino acid", description, ignore.case = TRUE) ~ "Amino acid biosynthesis",
    
    # NAD Metabolism
    grepl("NAD salvage|NAD de novo|nicotinamide|PNC", description, ignore.case = TRUE) ~ "NAD Metabolism",
    
    # Cofactor, carrier and vitamin biosynthesis
    grepl("tetrahydrofolate|thiamine|pantothenate|heme|glutamyl cycle", description, ignore.case = TRUE) ~ "Cofactor, carrier and vitamin biosynthesis",
    
    # Nucleotide or nucleoside biosynthesis
    grepl("purine nucleotide salvage|pyrimidine ribonucleoside|inosine|nucleotide|nucleoside", description, ignore.case = TRUE) ~ "Nucleotide or nucleoside biosynthesis",
    
    # Degradation, utilisation, assimilation
    grepl("degradation|nylon|inositol degradation", description, ignore.case = TRUE) ~ "Degradation, utilisation, assimilation",
    
    # Carbohydrate biosynthesis
    grepl("fructan biosynthesis|UDP-glucose|O-antigen building blocks biosynthesis|UDP-N-acetyl", description, ignore.case = TRUE) ~ "Carbohydrate biosynthesis",
    
    # Carbohydrate degradation
    grepl("lactose|galactose degradation|rhamnose degradation", description, ignore.case = TRUE) ~ "Carbohydrate degradation",
    
    # Fatty acid and lipid biosynthesis
    grepl("petroselinate|CDP-diacylglycerol|fatty acid|lipid", description, ignore.case = TRUE) ~ "Fatty acid and lipid biosynthesis",
    
    # Fermentation
    grepl("butanediol|homolactic fermentation|fermentation", description, ignore.case = TRUE) ~ "Fermentation",
    
    # Glycolysis
    grepl("glycolysis", description, ignore.case = TRUE) ~ "Glycolysis",
    
    # Photosynthesis
    grepl("photosynthetic carbon assimilation", description, ignore.case = TRUE) ~ "Photosynthesis",
    
    # Generation of precursor metabolites and energy
    grepl("PEPCK", description, ignore.case = TRUE) ~ "Generation of precursor metabolites and energy",

    # Secondary metabolite biosynthesis
    grepl("mevalonate|geranylgeranyldiphosphate", description, ignore.case = TRUE) ~ "Secondary metabolite biosynthesis",
    
    # Amines and polyamines biosynthesis
    grepl("UDP-N-acetyl-D-glucosamine", description, ignore.case = TRUE) ~ "Amines and polyamines biosynthesis",
    
    # Other
    TRUE ~ "Other"
  )
}

pathway_group_map <- sig_pathways_results %>%
  mutate(
    description = gsub(".*: ", "", Pathway),  # extract text after the colon
    Group = classify_pathway(description)
  ) %>%
  select(Pathway, Group)

# Check any that fell into "Other" - these need manual review
pathway_group_map %>% filter(Group == "Other") %>% select(Pathway)

heatmap_matrix <- pathabund %>%
  filter(rownames(.) %in% sig_pathway_names) %>%
  rownames_to_column("Pathway") %>%
  dplyr::left_join(pathway_group_map, by = "Pathway") %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Group")

# Rescale the rows to show how abundance changes from the overall average abundance of that pathway across samples
# similar concept to CLR, you're comparing how each feature in a sample deviates from the overall mean
# CLR uses geometric mean, scale() uses arithmetic means.
# To handle compositionality and differences in absolute abundances
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix))) 
# sample order matches the matrix columns
samples_ordered <- pfnd_metadata[colnames(heatmap_matrix_scaled), "Diagnosis", drop = FALSE]

annotation_col <- data.frame(Diagnosis = samples_ordered)
rownames(annotation_col) = colnames(heatmap_matrix_scaled)

# Create a row annotation data frame for the heatmap
# qval annotation
annotation_row <- sig_pathways_results %>%
  dplyr::left_join(pathway_group_map) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(min_Q = min(Qvalue, na.rm = TRUE)) %>% 
  dplyr::mutate(`Q < 0.25` = ifelse(min_Q < 0.25, 'Yes', 'No'),
                `Q < 0.05` = ifelse(min_Q < 0.05, 'Yes', 'No')) %>%
  dplyr::select(Group, `Q < 0.25`, `Q < 0.05`) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Group")

# Make sure the row order exactly matches your heatmap matrix rows
annotation_row <- annotation_row[rownames(heatmap_matrix_scaled), , drop = FALSE]

ann_colors = list(
  Diagnosis = c(Healthy = "#BC3C29FF", IPF = "#0072B5FF"), 
  `Q < 0.25` = c(Yes = "#E18727", No = "gray90"),
  `Q < 0.05` = c(Yes = "#20854E", No = "gray90"))

max_abs_value <- max(abs(min(heatmap_matrix_scaled)), abs(max(heatmap_matrix_scaled)))
breaks_list <- seq(-max_abs_value, max_abs_value, length.out = 101)
colors_list <- colorRampPalette(c("#003f5a", "white", "#BC3C29"))(100)

p <- as.ggplot(pheatmap(heatmap_matrix_scaled,
                        cluster_rows = TRUE, 
                        cluster_cols = FALSE,
                        annotation_col = annotation_col,
                        annotation_row = annotation_row,
                        annotation_colors = ann_colors,
                        show_colnames = FALSE,
                        fontsize = 12,
                        main = "Differential Pathway Abundance (Z-Score)",
                        color = colors_list,
                        breaks = breaks_list,
                        silent = TRUE,
                        legend = TRUE))

final_plot <- p + plot_annotation(
  caption = "
  Heatmap shows functional groups containing pathways with a nominal p < 0.05 for the Diagnosis effect. 
  A functional group is annotated as significant if at least one of its underlying pathways passed the respective FDR Q-value thresholds.") &
  theme(plot.caption = element_text(size = 12))

library(grid)
print(final_plot)
grid.text("Z-Score", x = 0.90, y = 0.7, gp = gpar(fontsize = 11, fontface = "bold"))

# library(Cairo)
# CairoSVG("Unstratified-heatmap-humann3.svg", width=15, height=10)
# print(final_plot)
# grid::grid.text("Z-Score", x = 0.90, y = 0.74, gp = gpar(fontsize = 11, fontface = "bold"))
# dev.off()

#### Stratified ####
pathabund = read_tsv("data/all_pathabundance_relabund_stratified.tsv")
pathway_regex <- paste(stringr::str_escape(sig_pathway_names), collapse = "|")

# Filter only on significant pathways names from previous
sig_ids <- gsub(":.*", "", sig_pathway_names)
id_regex <- paste(paste0("^", stringr::str_escape(sig_ids)), collapse = "|")
pathabund <- pathabund %>%
  dplyr::filter(str_detect(`# Pathway`, pattern = id_regex))

pathabund = pathabund %>%
  dplyr::filter(!grepl("UNMAPPED", `# Pathway`)) %>%
  dplyr::filter(!grepl("UNINTEGRATED", `# Pathway`)) %>%
  dplyr::filter(str_detect(`# Pathway`, pattern = id_regex)) %>%
  column_to_rownames("# Pathway")

colnames(pathabund) = str_split_i(colnames(pathabund), "_", 1)
pathabund = pathabund %>% dplyr::rename("PFND164" = "PFND164-amplified")

# --- Maaslin2 --- #
library(pheatmap)
fit_data <- Maaslin2(
  input_data = pathabund,
  input_metadata = pfnd_metadata,
  output = "data/Maaslin2-humann3-stratified",
  fixed_effects = c("Diagnosis"),
  reference=c("Diagnosis,Healthy"),
  standardize = FALSE, cores=1)

stratified_results = fit_data$results

maaslin2_names = make.names(rownames(pathabund))
original_names = rownames(pathabund)
name_map <- tibble(
  maaslin2_name = maaslin2_names,
  original_name = original_names
)

maaslin2_results = stratified_results %>%
  dplyr::left_join(name_map, by = c("feature" = "maaslin2_name"))

sig_contributions <- maaslin2_results %>%
  dplyr::filter(pval <= 0.05) %>%
  tidyr::separate(original_name, into = c("Pathway", "Species"), sep = "\\|", fill = "right", extra = "merge") %>%
  dplyr::filter(!is.na(Species), Species != "unclassified") %>%
  dplyr::mutate(
    Species = str_remove_all(Species, "[kpcofgs]__"),
    Species = sub("^[^.]+\\.(.*)$", "\\1", Species),
    Species = str_replace_all(Species, "\\|", " "),
    Species = trimws(Species),
    heatmap_metric = -log10(pval) * sign(coef)
  )

diag_contributions <- sig_contributions %>%
  dplyr::left_join(pathway_group_map, by = "Pathway")

# check no duplicates
diag_contributions %>% 
  count(Pathway, Species) %>% 
  filter(n > 1)

top_species <- diag_contributions %>% 
  dplyr::count(Species, wt = abs(heatmap_metric)) %>% 
  dplyr::top_n(30, n) %>% pull(Species)

plot_data <- diag_contributions %>%
  dplyr::filter(Species %in% top_species) %>%
  dplyr::filter(!Group == "Other")

species_order <- sort(unique(as.character(plot_data$Species)))

plot_data <- plot_data %>%
  dplyr::mutate(
    Species = factor(Species, levels = species_order),
    # Clean up long pathway names if needed for the plot display
    Pathway_Short = stringr::str_trunc(Pathway, 40) 
  )

contribution_heatmap <- ggplot(plot_data, 
                               aes(x = Species, y = Pathway_Short, fill = heatmap_metric)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#003f5a", 
    mid = "white",
    high = "#BC3C29",
    midpoint = 0,
    name = "(-log10(p-value) * sign(Coef))"
  ) +
  # Split the y-axis by functional group dynamically
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Species Contributions to Functional Pathways",
    subtitle = "Stratified Maaslin2 Analysis: IPF vs. Healthy",
    x = "Contributing Species",
    y = ""
  ) +
  theme_bw() + 
  theme(
    # Clean up X-axis: 45 degrees text alignment
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 9),
    
    # Clean up Facet Strips: Positioned horizontally on top of each group
    strip.text.x = element_text(size = 10, face = "bold", hjust = 0),
    strip.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    
    # Layout adjustments
    panel.spacing.y = unit(0.5, "lines"), # Slightly increased for clear visual separation
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 9),
    legend.position = "right"
  )

print(contribution_heatmap)

# ggsave("Stratified-heatmap-humann3.svg", dpi = 600,
#        width = 12, height = 10)
