 #### === PROFOUND BASELINE MS === ####
 # Descriptive analysis: alpha, beta, top ten, analysis with covariates.

# Data related to Figures 1, 2G-H, 3E-J, Supp Fig2, are not shared publicly, please ask. 
 
# Script to create the figures for the manuscript 

# Load libraries
library(renv)
setwd("~/Documents/GitHub/PROFOUND-Baseline/")
renv::load()
renv::status()
renv::restore()

library(dplyr)
library(tidyverse)
library(ggplot2)
library(qiime2R)
library(phyloseq)
library(decontam)
library(microbiome)
library(indicspecies)
library(ComplexHeatmap)
library(circlize)
library(effectsize)
library(ggpubr)
library(OTUtable)
library(glue)
library(patchwork)
library(survminer)
library("survival")
library(ecodist)
library(aplot)
library(ggplotify)
library(reshape2)
library(ggsci)
library(vegan)

#### Plot aesthetics ####
Palette_color = c("#BC3C29FF", "#0072B5FF")
Palette_fill = c("#BC3C29FF", "#0072B5FF")

Palette_IPF_color = c("#E18727FF","#20854EFF")
Palette_IPF_fill = c("#E18727FF", "#20854EFF")

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

Palette_diverging1 <- c("#003f5a", "#fec682", "#BC3C29")
Palette_diverging1 <- colorRampPalette(Palette_diverging1)
Palette_diverging1 <- Palette_diverging1(5)

##### === Phyloseq processing === ####
## not present in the Github repo, need to run it yourself ##
manifest_metadata <- read_tsv("data/BaselineStool_metadata.txt", trim_ws = T)
manifest_metadata <- column_to_rownames(manifest_metadata, var="sample-id")
physeq <- qza_to_phyloseq(features="data/BaselineStool_filtered-table-gg2.qza",
                          tree="data/BaselineStool_midrooted-tree.qza",
                          taxonomy = "data/BaselineStool_filtered-taxonomy-gg2.qza")
# 
# setdiff(rownames(manifest_metadata), colnames(otu_table(physeq)))
# 
mapping=manifest_metadata
sample_data(physeq)<-mapping
# 
# gplots::venn(list(rownames(sample_data(physeq)), rownames(mapping)))
# setdiff(rownames(mapping), rownames(sample_data(physeq)))
ps<-physeq
summarize_phyloseq(ps)

df_library <- data.frame(sample_data(ps))
df_library$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df_library <- df_library[order(df_library$LibrarySize),]
df_library$Index <- seq(nrow(df_library))
#
# # remove samples that have less than 500 reads.
# less_than_500 = df_library[df_library$LibrarySize<500 & !df_library$SampleType=="Sequencing Control",]
# remove_samples = rownames(less_than_500)
# manifest_metadata = manifest_metadata %>%
#   filter(!rownames(manifest_metadata) %in% remove_samples)
# sample_data(physeq)<-manifest_metadata
# ps <- physeq
# sum(taxa_sums(ps) == 0) # how many taxa aren't present in ANY samples
# summarize_phyloseq(ps)
# ps


tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla)

ps1 = subset_taxa(ps, !(!Phylum %in% filterPhyla))
ps1 # Only keep the Phyla in filterPhyla in the filtered reads dataset
summarize_phyloseq(ps1)

# Check if there are unique Phyla names
unique(as.data.frame(as(tax_table(ps1), "matrix"))$Phylum)

df_genus = aggregate_taxa(ps1, "Genus")
otu_genus = as.data.frame(otu_table(df_genus))

#Filter at 0.005%
minTotRelAbun = 0.00005
x = taxa_sums(ps1)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), ps1)

# Save OTU count data
otu_genus = aggregate_taxa(prunedSet, "Genus")
otu_genus = as.data.frame(otu_genus@otu_table)
otu_genus = rownames_to_column(otu_genus, "Genus")

## Now normalise the counts.
normalizeSample = function(x) {
  x/sum(x)
}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
Phylum_df<-as.data.frame(Controls_Phylum@otu_table)

Controls_Family <- aggregate_taxa(Controls_relative, 'Family')
Family_df<-as.data.frame(Controls_Family@otu_table)

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
Genus_df<-as.data.frame(Controls_Genus@otu_table)

Controls_species <- aggregate_taxa(Controls_relative, "Species")
Species_df <- as.data.frame(Controls_species@otu_table)

# Join reads to metadata
metadata <- rownames_to_column(manifest_metadata, var = "sample-id")

Phylum_df <- data.frame(t(Phylum_df))
Phylum_df <- rownames_to_column(Phylum_df, "sample-id")
Phylum_df <- inner_join(metadata, Phylum_df, by = "sample-id")

Family_df <- data.frame(t(Family_df))
Family_df <- rownames_to_column(Family_df, "sample-id")
Family_df <- inner_join(metadata, Family_df, by = "sample-id")

Genus_df <- data.frame(t(Genus_df))
Genus_df <- rownames_to_column(Genus_df, "sample-id")
Genus_df <- inner_join(metadata, Genus_df, by = "sample-id")

Species_df <- data.frame(t(Species_df))
Species_df <- rownames_to_column(Species_df, "sample-id")
Species_df <- inner_join(metadata, Species_df, by = "sample-id")

#### Taxa clean up ####
# Greengenes2 adds clades to the genus names so need to adjust for this
meta_table <- Genus_df %>%
  dplyr::select(`sample-id`:Diagnosis)

# Select sample-id and all taxanomic columns (columns 7 to end)
abund_table <- Genus_df %>% 
  dplyr::select(`sample-id`, 7:ncol(Genus_df)) %>%
  tibble::column_to_rownames(var = "sample-id")

# Transpose to Genus x Samples
abund_table <- as.data.frame(t(abund_table)) %>%
  tibble::rownames_to_column(var = "Genus")

# Clean Greengenes2 genus names (extract primary name before suffix/clade underscores)
abund_table$Genus <- str_extract(abund_table$Genus, "^[^_]+")

# Aggregate reads for duplicate cleaned genus names
abund_table_aggregated <- aggregate(. ~ Genus, abund_table, sum)

# Filter out "Unknown" / unassigned taxa BEFORE renormalisation
# (Adjust filter terms if Greengenes2 labels them differently, e.g., "g__" or "Unassigned")
abund_table_aggregated <- abund_table_aggregated %>%
  dplyr::filter(!is.na(Genus) & 
                  !Genus %in% c("Unknown", "unknown", "unassigned", "Unassigned", "g__", ""))

# Format back to Samples x Genus matrix
abund_table_aggregated <- abund_table_aggregated %>%
  tibble::column_to_rownames(var = "Genus")
abund_table_aggregated <- as.data.frame(t(abund_table_aggregated))

# Keep only samples with >0 total reads remaining
abund_table_aggregated <- abund_table_aggregated[rowSums(abund_table_aggregated) > 0, ]

# RENORMALISE to 100% after removing Unknowns
abund_table_aggregated <- (abund_table_aggregated / rowSums(abund_table_aggregated)) * 100

# Verify that all row sums equal 100
print(head(rowSums(abund_table_aggregated)))

# Join back with metadata
abund_table_aggregated <- tibble::rownames_to_column(abund_table_aggregated, var = "sample-id")
Genus_df <- dplyr::inner_join(meta_table, abund_table_aggregated, by = "sample-id")

write_csv(Genus_df, "data/Genus-aggregated-normalised-metadata.csv")

##### PERMANOVA #####
##Factors influencing the gut microbiota ##
df <- read_csv("data/Genus-aggregated-normalised-metadata.csv")
df_baseline <- df %>%
  filter(Diagnosis == "IPF" & Visit ==1 | Diagnosis =="Healthy" & Visit==1)
rownames(df_baseline) <- df_baseline$`sample-id`
abund_table <- df_baseline %>%
  dplyr::select(7:ncol(df_baseline))
rownames(abund_table) <- df_baseline$`sample-id`
meta_table <- df_baseline %>% 
  dplyr::select(`sample-id`, PatientID, Visit, Diagnosis)

more_metadata = read_csv("data/16SSequenced_complete-metadata.csv")
more_metadata = more_metadata %>%
  mutate(Smoking_history = ifelse(Smoking_history == "Current smoker" |
                                    Smoking_history == "Ex-smoker", "Ex-smoker", Smoking_history),
         Progression_status_1y = ifelse(Progression_status_1y == "Semi-progressive", "Stable", Progression_status_1y))
meta_table = dplyr::left_join(meta_table, more_metadata)

meta_table = meta_table %>%
  filter(!is.na(Sex) & !is.na(Smoking_history)) # lose 3 patients

abund_table = abund_table %>%
  filter(rownames(abund_table) %in% meta_table$`sample-id`)
meta_factors = c("Diagnosis", "Sex", "Reflux_treatment", "Smoking_history", "Age")

##### Healthy vs IPF #####
# Adonis2 to determine what metadata variables could potentially influence the gut microbiota
RNGversion("4.4.2")
set.seed(123)
perm_diagnosis = adonis2(abund_table ~ Diagnosis, permutations = 999, method = "robust.aitchison", data = meta_table)
perm_Sex = adonis2(abund_table ~ Sex, permutations = 999, method = "robust.aitchison", data = meta_table)
perm_refluxTreat = adonis2(abund_table ~ Reflux_treatment, permutations = 999, method = "robust.aitchison", data = meta_table)
perm_smoking = adonis2(abund_table ~ Smoking_history, permutations = 999, method = "robust.aitchison", data = meta_table)
perm_age = adonis2(abund_table ~ Age_at_recruitment, permutations = 999, method="robust.aitchison", data =meta_table)

permanova_p = data.frame(perm_diagnosis$`Pr(>F)`,
                         perm_Sex$`Pr(>F)`, 
                         perm_refluxTreat$`Pr(>F)`, perm_smoking$`Pr(>F)`,
                         perm_age$`Pr(>F)`) %>%
  drop_na()
colnames(permanova_p) = meta_factors
rownames(permanova_p) = "p_value"
permanova_p = as.data.frame(t(permanova_p))
permanova_p = rownames_to_column(permanova_p, "meta_factors")
permanova_p

####### PCoA plots for significant factors #######
###### Diagnosis ######
# PCoA plot- Diagnosis
abund_table_filtered_robust.aitchison = vegdist(abund_table, "robust.aitchison")
RA_pcoa = cmdscale(abund_table_filtered_robust.aitchison, k=2, eig=T)  

RA_pcoa_eig <- RA_pcoa$eig
RA_total_variance <- sum(RA_pcoa_eig[RA_pcoa_eig > 0])

percentage_explained_pco1 <- (RA_pcoa_eig[1] / RA_total_variance) * 100
percentage_explained_pco2 <- (RA_pcoa_eig[2] / RA_total_variance) * 100

RA_pcoa_coord = RA_pcoa$points

colnames(RA_pcoa_coord) = c("PCoA1", "PCoA2")

plot.data <- cbind(meta_table, RA_pcoa_coord)

fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Diagnosis, fill = Diagnosis),
             size = 3) + 
  stat_ellipse(aes(fill = Diagnosis, color = Diagnosis),
               alpha=0.2, level = 0.95, geom = "polygon") +
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_fill) + 
  scale_colour_manual(values = Palette_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Healthy vs. IPF patients",
       subtitle = glue("Permanova: p={round(permanova_p[1,2], 5)}\n"),
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

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA1,
                   fill = Diagnosis, alpha=0.2),
               show.legend = F) + coord_flip() +
  #scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(x = "none", y = "none") +
  scale_fill_manual(values = Palette_fill) 

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA2,
                   fill = Diagnosis, alpha=0.2),
               show.legend = F) +
  #scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  guides(x = "none", y = "none") +
  scale_fill_manual(values = Palette_fill)

fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Reflux treatment ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Reflux_treatment, fill = Reflux_treatment),
             size = 3) + 
  stat_ellipse(aes(fill = Reflux_treatment, color = Reflux_treatment),
               alpha=0.2, level = 0.95, geom = "polygon") +
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_fill) + 
  scale_colour_manual(values = Palette_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Reflux treatment vs. None",
       subtitle = glue("Permanova: p={round(permanova_p[3,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), 
    legend.title = element_text(size=18, face = "bold", colour = "firebrick"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Reflux_treatment, y=PCoA1,
                   fill = Reflux_treatment, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill) +
  guides(x = "none", y = "none")
  
fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Reflux_treatment, y=PCoA2,
                   fill = Reflux_treatment, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill) +
  guides(x = "none", y = "none")
  
fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Smoking status ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Smoking_history, fill = Smoking_history),
             size = 3) + 
  stat_ellipse(aes(fill = Smoking_history, color = Smoking_history),
               alpha=0.2, level = 0.95, geom = "polygon") +  # scale_"" is used to design the plot
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) + 
  scale_colour_manual(values = c("Ex-smoker" = "#B22222", "Never" = "black")) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Ex-smokers vs. Never: Beta Diversity (Robust Aitchison)",
       subtitle = glue("Permanova: p={round(permanova_p[4,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Smoking_history, y=PCoA1,
                   fill = Smoking_history, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) +
  guides(x = "none", y = "none")
  
fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Smoking_history, y=PCoA2,
                   fill = Smoking_history, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) +
  guides(x = "none", y = "none")
  
fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Sex ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Sex, fill = Sex),
             size = 3) + 
  stat_ellipse(aes(fill = Sex, color = Sex),
               alpha=0.2, level = 0.95, geom = "polygon") +  # scale_"" is used to design the plot
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_fill) + 
  scale_colour_manual(values = Palette_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Male vs. Female",
       subtitle = glue("Permanova: p={round(permanova_p[2,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Sex, y=PCoA1,
                   fill = Sex, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill) +
  guides(x = "none", y = "none")

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Sex, y=PCoA2,
                   fill = Sex, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill) +
  guides(x = "none", y = "none")
  
fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Age ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(colour = Age_at_recruitment),
             size = 3) + 
  scale_color_gradientn(colors = c(low = "#003f5a", medium="#fec682", high = "#BC3C29")) +
  labs(title = "Age",
       subtitle = glue("Permanova: p={round(permanova_p[5,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = ""),
       color = "Age", fill = "Age") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.title = element_text(size=18, face="bold", colour = "firebrick"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"),
    aspect.ratio = 1)

fig

##### IPF ONLY #####
df_IPF <- df %>%
  dplyr::filter(Diagnosis == "IPF" & Visit ==1)
rownames(df_IPF) <- df_IPF$`sample-id`
abund_IPF <- df_IPF %>%
  dplyr::select(7:ncol(df_IPF))
rownames(abund_IPF) <- df_IPF$`sample-id`
meta_IPF = meta_table %>%
  filter(`sample-id` %in% rownames(df_IPF))

abund_IPF = abund_IPF %>%
  filter(rownames(abund_IPF) %in% meta_IPF$`sample-id`)

meta_factors = c("Status", "Sex", "Reflux_treatment", "Smoking_history", "Age")

# Adonis2 to determine what metadata variables could potentially influence the gut microbiota
perm_diagnosis = adonis2(abund_IPF ~ Progression_status_1y, permutations = 999, method = "robust.aitchison", data = meta_IPF)
perm_Sex = adonis2(abund_IPF ~ Sex, permutations = 999, method = "robust.aitchison", data = meta_IPF)
perm_refluxTreat = adonis2(abund_IPF ~ Reflux_treatment, permutations = 999, method = "robust.aitchison", data = meta_IPF)
perm_smoking = adonis2(abund_IPF ~ Smoking_history, permutations = 999, method = "robust.aitchison", data = meta_IPF)
perm_age = adonis2(abund_IPF ~ Age_at_recruitment, permutations = 999, method="robust.aitchison", data =meta_IPF)

permanova_p = data.frame(perm_diagnosis$`Pr(>F)`,
                         perm_Sex$`Pr(>F)`, 
                         perm_refluxTreat$`Pr(>F)`, perm_smoking$`Pr(>F)`,
                         perm_age$`Pr(>F)`) %>%
  drop_na()
colnames(permanova_p) = meta_factors
rownames(permanova_p) = "p_value"
permanova_p = as.data.frame(t(permanova_p))
permanova_p = rownames_to_column(permanova_p, "meta_factors")

permanova_p

##### PCoA IPF data #####
abund_IPF_robust.aitchison = vegdist(abund_IPF, "robust.aitchison")
RA_pcoa_IPF = cmdscale(abund_IPF_robust.aitchison, k=2, eig=T)  

RA_pcoa_IPF_eig <- RA_pcoa_IPF$eig
RA_IPF_total_variance <- sum(RA_pcoa_IPF_eig[RA_pcoa_IPF_eig > 0])

percentage_explained_pco1 <- (RA_pcoa_IPF_eig[1] / RA_IPF_total_variance) * 100
percentage_explained_pco2 <- (RA_pcoa_IPF_eig[2] / RA_IPF_total_variance) * 100

RA_pcoa_IPF_coord = RA_pcoa_IPF$points

colnames(RA_pcoa_IPF_coord) = c("PCoA1", "PCoA2")

plot.data <- cbind(meta_IPF, RA_pcoa_IPF_coord)

###### Status ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Progression_status_1y, fill = Progression_status_1y),
             size = 3) + 
  stat_ellipse(aes(fill = Progression_status_1y, color = Progression_status_1y),
               alpha=0.2, level = 0.95, geom = "polygon") +  
  scale_fill_manual(values = Palette_IPF_fill) + 
  scale_colour_manual(values = Palette_IPF_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Stable vs. Progressive IPF",
       subtitle = glue("Permanova: p={round(permanova_p[1,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = ""),
       color = "Progression Status", fill = "Progression Status", shape = "Progression Status") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Progression_status_1y, y=PCoA1,
                   fill = Progression_status_1y, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Progression_status_1y, y=PCoA2,
                   fill = Progression_status_1y, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() +   
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Reflux treatment ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Reflux_treatment, fill = Reflux_treatment),
             size = 3) + 
  stat_ellipse(aes(fill = Reflux_treatment, color = Reflux_treatment),
               alpha=0.2, level = 0.95, geom = "polygon") +  # scale_"" is used to design the plot
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_IPF_fill) + 
  scale_colour_manual(values = Palette_IPF_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Reflux treatment vs. None",
       subtitle = glue("Permanova: p={round(permanova_p[3,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Reflux_treatment, y=PCoA1,
                   fill = Reflux_treatment, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Reflux_treatment, y=PCoA2,
                   fill = Reflux_treatment, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Smoking_history ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Smoking_history, fill = Smoking_history),
             size = 3) + 
  stat_ellipse(aes(fill = Smoking_history, color = Smoking_history),
               alpha=0.2, level = 0.95, geom = "polygon") +  # scale_"" is used to design the plot
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) + 
  scale_colour_manual(values = c("Ex-smoker" = "#B22222", "Never" = "black")) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "IPF Patients: Beta Diversity (Robust Aitchison)",
       subtitle = glue("Permanova: p={round(permanova_p[4,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Smoking_history, y=PCoA1,
                   fill = Smoking_history, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) +
  guides(x = "none", y = "none")

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Smoking_history, y=PCoA2,
                   fill = Smoking_history, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = c("Ex-smoker" = "#B22222", "Never" = "white")) +
  guides(x = "none", y = "none")

fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Sex ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(shape = Sex, fill = Sex),
             size = 2) + 
  stat_ellipse(aes(fill = Sex, color = Sex),
               alpha=0.2, level = 0.95, geom = "polygon") +  # scale_"" is used to design the plot  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_IPF_fill) + 
  scale_colour_manual(values = Palette_IPF_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "Male vs. Female",
       subtitle = glue("Permanova: p={round(permanova_p[2,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = c(0.2,0.1), legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"))

fig_pcoa1 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Sex, y=PCoA1,
                   fill = Sex, alpha=0.2),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig_pcoa2 = ggplot(plot.data) + 
  geom_boxplot(aes(x = Sex, y=PCoA2,
                   fill = Sex, alpha=0.2),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_IPF_fill) +
  guides(x="none", y="none")

fig = fig %>%
  aplot::insert_bottom(fig_pcoa1, height = 0.1) %>%
  aplot::insert_right(fig_pcoa2, width=0.1) %>%
  as.ggplot() + theme(aspect.ratio = 1)

fig

###### Age ######
fig = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(colour = Age_at_recruitment),
             size = 3) + 
  scale_color_gradientn(colors = c(low = "#18297b", medium="#eeb4cd", high = "#cd205a")) +
  labs(title = "Age",
       subtitle = glue("Permanova: p={round(permanova_p[5,2], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = ""),
       color = "Age", fill = "Age") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.title = element_text(size=18, colour = "black"),
    legend.text = element_text(size=16, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_rect(fill = "white"), 
    panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0,size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.2, linetype = "solid", colour = "grey"),
    aspect.ratio = 1)

fig

##### Alpha diversity - IPF and Healthy ####
df <- read_csv("Genus-aggregated-normalised-metadata.csv")
df <- df %>%
  filter(Diagnosis == "IPF" & Visit ==1 | Diagnosis =="Healthy" & Visit==1)
rownames(df) <- df$`sample-id`
abund_table <- df %>%
  select(7:ncol(df))
rownames(abund_table) <- df$`sample-id`
meta_table <- df %>% 
  select(`sample-id`, PatientID, Visit, Diagnosis)

rownames(meta_table) <- df$`sample-id`

meta_table = left_join(meta_table, more_metadata) 

meta_table = meta_table %>%
  mutate(Progression_status_1y = as.factor(Progression_status_1y),
         Sex = as.factor(Sex),
         Reflux = as.factor(Reflux))

###### Species number and diagnosis ####
species_number <- vegan::specnumber(abund_table)
species_number <- as.data.frame(species_number)
species_number_metadata <- cbind(meta_table, species_number)
wilcox = wilcox.test(species_number ~ Diagnosis, species_number_metadata)
p_value = wilcox$p.value

fig = ggplot(species_number_metadata, aes(x = Diagnosis, y = species_number),
                fill = Diagnosis, shape = Diagnosis) +
  geom_violin(aes(color = Diagnosis)) +
  geom_boxplot(aes(shape = Diagnosis, colour = Diagnosis), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis), 
             size = 2, stroke = .5, alpha=0.2) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_color)) +
  scale_fill_manual(values=Palette_fill) + 
  labs(title="", 
       caption = glue("Wilcoxon ranked sum p={round(p_value,3)}"),
       x = "", 
       y = "Species number") +
  theme(
    plot.title = element_text(size=20),
    plot.subtitle = element_text(size=18),
    plot.caption = element_text(size=14),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1, size=14),
    axis.title.x = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1)

fig

###### Shannon diversity and diagnosis######
shannon_diversity <- vegan::diversity(abund_table, index="shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(meta_table, shannon_diversity)
shannon_metadata$Diagnosis <- as.factor(shannon_metadata$Diagnosis)
wilcox = wilcox.test(shannon_diversity ~ Diagnosis, shannon_metadata) #p=0.5
p_value = wilcox$p.value

fig = ggplot(shannon_metadata, aes(x = Diagnosis, y = shannon_diversity),
                fill = Diagnosis, shape = Diagnosis) +
  geom_violin(aes(color=Diagnosis))+
  geom_boxplot(aes(shape = Diagnosis, colour = Diagnosis), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis), 
             size = 2, stroke = .5, alpha=0.2) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_color)) +
  scale_fill_manual(values=Palette_fill) +
  labs(title="", 
       caption = glue("Wilcoxon ranked sum p={round(p_value,3)}"),
       x = "",
       y = "Shannon diversity") +
  theme(
    plot.title = element_text(size=20),
    plot.subtitle = element_text(size=18),
    plot.caption = element_text(size=14),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1, size=14),
    axis.title.x = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1)

fig

## Create alpha diversity plot
Alpha_diversity_df = cbind(shannon_metadata, species_number)
Alpha_diversity_df = Alpha_diversity_df %>%
  select(`sample-id`, PatientID, Diagnosis, Age_at_recruitment, Sex, Reflux, Reflux_treatment, species_number, shannon_diversity) %>%
  rename("Observed" = "species_number") %>%
  rename("Shannon" = "shannon_diversity")
Alpha_diversity_df = melt(Alpha_diversity_df, id.vars = c("sample-id",
                                                          "PatientID",
                                                          "Diagnosis",
                                                          "Age_at_recruitment", "Sex",
                                                          "Reflux", "Reflux_treatment"),
                          variable.name = "Diversity index")

Alpha_plot = ggplot(Alpha_diversity_df, aes(x = Diagnosis, y = value),
       fill = Diagnosis, shape = Diagnosis) +
  geom_boxplot(aes(shape = Diagnosis, colour = Diagnosis), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis), 
             size = 3, stroke = .5, alpha=0.2) +
  geom_violin(alpha=0.1, aes(color=Diagnosis)) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_color)) +
  scale_fill_manual(values=c(Palette_fill)) +
  labs(title="", 
       x = "",
       y = "Diversity") + theme_classic() + 
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    legend.position = "none", panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "grey"),
    strip.text = element_text(size=12, face = "bold"),
    strip.background = element_rect(linewidth=0.8, colour = "grey")) + 
    facet_wrap(.~`Diversity index`, scales="free_y") 
Alpha_plot

##### Alpha diversity - IPF ####
df <- read_csv("Genus-aggregated-normalised-metadata.csv")
df_IPF <- df %>%
  filter(Diagnosis == "IPF" & Visit ==1)
rownames(df_IPF) <- df_IPF$`sample-id`
abund_IPF <- df_IPF %>%
  select(7:ncol(df_IPF))
rownames(abund_IPF) <- df_IPF$`sample-id`
meta_IPF <- df_IPF %>% 
  select(`sample-id`, PatientID, Visit, Diagnosis)

rownames(meta_IPF) <- df_IPF$`sample-id`

meta_IPF = left_join(meta_IPF, more_metadata) 

meta_IPF = meta_IPF %>%
  mutate(Progression_status_1y = as.factor(Progression_status_1y),
         Sex = as.factor(Sex),
         Reflux = as.factor(Reflux))

###### Species number by Status ####
species_number <- vegan::specnumber(abund_IPF)
species_number <- as.data.frame(species_number)
species_number_metadata <- cbind(meta_IPF, species_number)
wilcox = wilcox.test(species_number ~ Progression_status_1y, species_number_metadata)
p_value = wilcox$p.value

fig = ggplot(species_number_metadata, aes(x = Progression_status_1y, y = species_number),
             fill = Progression_status_1y, shape = Progression_status_1y) +
  geom_violin(aes(color = Progression_status_1y)) +
  geom_boxplot(aes(shape = Progression_status_1y, colour = Progression_status_1y), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Progression_status_1y, shape = Progression_status_1y, colour = Progression_status_1y), 
             size = 2, stroke = .5, alpha=0.2) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_IPF_color)) +
  scale_fill_manual(values=Palette_IPF_fill) + 
  labs(title="Stable vs. Progressive IPF", 
       subtitle = glue("Wilcoxon ranked sum p={round(p_value,2)}"),
       x = "", 
       y = "Species number") +
  theme(
    plot.title = element_text(size=20),
    plot.subtitle = element_text(size=18),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1, size=14),
    axis.title.x = element_text(size=14),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1)

fig

###### Shannon Diversity number Status ####
shannon_diversity <- vegan::diversity(abund_IPF, index="shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(meta_IPF, shannon_diversity)
wilcox = wilcox.test(shannon_diversity ~ Progression_status_1y, shannon_metadata) #p=0.5
p_value = wilcox$p.value

fig = ggplot(shannon_metadata, aes(x = Progression_status_1y, y = shannon_diversity),
             fill = Progression_status_1y, shape = Progression_status_1y) +
  geom_violin(aes(color=Progression_status_1y))+
  geom_boxplot(aes(shape = Progression_status_1y, colour = Progression_status_1y), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Progression_status_1y, shape = Progression_status_1y, colour = Progression_status_1y), 
             size = 2, stroke = .5, alpha=0.2) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_IPF_color)) +
  scale_fill_manual(values=Palette_IPF_fill) +
  labs(title="Stable vs. Progressive IPF", 
       subtitle = glue("Wilcoxon ranked sum p={round(p_value,2)}"),
       x = "",
       y = "Shannon diversity") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    legend.position = "none", panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=0, size=14),
    axis.title.x = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    aspect.ratio = 1)

fig

Alpha_diversity_df = cbind(shannon_metadata, species_number)
Alpha_diversity_df = Alpha_diversity_df %>%
  select(`sample-id`, PatientID, Diagnosis, Age_at_recruitment, Sex, Reflux, 
         Progression_status_1y, Reflux_treatment, species_number, shannon_diversity) %>%
  rename("Observed" = "species_number") %>%
  rename("Shannon" = "shannon_diversity")
Alpha_diversity_df = melt(Alpha_diversity_df, id.vars = c("sample-id",
                                                          "PatientID",
                                                          "Diagnosis",
                                                          "Age_at_recruitment", "Sex",
                                                          "Progression_status_1y",
                                                          "Reflux", "Reflux_treatment"),
                          variable.name = "Diversity index")

Alpha_plot = ggplot(Alpha_diversity_df, aes(x = Progression_status_1y, y = value),
                    fill = Progression_status_1y, shape = Progression_status_1y) +
  geom_boxplot(aes(shape = Progression_status_1y, colour = Progression_status_1y), width=0.2, outliers = F) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(fill = Progression_status_1y, shape = Progression_status_1y, colour = Progression_status_1y), 
             size = 3, stroke = .5, alpha=0.2) +
  geom_violin(alpha=0.1, aes(color=Progression_status_1y)) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_IPF_color)) +
  scale_fill_manual(values=c(Palette_IPF_fill)) +
  labs(title="Alpha diversity indices by 1y Progression status for IPF patients",
       subtitle="Progression: FVC decline of 10% or death within a year",
       x = "",
       y = "Diversity") + theme_classic() + 
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18, face = "italic"),
    legend.position = "none", panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "grey"),
    strip.text = element_text(size=12, face = "bold"),
    strip.background = element_rect(linewidth=0.8, colour = "grey")) + 
  facet_wrap(.~`Diversity index`, scales="free_y") 
Alpha_plot

##### Top ten genera #####
abund_table_taxa <- abund_table 
top <- abund_table_taxa[,order(colSums(abund_table_taxa),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
df_genus_top <- cbind(meta_table, top)

wilcox_results <- list()

for (tax in taxa_list) {
  form <- as.formula(paste0("`", tax, "` ~ Diagnosis"))
  
  sub_df <- df_genus_top %>% 
    dplyr::filter(!is.na(Diagnosis) & Diagnosis %in% c("IPF", "Healthy"))
  
  if (length(unique(sub_df$Diagnosis)) == 2) {
    
    test <- wilcox.test(form, data = sub_df, exact = FALSE)
    
    stats <- sub_df %>%
      dplyr::group_by(Diagnosis) %>%
      dplyr::summarise(
        median_abund = median(.data[[tax]], na.rm = TRUE),
        mean_abund   = mean(.data[[tax]], na.rm = TRUE),
        .groups      = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from  = Diagnosis, 
        values_from = c(median_abund, mean_abund)
      )
    
    wilcox_results[[tax]] <- data.frame(
      Taxon       = tax,
      W_statistic = unname(test$statistic),
      p_value     = test$p.value,
      stats,
      stringsAsFactors = FALSE
    )
  }
}

wilcox_summary <- dplyr::bind_rows(wilcox_results) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_value)

df_long <- reshape2::melt(df_genus_top, id.vars = c("Diagnosis", "Progression_status_1y"), 
                          measure.vars = c(taxa_list), 
                          variable.name = "Genus",
                          factorsAsStrings = TRUE, na.rm = TRUE)

df_long_summarised <- df_long %>%
  group_by(Diagnosis, Genus) %>%
  dplyr::summarise(mean_value=mean(value),
                   sd=sd(value),
                   median=median(value),
                   q1 = quantile(value, 0.25),  # 1st quartile
                   q3 = quantile(value, 0.75),
                   total = sum(value)) # 3rd quartile
df_long_summarised$Diagnosis <- factor(df_long_summarised$Diagnosis, levels = c("Healthy", "IPF"))

# Streptococcus p<0.001
dat_text <- data.frame(
  label=c("", "***"),
  Diagnosis=c("Healthy", "IPF"),
  x = c(10,10), 
  y = c(10,4.5)
)

fig <- ggplot(df_long_summarised, aes(x=Genus, y=median)) + 
  geom_bar(aes(y = median, x = Genus, fill = Genus),
           stat="identity", alpha=0.8) +
  geom_errorbar(aes(x=Genus, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5) 

fig <- fig + 
  scale_fill_manual(values = c(Palette_10)) +
  theme_classic() + facet_wrap(.~Diagnosis)+
  geom_text(data=dat_text, mapping= aes(x=x, y=y, label=label), size = 10) +
  labs(x = "", y = "Median relative abundance (%)",
       title = "Top ten most abundant genera: Healthy vs. IPF",
       caption = "Wilcoxon rank-sum test with \nBenjamini–Hochberg adjustment, ***p-adj < 0.001") + 
  theme(axis.text.x = element_blank(),
        legend.text = element_text(face="italic", size=16), 
        legend.title = element_text(size = 18),
        strip.text = element_text(size=16),
        plot.margin = unit(c(10,0,10,0), "pt"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
        axis.title = element_text(size=18), axis.text = element_text(size=16),
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"),    
        strip.background = element_rect(linewidth=0.8, colour = "grey"),
        plot.caption = element_text(size = 14)) +
  guides(fill = guide_legend(title = "Genus"))

fig

###### Strep abundance with smoking status ######
Strep_smoking = df_genus_top %>% select(`sample-id`,Streptococcus) %>%
  dplyr::left_join(meta_table[c("sample-id","Smoking_history","Diagnosis")]) %>%
  mutate(Smoking_history = ifelse(Smoking_history == "Current smoker", "Ex-smoker", Smoking_history)) %>% 
  drop_na(Smoking_history)
p_vals_df <- Strep_smoking %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(
    p_val = wilcox.test(Streptococcus ~ Smoking_history)$p.value) %>%
  mutate(label = glue("{Diagnosis} p={round(p_val, 3)}"))
plot_subtitle <- paste(p_vals_df$label, collapse = " | ")

ggplot(Strep_smoking, aes(x = Smoking_history, y = Streptococcus, color = Smoking_history)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) + # Jitter helps see density better than geom_point
  scale_color_manual(values = c("#6F99AD", "#E18727")) + 
  facet_wrap(~Diagnosis, scales = "free_y") + # scales="free_y" helps if abundance varies wildly
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.subtitle = element_text(size = 10, face = "italic", lineheight = 1.1)) +
  labs(x = "Smoking status", 
       y = "Relative abundance of Streptococcus (%)",
       subtitle = glue("Wilcoxon rank-sum test:\n{plot_subtitle}"))

#### MicrobiomeStat Volcano plot #####
#install.packages("rlang")
# If this mstat doesn't work restart
library(MicrobiomeStat)
library(microbiome)
library(tidyverse)
library(qiime2R)
library(phyloseq)
# remove failed samples
min.depth = 1000
sort(sample_sums(physeq))[1:10]
ps.filt = prune_samples(sample_sums(ps) >= min.depth, ps)

table(tax_table(ps.filt)[, "Genus"], useNA = "ifany")
ps.filt <- subset_taxa(ps.filt, !is.na(Genus))

metadata <- read_tsv("data/BaselineStool_metadata.txt", trim_ws = T)
metadata <- metadata %>% distinct() %>%
  filter(!duplicated(`sample-id`)) # just need to drop the duplicated sample ID
metadata <- column_to_rownames(metadata, var="sample-id")
metadata <- meta_table %>% 
  filter(Diagnosis == "IPF" & Visit ==1 | Diagnosis =="Healthy" & Visit==1) %>%
  column_to_rownames("sample-id") %>%
  filter(!is.na(Reflux_treatment))

sample_data(ps.filt)<-metadata

taxa_are_rows(ps.filt)
ps_genus = ps.filt
data.obj <- mStat_convert_phyloseq_to_data_obj(ps_genus)
str(data.obj) # assess structure

data.obj$meta.dat

### Adjust the factors
data.obj$meta.dat$PatientID <- factor(data.obj$meta.dat$PatientID)
data.obj$meta.dat$Visit <- factor(data.obj$meta.dat$Visit)
data.obj$meta.dat$Diagnosis <- factor(data.obj$meta.dat$Diagnosis, levels = c("Healthy",
                                                                              "IPF"))
data.obj$meta.dat$Reflux_treatment <- factor(data.obj$meta.dat$Reflux_treatment, levels = c("None", "Treatment"))
data.obj$meta.dat$Sex <- factor(data.obj$meta.dat$Sex, levels = c("Male", "Female"))
data.obj$meta.dat$Smoking_history <- factor(data.obj$meta.dat$Smoking_history, levels = c("Never", "Ex-smoker"))

## Create the DAA list
test.list <- generate_taxa_test_single(
  data.obj = data.obj,
  group.var = "Diagnosis",
  adj.vars = c("Reflux_treatment", "Age_at_recruitment", "Sex", "Smoking_history"),
  feature.dat.type = "count",
  feature.level = c("Genus"),
  prev.filter = 0.20, # 20% of samples
  abund.filter = 0.01 # minimum 10%
)

## Create the plot
fig <- generate_taxa_volcano_single(
  data.obj = data.obj,
  group.var = "Diagnosis",
  test.list = test.list,
  feature.sig.level = 0.05,
  feature.mt.method = "fdr",
  palette = Palette_diverging1,
  pdf = F)

fig = fig$Genus$`IPF vs Healthy (Reference)` +
  labs(title = "IPF vs Healthy (reference)",
       subtitle = "Prevalence of 20% and abundance of 10%",
       caption = "FDR: p-adjust < 0.05") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=14, hjust=0.5),
        panel.grid.minor = element_blank()
  )
        

fig

#### Reflux treatment #####
abund_streptococcus = abund_table %>% dplyr::select(Streptococcus)
rownames(abund_streptococcus) = rownames(abund_table)
abund_streptococcus = rownames_to_column(abund_streptococcus, "sample-id")

df_glm_strp = dplyr::left_join(meta_table, abund_streptococcus)
df_glm_strp = df_glm_strp %>%
  mutate(Streptococcus = Streptococcus/100,
         Diagnosis = factor(Diagnosis, levels =c("Healthy", "IPF")),
         Reflux_treatment = factor(Reflux_treatment, levels=c("None", "Treatment")))

# need to get the total read count to include weights
count_data = rownames_to_column(df_library, "sample-id") %>% dplyr::select(`sample-id`, LibrarySize)
df_glm_strp = dplyr::left_join(df_glm_strp, count_data)

lm_strep_reflux= glm(Streptococcus ~ Reflux_treatment*Diagnosis, data=df_glm_strp,
                     family = quasibinomial(link="logit"), weights = LibrarySize)
# microbiome data is overdispersed.
# greater variance with count due to complexity of communities. 
# We aggregated ASVs into genus, but individual ASVs can be different. 

# Quasibinomial accounts for some flexibility in the model. 
summary(lm_strep_reflux)
df_glm_strp_rel = df_glm_strp %>% mutate(Streptococcus = Streptococcus*100)
#par(mfrow = c(2, 2))
#plot(lm_strep_reflux, 1, main = "Residuals vs. Fitted")
#plot(lm_strep_reflux, 2, main = "Normal Q-Q Plot")
#plot(lm_strep_reflux, 3, main = "Scale-Location Plot")
#plot(lm_strep_reflux, 4, main = "Residuals vs. Leverage")
#par(mfrow = c(1, 1))
df_glm_strp_stats <- df_glm_strp_rel %>%
  group_by(Reflux_treatment, Diagnosis) %>%
  dplyr::summarise(median_value=median(Streptococcus),
                   sd=sd(Streptococcus),
                   q1 = quantile(Streptococcus, 0.25),  # 1st quartile
                   q3 = quantile(Streptococcus, 0.75)) # 3rd quartile
library(ggtext)
fig = ggplot(df_glm_strp_rel[!is.na(df_glm_strp_rel$Reflux_treatment),], aes(x=Reflux_treatment, y =Streptococcus, color=Reflux_treatment)) +
  geom_boxplot(aes(colour = Reflux_treatment), outliers = T) +
  geom_point(position="identity", aes(shape=Reflux_treatment, colour=Reflux_treatment)) +
  #stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.1, color = "black", size = 0.5) +
  facet_grid(.~Diagnosis) + theme_classic()+
  scale_shape_manual(values = c(21,21)) + scale_color_manual(values=c(Palette_fill)) + 
  labs(title = "*Streptococcus* relative abundance by diagnosis and reflux treatment status",
       subtitle = "*Streptococcus* ~ Reflux_treatment × Diagnosis",
       #subtitle = "GLM model: beta= 0.72, 0.89, -0.05, SE= 0.35, 0.27, 0.39, p= 0.04, 0.001, 0.89, 
       #for Reflux_treatment:Treatment, Diagnosis:IPF and RefluxTreatment;DiagnosisIPF respectively.",       
       x = "Reflux treatment", y = "*Streptococcus* abundance (relative abundance %)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.subtitle = element_markdown(size = 14, face = "italic"),
        plot.title = element_markdown(size=20),
        axis.title = element_text(size=16),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face="bold", size=16),
        axis.text = element_text(size = 14),
        axis.title.y = element_markdown(size = 14))
  
fig

library(marginaleffects)
pred_df <- predictions(
  lm_strep_reflux,
  newdata = datagrid(
    Diagnosis        = c("Healthy", "IPF"),
    Reflux_treatment = c("None", "Treatment"))) %>%
  as.data.frame() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~. * 100))

# Plot
ggplot(pred_df, aes(x = Diagnosis, y = estimate, 
                    colour = Reflux_treatment, 
                    group  = Reflux_treatment)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width    = 0.15,
                position = position_dodge(width = 0.3),
                linewidth = 0.8) +
  geom_line(position = position_dodge(width = 0.3),
            linewidth = 0.8, linetype = "dashed") +
  scale_colour_manual(
    values = c("None" = "#2166AC", "Treatment" = "#D6604D"),
    name   = "Reflux Treatment"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 16),
    plot.margin      = unit(c(10, 10, 10, 10), "pt"),
    plot.title       = element_text(size = 20),
    plot.subtitle    = element_text(size = 12, face = "italic"),
    axis.title       = element_text(size = 18),
    axis.text        = element_text(size = 16),
    panel.border     = element_rect(color = "grey", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "grey"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = expression(italic("Streptococcus") ~ "relative abundance by diagnosis and reflux treatment"),
    subtitle = "Predicted probabilities from quasibinomial GLM with interaction term",
    x        = "",
    y        = "Predicted relative abundance (%)",
    caption  = "Error bars represent 95% confidence intervals"
  )

##### Streptococcus and Age #####
breaks = c(0,40,50,60,70,80,90,100)
labels = c("<40", "<50", "<60", "<70", "<80", "90+", "90+")
lm_strep_age= glm(Streptococcus ~ Age_at_recruitment, data=df_glm_strp,
                     family = quasibinomial(link="logit"), weights = LibrarySize)

# Quasibinomial accounts for some flexibility in the model. 
summary(lm_strep_age)
df_strep_age = df_glm_strp_rel %>%
  filter(!is.na(Age_at_recruitment)) %>%
  mutate(AgeBin = cut(Age_at_recruitment, breaks = breaks, 
                      labels = labels, right = F, include.lowest = T)) %>%
  group_by(AgeBin, Diagnosis) %>%
  summarise(MeanAbundance = mean(Streptococcus,))

coefs <- coef(summary(lm_strep_age))
row   <- coefs["Age_at_recruitment", ]   # change to your actual term name

subtitle_txt <- sprintf(
  "GLM model: beta = %.3f, SE = %.3f, p-value = %s",
  row["Estimate"],
  row["Std. Error"],
  format.pval(row[ncol(coefs)], digits = 2, eps = 0.001)
)

fig = ggplot(df_strep_age, aes(x = AgeBin, y=MeanAbundance, fill = Diagnosis)) + 
  geom_col(position = position_dodge(preserve = "single")) +
  theme_classic() + scale_fill_manual(values = Palette_fill) +
  scale_color_manual(values= Palette_color) +
  labs(title = "Mean *Streptococcus* abundance with Age",
       subtitle = subtitle_txt,
       x = "Age bins", y = "Mean *Streptococcus* abundance (relative abundance %)") + 
  theme(plot.title = element_markdown(size = 20),
        plot.subtitle = element_text(size=14, face = "italic"),
        axis.title = element_text(size=16),
        axis.title.y = element_markdown(size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=16))
  
fig
