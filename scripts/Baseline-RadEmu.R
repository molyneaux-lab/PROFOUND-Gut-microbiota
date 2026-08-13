setwd("~/Documents/GitHub/PROFOUND-Baseline/")

library(radEmu)
library(magrittr)
library(dplyr)
library(ggplot2)
library(stringr)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(reshape2)

metadata = read_csv("data/16SSequenced_complete-metadata.csv")
manifest_metadata <- read_tsv("data/BaselineStool_metadata.txt", trim_ws = TRUE)
metadata <- left_join(manifest_metadata, metadata)
metadata = metadata %>%
  mutate(Progression_status_1y = ifelse(Progression_status_1y == "Semi-progressive", "Stable", Progression_status_1y),
         Smoking_history = ifelse(Smoking_history == "Current smoker", "Ex-smoker", Smoking_history))

#### Loading phyloseq object ####
physeq <- qza_to_phyloseq(
  features  = "data/BaselineStool_filtered-table-gg2.qza",
  tree      = "data/BaselineStool_midrooted-tree.qza",
  taxonomy  = "data/BaselineStool_filtered-taxonomy-gg2.qza"
)

mapping <- metadata %>% column_to_rownames("sample-id")
sample_data(physeq) <- mapping

# Sanity check sample overlap
gplots::venn(list(mapping = rownames(mapping), physeq = sample_names(physeq)))
setdiff(sample_names(physeq), rownames(mapping))

physeq <- prune_samples(colSums(otu_table(physeq)) > 0, physeq)

#### Restrict taxa ####
# Prevalence: present in >0 counts in at least 20% of samples
prevalence_threshold <- 0.20
n_samples <- nsamples(physeq)

prev_filter <- genefilter_sample(
  physeq,
  filterfun_sample(function(x) x > 0),
  A = prevalence_threshold * n_samples
)

# Relative abundance: mean relative abundance >= 0.1% across samples
rel_abund <- transform_sample_counts(physeq, function(x) x / sum(x))
mean_abund <- rowMeans(otu_table(rel_abund))
abund_filter <- mean_abund >= 0.001  

# Apply both filters
keep <- prev_filter & abund_filter
physeq_filtered <- prune_taxa(keep, physeq)
ntaxa(physeq_filtered)

#### Remove zero-sum taxa and samples ####
physeq_filtered = tax_glom(physeq_filtered, taxrank = "Genus")
taxa_names(physeq_filtered) <- tax_table(physeq_filtered)[, "Genus"]

# Check before cleaning
sum(rowSums(otu_table(physeq_filtered)) == 0)  # taxa with zero sum
sum(colSums(otu_table(physeq_filtered)) == 0)  # samples with zero sum

# Remove zero-sum taxa, then zero-sum samples
physeq_clean <- prune_taxa(
  rowSums(otu_table(physeq_filtered)) > 0, 
  physeq_filtered
)
physeq_clean <- prune_samples(
  colSums(otu_table(physeq_clean)) > 0, 
  physeq_clean
)

# Confirm clean
sum(rowSums(otu_table(physeq_clean)) == 0)  # should be 0
sum(colSums(otu_table(physeq_clean)) == 0)  # should be 0
ntaxa(physeq_clean)
nsamples(physeq_clean)

#### Set factor levels ####
sample_data(physeq_clean)$Diagnosis <- factor(
  sample_data(physeq_clean)$Diagnosis,
  levels = c("Healthy", "IPF")
)

sample_data(physeq_clean)$Progression_status_1y <- factor(
  sample_data(physeq_clean)$Progression_status_1y,
  levels = c("Healthy", "Stable", "Progression")
)

meta_clean <- as.data.frame(sample_data(physeq_clean))
colSums(is.na(meta_clean[, c("Diagnosis", "Reflux_treatment", "Sex", "Smoking_history", "Age_at_recruitment")]))
complete_samples <- complete.cases(meta_clean[, c("Diagnosis", "Reflux_treatment", "Sex", "Smoking_history", "Age_at_recruitment")])
physeq_clean <- prune_samples(complete_samples, physeq_clean)
nsamples(physeq_clean)  # confirm

## Design matrix ##
my_formula = ~ Diagnosis + Reflux_treatment + Sex + Age_at_recruitment + 
  Smoking_history
X <- radEmu::make_design_matrix(physeq_clean,
                                formula = my_formula)
print(X)
colnames(X)
n_taxa <- ntaxa(physeq_clean)

# Only diagnosis
## Do absolute abundances between IPF and Healthy significantly differ?
test_kj <- data.frame(
  k = c(2),  # k=2 for Diagnosis, k=3 for Reflux_treatment etc.
  j = 1:n_taxa
)

#### Pass 1 without score tests ####
ch_fit <- emuFit(formula = my_formula,
                 Y = physeq_clean,
                 data = as.data.frame(sample_data(physeq_clean)),
                 run_score_tests = F)


# Check suitability
ch_fit$estimation_converged # if TRUE == good
any(is.na(ch_fit$B)) # no NA present == good
any(is.infinite(ch_fit$B)) # no NA present == good

ch_fit$coef %>% # check for large deviations between groups, could suggest noise
  filter(abs(estimate) > 10) %>%
  arrange(desc(abs(estimate)))

plot(ch_fit)$plots # check is the spread is ok, around 0 

ch_fit$coef %>% # the mean estimates should be around 0
  # on average most taxa are similar, therefore should be 0. 
  group_by(covariate) %>%
  dplyr::summarize(mean_estimate = mean(estimate),
            sd_estimate = sd(estimate))

taxonomy <- as.data.frame(tax_table(physeq_clean)) %>%
  rownames_to_column("category")  # ASV codes become a column to match with ch_fit$coef

# Join taxonomy onto your coefficient table
coef_table <- ch_fit$coef %>%
  left_join(taxonomy, by = "category")

# Now you can view results with readable names
coef_table %>%
  filter(covariate == "DiagnosisIPF") %>%
  arrange(desc(abs(estimate))) %>%
  dplyr::select(covariate, Genus, Species, estimate, lower, upper) %>%
  head(20)

taxa_names <- taxonomy %>%
  mutate(cat_small = case_when(
    !is.na(Genus) ~ Genus,
    !is.na(Family) ~ paste("Unknown", Family),
    !is.na(Order) ~ paste("Unknown", Order),
    TRUE ~ category
  )) %>%
  dplyr::select(category, cat_small)

plot(ch_fit, taxon_names = taxa_names)$plots

#### Pass 2 - score tests for Diagnosis across all taxa ####
# Rebuild test_kj with correct n_taxa
test_kj <- data.frame(
  k = 2,
  j = 1:n_taxa
)
cores_to_use <- max(1, parallel::detectCores() - 2)

# This will take some time
ch_fit_tested <- emuFit(
  formula = my_formula,
  Y = physeq_clean,
  data = as.data.frame(sample_data(physeq_clean)),
  B = ch_fit$B,                # Pre-fit parameters speed up optimization
  run_score_tests = TRUE,
  test_kj = test_kj,
  n_cores = cores_to_use       # Parallelize score tests across CPU cores
)

# Join taxonomy onto tested results
coef_tested <- ch_fit_tested$coef %>%
  dplyr::left_join(taxonomy, by = "category") %>%
  mutate(label = case_when(
    !is.na(Genus) ~ Genus,
    !is.na(Family) ~ paste("Unknown", Family),
    !is.na(Order) ~ paste("Unknown", Order),
    TRUE ~ category
  ))

# Summary of p-values for Diagnosis
coef_tested %>%
  filter(covariate == "DiagnosisIPF") %>%
  dplyr::summarize(
    n_taxa = n(),
    n_nominal_sig = sum(pval < 0.05, na.rm = TRUE),      # before correction
    n_FDR_sig = sum(p.adjust(pval, method = "BH") < 0.05, na.rm = TRUE)  # after correction
  )

# Add FDR correction to the full table
coef_tested <- coef_tested %>%
  group_by(covariate) %>%
  mutate(FDR = p.adjust(pval, method = "holm")) %>%
  ungroup()

coef_tested %>%
  filter(covariate == "DiagnosisIPF") %>%
  mutate(significant = FDR < 0.2,
         label = reorder(label, estimate)) %>%
  ggplot(aes(x = estimate, y = label, colour = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "red3"),
                      labels = c("Not significant", "FDR < 0.20")) +
  labs(x = "Log fold-change (IPF vs Healthy)",
       y = NULL,
       colour = NULL,
       title = "Differential abundance: IPF vs Healthy",
       subtitle = "Adjusted for Reflux treatment, Age, Smoking history and Sex") +
  theme_bw() +
  theme(legend.position = "bottom")

coef_tested %>%
  filter(covariate == "DiagnosisIPF") %>%
  mutate(significant = FDR < 0.2) %>%
  ggplot(aes(x = estimate, y = -log10(pval), colour = significant)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text(aes(label = ifelse(significant, label, "")),
            hjust = -0.1, size = 3) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#BC3C29FF")) +
  labs(x = "Log fold-change (IPF vs Healthy)",
       y = "-log10(p-value)",
       colour = "FDR < 0.05",
       title = "Volcano plot: IPF vs Healthy") +
  theme_bw() +
  theme(legend.position = "bottom")

view(coef_tested %>%
  filter(covariate == "DiagnosisIPF") %>%
         #FDR < 0.05) %>%
  arrange(FDR) %>%
  dplyr::select(Genus, Species, estimate, lower, upper, pval, FDR) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))))

#### Decreases in IPF but Streptococcus ASVs so increase just not significant. 
## Robust consistent signal in IPF patients, higher Strep
## Could be a 'predisposition' that IPF patients have??
## Decrease in FB allows Strep to 'bloom'?
sd_df <- data.frame(sample_data(physeq_clean), check.names = FALSE)

library_sizes <- data.frame(
  Sample = sample_names(physeq_clean),
  library_size = sample_sums(physeq_clean)) %>%
  left_join(sd_df %>%
      rownames_to_column("Sample") %>%
      dplyr::select(Sample, Diagnosis),
    by = "Sample"
  )

library_sizes %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(
    mean_depth = mean(library_size),
    median_depth = median(library_size),
    sd_depth = sd(library_size)
  )
 
wilcox.test(library_size ~ Diagnosis, library_sizes)

ggplot(library_sizes, aes(x = Diagnosis, y = library_size, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = c("Healthy" = "#BC3C29FF", "IPF" = "#0072B5FF")) +
  labs(x = NULL, y = "Library size (log scale)",
       title = "Sequencing depth by diagnosis") +
  theme_bw() +
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

physeq_ra <- transform_sample_counts(physeq_clean, function(x) x / sum(x))

# Extract both genera
plot_data <- psmelt(physeq_ra) %>%
  filter(Genus %in% c("Streptococcus", "Coprococcus_A_121497")) %>%
  group_by(Sample, Diagnosis, Genus) %>%
  dplyr::summarize(abundance = sum(Abundance), .groups = "drop")

ggplot(plot_data, aes(x = Diagnosis, y = abundance, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  scale_fill_manual(values = c("Healthy" = "#0072B5FF", "IPF" = "#BC3C29FF")) +
  labs(x = NULL, y = "Relative abundance (log scale)",
       title = "Streptococcus vs Coprococcus by diagnosis") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        legend.position = "none",
        strip.text = element_text(size=16, face = "italic"),
        plot.margin = unit(c(10,0,10,0), "pt"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
        axis.title = element_text(size=18), axis.text = element_text(size=16),
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"),    
        strip.background = element_rect(linewidth=0.8, colour = "grey"))

# Relative abundance correlation - to show compositional relationship
cor_ra <- psmelt(physeq_ra) %>%
  filter(Genus %in% c("Streptococcus", "Coprococcus_A_121497")) %>%
  group_by(Sample, Genus, Diagnosis) %>%
  dplyr::summarize(abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = abundance)

cor_ra_results <- cor_ra %>%
  group_by(Diagnosis) %>%
  summarise(test = list(tidy(cor.test(Coprococcus_A_121497, Streptococcus, method = "spearman")))) %>%
  tidyr::unnest(test)

print(cor_ra_results)

# Plot both side by side for comparison
library(ggtext)
plot_ra <- ggplot(cor_ra[cor_ra$Diagnosis=="IPF",], 
                  aes(x = Coprococcus_A_121497, y = Streptococcus,
                                                        colour = Diagnosis)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("IPF" = "#0072B5FF")) +
  labs(title = "Relative abundance of *Streptococcus* and *Coprococcus_A_121497*",
       subtitle = "Spearman's rank: rho= -0.02, p-val= 0.745",
       x = "*Coprococcus_A_121497* (log(relative abundance (%)))",
       y = "*Streptococcus* (log(relative abundance (%)))") +
  theme_classic() +
  theme(legend.position = "none", 
        strip.text = element_text(size=16),
        plot.title = ggtext::element_markdown(size=20), 
        plot.subtitle = element_text(size=14, face = "italic"), 
        axis.title.x = ggtext::element_markdown(size=18), 
        axis.title.y = ggtext::element_markdown(size=18),
        axis.text = element_text(size=16),
        axis.text.x = element_text(size = 16), 
        panel.border = element_rect(color = "grey", fill = NA, linewidth=0.8),
        axis.line = element_line(color="grey"))

plot_ra
