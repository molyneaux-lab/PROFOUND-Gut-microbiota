## Streptococcus ASVs only ##
library(phyloseq)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(qiime2R)
library(microbiome)
library(janitor)
library(coin)
library(stringr)
library(patchwork)
library(survminer); library(survival)
set.seed(2026)
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

setwd("~/Documents/GitHub/PROFOUND-Baseline/")
ps <- qza_to_phyloseq(
  features = "data/BaselineStool_filtered-table-gg2.qza",
  taxonomy = "data/BaselineStool_filtered-taxonomy-gg2.qza",
  tree     = "data/BaselineStool_midrooted-tree.qza"
)

metadata = read_tsv("data/BaselineStool_metadata.txt",
                    comment = "#q2:types", show_col_types = FALSE) 
more_metadata = read_csv("data/16SSequenced_complete-metadata.csv",
                          show_col_types = FALSE)

dup <- setdiff(intersect(names(metadata), names(more_metadata)), "PatientID")
n0  <- nrow(metadata)
metadata <- dplyr::left_join(metadata, more_metadata, by = "PatientID", suffix = c("", ".dup"))
if (nrow(metadata) != n0) stop("join changed row count: key not unique")
if (length(dup)) {
  for (v in dup) {
    ne <- sum(!is.na(metadata[[v]]) & !is.na(metadata[[paste0(v, ".dup")]]) &
                as.character(metadata[[v]]) != as.character(metadata[[paste0(v, ".dup")]]))
    if (ne > 0) warning("'", v, "' disagrees between metadata files in ", ne, " samples")
  }
  metadata <- dplyr::select(metadata, -dplyr::all_of(paste0(dup, ".dup")))
}
metadata <- column_to_rownames(metadata, "sample-id")

absent <- setdiff(c("Diagnosis", "Age_at_recruitment", "Sex", "Reflux_treatment"),
                  names(metadata))
if (length(absent)) stop("missing from metadata: ", paste(absent, collapse = ", "))
sample_data(ps) <- metadata

#### What ASVs are contibuting to the Streptococcus genus ####
## Pull out Streptococcus ASVs
# Clean up the Streptococcal names 
tt <- as.data.frame(unclass(tax_table(ps)), stringsAsFactors = FALSE)

fix <- which(!is.na(tt$Genus) & tt$Genus == "Streptococcus" &
               grepl("^sp\\. [0-9]+$", tt$Species))
tt$Species[fix] <- paste("Streptococcus", tt$Species[fix])

tax_table(ps) <- tax_table(as.matrix(tt))
strep_ids <- rownames(tt)[!is.na(tt$Genus) & tt$Genus == "Streptococcus"]
strep_species    <- prune_taxa(strep_ids, ps)
strep_species_df <- as.data.frame(unclass(tax_table(strep_species))) %>%
  dplyr::select(Species) %>%
  rownames_to_column("ASV")

manual_names <- c(
  "RS-GCF-900104225.1-NZ-FNJK01000021.1"  = "Streptococcus mobilis",
  "GB-GCA-900767615.1-CAAFZU010000193.1"  = "Streptococcus salivarius #1",
  "MJ006-1-barcode37-umi151381bins-ubs-3" = "Streptococcus salivarius #2",
  "MJ006-2-barcode73-umi159901bins-ubs-3" = "Streptococcus salivarius #3",
  "RS-GCF-000339435.1-NZ-AHSQ01000040.1"  = "Streptococcus mutans #1",
  "RS-GCF-002157665.1-NZ-BDOS01000001.1"  = "Streptococcus mutans #2"
)

bad <- setdiff(names(manual_names), strep_species_df$ASV)
if (length(bad)) stop("manual_names IDs not in Streptococcus ASVs: ",
                      paste(bad, collapse = ", "))

strep_species_df <- strep_species_df %>%
  mutate(Species = dplyr::coalesce(manual_names[ASV], Species))

tt[strep_species_df$ASV, "Species"] <- strep_species_df$Species
tax_table(ps) <- tax_table(as.matrix(tt))

depth    <- sample_sums(ps)
sort(depth)[1:15] 
sum(depth < 5000)

# remove failed reads
ps_keep  <- prune_samples(depth >= 5000, ps)
ps_strep <- prune_taxa(strep_ids, ps_keep)

otu <- as(otu_table(ps_strep), "matrix")
if (taxa_are_rows(ps_strep)) otu <- t(otu)

## Streptococcus specific ## 
# Relative abundances of Streptococcus ASVs as a proportion of the entire Streptococcus genus
otu <- as(otu_table(ps_strep), "matrix")
if (!taxa_are_rows(ps_strep)) otu <- t(otu)

number_dupes <- function(x, sep = " #") {
  ave(x, x, FUN = function(g)
    if (length(g) > 1) paste0(g, sep, seq_along(g)) else g)
}

strep_df_counts <- as.data.frame(otu) %>%
  tibble::rownames_to_column("ASV") %>%
  mutate(Taxon = number_dupes(
    ifelse(is.na(tt[ASV, "Species"]),
           "Streptococcus unassigned",
           tt[ASV, "Species"])),
    .after = ASV)

#### LinDA ####
# Which Streptococcal ASV is more abundant in IPF vs. Healthy
cnt <- otu                                   # otu <- as(otu_table(ps_strep), "matrix"), Strep only.
cnt <- cnt[, colSums(cnt) >= 50]             # drop samples below 50 total reads of strep. 

md <- metadata[colnames(cnt), ]
md <- md[!is.na(md$Diagnosis) & md$Diagnosis %in% c("IPF", "Healthy"), ]
md$Diagnosis <- factor(md$Diagnosis, levels = c("Healthy", "IPF"))
cnt <- cnt[, rownames(md)]

md$Smoking_history <- factor(ifelse(as.character(md$Smoking_history) == "Current",
                                    "Ex-smoker", as.character(md$Smoking_history)))

library(MicrobiomeStat)
# If this doesn't work restart session and load MicrobiomeStat first then other packages
res_strep <- linda(
  feature.dat = as.data.frame(cnt),
  meta.dat    = md,
  formula     = "~ Diagnosis + Age_at_recruitment + Sex + Smoking_history + Reflux_treatment",
  feature.dat.type = "count",
  prev.filter = 0.1,
  is.winsor   = TRUE, outlier.pct = 0.03)

lin <- res_strep$output$DiagnosisIPF
lin$ASV <- rownames(lin)
lin <- lin %>%
  dplyr::left_join(strep_df_counts %>% dplyr::select(ASV, Taxon), by = "ASV") %>%
  arrange(pvalue)

fp <- lin %>%                                   # LinDA output, already joined to Taxon
  mutate(
    label = sub("^Streptococcus ", "S. ", Taxon),
    label = sub("_[A-Z](_[0-9]+)?$", "", label),
    lo    = log2FoldChange - 1.96 * lfcSE,
    hi    = log2FoldChange + 1.96 * lfcSE,
    sig   = padj < 0.05
  ) %>%
  arrange(log2FoldChange)

ggplot(fp, aes(x = log2FoldChange, y = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2,
                 colour = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = sig), shape = 21, size = 5,
             colour = "grey20", stroke = 0.4) +
  scale_fill_manual(values = c(`FALSE` = "gray", `TRUE` = "#BC3C29"),
                    guide = "none") +
  labs(x = expression(log[2]~"fold change (IPF vs healthy)"), y = NULL,
       caption = "Within-genus composition") +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=18),
        plot.caption = element_text(size=14)
  )

#### Streptococcus ASV contribution to the genus % of the cohort ####
# Normalise reads 
strep_rel = strep_df_counts %>%
  mutate(across(-c(ASV, Taxon), ~ if (sum(.x) > 0) .x / sum(.x) else 0)) %>%
  column_to_rownames("Taxon") %>% dplyr::select(-ASV)
strep_rel <- data.frame(t(strep_rel), check.names = FALSE)

strep_long <- strep_rel %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(-SampleID, names_to = "Taxon", values_to = "RelAbund") %>%
  mutate(PctStrep = 100 * RelAbund)

strep_long %>% group_by(Taxon) %>%
  summarise(mean = mean(RelAbund),
            median = median(RelAbund))

strep_wide <- strep_long %>%
  dplyr::select(SampleID, Taxon, PctStrep) %>%
  tidyr::pivot_wider(names_from = Taxon, values_from = PctStrep) %>% janitor::clean_names(.)
strep_wide <- strep_wide %>%
  mutate(strep_depth = colSums(otu)[sample_id],
         PatientID = str_extract(sample_id, "^[^-]+")) %>%
  left_join(metadata) %>% dplyr::filter(Diagnosis == "IPF") %>%
  dplyr::mutate(Death = ifelse(Death == "Yes", 1, 0),
                Progression_status_1y = ifelse(Progression_status_1y == "Semi-progressive", "Stable", Progression_status_1y),
                PFS_event = ifelse(Death == 1 & FVC_10 == "Yes", 1, 0),
                Death_1y = ifelse(Death == 1 & Survival_death < 366, 1, 0),
                Smoking_history = ifelse(Smoking_history == "Current smoker", "Ex-smoker", Smoking_history),
                Sex = factor(Sex, levels = c("Male", "Female")))

pt <- strep_long %>%
  dplyr::filter(SampleID %in% rownames(strep_rel)) %>%
  dplyr::mutate(Diagnosis = dplyr::if_else(stringr::str_detect(SampleID, "PFND"), "IPF", "Healthy")) %>%
  dplyr::mutate(
    label = sub("^Streptococcus ", "S. ", Taxon),
    label = sub("_[A-Z](_[0-9]+)?$", "", label)
  ) %>%
  dplyr::group_by(label, Diagnosis) %>%
  dplyr::summarise(
    pct_strep = mean(PctStrep),
    prev      = mean(PctStrep > 0), 
    .groups   = "drop"
  ) %>%
  dplyr::filter(pct_strep > 0) %>%
  dplyr::arrange(pct_strep) %>%
  dplyr::mutate(label = factor(label, levels = unique(label)))

lv <- levels(pt$label)

strep_wilcox_df <- strep_long %>%
  # Match metadata using SampleID to prevent misclassification
  dplyr::left_join(
    metadata %>% tibble::rownames_to_column("SampleID"), 
    by = "SampleID"
  ) %>%
  dplyr::filter(Diagnosis %in% c("IPF", "Healthy")) %>%
  dplyr::mutate(
    Diagnosis = factor(Diagnosis, levels = c("Healthy", "IPF")),
    label = sub("^Streptococcus ", "S. ", Taxon),
    label = sub("_[A-Z](_[0-9]+)?$", "", label)
  )

taxa_list <- unique(strep_wilcox_df$label)

wilcox_results <- list()

for (tax in taxa_list) {
  sub_df <- strep_wilcox_df %>% dplyr::filter(label == tax)
  group_counts <- table(sub_df$Diagnosis[sub_df$PctStrep > 0])
  
  if (length(unique(sub_df$Diagnosis)) == 2) {
    
    test <- wilcox.test(PctStrep ~ Diagnosis, data = sub_df, exact = FALSE)
    
    stats <- sub_df %>%
      dplyr::group_by(Diagnosis) %>%
      dplyr::summarise(
        median_pct = median(PctStrep, na.rm = TRUE),
        mean_pct   = mean(PctStrep, na.rm = TRUE),
        prev       = mean(PctStrep > 0, na.rm = TRUE),
        .groups    = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = Diagnosis, 
        values_from = c(median_pct, mean_pct, prev)
      )
    
    wilcox_results[[tax]] <- data.frame(
      Taxon = tax,
      W_statistic = unname(test$statistic),
      p_value = test$p.value,
      stats,
      stringsAsFactors = FALSE
    )
  }
}

wilcox_summary <- dplyr::bind_rows(wilcox_results) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_value)
sig_taxa <- c("S. parasanguinis", "S. mutans #2", "S. mobilis")
pt <- pt %>%
  mutate(is_sig = label %in% sig_taxa)

ggplot(pt, aes(x = pct_strep, y = label, fill = is_sig)) + 
  facet_wrap(~ Diagnosis) +
  geom_col(width = 0.65) +
  scale_fill_manual(
    values = c(`FALSE` = "#20854E", `TRUE` = "#BC3C29"), 
    guide = "none" # hides the fill legend
  ) +
  scale_y_discrete(limits = lv, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = expression(paste("% of gut ", italic("Streptococcus"))), 
    y = NULL,
    caption = "Within-genus composition \nIn red are taxa where p-adj<0.05"
  ) + 
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "italic"),
    axis.text          = element_text(size = 14),
    axis.title         = element_text(size = 16),
    plot.title         = element_text(size = 20),
    plot.subtitle      = element_text(size = 18),
    plot.caption       = element_text(size = 14),
    strip.text         = element_text(face = "bold", size = 16)
  )

#### Survival analysis ####
# Mortality
strep_wide_ipf = strep_wide 
coxph(Surv(Survival_death, Death) ~ streptococcus_parasanguinis, data = strep_wide)

# Composite: Death or a FVC decline of 10%, death trumps FVC decline
coxph(Surv(PFS_FVC_death, PFS_event) ~ streptococcus_parasanguinis + Age_at_recruitment + 
        Sex + Smoking_history + Baseline_ppFVC, data = strep_wide)

strep_res= surv_cutpoint(strep_wide,
                         time = "PFS_FVC_death", 
                         event = "PFS_event", minprop = 0.1, 
                         variable = c("streptococcus_parasanguinis"))

summary(strep_res)
plot(strep_res)

strep_cat = surv_categorize(strep_res)
head(strep_cat)

plot(strep_res, c("streptococcus_parasanguinis"), palette = "npg")

fit1 = survfit(Surv(PFS_FVC_death, PFS_event) ~ streptococcus_parasanguinis, data = strep_cat)
names(fit1$strata)

fig = ggsurvplot(fit1, data = strep_cat, risk.table = TRUE, conf.int = FALSE,
                 pval = T,
                 palette = c("#BC3C29FF", "#0072B5FF"),
                 xlab = "Time (years)",
                 ylab = "PFS probability",
                 risk.table.height = 0.3,
                 xlim = c(0,1096), xscale=365.25, break.time.by = 365.25,
                 legend.labs = c("High", "Low"),
                 legend.title = expression(paste(italic("S. parasanguinis")," group")))

fig_plot = fig$plot
fig_table = fig$table

fig = (fig_plot) / (fig_table) + plot_layout(heights = c(4,1))
fig
