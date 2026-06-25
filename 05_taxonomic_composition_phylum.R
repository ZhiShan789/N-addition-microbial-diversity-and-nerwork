# Dominant phylum composition analyses
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 4
# Generated on: 2026-06-25


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/summarize_plot_phylum_relative_abundance.R
# ================================================================

source("R/utils.R")
source("R/00_config.R")

need_package("ggplot2")
need_package("dplyr")
need_package("tidyr")
library(ggplot2)
library(dplyr)
library(tidyr)

cfg <- load_config()
out_dir <- file.path(cfg$output_dir, "05_significance_summary")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

clean_phylum <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^[a-z]__", "", x)
  x[x == "" | is.na(x) | x %in% c("NA", "na", "Unassigned", "unassigned", "unknown", "Unknown")] <- "Unassigned"
  dplyr::recode(
    x,
    "Pseudomonadota" = "Proteobacteria",
    "Actinomycetota" = "Actinobacteria",
    "Bacillota" = "Firmicutes",
    "Bacteroidota" = "Bacteroidetes",
    "Acidobacteriota" = "Acidobacteria",
    "Planctomycetota" = "Planctomycetes",
    "Verrucomicrobiota" = "Verrucomicrobia",
    "Chloroflexota" = "Chloroflexi",
    "Gemmatimonadota" = "Gemmatimonadetes",
    "Cyanobacteriota" = "Cyanobacteria",
    .default = x
  )
}

read_taxonomy <- function(domain_row) {
  tax_path <- file.path(domain_row$ResultPath, "taxonomy.txt")
  tax <- safe_read_tsv(tax_path)
  if (is.null(tax) || nrow(tax) == 0 || !("Phylum" %in% names(tax))) return(NULL)
  names(tax)[1] <- "OTUID"
  tax$OTUID <- trimws(as.character(tax$OTUID))
  tax$PhylumRaw <- trimws(as.character(tax$Phylum))
  tax$Phylum <- clean_phylum(tax$PhylumRaw)
  tax[, c("OTUID", "PhylumRaw", "Phylum"), drop = FALSE]
}

domain_phylum_abundance <- function(domain_row) {
  md <- read_domain_metadata(domain_row)
  otu <- read_otutab_rare(domain_row)
  tax <- read_taxonomy(domain_row)
  if (is.null(md) || is.null(otu) || is.null(tax)) return(NULL)
  if (!("Group" %in% names(md))) return(NULL)

  md$GroupRole <- ifelse(
    is_control_group(md$Group, cfg), "CK",
    ifelse(is_treatment_group(md$Group, cfg), "N addition", NA_character_)
  )
  md <- md[!is.na(md$GroupRole), , drop = FALSE]
  if (nrow(md) == 0) return(NULL)

  sample_ids <- intersect(colnames(otu), md$SampleID)
  if (length(sample_ids) < 2) return(NULL)

  otu <- otu[, sample_ids, drop = FALSE]
  otu$OTUID <- rownames(otu)
  otu_tax <- merge(otu, tax, by = "OTUID", all.x = TRUE, sort = FALSE)
  otu_tax$Phylum[is.na(otu_tax$Phylum)] <- "Unassigned"
  otu_tax$PhylumRaw[is.na(otu_tax$PhylumRaw)] <- "Unassigned"

  sample_totals <- colSums(otu_tax[, sample_ids, drop = FALSE], na.rm = TRUE)
  sample_totals[sample_totals <= 0 | !is.finite(sample_totals)] <- NA_real_

  phylum_counts <- stats::aggregate(
    otu_tax[, sample_ids, drop = FALSE],
    by = list(Phylum = otu_tax$Phylum),
    FUN = sum,
    na.rm = TRUE
  )
  phylum_counts <- phylum_counts[phylum_counts$Phylum != "", , drop = FALSE]
  count_mat <- as.matrix(phylum_counts[, sample_ids, drop = FALSE])
  rel_mat <- sweep(count_mat, 2, sample_totals, "/") * 100
  rel_mat[!is.finite(rel_mat)] <- 0

  rel_df <- data.frame(Phylum = phylum_counts$Phylum, rel_mat, check.names = FALSE)
  long <- tidyr::pivot_longer(
    rel_df,
    cols = dplyr::all_of(sample_ids),
    names_to = "SampleID",
    values_to = "RelativeAbundance"
  )
  md_keep <- md[, c("SampleID", "Group", "GroupRole"), drop = FALSE]
  long <- merge(long, md_keep, by = "SampleID", all.x = TRUE, sort = FALSE)
  long$StudyID <- domain_row$StudyID
  long$Domain <- domain_row$Domain
  long$DomainPath <- domain_row$DomainPath
  long
}

domain_rows <- scan_domain_dirs(cfg)
domain_rows <- domain_rows[domain_rows$HasOtuRare & domain_rows$HasTaxonomy & domain_rows$HasMetadata, , drop = FALSE]

message(sprintf("Scanning %d domains with rarefied OTU tables and taxonomy...", nrow(domain_rows)))
sample_rows <- lapply(seq_len(nrow(domain_rows)), function(i) {
  if (i %% 25 == 0) message(sprintf("  processed %d / %d", i, nrow(domain_rows)))
  domain_phylum_abundance(domain_rows[i, ])
})
sample_level <- dplyr::bind_rows(sample_rows)
if (nrow(sample_level) == 0) {
  stop("No phylum abundance data were generated.", call. = FALSE)
}

sample_level <- sample_level %>%
  dplyr::select(StudyID, Domain, DomainPath, SampleID, Group, GroupRole, Phylum, RelativeAbundance)

domain_level <- sample_level %>%
  dplyr::group_by(StudyID, Domain, DomainPath, GroupRole, Phylum) %>%
  dplyr::summarise(
    SampleN = dplyr::n_distinct(SampleID),
    DomainMeanRelativeAbundance = mean(RelativeAbundance, na.rm = TRUE),
    .groups = "drop"
  )

domain_group_keys <- domain_level %>%
  dplyr::distinct(StudyID, Domain, DomainPath, GroupRole)
phylum_keys <- domain_level %>%
  dplyr::distinct(Domain, Phylum)
domain_level_complete <- merge(domain_group_keys, phylum_keys, by = "Domain", all = TRUE) %>%
  dplyr::left_join(
    domain_level,
    by = c("StudyID", "Domain", "DomainPath", "GroupRole", "Phylum")
  ) %>%
  dplyr::mutate(
    SampleN = dplyr::coalesce(SampleN, 0L),
    DomainMeanRelativeAbundance = dplyr::coalesce(DomainMeanRelativeAbundance, 0)
  )

group_summary <- domain_level_complete %>%
  dplyr::group_by(Domain, GroupRole, Phylum) %>%
  dplyr::summarise(
    DomainN = dplyr::n(),
    SampleN = sum(SampleN, na.rm = TRUE),
    MeanRelativeAbundance = mean(DomainMeanRelativeAbundance, na.rm = TRUE),
    MedianRelativeAbundance = stats::median(DomainMeanRelativeAbundance, na.rm = TRUE),
    SERelativeAbundance = stats::sd(DomainMeanRelativeAbundance, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

top_phyla <- group_summary %>%
  dplyr::filter(Phylum != "Unassigned") %>%
  dplyr::group_by(Domain, Phylum) %>%
  dplyr::summarise(MeanOverall = mean(MeanRelativeAbundance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Domain) %>%
  dplyr::slice_max(MeanOverall, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

plot_domain_level <- domain_level_complete %>%
  dplyr::left_join(top_phyla %>% dplyr::mutate(IsTop = TRUE), by = c("Domain", "Phylum")) %>%
  dplyr::mutate(PhylumPlot = ifelse(is.na(IsTop), "Other/Unassigned", Phylum)) %>%
  dplyr::group_by(StudyID, Domain, DomainPath, GroupRole, PhylumPlot) %>%
  dplyr::summarise(
    SampleN = sum(SampleN, na.rm = TRUE),
    DomainMeanRelativeAbundance = sum(DomainMeanRelativeAbundance, na.rm = TRUE),
    .groups = "drop"
  )

plot_summary <- plot_domain_level %>%
  dplyr::group_by(Domain, GroupRole, Phylum = PhylumPlot) %>%
  dplyr::summarise(
    DomainN = dplyr::n(),
    SampleN = sum(SampleN, na.rm = TRUE),
    MeanRelativeAbundance = mean(DomainMeanRelativeAbundance, na.rm = TRUE),
    SERelativeAbundance = stats::sd(DomainMeanRelativeAbundance, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

paired_top10 <- plot_summary %>%
  dplyr::filter(Phylum != "Other/Unassigned") %>%
  dplyr::select(Domain, GroupRole, Phylum, MeanRelativeAbundance) %>%
  tidyr::pivot_wider(names_from = GroupRole, values_from = MeanRelativeAbundance) %>%
  dplyr::mutate(
    CK = dplyr::coalesce(CK, 0),
    `N addition` = dplyr::coalesce(`N addition`, 0),
    Delta_N_minus_CK = `N addition` - CK,
    lnRR_N_vs_CK = log((`N addition` + 1e-6) / (CK + 1e-6)),
    MeanOverall = (`N addition` + CK) / 2
  ) %>%
  dplyr::arrange(Domain, dplyr::desc(MeanOverall))

write_csv_utf8(sample_level, file.path(out_dir, "phylum_relative_abundance_sample_level.csv"))
write_csv_utf8(domain_level, file.path(out_dir, "phylum_relative_abundance_domain_level.csv"))
write_csv_utf8(group_summary, file.path(out_dir, "phylum_relative_abundance_group_summary.csv"))
write_csv_utf8(plot_summary, file.path(out_dir, "phylum_relative_abundance_top10_plot_summary.csv"))
write_csv_utf8(paired_top10, file.path(out_dir, "phylum_relative_abundance_top10_ck_n.csv"))

palette_values <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#B07AA1",
  "#EDC948", "#9C755F", "#FF9DA7", "#86BCB6", "#BAB0AC", "#6B7280",
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#F1CE63", "#D37295", "#499894",
  "#B6992D", "#D4A6C8", "#7F7F7F", "#C7C7C7"
)

plot_summary$GroupRole <- factor(plot_summary$GroupRole, levels = c("CK", "N addition"))
plot_summary$Phylum <- factor(plot_summary$Phylum)

stacked <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(x = GroupRole, y = MeanRelativeAbundance, fill = Phylum)
) +
  ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~Domain, nrow = 1) +
  ggplot2::scale_fill_manual(values = rep(palette_values, length.out = length(unique(plot_summary$Phylum)))) +
  ggplot2::labs(x = NULL, y = "Mean relative abundance (%)", fill = "Phylum") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", color = "grey75"),
    legend.position = "right",
    legend.key.size = grid::unit(0.42, "cm")
  )

ggplot2::ggsave(file.path(fig_dir, "Figure_phylum_relative_abundance_top10_stacked.png"), stacked, width = 10.5, height = 5.6, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "Figure_phylum_relative_abundance_top10_stacked.pdf"), stacked, width = 10.5, height = 5.6)

paired_plot_df <- plot_summary %>%
  dplyr::filter(Phylum != "Other/Unassigned") %>%
  dplyr::left_join(
    paired_top10 %>% dplyr::select(Domain, Phylum, MeanOverall),
    by = c("Domain", "Phylum")
  ) %>%
  dplyr::mutate(PhylumLabel = reorder(paste(Domain, Phylum, sep = "__"), MeanOverall))

label_map <- setNames(sub("^.*__", "", levels(paired_plot_df$PhylumLabel)), levels(paired_plot_df$PhylumLabel))

dot <- ggplot2::ggplot(
  paired_plot_df,
  ggplot2::aes(x = MeanRelativeAbundance, y = PhylumLabel, color = GroupRole)
) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = pmax(0, MeanRelativeAbundance - SERelativeAbundance), xmax = MeanRelativeAbundance + SERelativeAbundance),
    height = 0.18,
    alpha = 0.55,
    linewidth = 0.45,
    position = ggplot2::position_dodge(width = 0.55)
  ) +
  ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::facet_wrap(~Domain, scales = "free_y", nrow = 1) +
  ggplot2::scale_y_discrete(labels = label_map) +
  ggplot2::scale_color_manual(values = c("CK" = "#F2B43F", "N addition" = "#4E9BD3")) +
  ggplot2::labs(x = "Mean relative abundance (%)", y = NULL, color = NULL) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", color = "grey75"),
    legend.position = "top"
  )

ggplot2::ggsave(file.path(fig_dir, "Figure_phylum_relative_abundance_top10_ck_vs_n.png"), dot, width = 9.5, height = 5.8, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "Figure_phylum_relative_abundance_top10_ck_vs_n.pdf"), dot, width = 9.5, height = 5.8)

message("Wrote phylum relative abundance tables and figures:")
message("  ", file.path(out_dir, "phylum_relative_abundance_top10_ck_n.csv"))
message("  ", file.path(fig_dir, "Figure_phylum_relative_abundance_top10_stacked.png"))
message("  ", file.path(fig_dir, "Figure_phylum_relative_abundance_top10_ck_vs_n.png"))


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/plot_phylum_relative_abundance_forest.R
# ================================================================

source("R/utils.R")
source("R/00_config.R")

need_package("ggplot2")
need_package("dplyr")
need_package("tidyr")
library(ggplot2)
library(dplyr)
library(tidyr)

cfg <- load_config()
summary_dir <- file.path(cfg$output_dir, "05_significance_summary")
fig_dir <- file.path(summary_dir, "figures")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sample_path <- file.path(summary_dir, "phylum_relative_abundance_sample_level.csv")
top_path <- file.path(summary_dir, "phylum_relative_abundance_top10_ck_n.csv")
if (!file.exists(sample_path) || !file.exists(top_path)) {
  stop("Run scripts/summarize_plot_phylum_relative_abundance.R first.", call. = FALSE)
}

numify <- function(x) suppressWarnings(as.numeric(as.character(x)))

random_effect_dl <- function(yi, vi) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]
  vi <- vi[keep]
  k <- length(yi)
  if (k == 0) {
    return(data.frame(k = 0, estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
                      ci_high = NA_real_, p_value = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  }
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (k > 1 && c_val > 0) max(0, (q - (k - 1)) / c_val) else 0
  w_re <- 1 / (vi + tau2)
  est <- sum(w_re * yi) / sum(w_re)
  se <- sqrt(1 / sum(w_re))
  ci <- est + c(-1, 1) * 1.96 * se
  z <- est / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  i2 <- if (k > 1 && q > 0) max(0, (q - (k - 1)) / q) * 100 else 0
  data.frame(k = k, estimate = est, se = se, ci_low = ci[1], ci_high = ci[2],
             p_value = p, tau2 = tau2, I2 = i2)
}

effect_for_domain_phylum <- function(d) {
  control <- d$RelativeAbundance[d$GroupRole == "CK"]
  treat <- d$RelativeAbundance[d$GroupRole == "N addition"]
  e <- lnrr_effect(treat, control, min_n = cfg$min_group_replicates)
  data.frame(
    StudyID = unique(d$StudyID)[1],
    Domain = unique(d$Domain)[1],
    DomainPath = unique(d$DomainPath)[1],
    Phylum = unique(d$Phylum)[1],
    e,
    stringsAsFactors = FALSE
  )
}

sample_level <- read.csv(sample_path, check.names = FALSE, stringsAsFactors = FALSE)
top10 <- read.csv(top_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_level$RelativeAbundance <- numify(sample_level$RelativeAbundance)

top_keys <- top10 %>%
  select(Domain, Phylum, MeanOverall) %>%
  distinct()

forest_input <- sample_level %>%
  inner_join(top_keys, by = c("Domain", "Phylum")) %>%
  filter(GroupRole %in% c("CK", "N addition"))

split_key <- paste(forest_input$DomainPath, forest_input$Phylum, sep = "||")
effects <- bind_rows(lapply(split(forest_input, split_key), effect_for_domain_phylum)) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

summary <- bind_rows(lapply(split(effects, list(effects$Domain, effects$Phylum), drop = TRUE), function(d) {
  re <- random_effect_dl(d$yi, d$vi)
  data.frame(
    Domain = unique(d$Domain)[1],
    Phylum = unique(d$Phylum)[1],
    re,
    stringsAsFactors = FALSE
  )
})) %>%
  left_join(top_keys, by = c("Domain", "Phylum")) %>%
  filter(k > 0, is.finite(estimate)) %>%
  group_by(Domain) %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH"),
    SigClass = case_when(
      FDR_BH < 0.05 & estimate > 0 ~ "Significant increase",
      FDR_BH < 0.05 & estimate < 0 ~ "Significant decrease",
      TRUE ~ "Not significant"
    ),
    Label = ifelse(
      is.finite(FDR_BH),
      paste0("k=", k, "; FDR=", format.pval(FDR_BH, digits = 2, eps = 0.001)),
      paste0("k=", k)
    )
  ) %>%
  ungroup()

write_csv_utf8(effects, file.path(summary_dir, "phylum_relative_abundance_lnRR_effect_sizes.csv"))
write_csv_utf8(summary, file.path(summary_dir, "phylum_relative_abundance_lnRR_forest_summary.csv"))

summary <- summary %>%
  arrange(Domain, desc(MeanOverall)) %>%
  mutate(PhylumPanel = factor(paste(Domain, Phylum, sep = "__"), levels = rev(paste(Domain, Phylum, sep = "__"))))

label_map <- setNames(sub("^.*__", "", levels(summary$PhylumPanel)), levels(summary$PhylumPanel))
cols <- c(
  "Significant decrease" = "#0072B2",
  "Not significant" = "grey55",
  "Significant increase" = "#D55E00"
)

x_min <- min(summary$ci_low, na.rm = TRUE)
x_max <- max(summary$ci_high, na.rm = TRUE)
x_pad <- max(0.12, (x_max - x_min) * 0.22)

p <- ggplot(summary, aes(x = estimate, y = PhylumPanel, color = SigClass)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey35") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.16, linewidth = 0.48) +
  geom_point(aes(size = k), alpha = 0.96) +
  geom_text(
    aes(label = Label),
    x = x_max + x_pad * 0.18,
    hjust = 0,
    size = 2.45,
    color = "black"
  ) +
  facet_wrap(~Domain, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = label_map) +
  scale_color_manual(
    values = cols,
    labels = c(
      "Significant decrease" = "Decrease",
      "Not significant" = "NS",
      "Significant increase" = "Increase"
    ),
    drop = FALSE
  ) +
  scale_size_continuous(range = c(1.8, 4.5)) +
  scale_x_continuous(limits = c(x_min - x_pad, x_max + x_pad), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Major phylum responses to nitrogen addition",
    subtitle = "Random-effects meta-analysis; effect size is lnRR of relative abundance (N addition / CK)",
    x = "lnRR of relative abundance",
    y = NULL,
    color = NULL,
    size = "Study-domain effects"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    panel.spacing.y = grid::unit(0.55, "cm")
  ) +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

png_path <- file.path(fig_dir, "Figure_phylum_relative_abundance_lnRR_forest.png")
pdf_path <- file.path(fig_dir, "Figure_phylum_relative_abundance_lnRR_forest.pdf")
ggsave(png_path, p, width = 7.6, height = 7.2, dpi = 300)
ggsave(pdf_path, p, width = 7.6, height = 7.2)

message("Wrote phylum lnRR effects: ", file.path(summary_dir, "phylum_relative_abundance_lnRR_effect_sizes.csv"))
message("Wrote phylum forest summary: ", file.path(summary_dir, "phylum_relative_abundance_lnRR_forest_summary.csv"))
message("Wrote phylum forest figure: ", png_path)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/plot_phylum_composition_effect_sizes.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

pipeline_dir <- Sys.getenv(
  "PIPELINE_DIR",
  unset = "E:/BaiduSyncdisk/N_deposition1/Data/Global_Nitrogen_Pipeline"
)

summary_dir <- file.path(pipeline_dir, "output", "05_significance_summary")
fig_dir <- file.path(summary_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(summary_dir, "phylum_composition_random_effect_lnRR_summary.csv")
if (!file.exists(input_path)) {
  stop("Missing phylum composition summary: ", input_path, call. = FALSE)
}

x <- read.csv(input_path, check.names = FALSE)
num_cols <- c("MeanControl", "MeanTreatment", "k", "estimate", "se", "ci_low", "ci_high", "p_value", "FDR_BH", "I2")
for (col in intersect(num_cols, names(x))) {
  x[[col]] <- suppressWarnings(as.numeric(x[[col]]))
}

x <- x %>%
  mutate(
    Domain = factor(Domain, levels = c("Bacterial", "Fungi")),
    MeanRelativeAbundance = rowMeans(cbind(MeanControl, MeanTreatment), na.rm = TRUE),
    SigClass = case_when(
      FDR_BH < 0.05 & estimate > 0 ~ "Significant increase",
      FDR_BH < 0.05 & estimate < 0 ~ "Significant decrease",
      TRUE ~ "Not significant"
    ),
    SigClass = factor(SigClass, levels = c("Significant decrease", "Not significant", "Significant increase")),
    Label = ifelse(FDR_BH < 0.001, "FDR < 0.001", paste0("FDR = ", sprintf("%.3f", FDR_BH)))
  )

# Keep common phyla for the overview; always keep significant taxa.
overview <- x %>%
  filter(k >= 20 | FDR_BH < 0.05) %>%
  group_by(Domain) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(TaxonPlot = factor(Taxon, levels = unique(Taxon))) %>%
  ungroup()

sig_only <- x %>%
  filter(FDR_BH < 0.05) %>%
  group_by(Domain) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(TaxonPlot = factor(Taxon, levels = unique(Taxon))) %>%
  ungroup()

top10 <- x %>%
  group_by(Domain) %>%
  slice_max(order_by = MeanRelativeAbundance, n = 10, with_ties = FALSE) %>%
  arrange(Domain, estimate) %>%
  mutate(
    TaxonPlot = factor(Taxon, levels = unique(Taxon)),
    SigLabel = ifelse(
      FDR_BH < 0.05,
      ifelse(FDR_BH < 0.001, "FDR < 0.001", paste0("FDR = ", sprintf("%.3f", FDR_BH))),
      ""
    )
  ) %>%
  ungroup()

cols <- c(
  "Significant decrease" = "#0072B2",
  "Not significant" = "grey65",
  "Significant increase" = "#D55E00"
)

base_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey25")
  )

save_plot <- function(plot, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png(filename = png_path, width = width, height = height, units = "in", res = 450)
  print(plot)
  dev.off()
  pdf(file = pdf_path, width = width, height = height, family = "sans", useDingbats = FALSE)
  print(plot)
  dev.off()
  c(png = png_path, pdf = pdf_path)
}

p_overview <- ggplot(overview, aes(x = estimate, y = TaxonPlot, color = SigClass)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, linewidth = 0.35) +
  geom_point(aes(size = k), alpha = 0.92) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_size_continuous(range = c(1.3, 3.6), breaks = c(20, 60, 100, 140)) +
  labs(
    title = "Phylum-level composition responses to nitrogen addition",
    subtitle = "Effect size is ln response ratio of relative abundance percentage (N addition / CK)",
    x = "lnRR of relative abundance",
    y = NULL,
    size = "Effects"
  ) +
  base_theme

p_sig <- ggplot(sig_only, aes(x = estimate, y = TaxonPlot, color = SigClass)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.17, linewidth = 0.42) +
  geom_point(aes(size = k), alpha = 0.96) +
  geom_text(aes(label = Label), hjust = ifelse(sig_only$estimate >= 0, -0.05, 1.05), size = 2.3, color = "black") +
  facet_wrap(~ Domain, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_size_continuous(range = c(1.6, 3.8), breaks = c(20, 60, 100, 140)) +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.28))) +
  labs(
    title = "Significant phylum-level composition responses",
    subtitle = "Only phyla with BH-FDR < 0.05 are shown",
    x = "lnRR of relative abundance",
    y = NULL,
    size = "Effects"
  ) +
  base_theme

p_top10 <- ggplot(top10, aes(x = estimate, y = TaxonPlot, color = SigClass)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.17, linewidth = 0.42) +
  geom_point(aes(size = MeanRelativeAbundance), alpha = 0.96) +
  geom_text(
    aes(label = SigLabel),
    hjust = ifelse(top10$estimate >= 0, -0.05, 1.05),
    size = 2.2,
    color = "black"
  ) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 1) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_size_continuous(range = c(1.8, 4.2)) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.30))) +
  labs(
    title = "Top 10 abundant phyla: composition effect sizes",
    subtitle = "Top phyla are ranked by mean relative abundance across CK and N-addition samples",
    x = "lnRR of relative abundance",
    y = NULL,
    size = "Mean relative abundance"
  ) +
  base_theme

paths_overview <- save_plot(p_overview, "Figure_phylum_composition_effect_sizes_overview", 7.0, 7.2)
paths_sig <- save_plot(p_sig, "Figure_phylum_composition_effect_sizes_significant_only", 7.0, 4.8)
paths_top10 <- save_plot(p_top10, "Figure_phylum_composition_effect_sizes_top10_abundant", 7.2, 5.8)

sig_table <- sig_only %>%
  arrange(Domain, estimate) %>%
  select(Domain, Taxon, k, MeanControl, MeanTreatment, estimate, ci_low, ci_high, p_value, FDR_BH, SigClass)
write.csv(sig_table, file.path(summary_dir, "phylum_composition_significant_effects_only.csv"), row.names = FALSE, fileEncoding = "UTF-8")

top10_table <- top10 %>%
  arrange(Domain, desc(MeanRelativeAbundance)) %>%
  select(Domain, Taxon, k, MeanControl, MeanTreatment, MeanRelativeAbundance, estimate, ci_low, ci_high, p_value, FDR_BH, SigClass)
write.csv(top10_table, file.path(summary_dir, "phylum_composition_top10_abundant_effects.csv"), row.names = FALSE, fileEncoding = "UTF-8")

index <- data.frame(
  Figure = c("Overview", "Significant only", "Top 10 abundant"),
  PNG = c(paths_overview[["png"]], paths_sig[["png"]], paths_top10[["png"]]),
  PDF = c(paths_overview[["pdf"]], paths_sig[["pdf"]], paths_top10[["pdf"]]),
  stringsAsFactors = FALSE
)
write.csv(index, file.path(fig_dir, "phylum_composition_effect_figure_index.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("Wrote overview: ", paths_overview[["png"]])
message("Wrote significant-only: ", paths_sig[["png"]])
message("Wrote top 10 abundant: ", paths_top10[["png"]])
message("Significant phyla: ", nrow(sig_table))


# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/analyze_phylum_composition_significance.R
# ================================================================

library(dplyr)
library(ggplot2)
library(metafor)
library(openxlsx)

setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

out_dir <- "outputs/plots/phylum_composition"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

effect_file <- "outputs/meta_effects_all.csv"

out_effect_rows <- file.path(out_dir, "phylum_meta_ready_effect_rows.csv")
out_sig_all <- file.path(out_dir, "phylum_lnRR_effect_significance_all.csv")
out_sig_top10 <- file.path(out_dir, "phylum_lnRR_effect_significance_top10.csv")
out_sig_all_xlsx <- file.path(out_dir, "phylum_lnRR_effect_significance_all.xlsx")
out_sig_top10_xlsx <- file.path(out_dir, "phylum_lnRR_effect_significance_top10.xlsx")
out_fig_png <- file.path(out_dir, "phylum_lnRR_effect_size_with_significance.png")
out_fig_pdf <- file.path(out_dir, "phylum_lnRR_effect_size_with_significance.pdf")

domain_order <- c("Bacterial", "Fungi")
effect_colors <- c("N lower" = "#F2B84B", "N higher" = "#4E9BD3", "No change" = "#9A9A9A")

non_phylum_metrics <- c(
  "richness", "shannon", "simpson",
  "RRBeta_ln_Dt_over_Dc",
  "RRStructure_ln_Db_over_meanDcDt"
)

p_text <- function(p) {
  ifelse(is.na(p), "P = NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
}

p_star <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}

clean_effect_rows <- function(effects) {
  if (!("RR" %in% names(effects))) effects$RR <- effects$yi
  if (!("Vi" %in% names(effects))) effects$Vi <- effects$vi
  if (!("StudyID" %in% names(effects))) effects$StudyID <- effects$ID

  effects %>%
    mutate(
      RR = suppressWarnings(as.numeric(RR)),
      Vi = suppressWarnings(as.numeric(Vi)),
      StudyID = as.character(StudyID)
    ) %>%
    filter(
      Domain %in% domain_order,
      !(Metric %in% non_phylum_metrics),
      !is.na(RR),
      !is.na(Vi),
      is.finite(RR),
      is.finite(Vi),
      Vi > 0
    )
}

fit_one_phylum <- function(dat) {
  if (nrow(dat) < 3 || length(unique(dat$StudyID)) < 2) {
    return(data.frame(
      Estimate = mean(dat$RR, na.rm = TRUE),
      CI_low = NA_real_,
      CI_high = NA_real_,
      P = NA_real_
    ))
  }

  fit <- tryCatch(
    suppressWarnings(
      rma.mv(
        yi = RR,
        V = Vi,
        random = ~ 1 | StudyID,
        method = "REML",
        data = dat
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      Estimate = mean(dat$RR, na.rm = TRUE),
      CI_low = NA_real_,
      CI_high = NA_real_,
      P = NA_real_
    ))
  }

  data.frame(
    Estimate = as.numeric(fit$beta[1]),
    CI_low = as.numeric(fit$ci.lb),
    CI_high = as.numeric(fit$ci.ub),
    P = as.numeric(fit$pval[1])
  )
}

summarise_phylum_effects <- function(effect_rows) {
  split_rows <- split(
    effect_rows,
    list(effect_rows$Domain, effect_rows$Metric),
    drop = TRUE
  )

  bind_rows(lapply(split_rows, function(dat) {
    model_row <- fit_one_phylum(dat)

    data.frame(
      Domain = unique(dat$Domain)[1],
      Phylum = unique(dat$Metric)[1],
      PairN = nrow(dat),
      ArticleN = length(unique(dat$StudyID)),
      Mean_lnRR = mean(dat$RR, na.rm = TRUE),
      SD_lnRR = sd(dat$RR, na.rm = TRUE),
      model_row,
      stringsAsFactors = FALSE
    )
  })) %>%
    group_by(Domain) %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%
    ungroup() %>%
    mutate(
      Direction = case_when(
        Estimate > 0 ~ "N higher",
        Estimate < 0 ~ "N lower",
        TRUE ~ "No change"
      ),
      FDR_label = paste0(
        p_text(FDR),
        ifelse(p_star(FDR) == "", "", paste0(" ", p_star(FDR)))
      ),
      Count_label = paste0(ArticleN, "(", PairN, ")")
    ) %>%
    arrange(Domain, desc(abs(Estimate)))
}

get_top_phyla <- function(effect_df, n = 10) {
  effect_df %>%
    group_by(Domain) %>%
    slice_max(order_by = PairN, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Domain, desc(PairN))
}

draw_effect_plot <- function(effect_df) {
  top_phyla <- get_top_phyla(effect_df, 10)
  plot_df <- effect_df %>%
    semi_join(top_phyla %>% select(Domain, Phylum), by = c("Domain", "Phylum"))

  x_effect_min <- min(plot_df$CI_low, plot_df$Estimate, na.rm = TRUE)
  x_effect_max <- max(plot_df$CI_high, plot_df$Estimate, na.rm = TRUE)
  effect_abs <- max(abs(c(x_effect_min, x_effect_max)), na.rm = TRUE)
  if (!is.finite(effect_abs) || effect_abs == 0) effect_abs <- 0.2
  effect_xlim <- c(-effect_abs * 1.15, effect_abs * 1.15)
  x_p <- effect_xlim[2] + effect_abs * 0.48
  x_count <- effect_xlim[2] + effect_abs * 1.40
  x_max <- effect_xlim[2] + effect_abs * 2.02

  draw_one_device <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(0.38, 4.2))
    par(family = "sans", mar = c(0, 0, 0, 0), oma = c(2.4, 0.2, 0.2, 0.2), xpd = FALSE)

    plot.new()
    legend(
      "center",
      legend = c("N lower", "N higher"),
      col = effect_colors[c("N lower", "N higher")],
      pch = 16,
      lty = 1,
      bty = "n",
      horiz = TRUE,
      cex = 1.08,
      x.intersp = 0.9
    )

    par(mar = c(4.5, 8.0, 2.4, 1.0))

    for (domain_i in domain_order) {
      phyla <- top_phyla %>%
        filter(Domain == domain_i) %>%
        arrange(desc(PairN)) %>%
        pull(Phylum)

      dat_domain <- plot_df %>%
        filter(Domain == domain_i, Phylum %in% phyla)

      y_pos <- rev(seq_along(phyla))
      names(y_pos) <- phyla
      n_y <- length(phyla)

      plot(
        NA,
        xlim = c(effect_xlim[1], x_max),
        ylim = c(0.4, n_y + 0.85),
        xlab = "",
        ylab = "",
        yaxt = "n",
        bty = "l",
        axes = FALSE
      )
      abline(v = pretty(effect_xlim, n = 5), col = "#E7E7E7", lwd = 1)
      abline(v = 0, col = "#555555", lty = 2, lwd = 1.1)
      axis(1, las = 1, cex.axis = 0.92, lwd = 1)
      axis(2, at = y_pos, labels = phyla, las = 1, cex.axis = 0.92, lwd = 1)
      box(col = "#222222", lwd = 1)

      usr <- par("usr")
      rect(usr[1], n_y + 0.45, usr[2], n_y + 0.85, col = "#EBEBEB", border = "#BDBDBD")
      text(mean(c(usr[1], x_p - effect_abs * 0.12)), n_y + 0.65, domain_i, cex = 1.1)
      text(x_p, n_y + 0.65, "FDR P", adj = 0, cex = 0.84)
      text(x_count, n_y + 0.65, "study(pair)", adj = 0, cex = 0.84)

      for (i in seq_len(nrow(dat_domain))) {
        y_i <- y_pos[[dat_domain$Phylum[i]]]
        col_i <- effect_colors[[dat_domain$Direction[i]]]

        ci_low <- dat_domain$CI_low[i]
        ci_high <- dat_domain$CI_high[i]
        if (is.na(ci_low) || is.na(ci_high)) {
          ci_low <- dat_domain$Estimate[i]
          ci_high <- dat_domain$Estimate[i]
        }

        segments(ci_low, y_i, ci_high, y_i, col = col_i, lwd = 1.8)
        segments(ci_low, y_i - 0.07, ci_low, y_i + 0.07, col = col_i, lwd = 1.4)
        segments(ci_high, y_i - 0.07, ci_high, y_i + 0.07, col = col_i, lwd = 1.4)
        points(dat_domain$Estimate[i], y_i, pch = 16, cex = 1.08, col = col_i)
        text(x_p, y_i, dat_domain$FDR_label[i], adj = 0, cex = 0.82, col = "#222222")
        text(x_count, y_i, dat_domain$Count_label[i], adj = 0, cex = 0.82, col = "#222222")
      }
    }

    mtext("Meta-analytic effect size of relative abundance (lnRR)",
      side = 1, outer = TRUE, line = 0.7, cex = 1.05
    )
  }

  png(out_fig_png, width = 11.3, height = 6.2, units = "in", res = 300)
  draw_one_device()
  dev.off()

  pdf(out_fig_pdf, width = 11.3, height = 6.2, useDingbats = FALSE)
  draw_one_device()
  dev.off()
}

main <- function() {
  effects <- read.csv(effect_file, check.names = FALSE, stringsAsFactors = FALSE)
  effect_rows <- clean_effect_rows(effects)
  effect_df <- summarise_phylum_effects(effect_rows)
  top_phyla <- get_top_phyla(effect_df, 10)
  effect_top <- effect_df %>%
    semi_join(top_phyla %>% select(Domain, Phylum), by = c("Domain", "Phylum"))

  write.csv(effect_rows, out_effect_rows, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(effect_df, out_sig_all, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(effect_top, out_sig_top10, row.names = FALSE, fileEncoding = "UTF-8")
  openxlsx::write.xlsx(effect_df, out_sig_all_xlsx, overwrite = TRUE)
  openxlsx::write.xlsx(effect_top, out_sig_top10_xlsx, overwrite = TRUE)
  draw_effect_plot(effect_df)

  message("Meta-ready phylum effect rows: ", out_effect_rows)
  message("All phylum effects: ", out_sig_all_xlsx)
  message("Top phylum effects: ", out_sig_top10_xlsx)
  message("Effect-size figure: ", out_fig_png)
}

main()

