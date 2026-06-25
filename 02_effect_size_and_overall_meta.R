# Effect-size tables and overall meta-analysis
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 3
# Generated on: 2026-06-25


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/summarize_alpha_beta_tax_significance.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(vegan)
  library(metafor)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "scripts/summarize_alpha_beta_tax_significance.R"
root_dir <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
setwd(root_dir)

source("R/utils.R")
source("R/00_config.R")
cfg <- load_config()

out_dir <- file.path(cfg$output_dir, "05_significance_summary")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

alpha_path <- file.path(cfg$output_dir, "03_effect_size", "alpha_lnrr_effect_sizes.csv")
beta_existing_path <- file.path(cfg$output_dir, "02_collected", "beta_domain_summary.csv")

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

fit_random_effect <- function(df) {
  df <- df[is.finite(df$yi) & is.finite(df$vi) & df$vi > 0, , drop = FALSE]
  if (nrow(df) < 3) {
    return(data.frame(
      k = nrow(df), estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, p_value = NA_real_, tau2 = NA_real_, I2 = NA_real_
    ))
  }
  fit <- tryCatch(metafor::rma.uni(yi = yi, vi = vi, data = df, method = "REML"), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(
      k = nrow(df), estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, p_value = NA_real_, tau2 = NA_real_, I2 = NA_real_
    ))
  }
  data.frame(
    k = fit$k,
    estimate = as.numeric(fit$b[1]),
    se = fit$se,
    ci_low = fit$ci.lb,
    ci_high = fit$ci.ub,
    p_value = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2
  )
}

write_csv <- function(x, name) {
  path <- file.path(out_dir, name)
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

# Alpha diversity: random-effects meta-analysis of lnRR.
alpha <- read.csv(alpha_path, check.names = FALSE)
alpha$yi <- suppressWarnings(as.numeric(alpha$yi))
alpha$vi <- suppressWarnings(as.numeric(alpha$vi))
alpha$Domain <- factor(alpha$Domain, levels = c("Bacterial", "Fungi"))

alpha_summary <- alpha %>%
  group_by(Domain, Metric) %>%
  group_modify(~ fit_random_effect(.x)) %>%
  ungroup() %>%
  group_by(Domain) %>%
  mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Direction = case_when(
      is.na(estimate) ~ "NA",
      FDR_BH < 0.05 & estimate > 0 ~ "Significant increase",
      FDR_BH < 0.05 & estimate < 0 ~ "Significant decrease",
      estimate > 0 ~ "Non-significant increase",
      estimate < 0 ~ "Non-significant decrease",
      TRUE ~ "No change"
    )
  )

write_csv(alpha_summary, "alpha_random_effect_lnRR_summary.csv")

# Beta diversity: recompute PERMANOVA/ANOSIM/MRPP from Bray-Curtis matrices using CK vs N role.
read_bray <- function(path) {
  x <- tryCatch(read.table(path, header = TRUE, sep = "\t", row.names = 1, quote = "", comment.char = "", check.names = FALSE), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- NA_real_
  m
}

beta_rows <- list()
domains <- scan_domain_dirs(cfg)
domains <- domains[domains$HasMetadata, , drop = FALSE]
for (i in seq_len(nrow(domains))) {
  dr <- domains[i, ]
  bray_path <- file.path(dr$ResultPath, "beta", "bray_curtis.txt")
  if (!file.exists(bray_path)) next
  md <- read_domain_metadata(dr)
  m <- read_bray(bray_path)
  if (is.null(md) || is.null(m)) next
  common <- intersect(rownames(m), md$SampleID)
  if (length(common) < 4) next
  md <- md[match(common, md$SampleID), , drop = FALSE]
  m <- m[common, common, drop = FALSE]
  role <- ifelse(is_control_group(md$Group, cfg), "Control",
                 ifelse(is_treatment_group(md$Group, cfg), "Treatment", NA_character_))
  keep <- !is.na(role)
  if (sum(keep) < 4) next
  md <- md[keep, , drop = FALSE]
  role <- role[keep]
  m <- m[keep, keep, drop = FALSE]
  n_c <- sum(role == "Control")
  n_t <- sum(role == "Treatment")
  if (n_c < 2 || n_t < 2) next

  dist_obj <- as.dist(m)
  meta <- data.frame(GroupRole = factor(role, levels = c("Control", "Treatment")))

  ad <- tryCatch(vegan::adonis2(dist_obj ~ GroupRole, data = meta, permutations = 999), error = function(e) NULL)
  an <- tryCatch(vegan::anosim(dist_obj, grouping = meta$GroupRole, permutations = 999), error = function(e) NULL)
  mr <- tryCatch(vegan::mrpp(as.matrix(dist_obj), grouping = meta$GroupRole, distance = "euclidean", permutations = 999), error = function(e) NULL)
  bd <- tryCatch({
    b <- vegan::betadisper(dist_obj, meta$GroupRole)
    p <- permutest(b, permutations = 999)
    p$tab$`Pr(>F)`[1]
  }, error = function(e) NA_real_)

  lower <- m
  lower[upper.tri(lower, diag = TRUE)] <- NA_real_
  within_c <- mean_or_na(as.numeric(m[role == "Control", role == "Control"][lower.tri(m[role == "Control", role == "Control"], diag = FALSE)]))
  within_t <- mean_or_na(as.numeric(m[role == "Treatment", role == "Treatment"][lower.tri(m[role == "Treatment", role == "Treatment"], diag = FALSE)]))
  between <- mean_or_na(as.numeric(m[role == "Control", role == "Treatment"]))

  beta_rows[[length(beta_rows) + 1]] <- data.frame(
    StudyID = dr$StudyID,
    Domain = dr$Domain,
    Samples = length(role),
    ControlSamples = n_c,
    TreatmentSamples = n_t,
    BrayWithinControl = within_c,
    BrayWithinTreatment = within_t,
    BrayBetweenControlTreatment = between,
    BetweenMinusMeanWithin = between - mean(c(within_c, within_t), na.rm = TRUE),
    PERMANOVA_R2 = if (!is.null(ad)) ad$R2[1] else NA_real_,
    PERMANOVA_P = if (!is.null(ad)) ad$`Pr(>F)`[1] else NA_real_,
    ANOSIM_R = if (!is.null(an)) unname(an$statistic) else NA_real_,
    ANOSIM_P = if (!is.null(an)) an$signif else NA_real_,
    MRPP_A = if (!is.null(mr)) mr$A else NA_real_,
    MRPP_P = if (!is.null(mr)) mr$Pvalue else NA_real_,
    BetaDispersion_P = bd,
    stringsAsFactors = FALSE
  )
}

beta_tests <- if (length(beta_rows) == 0) data.frame() else bind_rows(beta_rows)
if (nrow(beta_tests) > 0) {
  beta_tests <- beta_tests %>%
    group_by(Domain) %>%
    mutate(
      PERMANOVA_FDR_BH = p.adjust(PERMANOVA_P, method = "BH"),
      ANOSIM_FDR_BH = p.adjust(ANOSIM_P, method = "BH"),
      MRPP_FDR_BH = p.adjust(MRPP_P, method = "BH"),
      BetaDispersion_FDR_BH = p.adjust(BetaDispersion_P, method = "BH")
    ) %>%
    ungroup()
}
write_csv(beta_tests, "beta_recomputed_permanova_anosim_mrpp.csv")

beta_summary <- beta_tests %>%
  group_by(Domain) %>%
  summarise(
    DomainsTested = n(),
    PERMANOVA_P_lt_0.05 = sum(PERMANOVA_P < 0.05, na.rm = TRUE),
    PERMANOVA_FDR_lt_0.05 = sum(PERMANOVA_FDR_BH < 0.05, na.rm = TRUE),
    Median_PERMANOVA_R2 = median(PERMANOVA_R2, na.rm = TRUE),
    ANOSIM_P_lt_0.05 = sum(ANOSIM_P < 0.05, na.rm = TRUE),
    ANOSIM_FDR_lt_0.05 = sum(ANOSIM_FDR_BH < 0.05, na.rm = TRUE),
    MRPP_P_lt_0.05 = sum(MRPP_P < 0.05, na.rm = TRUE),
    MRPP_FDR_lt_0.05 = sum(MRPP_FDR_BH < 0.05, na.rm = TRUE),
    BetaDispersion_P_lt_0.05 = sum(BetaDispersion_P < 0.05, na.rm = TRUE),
    BetaDispersion_FDR_lt_0.05 = sum(BetaDispersion_FDR_BH < 0.05, na.rm = TRUE),
    Median_BetweenMinusMeanWithin = median(BetweenMinusMeanWithin, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(beta_summary, "beta_significance_summary.csv")

# Taxonomic composition: phylum-level CK vs N lnRR across domains.
composition_rows <- list()
tax_files <- list.files(cfg$base_data_dir, pattern = "sum_p\\.txt$", recursive = TRUE, full.names = TRUE)
tax_files <- tax_files[grepl("[/\\\\]result[/\\\\]tax[/\\\\]sum_p\\.txt$", tax_files)]

for (tax_path in tax_files) {
  domain_dir <- dirname(dirname(dirname(tax_path)))
  study_dir <- dirname(domain_dir)
  domain <- basename(domain_dir)
  study_id <- basename(study_dir)
  if (!(domain %in% cfg$valid_domains)) next
  md_path <- file.path(domain_dir, "result", "metadata.txt")
  raw_path <- file.path(domain_dir, "result", "metadata_raw.txt")
  md <- safe_read_tsv(if (file.exists(md_path)) md_path else raw_path)
  if (is.null(md) || !all(c("SampleID", "Group") %in% names(md))) next

  tax <- tryCatch(read.table(tax_path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE), error = function(e) NULL)
  if (is.null(tax) || nrow(tax) == 0) next
  rank_col <- names(tax)[1]
  sample_cols <- intersect(names(tax), md$SampleID)
  if (length(sample_cols) < 4) next
  mat <- as.matrix(tax[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- tax[[rank_col]]

  md2 <- md[match(sample_cols, md$SampleID), , drop = FALSE]
  role <- ifelse(is_control_group(md2$Group, cfg), "Control",
                 ifelse(is_treatment_group(md2$Group, cfg), "Treatment", NA_character_))
  if (sum(role == "Control", na.rm = TRUE) < 2 || sum(role == "Treatment", na.rm = TRUE) < 2) next

  for (taxon in rownames(mat)) {
    vals <- mat[taxon, ]
    cvals <- vals[role == "Control"]
    tvals <- vals[role == "Treatment"]
    cvals <- cvals[is.finite(cvals)]
    tvals <- tvals[is.finite(tvals)]
    if (length(cvals) < 2 || length(tvals) < 2) next
    eps <- min(c(cvals[cvals > 0], tvals[tvals > 0]), na.rm = TRUE) / 2
    if (!is.finite(eps)) eps <- 1e-6
    mt <- mean(tvals + eps)
    mc <- mean(cvals + eps)
    sdt <- sd(tvals + eps)
    sdc <- sd(cvals + eps)
    yi <- log(mt / mc)
    vi <- (sdt^2 / (length(tvals) * mt^2)) + (sdc^2 / (length(cvals) * mc^2))
    p <- tryCatch(wilcox.test(tvals, cvals, exact = FALSE)$p.value, error = function(e) NA_real_)
    composition_rows[[length(composition_rows) + 1]] <- data.frame(
      StudyID = study_id,
      Domain = domain,
      Rank = rank_col,
      Taxon = taxon,
      yi = yi,
      vi = vi,
      n_t = length(tvals),
      n_c = length(cvals),
      MeanTreatment = mean(tvals),
      MeanControl = mean(cvals),
      WilcoxP = p,
      stringsAsFactors = FALSE
    )
  }
}

composition_effects <- if (length(composition_rows) == 0) data.frame() else bind_rows(composition_rows)
write_csv(composition_effects, "phylum_composition_lnRR_effects.csv")

composition_summary_rows <- list()
if (nrow(composition_effects) > 0) {
  top_taxa <- composition_effects %>%
    group_by(Domain, Taxon) %>%
    summarise(k = n(), mean_control = mean(MeanControl, na.rm = TRUE), mean_treatment = mean(MeanTreatment, na.rm = TRUE), .groups = "drop") %>%
    filter(k >= 5) %>%
    group_by(Domain) %>%
    arrange(desc(mean_control + mean_treatment), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
  for (i in seq_len(nrow(top_taxa))) {
    domain <- top_taxa$Domain[i]
    taxon <- top_taxa$Taxon[i]
    d <- composition_effects[composition_effects$Domain == domain & composition_effects$Taxon == taxon, , drop = FALSE]
    fit <- fit_random_effect(d)
    composition_summary_rows[[length(composition_summary_rows) + 1]] <- cbind(
      data.frame(Domain = domain, Taxon = taxon, MeanControl = top_taxa$mean_control[i], MeanTreatment = top_taxa$mean_treatment[i]),
      fit
    )
  }
}
composition_summary <- if (length(composition_summary_rows) == 0) data.frame() else bind_rows(composition_summary_rows)
if (nrow(composition_summary) > 0) {
  composition_summary <- composition_summary %>%
    group_by(Domain) %>%
    mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(Direction = case_when(
      is.na(estimate) ~ "NA",
      FDR_BH < 0.05 & estimate > 0 ~ "Significant increase",
      FDR_BH < 0.05 & estimate < 0 ~ "Significant decrease",
      estimate > 0 ~ "Non-significant increase",
      estimate < 0 ~ "Non-significant decrease",
      TRUE ~ "No change"
    ))
}
write_csv(composition_summary, "phylum_composition_random_effect_lnRR_summary.csv")

# Simple figures for quick inspection.
alpha_plot <- alpha_summary %>%
  mutate(Metric = factor(Metric, levels = unique(Metric))) %>%
  ggplot(aes(Metric, estimate, color = Domain)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18, position = position_dodge(width = 0.45)) +
  geom_point(position = position_dodge(width = 0.45), size = 2) +
  coord_flip() +
  labs(title = "Alpha diversity lnRR summary", x = NULL, y = "ln response ratio (N / CK)") +
  theme_classic(base_size = 8)
ggsave(file.path(fig_dir, "alpha_lnRR_summary.png"), alpha_plot, width = 5.8, height = 3.8, dpi = 450)
ggsave(file.path(fig_dir, "alpha_lnRR_summary.pdf"), alpha_plot, width = 5.8, height = 3.8)

if (nrow(composition_summary) > 0) {
  comp_plot <- composition_summary %>%
    group_by(Domain) %>%
    arrange(p_value, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(Taxon = factor(Taxon, levels = rev(unique(Taxon)))) %>%
    ggplot(aes(Taxon, estimate, color = Domain)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18) +
    geom_point(size = 2) +
    facet_wrap(~ Domain, scales = "free_y") +
    coord_flip() +
    labs(title = "Phylum composition lnRR summary", x = NULL, y = "ln response ratio (N / CK)") +
    theme_classic(base_size = 8)
  ggsave(file.path(fig_dir, "phylum_composition_lnRR_summary.png"), comp_plot, width = 7.2, height = 4.5, dpi = 450)
  ggsave(file.path(fig_dir, "phylum_composition_lnRR_summary.pdf"), comp_plot, width = 7.2, height = 4.5)
}

readme <- c(
  "Alpha, beta, and phylum composition significance summary",
  "",
  "Alpha:",
  "- Input: output/03_effect_size/alpha_lnrr_effect_sizes.csv.",
  "- Model: random-effects meta-analysis of lnRR (N addition / CK) using metafor.",
  "",
  "Beta:",
  "- Input: result/beta/bray_curtis.txt and result/metadata.txt for each domain.",
  "- Tests recomputed using CK-vs-treatment role: PERMANOVA/adonis2, ANOSIM, MRPP, beta dispersion.",
  "- P values are BH-FDR adjusted within each microbial Domain.",
  "",
  "Composition:",
  "- Input: result/tax/sum_p.txt for each domain.",
  "- Rank: phylum.",
  "- Effect size: lnRR of phylum relative abundance percentage (N addition / CK).",
  "- Summary: random-effects meta-analysis for common phyla with >=5 study-domain effects.",
  "",
  "Outputs are in this folder:",
  out_dir
)
writeLines(readme, file.path(out_dir, "README.txt"), useBytes = TRUE)

message("Alpha summary rows: ", nrow(alpha_summary))
message("Beta test rows: ", nrow(beta_tests))
message("Composition effect rows: ", nrow(composition_effects))
message("Composition summary rows: ", nrow(composition_summary))
message("Output dir: ", out_dir)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/run_effect_size_moderator_models.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(lme4)
  library(lmerTest)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "scripts/run_effect_size_moderator_models.R"
root_dir <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
setwd(root_dir)

source("R/utils.R")
source("R/00_config.R")
cfg <- load_config()

out_dir <- file.path(cfg$output_dir, "06_moderator_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

beta_path <- file.path(cfg$output_dir, "05_significance_summary", "beta_recomputed_permanova_anosim_mrpp.csv")
alpha_path <- file.path(cfg$output_dir, "03_effect_size", "alpha_lnrr_effect_sizes.csv")

safe_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))

write_csv <- function(x, name) {
  path <- file.path(out_dir, name)
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

master <- read_master_metadata(cfg)
master_domain <- master %>%
  group_by(StudyID, Domain) %>%
  summarise(
    Title = first(Title),
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    MAP = mean(MAP, na.rm = TRUE),
    MAT = mean(MAT, na.rm = TRUE),
    pH = mean(pH, na.rm = TRUE),
    SOC = mean(SOC, na.rm = TRUE),
    Nitrogen = mean(Nitrogen, na.rm = TRUE),
    Ecosystem = first(na.omit(Ecosystem)),
    N_form = first(na.omit(N_form)),
    N_dose = mean(N_dose, na.rm = TRUE),
    Duration = mean(Duration, na.rm = TRUE),
    Soil_layer = first(na.omit(Soil_layer)),
    Platform = first(na.omit(Platform)),
    Target_Region = first(na.omit(Target_Region)),
    .groups = "drop"
  )

for (col in c("Latitude", "Longitude", "MAP", "MAT", "pH", "SOC", "Nitrogen", "N_dose", "Duration")) {
  master_domain[[col]][!is.finite(master_domain[[col]])] <- NA_real_
}

beta <- read.csv(beta_path, check.names = FALSE)
beta <- beta %>%
  mutate(
    BetaEffect_BetweenMinusMeanWithin = safe_num(BetweenMinusMeanWithin),
    BetaEffect_PERMANOVA_R2 = safe_num(PERMANOVA_R2),
    BetaEffect_ANOSIM_R = safe_num(ANOSIM_R),
    BetaEffect_MRPP_A = safe_num(MRPP_A)
  ) %>%
  left_join(master_domain, by = c("StudyID", "Domain"))

write_csv(beta, "beta_effect_sizes_with_moderators.csv")

alpha <- read.csv(alpha_path, check.names = FALSE) %>%
  mutate(
    AlphaEffect_lnRR = safe_num(yi),
    AlphaVariance = safe_num(vi),
    N_dose = safe_num(N_dose_mean),
    Duration = safe_num(Duration_mean)
  )
write_csv(alpha, "alpha_effect_sizes_with_moderators.csv")

moderators_numeric <- c("N_dose", "Duration", "pH", "SOC", "Nitrogen", "MAT", "MAP", "Latitude", "Longitude")
moderators_categorical <- c("Ecosystem", "N_form", "Platform", "Target_Region")

scale_numeric <- function(x) {
  x <- safe_num(x)
  if (sum(is.finite(x)) < 3 || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

fit_single_lmer <- function(data, response, moderator, domain = NULL, metric = NULL) {
  cols <- c("StudyID", response, moderator)
  d <- data[, intersect(cols, names(data)), drop = FALSE]
  if (!all(cols %in% names(d))) return(NULL)
  d <- d[!is.na(d[[response]]) & !is.na(d[[moderator]]) & !is.na(d$StudyID), , drop = FALSE]
  if (nrow(d) < 8 || length(unique(d$StudyID)) < 3) return(NULL)

  if (is.numeric(d[[moderator]])) {
    if (sd(d[[moderator]], na.rm = TRUE) == 0) return(NULL)
    d$ModeratorValue <- scale_numeric(d[[moderator]])
    if (all(is.na(d$ModeratorValue))) return(NULL)
  } else {
    d[[moderator]] <- trimws(as.character(d[[moderator]]))
    d <- d[d[[moderator]] != "", , drop = FALSE]
    keep_levels <- names(which(table(d[[moderator]]) >= 3))
    d <- d[d[[moderator]] %in% keep_levels, , drop = FALSE]
    if (nrow(d) < 8 || length(unique(d[[moderator]])) < 2) return(NULL)
    d$ModeratorValue <- factor(d[[moderator]])
  }

  form <- stats::as.formula(paste(response, "~ ModeratorValue + (1|StudyID)"))
  fit <- tryCatch(lmer(form, data = d, REML = TRUE), error = function(e) NULL)
  if (is.null(fit)) {
    fit <- tryCatch(lm(stats::as.formula(paste(response, "~ ModeratorValue")), data = d), error = function(e) NULL)
    model_type <- "lm_fallback"
  } else {
    model_type <- "lmer"
  }
  if (is.null(fit)) return(NULL)
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(co) || nrow(co) < 2) return(NULL)
  rows <- rownames(co)[rownames(co) != "(Intercept)"]
  out <- lapply(rows, function(term) {
    data.frame(
      Response = response,
      Domain = domain %||% NA_character_,
      Metric = metric %||% NA_character_,
      Moderator = moderator,
      Term = term,
      Estimate = co[term, "Estimate"],
      SE = co[term, "Std. Error"],
      DF = if ("df" %in% colnames(co)) co[term, "df"] else NA_real_,
      Tvalue = if ("t value" %in% colnames(co)) co[term, "t value"] else NA_real_,
      Pvalue = if ("Pr(>|t|)" %in% colnames(co)) co[term, "Pr(>|t|)"] else NA_real_,
      N = nrow(d),
      Studies = length(unique(d$StudyID)),
      Model = model_type,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

fit_set <- function(data, response, domain = NULL, metric = NULL) {
  rows <- list()
  for (mod in moderators_numeric) {
    if (!(mod %in% names(data))) next
    dd <- data
    dd[[mod]] <- safe_num(dd[[mod]])
    res <- fit_single_lmer(dd, response, mod, domain, metric)
    if (!is.null(res)) rows[[length(rows) + 1]] <- res
  }
  for (mod in moderators_categorical) {
    if (!(mod %in% names(data))) next
    res <- fit_single_lmer(data, response, mod, domain, metric)
    if (!is.null(res)) rows[[length(rows) + 1]] <- res
  }
  if (length(rows) == 0) data.frame() else bind_rows(rows)
}

beta_model_rows <- list()
for (domain in sort(unique(beta$Domain))) {
  bd <- beta[beta$Domain == domain, , drop = FALSE]
  beta_model_rows[[length(beta_model_rows) + 1]] <- fit_set(bd, "BetaEffect_BetweenMinusMeanWithin", domain, "BetweenMinusMeanWithin")
  beta_model_rows[[length(beta_model_rows) + 1]] <- fit_set(bd, "BetaEffect_PERMANOVA_R2", domain, "PERMANOVA_R2")
}
beta_models <- bind_rows(beta_model_rows)
if (nrow(beta_models) > 0) {
  beta_models <- beta_models %>%
    group_by(Response, Domain) %>%
    mutate(FDR_BH = p.adjust(Pvalue, method = "BH")) %>%
    ungroup()
}
write_csv(beta_models, "beta_moderator_lmer_results.csv")

alpha_model_rows <- list()
for (domain in sort(unique(alpha$Domain))) {
  for (metric in sort(unique(alpha$Metric))) {
    ad <- alpha[alpha$Domain == domain & alpha$Metric == metric, , drop = FALSE]
    alpha_model_rows[[length(alpha_model_rows) + 1]] <- fit_set(ad, "AlphaEffect_lnRR", domain, metric)
  }
}
alpha_models <- bind_rows(alpha_model_rows)
if (nrow(alpha_models) > 0) {
  alpha_models <- alpha_models %>%
    group_by(Domain, Metric) %>%
    mutate(FDR_BH = p.adjust(Pvalue, method = "BH")) %>%
    ungroup()
}
write_csv(alpha_models, "alpha_moderator_lmer_results.csv")

summarise_sig <- function(x) {
  if (nrow(x) == 0) return(data.frame())
  x %>%
    mutate(SignificantFDR = FDR_BH < 0.05) %>%
    arrange(FDR_BH, Pvalue)
}

write_csv(summarise_sig(beta_models), "beta_moderator_lmer_results_sorted.csv")
write_csv(summarise_sig(alpha_models), "alpha_moderator_lmer_results_sorted.csv")

readme <- c(
  "Effect size moderator models",
  "",
  "Beta diversity calculation used here:",
  "BetaEffect_BetweenMinusMeanWithin = mean Bray-Curtis distance between CK and N samples - mean within-group Bray-Curtis distance.",
  "Positive values mean CK and N are more compositionally separated than expected from within-group heterogeneity.",
  "",
  "Additional beta effect metrics included:",
  "- PERMANOVA R2.",
  "- ANOSIM R.",
  "- MRPP A.",
  "",
  "Moderator models:",
  "EffectSize ~ Moderator + (1 | StudyID), fitted with lmerTest/lme4 using REML.",
  "Numeric moderators are standardized before modeling.",
  "Categorical moderators are modeled as factors when at least two levels have >=3 observations.",
  "",
  "Numeric moderators tested:",
  paste(moderators_numeric, collapse = ", "),
  "",
  "Categorical moderators tested:",
  paste(moderators_categorical, collapse = ", "),
  "",
  "Outputs:",
  file.path(out_dir, "beta_effect_sizes_with_moderators.csv"),
  file.path(out_dir, "beta_moderator_lmer_results.csv"),
  file.path(out_dir, "beta_moderator_lmer_results_sorted.csv"),
  file.path(out_dir, "alpha_effect_sizes_with_moderators.csv"),
  file.path(out_dir, "alpha_moderator_lmer_results.csv")
)
writeLines(readme, file.path(out_dir, "README.txt"), useBytes = TRUE)

message("Beta effect rows: ", nrow(beta))
message("Beta moderator coefficient rows: ", nrow(beta_models))
message("Alpha effect rows: ", nrow(alpha))
message("Alpha moderator coefficient rows: ", nrow(alpha_models))
message("Output dir: ", out_dir)


# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/Meta_analysis.R
# ================================================================

rm(list = ls())
setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

# Meta-analysis workflow for microbial diversity and community composition.
# Input:  Metadata_new.xlsx
# Output: outputs/meta_effects_all.csv and overall summary/forest plots.

required_packages <- c("readxl", "dplyr", "metafor", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Please install them before running this script.")
}

library(readxl)
library(dplyr)
library(metafor)
library(ggplot2)

metadata_file <- if (file.exists("Metadata_new_Nform_reclassified.xlsx")) {
  "Metadata_new_Nform_reclassified.xlsx"
} else {
  "Metadata_new.xlsx"
}
out_dir <- "outputs"
plot_dir <- file.path(out_dir, "plots")
effect_table_dir <- file.path(out_dir, "effect_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(effect_table_dir, showWarnings = FALSE, recursive = TRUE)

zero_correction <- 0.001
min_rows_for_model <- 3

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | x == "NA" | x == "NaN"] <- NA_character_
  x
}

safe_factor_text <- function(x) {
  x <- clean_text(x)
  x[is.na(x)] <- "Unknown"
  x
}

make_n_dose_class <- function(x) {
  x <- to_num(x)
  cut(
    x,
    breaks = c(-Inf, 50, 100, 200, Inf),
    labels = c("<=50", "50-100", "100-200", ">200"),
    right = TRUE
  )
}

make_duration_class <- function(x) {
  x <- to_num(x)
  cut(
    x,
    breaks = c(-Inf, 2, 5, 10, Inf),
    labels = c("<=2", "2-5", "5-10", ">10"),
    right = TRUE
  )
}

make_soil_layer_class <- function(x) {
  x <- to_num(x)
  cut(
    x,
    breaks = c(-Inf, 5, 15, 30, Inf),
    labels = c("<=5", "5-15", "15-30", ">30"),
    right = TRUE
  )
}

make_row_id <- function(dat) {
  paste(dat$Domain, dat$ID, dat$Sample, dat$N_dose, dat$Soil_layer, sep = "__")
}

add_common_columns <- function(es, source_dat, response_group, metric, measure, note) {
  keep_cols <- c(
    "Domain", "Sample", "Note", "ID", "Title", "Site", "Latitude", "Longitude",
    "Elevation", "MAP", "MAT", "pH", "SOC", "Nitrogen", "Bdod", "Clay", "Sand",
    "Silt", "Ecosystem", "Ecosystem_Origin", "N_form", "N_form_Origin",
    "N_form_reclassified",
    "N_context", "N_dose", "N_dose_Origin", "N_addition_times_per_year",
    "N_addition_position", "Duration", "Soil_layer", "Soil_layer_origin",
    "C_Replication", "T_Replication", "Platform", "Accession_number",
    "Target_Region", "Primer_Type"
  )
  keep_cols <- intersect(keep_cols, names(source_dat))
  out <- cbind(source_dat[, keep_cols, drop = FALSE], es[, c("yi", "vi"), drop = FALSE])
  out$ResponseGroup <- response_group
  out$Metric <- metric
  out$Measure <- measure
  out$EffectNote <- note
  out$RowID <- make_row_id(out)
  out$Domain <- safe_factor_text(out$Domain)
  out$ID <- safe_factor_text(out$ID)
  out$Ecosystem <- safe_factor_text(out$Ecosystem)
  out$N_form <- safe_factor_text(out$N_form)
  if ("N_form_reclassified" %in% names(out)) {
    out$N_form_reclassified <- safe_factor_text(out$N_form_reclassified)
  }
  out$N_context <- safe_factor_text(out$N_context)
  out$N_addition_position <- safe_factor_text(out$N_addition_position)
  out$Platform <- safe_factor_text(out$Platform)
  out$Target_Region <- safe_factor_text(out$Target_Region)
  out$Primer_Type <- safe_factor_text(out$Primer_Type)
  out$N_dose_class <- safe_factor_text(make_n_dose_class(out$N_dose))
  out$Duration_class <- safe_factor_text(make_duration_class(out$Duration))
  out$Soil_layer_class <- safe_factor_text(make_soil_layer_class(out$Soil_layer))
  out
}

write_reference_style_effect_tables <- function(meta_effects) {
  x <- meta_effects
  x$StudyID <- x$ID
  x$RR <- x$yi
  x$Vi <- x$vi

  metric_name <- function(metric) {
    if (metric == "RRBeta_ln_Dt_over_Dc") return("Beta")
    if (metric == "RRStructure_ln_Db_over_meanDcDt") return("Structure")
    paste0(toupper(substr(metric, 1, 1)), substring(metric, 2))
  }

  for (domain in sort(unique(x$Domain))) {
    domain_prefix <- ifelse(domain == "Bacterial", "Bacterial",
                            ifelse(domain == "Fungi", "Fungal", domain))
    domain_dat <- x[x$Domain == domain, , drop = FALSE]
    for (metric in sort(unique(domain_dat$Metric))) {
      one <- domain_dat[domain_dat$Metric == metric, , drop = FALSE]
      if (nrow(one) == 0) next
      file_name <- paste0(domain_prefix, metric_name(metric), ".csv")
      keep <- c(
        "StudyID", "RR", "Vi", "Domain", "ResponseGroup", "Metric",
        "Sample", "ID", "Title", "Ecosystem", "N_form", "N_form_reclassified",
        "N_addition_position",
        "N_context", "N_dose", "Duration", "Soil_layer", "C_Replication",
        "T_Replication", "Platform", "Target_Region", "Primer_Type"
      )
      keep <- intersect(keep, names(one))
      write.csv(one[, keep, drop = FALSE], file.path(effect_table_dir, file_name),
                row.names = FALSE, fileEncoding = "UTF-8")
    }
  }
}

build_rom_effect <- function(dat, metric, response_group, prefix_ck = "ck", prefix_n = "N",
                             zero_correct = FALSE) {
  ck_mean <- paste0(prefix_ck, "_", metric, "_mean")
  ck_sd <- paste0(prefix_ck, "_", metric, "_sd")
  n_mean <- paste0(prefix_n, "_", metric, "_mean")
  n_sd <- paste0(prefix_n, "_", metric, "_sd")
  needed <- c(ck_mean, ck_sd, n_mean, n_sd, "C_Replication", "T_Replication")
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0) {
    warning("Skip ", metric, ": missing columns: ", paste(missing_cols, collapse = ", "))
    return(data.frame())
  }

  x <- dat
  x$m_control <- to_num(x[[ck_mean]])
  x$sd_control <- to_num(x[[ck_sd]])
  x$m_treatment <- to_num(x[[n_mean]])
  x$sd_treatment <- to_num(x[[n_sd]])
  x$n_control <- to_num(x$C_Replication)
  x$n_treatment <- to_num(x$T_Replication)

  if (zero_correct) {
    x$m_control <- ifelse(!is.na(x$m_control) & x$m_control <= 0, zero_correction, x$m_control)
    x$m_treatment <- ifelse(!is.na(x$m_treatment) & x$m_treatment <= 0, zero_correction, x$m_treatment)
  }

  ok <- !is.na(x$m_control) & !is.na(x$m_treatment) &
    !is.na(x$sd_control) & !is.na(x$sd_treatment) &
    !is.na(x$n_control) & !is.na(x$n_treatment) &
    x$m_control > 0 & x$m_treatment > 0 &
    x$sd_control >= 0 & x$sd_treatment >= 0 &
    x$n_control >= 2 & x$n_treatment >= 2

  x <- x[ok, , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())

  es <- metafor::escalc(
    measure = "ROM",
    m1i = m_treatment, sd1i = sd_treatment, n1i = n_treatment,
    m2i = m_control, sd2i = sd_control, n2i = n_control,
    data = x
  )
  add_common_columns(
    es = es,
    source_dat = x,
    response_group = response_group,
    metric = metric,
    measure = "lnRR",
    note = ifelse(zero_correct, paste0("ROM with zero correction ", zero_correction), "ROM")
  )
}

pair_count <- function(n) {
  n <- to_num(n)
  ifelse(!is.na(n) & n >= 2, n * (n - 1) / 2, NA_real_)
}

build_beta_effect <- function(dat, metric_col, metric_label) {
  needed <- c(
    metric_col,
    "Dc_BrayWithinCK_mean", "Dc_BrayWithinCK_sd",
    "Dt_BrayWithinN_mean", "Dt_BrayWithinN_sd",
    "Db_BrayBetweenCKN_mean", "Db_BrayBetweenCKN_sd",
    "Db_pairwise_n", "C_Replication", "T_Replication"
  )
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0) {
    warning("Skip beta metric ", metric_label, ": missing columns: ", paste(missing_cols, collapse = ", "))
    return(data.frame())
  }

  x <- dat
  x$yi <- to_num(x[[metric_col]])
  x$dc <- to_num(x$Dc_BrayWithinCK_mean)
  x$dc_sd <- to_num(x$Dc_BrayWithinCK_sd)
  x$dt <- to_num(x$Dt_BrayWithinN_mean)
  x$dt_sd <- to_num(x$Dt_BrayWithinN_sd)
  x$db <- to_num(x$Db_BrayBetweenCKN_mean)
  x$db_sd <- to_num(x$Db_BrayBetweenCKN_sd)
  x$n_dc <- pair_count(x$C_Replication)
  x$n_dt <- pair_count(x$T_Replication)
  x$n_db <- to_num(x$Db_pairwise_n)
  fallback_n_db <- to_num(x$C_Replication) * to_num(x$T_Replication)
  x$n_db <- ifelse(is.na(x$n_db) | x$n_db <= 0, fallback_n_db, x$n_db)

  if (metric_col == "beta_RRBeta_ln_Dt_over_Dc") {
    x$vi <- (x$dt_sd^2) / (x$n_dt * x$dt^2) + (x$dc_sd^2) / (x$n_dc * x$dc^2)
  } else {
    x$vi <- (x$db_sd^2) / (x$n_db * x$db^2) +
      ((x$dc_sd^2 / x$n_dc) + (x$dt_sd^2 / x$n_dt)) / ((x$dc + x$dt)^2)
  }

  ok <- !is.na(x$yi) & !is.na(x$vi) & is.finite(x$yi) & is.finite(x$vi) & x$vi > 0
  x <- x[ok, , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())

  add_common_columns(
    es = x[, c("yi", "vi"), drop = FALSE],
    source_dat = x,
    response_group = "Beta diversity",
    metric = metric_label,
    measure = "ln distance ratio",
    note = "Beta yi from workbook; vi approximated by delta method from Bray distance means/SDs"
  )
}

fit_meta_model <- function(dat) {
  dat <- dat[!is.na(dat$yi) & !is.na(dat$vi) & dat$vi > 0 & is.finite(dat$yi) & is.finite(dat$vi), ]
  if (nrow(dat) < min_rows_for_model) return(NULL)
  dat$ID <- safe_factor_text(dat$ID)
  dat$RowID <- safe_factor_text(dat$RowID)
  tryCatch(
    metafor::rma.mv(yi = yi, V = vi, random = ~ 1 | ID, data = dat, method = "REML"),
    error = function(e) {
      tryCatch(metafor::rma(yi = yi, vi = vi, data = dat, method = "REML"), error = function(e2) NULL)
    }
  )
}

summarize_model <- function(dat) {
  dat <- dat[!is.na(dat$yi) & !is.na(dat$vi) & dat$vi > 0, ]
  if (nrow(dat) < min_rows_for_model) return(NULL)
  res <- fit_meta_model(dat)
  if (is.null(res)) return(NULL)
  data.frame(
    ResponseGroup = unique(dat$ResponseGroup)[1],
    Domain = unique(dat$Domain)[1],
    Metric = unique(dat$Metric)[1],
    k = nrow(dat),
    studies = length(unique(dat$ID)),
    estimate = as.numeric(res$beta[1]),
    se = as.numeric(res$se[1]),
    ci_low = as.numeric(res$ci.lb[1]),
    ci_high = as.numeric(res$ci.ub[1]),
    p_value = as.numeric(res$pval[1]),
    QE = ifelse(!is.null(res$QE), as.numeric(res$QE), NA_real_),
    QEp = ifelse(!is.null(res$QEp), as.numeric(res$QEp), NA_real_),
    stringsAsFactors = FALSE
  )
}

summarize_bias_diagnostics <- function(dat) {
  dat <- dat[!is.na(dat$yi) & !is.na(dat$vi) & dat$vi > 0 &
               is.finite(dat$yi) & is.finite(dat$vi), , drop = FALSE]
  if (nrow(dat) < min_rows_for_model) return(NULL)
  simple_model <- tryCatch(
    metafor::rma(yi = yi, vi = vi, data = dat, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(simple_model)) return(NULL)

  egger <- tryCatch(metafor::regtest(simple_model), error = function(e) NULL)
  fsn_res <- tryCatch(metafor::fsn(x = simple_model, type = "Rosenthal"), error = function(e) NULL)

  data.frame(
    ResponseGroup = unique(dat$ResponseGroup)[1],
    Domain = unique(dat$Domain)[1],
    Metric = unique(dat$Metric)[1],
    k = nrow(dat),
    studies = length(unique(dat$ID)),
    egger_z = ifelse(is.null(egger), NA_real_, as.numeric(egger$zval)),
    egger_p = ifelse(is.null(egger), NA_real_, as.numeric(egger$pval)),
    fsn_rosenthal = ifelse(is.null(fsn_res), NA_real_, as.numeric(fsn_res$fsnum)),
    stringsAsFactors = FALSE
  )
}

plot_overall_forest <- function(summary_dat, raw_dat, response_group, filename) {
  s <- summary_dat[summary_dat$ResponseGroup == response_group, , drop = FALSE]
  r <- raw_dat[raw_dat$ResponseGroup == response_group, , drop = FALSE]
  if (nrow(s) == 0) return(invisible(NULL))

  s$Label <- paste(s$Domain, s$Metric, sep = " - ")
  s$Label <- factor(s$Label, levels = rev(s$Label[order(s$Domain, s$Metric)]))
  r$Label <- paste(r$Domain, r$Metric, sep = " - ")
  r$Label <- factor(r$Label, levels = levels(s$Label))
  s$stars <- ifelse(s$p_value < 0.001, "***",
                    ifelse(s$p_value < 0.01, "**",
                           ifelse(s$p_value < 0.05, "*", "")))
  s$n_label <- paste0("k=", s$k, ", study=", s$studies, s$stars)
  text_x <- max(s$ci_high, na.rm = TRUE)

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_jitter(data = r, aes(x = yi, y = Label, color = Domain),
                height = 0.16, width = 0, alpha = 0.25, size = 1.6) +
    geom_errorbar(data = s, aes(y = Label, xmin = ci_low, xmax = ci_high),
                  orientation = "y", width = 0.18, linewidth = 0.8, color = "black") +
    geom_point(data = s, aes(x = estimate, y = Label, fill = Domain),
               shape = 21, size = 3.5, color = "black") +
    geom_text(data = s, aes(y = Label, label = n_label), x = text_x,
              hjust = -0.05, size = 3.2) +
    scale_color_manual(values = c(Bacterial = "#4E79A7", Fungi = "#59A14F", Unknown = "grey50")) +
    scale_fill_manual(values = c(Bacterial = "#4E79A7", Fungi = "#59A14F", Unknown = "grey50")) +
    labs(x = "Effect size", y = NULL, title = response_group) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text = element_text(color = "black")
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 60, 10, 10))

  ggsave(file.path(plot_dir, filename), p, width = 9, height = max(4, 0.35 * nrow(s) + 1.5))
  invisible(p)
}

d0 <- readxl::read_excel(metadata_file)

numeric_cols <- c(
  "Latitude", "Longitude", "Elevation", "MAP", "MAT", "pH", "SOC", "Nitrogen",
  "Bdod", "Clay", "Sand", "Silt", "N_dose", "N_addition_times_per_year",
  "Duration", "Soil_layer", "C_Replication", "T_Replication"
)
for (col in intersect(numeric_cols, names(d0))) d0[[col]] <- to_num(d0[[col]])

alpha_metrics <- c("richness", "shannon", "simpson")
bacterial_phyla <- c("Proteobacteria", "Actinobacteria", "Acidobacteria", "Firmicutes", "Chloroflexi")
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Mortierellomycota", "Chytridiomycota", "Rozellomycota")

effects <- list()

for (m in alpha_metrics) {
  effects[[paste0("alpha_", m)]] <- build_rom_effect(d0, m, "Alpha diversity", zero_correct = FALSE)
}

effects[["beta_RRBeta"]] <- build_beta_effect(
  d0,
  metric_col = "beta_RRBeta_ln_Dt_over_Dc",
  metric_label = "RRBeta_ln_Dt_over_Dc"
)
effects[["beta_RRStructure"]] <- build_beta_effect(
  d0,
  metric_col = "beta_RRStructure_ln_Db_over_meanDcDt",
  metric_label = "RRStructure_ln_Db_over_meanDcDt"
)

for (m in bacterial_phyla) {
  effects[[paste0("composition_", m)]] <- build_rom_effect(d0, m, "Phylum composition", zero_correct = TRUE)
}
for (m in fungal_phyla) {
  effects[[paste0("composition_", m)]] <- build_rom_effect(d0, m, "Phylum composition", zero_correct = TRUE)
}

meta_effects <- bind_rows(effects)
meta_effects <- meta_effects[!is.na(meta_effects$yi) & !is.na(meta_effects$vi) &
                               is.finite(meta_effects$yi) & is.finite(meta_effects$vi) &
                               meta_effects$vi > 0, , drop = FALSE]
meta_effects$StudyID <- meta_effects$ID
meta_effects$RR <- meta_effects$yi
meta_effects$Vi <- meta_effects$vi

write.csv(meta_effects, file.path(out_dir, "meta_effects_all.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write_reference_style_effect_tables(meta_effects)

overall_summary <- bind_rows(lapply(
  split(meta_effects, list(meta_effects$ResponseGroup, meta_effects$Domain, meta_effects$Metric), drop = TRUE),
  summarize_model
))

if (nrow(overall_summary) > 0) {
  overall_summary$significance <- ifelse(overall_summary$p_value < 0.001, "***",
                                         ifelse(overall_summary$p_value < 0.01, "**",
                                                ifelse(overall_summary$p_value < 0.05, "*", "")))
}

write.csv(overall_summary, file.path(out_dir, "meta_overall_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

bias_diagnostics <- bind_rows(lapply(
  split(meta_effects, list(meta_effects$ResponseGroup, meta_effects$Domain, meta_effects$Metric), drop = TRUE),
  summarize_bias_diagnostics
))
write.csv(bias_diagnostics, file.path(out_dir, "meta_bias_diagnostics_NC_style.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

plot_overall_forest(overall_summary, meta_effects, "Alpha diversity", "forest_alpha_diversity.pdf")
plot_overall_forest(overall_summary, meta_effects, "Beta diversity", "forest_beta_diversity.pdf")
plot_overall_forest(overall_summary, meta_effects, "Phylum composition", "forest_phylum_composition.pdf")

cat("Meta-analysis data prepared.\n")
cat("Effect-size rows:", nrow(meta_effects), "\n")
cat("Overall summaries:", nrow(overall_summary), "\n")
cat("Outputs written to:", normalizePath(out_dir, winslash = "/"), "\n")

