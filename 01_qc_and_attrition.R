# QC, sample retention and attrition summaries
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 3
# Generated on: 2026-06-25


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/audit_sample_retention_and_analysis_counts.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "scripts/audit_sample_retention_and_analysis_counts.R"
}
root_dir <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
setwd(root_dir)

source("R/utils.R")
source("R/00_config.R")
cfg <- load_config()

out_dir <- file.path(cfg$output_dir, "05_significance_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_id_from_fastq <- function(files) {
  x <- basename(files)
  x <- sub("\\.fastq\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.fq\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("(_R?[12])$", "", x, ignore.case = TRUE)
  unique(x[nzchar(x)])
}

read_alpha_samples <- function(domain_row) {
  path <- file.path(domain_row$ResultPath, "alpha", "vegan.txt")
  x <- safe_read_tsv(path)
  if (is.null(x) || !("SampleID" %in% names(x))) character(0) else unique(trimws(as.character(x$SampleID)))
}

read_beta_samples <- function(domain_row) {
  path <- file.path(domain_row$ResultPath, "beta", "bray_curtis.txt")
  x <- safe_read_tsv(path)
  if (is.null(x) || ncol(x) < 2) return(character(0))
  unique(trimws(as.character(x[[1]])))
}

read_otutab_rare_samples <- function(domain_row) {
  path <- file.path(domain_row$ResultPath, "otutab_rare.txt")
  x <- safe_read_tsv(path)
  if (is.null(x) || ncol(x) < 2) return(character(0))
  unique(trimws(names(x)[-1]))
}

map_raw_seq_ids_to_internal <- function(seq_ids, domain_row) {
  map_path <- file.path(domain_row$ResultPath, "sample_id_map.tsv")
  x <- safe_read_tsv(map_path)
  if (is.null(x) || !all(c("SampleID", "InternalSampleID") %in% names(x))) {
    return(seq_ids)
  }
  raw <- trimws(as.character(x$SampleID))
  internal <- trimws(as.character(x$InternalSampleID))
  keep <- nzchar(raw) & nzchar(internal)
  raw <- raw[keep]
  internal <- internal[keep]
  mapped <- internal[match(seq_ids, raw)]
  ifelse(is.na(mapped) | !nzchar(mapped), seq_ids, mapped)
}

read_rarefaction_report <- function(domain_row) {
  path <- file.path(domain_row$ResultPath, "rarefaction_sample_report.tsv")
  x <- safe_read_tsv(path)
  if (is.null(x) || !("SampleID" %in% names(x))) return(NULL)
  x$SampleID <- trimws(as.character(x$SampleID))
  if (!("Status" %in% names(x))) x$Status <- NA_character_
  if (!("Reason" %in% names(x))) x$Reason <- NA_character_
  x
}

count_roles <- function(ids, md, cfg) {
  ids <- unique(ids[nzchar(ids)])
  if (length(ids) == 0 || is.null(md) || nrow(md) == 0) {
    return(c(total = length(ids), control = 0, treatment = 0, other = 0, missing_group = length(ids), no_metadata = length(ids)))
  }
  mm <- md[match(ids, md$SampleID), , drop = FALSE]
  has_md <- !is.na(mm$SampleID)
  grp <- mm$Group
  is_c <- has_md & is_control_group(grp, cfg)
  is_t <- has_md & is_treatment_group(grp, cfg)
  c(
    total = length(ids),
    control = sum(is_c, na.rm = TRUE),
    treatment = sum(is_t, na.rm = TRUE),
    other = sum(has_md & !is_c & !is_t, na.rm = TRUE),
    missing_group = sum(has_md & (is.na(grp) | trimws(as.character(grp)) == ""), na.rm = TRUE),
    no_metadata = sum(!has_md, na.rm = TRUE)
  )
}

add_issue <- function(issues, study_id, domain, issue_type, sample_ids, detail) {
  sample_ids <- unique(sample_ids[!is.na(sample_ids) & nzchar(sample_ids)])
  if (length(sample_ids) == 0) return(issues)
  issues[[length(issues) + 1]] <- data.frame(
    StudyID = study_id,
    Domain = domain,
    IssueType = issue_type,
    SampleID = sample_ids,
    Detail = detail,
    stringsAsFactors = FALSE
  )
  issues
}

domains <- scan_domain_dirs(cfg)
domains <- domains[domains$HasMetadata, , drop = FALSE]

rows <- list()
issues <- list()

for (i in seq_len(nrow(domains))) {
  dr <- domains[i, ]
  md <- read_domain_metadata(dr)
  if (is.null(md)) next
  md$SampleID <- trimws(as.character(md$SampleID))
  md <- md[!is.na(md$SampleID) & md$SampleID != "", , drop = FALSE]
  md_unique <- md[!duplicated(md$SampleID), , drop = FALSE]

  meta_ids <- unique(md_unique$SampleID)
  seq_files <- if (dir.exists(dr$SeqPath)) {
    list.files(dr$SeqPath, pattern = "\\.(fastq|fq)\\.gz$", full.names = TRUE, ignore.case = TRUE)
  } else {
    character(0)
  }
  seq_ids <- sample_id_from_fastq(seq_files)
  seq_ids <- unique(map_raw_seq_ids_to_internal(seq_ids, dr))

  rare <- read_rarefaction_report(dr)
  has_rare <- !is.null(rare)
  rare_ids <- if (is.null(rare)) character(0) else unique(rare$SampleID)
  retained_ids <- if (is.null(rare)) character(0) else unique(rare$SampleID[tolower(rare$Status) == "retained"])
  low_ids <- if (is.null(rare)) character(0) else unique(rare$SampleID[tolower(rare$Status) != "retained"])

  otu_ids <- read_otutab_rare_samples(dr)
  alpha_ids <- read_alpha_samples(dr)
  beta_ids <- read_beta_samples(dr)

  meta_roles <- count_roles(meta_ids, md_unique, cfg)
  seq_roles <- count_roles(intersect(seq_ids, meta_ids), md_unique, cfg)
  retained_roles <- count_roles(intersect(retained_ids, meta_ids), md_unique, cfg)
  otu_roles <- count_roles(intersect(otu_ids, meta_ids), md_unique, cfg)
  alpha_roles <- count_roles(intersect(alpha_ids, meta_ids), md_unique, cfg)
  beta_roles <- count_roles(intersect(beta_ids, meta_ids), md_unique, cfg)

  alpha_valid_for_lnrr <- alpha_roles[["control"]] >= cfg$min_group_replicates &&
    alpha_roles[["treatment"]] >= cfg$min_group_replicates
  beta_valid_for_tests <- beta_roles[["control"]] >= 2 &&
    beta_roles[["treatment"]] >= 2 &&
    beta_roles[["total"]] >= 4

  reasons <- character(0)
  if (length(setdiff(meta_ids, seq_ids)) > 0) reasons <- c(reasons, "Metadata samples missing FASTQ in seq")
  if (length(setdiff(seq_ids, meta_ids)) > 0) reasons <- c(reasons, "Extra FASTQ samples not in metadata")
  if (!has_rare) reasons <- c(reasons, "Missing rarefaction sample report")
  if (has_rare && length(setdiff(meta_ids, retained_ids)) > 0) reasons <- c(reasons, "Metadata samples not retained after rarefaction")
  if (has_rare && length(setdiff(retained_ids, otu_ids)) > 0) reasons <- c(reasons, "Retained samples missing from otutab_rare")
  if (has_rare && length(setdiff(retained_ids, alpha_ids)) > 0) reasons <- c(reasons, "Retained samples missing from alpha table")
  if (length(setdiff(alpha_ids, meta_ids)) > 0) reasons <- c(reasons, "Alpha samples absent from metadata")
  if (length(setdiff(beta_ids, meta_ids)) > 0) reasons <- c(reasons, "Beta samples absent from metadata")
  if (!alpha_valid_for_lnrr) reasons <- c(reasons, "Alpha CK/N repeats after processing are insufficient")
  if (!beta_valid_for_tests) reasons <- c(reasons, "Beta CK/N repeats after processing are insufficient")

  rows[[length(rows) + 1]] <- data.frame(
    StudyID = dr$StudyID,
    Domain = dr$Domain,
    MetadataSamples = length(meta_ids),
    MetadataControl = meta_roles[["control"]],
    MetadataTreatment = meta_roles[["treatment"]],
    MetadataOther = meta_roles[["other"]],
    MetadataMissingGroup = meta_roles[["missing_group"]],
    SeqFastqFiles = length(seq_files),
    SeqUniqueSamples = length(seq_ids),
    SeqMatchedMetadataSamples = length(intersect(seq_ids, meta_ids)),
    HasRarefactionReport = has_rare,
    RarefactionReportSamples = length(rare_ids),
    RarefactionRetainedSamples = length(retained_ids),
    RarefactionDroppedSamples = length(low_ids),
    RetainedControl = retained_roles[["control"]],
    RetainedTreatment = retained_roles[["treatment"]],
    RetainedOther = retained_roles[["other"]],
    OtuRareSamples = length(otu_ids),
    OtuRareMatchedMetadataSamples = length(intersect(otu_ids, meta_ids)),
    AlphaSamples = length(alpha_ids),
    AlphaMatchedMetadataSamples = length(intersect(alpha_ids, meta_ids)),
    AlphaControl = alpha_roles[["control"]],
    AlphaTreatment = alpha_roles[["treatment"]],
    AlphaOther = alpha_roles[["other"]],
    AlphaNoMetadata = length(setdiff(alpha_ids, meta_ids)),
    BetaSamples = length(beta_ids),
    BetaMatchedMetadataSamples = length(intersect(beta_ids, meta_ids)),
    BetaControl = beta_roles[["control"]],
    BetaTreatment = beta_roles[["treatment"]],
    BetaOther = beta_roles[["other"]],
    BetaNoMetadata = length(setdiff(beta_ids, meta_ids)),
    AlphaEligibleStrict = alpha_valid_for_lnrr,
    BetaEligibleStrict = beta_valid_for_tests,
    MetadataToRetainedDelta = if (has_rare) length(meta_ids) - length(intersect(meta_ids, retained_ids)) else NA_integer_,
    RetainedToAlphaDelta = if (has_rare) length(retained_ids) - length(intersect(retained_ids, alpha_ids)) else NA_integer_,
    RetainedToBetaDelta = if (has_rare) length(retained_ids) - length(intersect(retained_ids, beta_ids)) else NA_integer_,
    ProblemFlag = length(reasons) > 0,
    ProblemSummary = paste(unique(reasons), collapse = "; "),
    DomainPath = dr$DomainPath,
    stringsAsFactors = FALSE
  )

  issues <- add_issue(issues, dr$StudyID, dr$Domain, "MetadataMissingFASTQ", setdiff(meta_ids, seq_ids), "SampleID exists in metadata but no matching .fastq.gz/.fq.gz basename was found in seq.")
  issues <- add_issue(issues, dr$StudyID, dr$Domain, "ExtraFASTQNotInMetadata", setdiff(seq_ids, meta_ids), "FASTQ basename exists in seq but is absent from metadata.")
  if (!has_rare) {
    issues <- add_issue(issues, dr$StudyID, dr$Domain, "MissingRarefactionReport", meta_ids[1], "rarefaction_sample_report.tsv is missing, so metadata-to-retained sample loss cannot be audited directly.")
  }
  if (has_rare) {
    issues <- add_issue(issues, dr$StudyID, dr$Domain, "MetadataNotRetainedAfterRarefaction", setdiff(meta_ids, retained_ids), "Metadata sample is not retained in rarefaction_sample_report.tsv.")
    if (length(retained_ids) >= 2) {
      issues <- add_issue(issues, dr$StudyID, dr$Domain, "RetainedMissingFromOtuRare", setdiff(retained_ids, otu_ids), "Sample retained by rarefaction report but absent from otutab_rare.txt.")
      issues <- add_issue(issues, dr$StudyID, dr$Domain, "RetainedMissingFromAlpha", setdiff(retained_ids, alpha_ids), "Sample retained by rarefaction report but absent from alpha/vegan.txt.")
    }
  }
  issues <- add_issue(issues, dr$StudyID, dr$Domain, "AlphaWithoutMetadata", setdiff(alpha_ids, meta_ids), "Alpha sample is absent from current metadata.")
  issues <- add_issue(issues, dr$StudyID, dr$Domain, "BetaWithoutMetadata", setdiff(beta_ids, meta_ids), "Beta matrix sample is absent from current metadata.")
}

audit <- bind_rows(rows)
issues_df <- if (length(issues) == 0) data.frame() else bind_rows(issues)

audit_path <- file.path(out_dir, "global_sample_retention_analysis_audit.csv")
issue_path <- file.path(out_dir, "global_sample_retention_analysis_issues_long.csv")
summary_path <- file.path(out_dir, "global_sample_retention_analysis_issue_summary.csv")

write.csv(audit, audit_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(issues_df, issue_path, row.names = FALSE, fileEncoding = "UTF-8")

summary <- audit %>%
  summarise(
    DomainUnits = n(),
    ProblemDomainUnits = sum(ProblemFlag, na.rm = TRUE),
    AlphaEligibleStrict = sum(AlphaEligibleStrict, na.rm = TRUE),
    BetaEligibleStrict = sum(BetaEligibleStrict, na.rm = TRUE),
    MetadataSamples = sum(MetadataSamples, na.rm = TRUE),
    DomainUnitsWithRarefactionReport = sum(HasRarefactionReport, na.rm = TRUE),
    RarefactionRetainedSamples = sum(RarefactionRetainedSamples, na.rm = TRUE),
    AlphaSamples = sum(AlphaSamples, na.rm = TRUE),
    BetaSamples = sum(BetaSamples, na.rm = TRUE)
  )

by_domain <- audit %>%
  group_by(Domain) %>%
  summarise(
    DomainUnits = n(),
    ProblemDomainUnits = sum(ProblemFlag, na.rm = TRUE),
    AlphaEligibleStrict = sum(AlphaEligibleStrict, na.rm = TRUE),
    BetaEligibleStrict = sum(BetaEligibleStrict, na.rm = TRUE),
    MetadataSamples = sum(MetadataSamples, na.rm = TRUE),
    DomainUnitsWithRarefactionReport = sum(HasRarefactionReport, na.rm = TRUE),
    RarefactionRetainedSamples = sum(RarefactionRetainedSamples, na.rm = TRUE),
    AlphaSamples = sum(AlphaSamples, na.rm = TRUE),
    BetaSamples = sum(BetaSamples, na.rm = TRUE),
    .groups = "drop"
  )

issue_summary <- audit %>%
  mutate(
    HasMetadataMissingFASTQ = grepl("Metadata samples missing FASTQ", ProblemSummary, fixed = TRUE),
    HasExtraFASTQ = grepl("Extra FASTQ samples", ProblemSummary, fixed = TRUE),
    HasMissingRarefactionReport = grepl("Missing rarefaction sample report", ProblemSummary, fixed = TRUE),
    HasRarefactionLoss = grepl("Metadata samples not retained", ProblemSummary, fixed = TRUE),
    HasRetainedAlphaMismatch = grepl("Retained samples missing from alpha", ProblemSummary, fixed = TRUE),
    HasAlphaNoMetadata = grepl("Alpha samples absent", ProblemSummary, fixed = TRUE),
    HasBetaNoMetadata = grepl("Beta samples absent", ProblemSummary, fixed = TRUE),
    HasAlphaInsufficientRepeats = grepl("Alpha CK/N repeats", ProblemSummary, fixed = TRUE),
    HasBetaInsufficientRepeats = grepl("Beta CK/N repeats", ProblemSummary, fixed = TRUE)
  ) %>%
  group_by(Domain) %>%
  summarise(
    DomainUnits = n(),
    MetadataMissingFASTQ = sum(HasMetadataMissingFASTQ, na.rm = TRUE),
    ExtraFASTQ = sum(HasExtraFASTQ, na.rm = TRUE),
    MissingRarefactionReport = sum(HasMissingRarefactionReport, na.rm = TRUE),
    MetadataNotRetainedAfterRarefaction = sum(HasRarefactionLoss, na.rm = TRUE),
    RetainedMissingFromAlpha = sum(HasRetainedAlphaMismatch, na.rm = TRUE),
    AlphaWithoutMetadata = sum(HasAlphaNoMetadata, na.rm = TRUE),
    BetaWithoutMetadata = sum(HasBetaNoMetadata, na.rm = TRUE),
    AlphaInsufficientRepeats = sum(HasAlphaInsufficientRepeats, na.rm = TRUE),
    BetaInsufficientRepeats = sum(HasBetaInsufficientRepeats, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  bind_rows(
    cbind(Level = "All", summary),
    cbind(Level = "ByDomain", by_domain)
  ),
  summary_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(issue_summary, sub("\\.csv$", "_by_issue_type.csv", summary_path), row.names = FALSE, fileEncoding = "UTF-8")

message("Wrote audit: ", audit_path)
message("Wrote issues: ", issue_path)
message("Wrote summary: ", summary_path)
print(summary)
print(by_domain)
print(issue_summary)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/summarize_overall_qc.R
# ================================================================

options(stringsAsFactors = FALSE)
root <- "E:/BaiduSyncdisk/N_deposition1/Data/Global_Nitrogen_Pipeline"
out_dir <- file.path(root, "output", "08_overall_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) if (file.exists(path)) read.csv(path, check.names = FALSE) else data.frame()
read_tsv <- function(path) if (file.exists(path)) read.delim(path, check.names = FALSE) else data.frame()
inv <- read_csv(file.path(root, "output/01_inventory/domain_inventory.csv"))
master <- read_csv(file.path(root, "output/01_inventory/master_metadata_clean.csv"))
alpha <- read_csv(file.path(root, "output/02_collected/alpha_sample_table.csv"))
beta <- read_csv(file.path(root, "output/02_collected/beta_domain_summary.csv"))
alpha_eff <- read_csv(file.path(root, "output/03_effect_size/alpha_lnrr_effect_sizes.csv"))
rare <- read_csv(file.path(root, "output/02_collected/rarefaction_sample_reports.csv"))
audit <- read_csv(file.path(root, "output/05_significance_summary/global_sample_retention_analysis_audit.csv"))
issues <- read_csv(file.path(root, "output/05_significance_summary/global_sample_retention_analysis_issues_long.csv"))
rerun <- read_csv(file.path(root, "output/00_amplicon_logs/rerun_primer_risk_domains_final_20260602_060200.csv"))
key <- function(df) paste(df$StudyID, df$Domain, sep = "||")
if (nrow(master)) {
  master$DomainKey <- paste(master$StudyID, master$Domain, sep = "||")
  master_domain <- unique(master[, intersect(c("DomainKey","StudyID","Domain","ID","Title","Accession_number","Target_Region","Primer_Type","Platform","Ecosystem","N_dose","N_time","Latitude","Longitude"), names(master))])
} else master_domain <- data.frame()
if (nrow(inv)) inv$DomainKey <- key(inv)
if (nrow(alpha)) alpha$DomainKey <- key(alpha)
if (nrow(beta)) beta$DomainKey <- key(beta)
if (nrow(alpha_eff)) alpha_eff$DomainKey <- key(alpha_eff)
if (nrow(rare)) rare$DomainKey <- key(rare)
if (nrow(audit)) audit$DomainKey <- key(audit)
# domain-level summary
all_keys <- unique(c(inv$DomainKey, audit$DomainKey, alpha$DomainKey, beta$DomainKey, alpha_eff$DomainKey, rare$DomainKey))
dom <- data.frame(DomainKey = all_keys)
parts <- strsplit(dom$DomainKey, "\\|\\|")
dom$StudyID <- vapply(parts, function(x) x[1], "")
dom$Domain <- vapply(parts, function(x) ifelse(length(x)>1, x[2], ""), "")
add_flag <- function(keys) dom$DomainKey %in% unique(keys)
dom$HasDomainFolder <- add_flag(inv$DomainKey)
dom$HasMasterMetadata <- add_flag(master_domain$DomainKey)
dom$HasRarefactionReport <- add_flag(rare$DomainKey)
dom$OutputHasAlpha <- add_flag(alpha$DomainKey)
dom$OutputHasBeta <- add_flag(beta$DomainKey)
dom$HasAlphaEffectSize <- add_flag(alpha_eff$DomainKey)
# sample counts from inventory/audit/rare/alpha
inv2 <- inv[, intersect(c("DomainKey","MetadataRows","FastqGzFiles","HasOtuRare","HasAlpha","HasBeta","HasTaxonomy"), names(inv))]
dom <- merge(dom, inv2, by="DomainKey", all.x=TRUE)
if (nrow(rare)) {
  rare_counts <- aggregate(SampleID ~ DomainKey + Status, rare, length)
  rare_wide <- reshape(rare_counts, idvar="DomainKey", timevar="Status", direction="wide")
  names(rare_wide) <- sub("SampleID\\.", "Rare_", names(rare_wide))
  dom <- merge(dom, rare_wide, by="DomainKey", all.x=TRUE)
}
if (nrow(alpha)) {
  ac <- aggregate(SampleID ~ DomainKey, alpha, function(x) length(unique(x)))
  names(ac)[2] <- "AlphaSampleN"
  dom <- merge(dom, ac, by="DomainKey", all.x=TRUE)
}
if (nrow(beta)) {
  bc <- aggregate(Domain ~ DomainKey, beta, length)
  names(bc)[2] <- "BetaRows"
  dom <- merge(dom, bc, by="DomainKey", all.x=TRUE)
}
if (nrow(alpha_eff)) {
  ec <- aggregate(Metric ~ DomainKey, alpha_eff, length)
  names(ec)[2] <- "AlphaEffectSizeN"
  dom <- merge(dom, ec, by="DomainKey", all.x=TRUE)
}
# add title; many-to-one by DomainKey may have several rows due treatments, keep first nonblank title
if (nrow(master_domain)) {
  md <- master_domain[!duplicated(master_domain$DomainKey), ]
  dom <- merge(dom, md[, intersect(c("DomainKey","Title","Accession_number","Target_Region","Primer_Type","Platform","Ecosystem"), names(md))], by="DomainKey", all.x=TRUE)
}
num_cols <- c("MetadataRows","FastqGzFiles","Rare_Retained","Rare_Discarded","AlphaSampleN","BetaRows","AlphaEffectSizeN")
for (c in intersect(num_cols, names(dom))) dom[[c]][is.na(dom[[c]])] <- 0
# issue counts by domain
if (nrow(issues)) {
  issues$DomainKey <- paste(issues$StudyID, issues$Domain, sep="||")
  ic <- aggregate(SampleID ~ DomainKey + IssueType, issues, length)
  names(ic)[3] <- "Count"
  issue_wide <- reshape(ic, idvar="DomainKey", timevar="IssueType", direction="wide")
  names(issue_wide) <- sub("Count\\.", "Issue_", names(issue_wide))
  dom <- merge(dom, issue_wide, by="DomainKey", all.x=TRUE)
}
# status category
has_issue_col <- grepl("^Issue_", names(dom))
for (c in names(dom)[has_issue_col]) dom[[c]][is.na(dom[[c]])] <- 0
dom$QCStatus <- "OK_or_usable"
dom$QCStatus[!dom$HasRarefactionReport] <- "NoRarefactionReport"
dom$QCStatus[dom$HasRarefactionReport & dom$Rare_Retained == 0] <- "ZeroRetainedAfterRarefaction"
dom$QCStatus[dom$HasRarefactionReport & dom$Rare_Retained > 0 & dom$Rare_Retained < 2] <- "OnlyOneRetainedAfterRarefaction"
dom$QCStatus[dom$HasRarefactionReport & dom$Rare_Retained >= 2 & !dom$OutputHasAlpha] <- "RetainedButNoAlpha"
# summaries
summ <- list()
add <- function(name, value) summ[[length(summ)+1]] <<- data.frame(Item=name, Value=as.character(value))
add("Domain folders in inventory", nrow(inv))
add("Domain units in audit", nrow(audit))
add("Bacterial domain units", sum(dom$Domain=="Bacterial", na.rm=TRUE))
add("Fungi domain units", sum(dom$Domain=="Fungi", na.rm=TRUE))
add("Metadata samples across audited domain units", sum(audit$MetadataSamples, na.rm=TRUE))
add("Domain units with rarefaction report", length(unique(rare$DomainKey)))
add("Samples with rarefaction report", nrow(rare))
add("Samples retained after rarefaction", sum(rare$Status=="Retained", na.rm=TRUE))
add("Samples discarded after rarefaction", sum(rare$Status=="Discarded", na.rm=TRUE))
add("Discarded percentage among reported samples", sprintf("%.1f%%", 100*sum(rare$Status=="Discarded", na.rm=TRUE)/max(1,nrow(rare))))
add("Alpha sample rows", nrow(alpha))
add("Beta sample rows", if("SampleID" %in% names(beta)) length(unique(beta$SampleID)) else NA)
add("Beta domain summary rows", nrow(beta))
add("Alpha lnRR effect sizes", nrow(alpha_eff))
add("Domain units with alpha", length(unique(alpha$DomainKey)))
add("Domain units with beta", length(unique(beta$DomainKey)))
add("Domain units with alpha effect sizes", length(unique(alpha_eff$DomainKey)))
add("Unique StudyID with alpha", length(unique(alpha$StudyID)))
add("Unique StudyID with beta", length(unique(beta$StudyID)))
add("Unique StudyID with alpha effect sizes", length(unique(alpha_eff$StudyID)))
if (nrow(master)) {
  add("Unique StudyID in master metadata", length(unique(master$StudyID)))
  add("Unique titles in master metadata", length(unique(master$Title[!is.na(master$Title) & master$Title!=""])))
  alpha_titles <- unique(merge(unique(alpha[,c("StudyID","Domain")]), master[,c("StudyID","Domain","Title")], by=c("StudyID","Domain"), all.x=FALSE)$Title)
  beta_titles <- unique(merge(unique(beta[,c("StudyID","Domain")]), master[,c("StudyID","Domain","Title")], by=c("StudyID","Domain"), all.x=FALSE)$Title)
  eff_titles <- unique(merge(unique(alpha_eff[,c("StudyID","Domain")]), master[,c("StudyID","Domain","Title")], by=c("StudyID","Domain"), all.x=FALSE)$Title)
  add("Unique titles represented in alpha", length(unique(alpha_titles[!is.na(alpha_titles) & alpha_titles!=""])))
  add("Unique titles represented in beta", length(unique(beta_titles[!is.na(beta_titles) & beta_titles!=""])))
  add("Unique titles represented in alpha effect sizes", length(unique(eff_titles[!is.na(eff_titles) & eff_titles!=""])))
}
summary_df <- do.call(rbind, summ)
# QC status and issue summaries
qc_status <- aggregate(DomainKey ~ QCStatus, dom, length); names(qc_status)[2] <- "DomainUnits"; qc_status <- qc_status[order(-qc_status$DomainUnits),]
issue_summary <- if(nrow(issues)) aggregate(SampleID ~ IssueType, issues, length) else data.frame(); if(nrow(issue_summary)){names(issue_summary)[2] <- "Count"; issue_summary <- issue_summary[order(-issue_summary$Count),]}
rare_by_domain <- if(nrow(rare)) aggregate(SampleID ~ Domain + Status, rare, length) else data.frame(); if(nrow(rare_by_domain)) names(rare_by_domain)[3] <- "Samples"
write.csv(summary_df, file.path(out_dir, "overall_qc_key_numbers.csv"), row.names=FALSE)
write.csv(qc_status, file.path(out_dir, "overall_qc_domain_status_summary.csv"), row.names=FALSE)
write.csv(issue_summary, file.path(out_dir, "overall_qc_issue_type_summary.csv"), row.names=FALSE)
write.csv(rare_by_domain, file.path(out_dir, "overall_qc_rarefaction_by_domain.csv"), row.names=FALSE)
write.csv(dom[order(dom$QCStatus, dom$StudyID, dom$Domain),], file.path(out_dir, "overall_qc_domain_level_detail.csv"), row.names=FALSE)
# literature/title list represented
if(nrow(master)){
  lit_dom <- unique(dom[dom$OutputHasAlpha | dom$OutputHasBeta | dom$HasAlphaEffectSize, c("StudyID","Domain","Title","OutputHasAlpha","OutputHasBeta","HasAlphaEffectSize","HasMasterMetadata")])
  write.csv(lit_dom, file.path(out_dir, "overall_qc_represented_literature_domains.csv"), row.names=FALSE)
}
cat("Wrote overall QC outputs to", out_dir, "\n")
print(summary_df)
print(qc_status)
print(issue_summary)



# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/make_meta_analysis_attrition_table.R
# ================================================================

root <- "E:/BaiduSyncdisk/N_deposition1/Data/Global_Nitrogen_Pipeline"
out_dir <- file.path(root, "output", "08_overall_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(rel) {
  read.csv(file.path(root, rel), check.names = FALSE, stringsAsFactors = FALSE)
}

audit <- read_csv("output/05_significance_summary/global_sample_retention_analysis_audit.csv")
rare <- read_csv("output/02_collected/rarefaction_sample_reports.csv")
alpha <- read_csv("output/02_collected/alpha_sample_table.csv")
beta <- read_csv("output/02_collected/beta_domain_summary.csv")
eff <- read_csv("output/03_effect_size/alpha_lnrr_effect_sizes.csv")
issues <- read_csv("output/05_significance_summary/global_sample_retention_analysis_issues_long.csv")
master <- read_csv("output/01_inventory/master_metadata_clean.csv")

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) == "TRUE"
}

domain_key <- function(df) paste(df$StudyID, df$Domain, sep = "||")

metadata_total <- sum(audit$MetadataSamples, na.rm = TRUE)
seq_matched <- sum(audit$SeqMatchedMetadataSamples, na.rm = TRUE)
missing_fastq <- metadata_total - seq_matched

has_rare_report <- as_bool(audit$HasRarefactionReport)
rare_report_samples <- sum(audit$MetadataSamples[has_rare_report], na.rm = TRUE)
not_auditable_by_rarefaction <- seq_matched - rare_report_samples
rare_retained <- sum(audit$RarefactionRetainedSamples, na.rm = TRUE)
rare_dropped <- rare_report_samples - rare_retained

otu_matched <- sum(audit$OtuRareMatchedMetadataSamples, na.rm = TRUE)
alpha_matched <- sum(audit$AlphaMatchedMetadataSamples, na.rm = TRUE)
beta_matched <- sum(audit$BetaMatchedMetadataSamples, na.rm = TRUE)

attrition <- data.frame(
  Step = c(
    "1_Metadata_samples",
    "2_FASTQ_matched_to_metadata",
    "3_Rarefaction_report_available",
    "4_Fixed_depth_rarefaction_retained",
    "5_OTU_rare_table_matched_to_metadata",
    "6_Alpha_table_matched_to_metadata",
    "7_Beta_table_matched_to_metadata",
    "8_Alpha_lnRR_effect_sizes"
  ),
  Kept = c(
    metadata_total,
    seq_matched,
    rare_report_samples,
    rare_retained,
    otu_matched,
    alpha_matched,
    beta_matched,
    nrow(eff)
  ),
  RemovedFromPreviousStep = c(
    NA,
    missing_fastq,
    not_auditable_by_rarefaction,
    rare_dropped,
    rare_retained - otu_matched,
    otu_matched - alpha_matched,
    otu_matched - beta_matched,
    NA
  ),
  MainReason = c(
    "Starting metadata samples across audited Bacterial/Fungi domain units.",
    "Metadata samples without matched FASTQ were removed or could not enter sequence processing.",
    "Domain units without current rarefaction_sample_report.tsv cannot be audited at sample level.",
    "Samples below fixed rarefaction depth were discarded: Bacterial <8000, Fungi <6000.",
    "Retained samples not entering otutab_rare, mainly because their domain had fewer than 2 retained samples.",
    "Retained OTU samples missing from alpha output or not matched to metadata.",
    "Retained OTU samples missing from beta output or not matched to metadata.",
    "Effect-size records are metric-level CK-vs-N lnRR observations, not sample counts."
  ),
  stringsAsFactors = FALSE
)

attrition$RemovedPercentOfPrevious <- c(
  NA,
  sprintf("%.1f%%", 100 * missing_fastq / metadata_total),
  sprintf("%.1f%%", 100 * not_auditable_by_rarefaction / seq_matched),
  sprintf("%.1f%%", 100 * rare_dropped / rare_report_samples),
  sprintf("%.1f%%", 100 * (rare_retained - otu_matched) / rare_retained),
  sprintf("%.1f%%", 100 * (otu_matched - alpha_matched) / otu_matched),
  sprintf("%.1f%%", 100 * (otu_matched - beta_matched) / otu_matched),
  NA
)

domain_summary <- data.frame(
  Metric = c(
    "Domain units audited",
    "Domain units with rarefaction report",
    "Domain units without rarefaction report",
    "Strict alpha-eligible domain units",
    "Strict beta-eligible domain units",
    "Problem domain units",
    "Bacterial strict alpha/beta eligible",
    "Fungi strict alpha/beta eligible"
  ),
  Value = c(
    nrow(audit),
    sum(as_bool(audit$HasRarefactionReport), na.rm = TRUE),
    sum(!as_bool(audit$HasRarefactionReport), na.rm = TRUE),
    sum(as_bool(audit$AlphaEligibleStrict), na.rm = TRUE),
    sum(as_bool(audit$BetaEligibleStrict), na.rm = TRUE),
    sum(as_bool(audit$ProblemFlag), na.rm = TRUE),
    sum(as_bool(audit$AlphaEligibleStrict) & audit$Domain == "Bacterial", na.rm = TRUE),
    sum(as_bool(audit$AlphaEligibleStrict) & audit$Domain == "Fungi", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

by_domain <- aggregate(
  audit[, c(
    "MetadataSamples",
    "SeqMatchedMetadataSamples",
    "RarefactionReportSamples",
    "RarefactionRetainedSamples",
    "RarefactionDroppedSamples",
    "AlphaMatchedMetadataSamples",
    "BetaMatchedMetadataSamples"
  )],
  list(Domain = audit$Domain),
  sum,
  na.rm = TRUE
)

issue_tab <- sort(table(issues$IssueType), decreasing = TRUE)
issue_summary <- data.frame(
  IssueType = names(issue_tab),
  Count = as.integer(issue_tab),
  stringsAsFactors = FALSE
)

rare_status <- as.data.frame(table(rare$Domain, rare$Status))
names(rare_status) <- c("Domain", "Status", "Samples")

alpha_keys <- unique(domain_key(alpha))
beta_keys <- unique(domain_key(beta))
eff_keys <- unique(domain_key(eff))

master_titles <- unique(master[, c("StudyID", "Domain", "Title")])
count_titles <- function(df) {
  merged <- merge(unique(df[, c("StudyID", "Domain")]), master_titles, by = c("StudyID", "Domain"))
  length(unique(merged$Title[!is.na(merged$Title) & merged$Title != ""]))
}

literature_summary <- data.frame(
  Metric = c(
    "Unique StudyID in master metadata",
    "Unique article titles in master metadata",
    "Unique StudyID represented in alpha outputs",
    "Unique StudyID represented in beta outputs",
    "Unique StudyID represented in alpha effect sizes",
    "Unique article titles matched in alpha outputs",
    "Unique article titles matched in beta outputs",
    "Unique article titles matched in alpha effect sizes"
  ),
  Value = c(
    length(unique(master$StudyID)),
    length(unique(master$Title[!is.na(master$Title) & master$Title != ""])),
    length(unique(alpha$StudyID)),
    length(unique(beta$StudyID)),
    length(unique(eff$StudyID)),
    count_titles(alpha),
    count_titles(beta),
    count_titles(eff)
  ),
  stringsAsFactors = FALSE
)

write.csv(attrition, file.path(out_dir, "meta_analysis_stepwise_sample_attrition.csv"), row.names = FALSE)
write.csv(domain_summary, file.path(out_dir, "meta_analysis_stepwise_domain_summary.csv"), row.names = FALSE)
write.csv(by_domain, file.path(out_dir, "meta_analysis_stepwise_by_domain.csv"), row.names = FALSE)
write.csv(issue_summary, file.path(out_dir, "meta_analysis_stepwise_issue_summary.csv"), row.names = FALSE)
write.csv(rare_status, file.path(out_dir, "meta_analysis_stepwise_rarefaction_status.csv"), row.names = FALSE)
write.csv(literature_summary, file.path(out_dir, "meta_analysis_stepwise_literature_summary.csv"), row.names = FALSE)

print(attrition)
print(domain_summary)
print(by_domain)
print(issue_summary)
print(literature_summary)

