# Network metric summaries and network figures
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 7
# Generated on: 2026-06-25


# ================================================================
# Source: network_analysis_pipeline/plot_network_stability_results.R
# ================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_dir <- "E:/BaiduSyncdisk/N_deposition1"
out_dir <- file.path(
  base_dir,
  "Data/Global_Nitrogen_Pipeline/output/11_network_stability_20260605"
)
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

delta_path <- file.path(table_dir, "network_stability_pair_delta_all.csv")
delta <- read_csv(delta_path, show_col_types = FALSE)

metrics <- c(
  natural_connectivity = "Natural connectivity",
  random_robustness_50 = "Random robustness (50% removed)",
  targeted_robustness_4 = "Targeted robustness (4 hubs)",
  neg_pos_cohesion_ratio = "|Negative| / positive cohesion"
)

paired_long <- list()
for (metric in names(metrics)) {
  ck_col <- paste0("ck_", metric)
  n_col <- paste0("N_", metric)
  if (!all(c(ck_col, n_col) %in% names(delta))) {
    next
  }
  tmp <- delta %>%
    select(study_id, domain, control_group, treatment_group, all_of(ck_col), all_of(n_col)) %>%
    rename(CK = all_of(ck_col), N = all_of(n_col)) %>%
    mutate(
      Metric = metrics[[metric]],
      PairID = paste(study_id, domain, treatment_group, sep = "|")
    ) %>%
    pivot_longer(c(CK, N), names_to = "Network", values_to = "Value")
  paired_long[[metric]] <- tmp
}

paired_long <- bind_rows(paired_long) %>%
  filter(!is.na(Value))

write_csv(
  paired_long,
  file.path(table_dir, "network_stability_paired_metrics_long.csv")
)

wilcox_summary <- delta %>%
  group_by(domain) %>%
  group_modify(function(df, key) {
    rows <- list()
    for (metric in names(metrics)) {
      ck_col <- paste0("ck_", metric)
      n_col <- paste0("N_", metric)
      if (!all(c(ck_col, n_col) %in% names(df))) {
        next
      }
      x <- df[[ck_col]]
      y <- df[[n_col]]
      keep <- !is.na(x) & !is.na(y)
      if (sum(keep) < 3) {
        p <- NA_real_
        statistic <- NA_real_
      } else {
        test <- wilcox.test(y[keep], x[keep], paired = TRUE, exact = FALSE)
        p <- test$p.value
        statistic <- unname(test$statistic)
      }
      rows[[length(rows) + 1]] <- tibble(
        Metric = metrics[[metric]],
        n_pairs = sum(keep),
        CK_mean = mean(x[keep], na.rm = TRUE),
        N_mean = mean(y[keep], na.rm = TRUE),
        Delta_mean = mean(y[keep] - x[keep], na.rm = TRUE),
        CK_median = median(x[keep], na.rm = TRUE),
        N_median = median(y[keep], na.rm = TRUE),
        Delta_median = median(y[keep] - x[keep], na.rm = TRUE),
        Wilcoxon_W = statistic,
        p_value = p
      )
    }
    bind_rows(rows)
  }) %>%
  ungroup() %>%
  group_by(domain) %>%
  mutate(FDR_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(
  wilcox_summary,
  file.path(table_dir, "network_stability_wilcoxon_summary.csv")
)

plot_domain <- function(domain_name) {
  df <- paired_long %>% filter(domain == domain_name)
  if (nrow(df) == 0) {
    return(NULL)
  }
  p <- ggplot(df, aes(x = Network, y = Value, group = PairID)) +
    geom_line(color = "grey72", linewidth = 0.35, alpha = 0.65) +
    geom_point(aes(color = Network), size = 1.8, alpha = 0.85) +
    stat_summary(
      aes(group = Network),
      fun = mean,
      geom = "point",
      shape = 21,
      size = 3.2,
      color = "black",
      fill = "white"
    ) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(CK = "#F2B84B", N = "#5AA0D8")) +
    labs(
      title = paste(domain_name, "network stability"),
      x = NULL,
      y = "Network metric"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      plot.title = element_text(face = "bold")
    )
  ggsave(
    file.path(fig_dir, paste0("Figure_network_stability_paired_", domain_name, ".pdf")),
    p,
    width = 8,
    height = 6
  )
  ggsave(
    file.path(fig_dir, paste0("Figure_network_stability_paired_", domain_name, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  p
}

plot_delta <- function(domain_name) {
  df <- delta %>% filter(domain == domain_name)
  if (nrow(df) == 0) {
    return(NULL)
  }
  cols <- paste0("delta_N_minus_CK_", names(metrics))
  use_cols <- cols[cols %in% names(df)]
  long <- df %>%
    select(study_id, domain, treatment_group, all_of(use_cols)) %>%
    pivot_longer(all_of(use_cols), names_to = "MetricRaw", values_to = "Delta") %>%
    mutate(
      MetricKey = sub("^delta_N_minus_CK_", "", MetricRaw),
      Metric = unname(metrics[MetricKey])
    ) %>%
    filter(!is.na(Delta))
  p <- ggplot(long, aes(x = Metric, y = Delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey92", color = "grey35") +
    geom_jitter(width = 0.16, height = 0, size = 1.2, alpha = 0.55, color = "#3E6C9A") +
    labs(
      title = paste(domain_name, "N - CK network stability effects"),
      x = NULL,
      y = "N minus CK"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(
    file.path(fig_dir, paste0("Figure_network_stability_delta_", domain_name, ".pdf")),
    p,
    width = 8,
    height = 5
  )
  ggsave(
    file.path(fig_dir, paste0("Figure_network_stability_delta_", domain_name, ".png")),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  p
}

for (domain_name in sort(unique(paired_long$domain))) {
  plot_domain(domain_name)
  plot_delta(domain_name)
}

message("Saved network stability plots and Wilcoxon summaries to: ", out_dir)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/run_group_network_analysis.R
# ================================================================

options(stringsAsFactors = FALSE)

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "scripts/run_group_network_analysis.R"
root_dir <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
setwd(root_dir)

source("R/utils.R")
source("R/00_config.R")
source("R/04_network_stability.R")

cfg <- load_config()
out_dir <- file.path(cfg$output_dir, "04_network_groups")
edge_dir <- file.path(out_dir, "edges")
node_dir <- file.path(out_dir, "nodes")
dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)

min_samples <- as.integer(Sys.getenv("GROUP_NETWORK_MIN_SAMPLES", unset = "6"))
top_taxa <- as.integer(Sys.getenv("GROUP_NETWORK_TOP_TAXA", unset = "80"))
min_taxa <- as.integer(Sys.getenv("GROUP_NETWORK_MIN_TAXA", unset = "20"))
cor_cutoff <- as.numeric(Sys.getenv("GROUP_NETWORK_COR_CUTOFF", unset = as.character(cfg$network_cor_cutoff)))
p_cutoff <- as.numeric(Sys.getenv("GROUP_NETWORK_P_CUTOFF", unset = as.character(cfg$network_p_cutoff)))

safe_file <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

write_process_note <- function(path) {
  txt <- c(
    "Group-level co-occurrence network workflow",
    "",
    "Input:",
    "- Each completed StudyID/Domain result/otutab_rare.txt.",
    "- Matching result/metadata.txt or metadata_raw.txt with SampleID and Group.",
    "",
    "Eligibility:",
    sprintf("- Bacterial and Fungi are analyzed separately."),
    sprintf("- A StudyID/Domain/Group subset is analyzed only when sample count >= %d (>5).", min_samples),
    sprintf("- Taxa are filtered to prevalence >= 2 samples and the top %d taxa by total abundance.", top_taxa),
    sprintf("- Subsets with fewer than %d retained taxa after filtering are reported but not networked.", min_taxa),
    "",
    "Network inference:",
    "- Relative abundance is transformed as log1p(relative_abundance * 10000).",
    "- Pairwise taxon associations use Spearman correlation.",
    sprintf("- Edges are retained when |rho| >= %.3f and FDR-adjusted P <= %.3f.", cor_cutoff, p_cutoff),
    "- Positive and negative edges are retained in the same network and marked by sign.",
    "",
    "Outputs:",
    "- eligible_groups_over5.csv: all StudyID/Domain/Group counts and eligibility.",
    "- network_metrics_over5.csv: topology and stability metrics for eligible subsets.",
    "- edges/*.csv and nodes/*.csv: network edge and node tables for subsets with edges.",
    "",
    "Interpretation note:",
    "- These are within-domain networks. Combined bacteria-fungi networks should only be built for studies where the same plots/samples have both marker datasets and compatible pairing metadata."
  )
  writeLines(txt, path, useBytes = TRUE)
}

build_edges_nodes <- function(sample_by_taxa, cfg_local) {
  filtered <- filter_top_taxa(sample_by_taxa, top_n = cfg_local$top_taxa, min_prevalence = 2)
  if (is.null(filtered)) {
    return(list(status = "TooFewTaxa", filtered = NULL, adj = NULL, edges = data.frame(), nodes = data.frame()))
  }
  if (nrow(filtered) < cfg_local$min_samples || ncol(filtered) < cfg_local$min_taxa) {
    return(list(status = "TooFewSamplesOrTaxa", filtered = filtered, adj = NULL, edges = data.frame(), nodes = data.frame()))
  }

  rel <- relative_abundance(filtered)
  rel <- log1p(rel * 10000)
  cmat <- suppressWarnings(stats::cor(rel, method = "spearman", use = "pairwise.complete.obs"))
  pmat <- cor_p_matrix(rel, method = "spearman")
  padj <- matrix(stats::p.adjust(as.vector(pmat), method = "fdr"), nrow = nrow(pmat), dimnames = dimnames(pmat))

  adj <- cmat
  adj[is.na(adj)] <- 0
  adj[abs(adj) < cfg_local$cor_cutoff | padj > cfg_local$p_cutoff] <- 0
  diag(adj) <- 0

  keep <- rowSums(adj != 0) > 0
  adj <- adj[keep, keep, drop = FALSE]
  if (nrow(adj) < 2 || sum(adj[upper.tri(adj)] != 0) == 0) {
    nodes <- data.frame(Taxon = colnames(filtered), TotalAbundance = colSums(filtered), Prevalence = colSums(filtered > 0), stringsAsFactors = FALSE)
    return(list(status = "NoEdges", filtered = filtered, adj = adj, edges = data.frame(), nodes = nodes))
  }

  taxa <- colnames(adj)
  idx <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  edges <- data.frame(
    Source = taxa[idx[, 1]],
    Target = taxa[idx[, 2]],
    Rho = adj[idx],
    Sign = ifelse(adj[idx] > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(
    Taxon = taxa,
    TotalAbundance = colSums(filtered[, taxa, drop = FALSE]),
    Prevalence = colSums(filtered[, taxa, drop = FALSE] > 0),
    Degree = rowSums(adj != 0),
    stringsAsFactors = FALSE
  )
  list(status = "OK", filtered = filtered, adj = adj, edges = edges, nodes = nodes)
}

metrics_from_network <- function(net, sample_n, taxa_n) {
  if (is.null(net$adj) || nrow(net$adj) < 2 || nrow(net$edges) == 0) {
    return(data.frame(
      Samples = sample_n, Taxa = taxa_n, Nodes = if (is.null(net$adj)) 0 else nrow(net$adj),
      Edges = 0, PositiveEdges = 0, NegativeEdges = 0,
      Connectance = NA_real_, MeanDegree = NA_real_, Modularity = NA_real_,
      NaturalConnectivity = NA_real_, RobustnessSlope = NA_real_,
      Status = net$status
    ))
  }
  adj <- net$adj
  upper <- adj[upper.tri(adj)]
  edges <- sum(upper != 0)
  deg <- rowSums(adj != 0)
  modularity <- NA_real_
  if (optional_package("igraph") && edges > 0) {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    cl <- tryCatch(igraph::cluster_fast_greedy(g), error = function(e) NULL)
    if (!is.null(cl)) modularity <- igraph::modularity(cl)
  }
  data.frame(
    Samples = sample_n,
    Taxa = taxa_n,
    Nodes = nrow(adj),
    Edges = edges,
    PositiveEdges = sum(upper > 0),
    NegativeEdges = sum(upper < 0),
    Connectance = edges / (nrow(adj) * (nrow(adj) - 1) / 2),
    MeanDegree = mean(deg),
    Modularity = modularity,
    NaturalConnectivity = natural_connectivity((adj != 0) * 1),
    RobustnessSlope = robustness_slope((adj != 0) * 1),
    Status = net$status
  )
}

domains <- scan_domain_dirs(cfg)
domains <- domains[domains$HasOtuRare & domains$HasMetadata & domains$HasAlpha & domains$HasBeta, , drop = FALSE]

elig_rows <- list()
metric_rows <- list()
cfg_local <- list(
  min_samples = min_samples,
  top_taxa = top_taxa,
  min_taxa = min_taxa,
  cor_cutoff = cor_cutoff,
  p_cutoff = p_cutoff
)

for (i in seq_len(nrow(domains))) {
  dr <- domains[i, ]
  otu <- read_otutab_rare(dr)
  md <- read_domain_metadata(dr)
  if (is.null(otu) || is.null(md) || !("Group" %in% names(md))) next

  sample_by_taxa <- otu_samples_by_rows(otu)
  common <- intersect(rownames(sample_by_taxa), md$SampleID)
  if (length(common) == 0) next
  sample_by_taxa <- sample_by_taxa[common, , drop = FALSE]
  md <- md[match(common, md$SampleID), , drop = FALSE]
  groups <- sort(unique(md$Group[!is.na(md$Group) & md$Group != ""]))
  if (length(groups) == 0) next

  for (group in groups) {
    sample_ids <- md$SampleID[md$Group == group]
    sample_n <- length(sample_ids)
    role <- if (is_control_group(group, cfg)) "Control" else if (is_treatment_group(group, cfg)) "Treatment" else "Other"
    eligible <- sample_n >= min_samples
    elig_rows[[length(elig_rows) + 1]] <- data.frame(
      StudyID = dr$StudyID,
      Domain = dr$Domain,
      Group = group,
      GroupRole = role,
      SampleN = sample_n,
      EligibleOver5 = eligible,
      stringsAsFactors = FALSE
    )
    if (!eligible) next

    subset_otu <- sample_by_taxa[sample_ids, , drop = FALSE]
    net <- build_edges_nodes(subset_otu, cfg_local)
    taxa_n <- if (is.null(net$filtered)) 0 else ncol(net$filtered)
    metrics <- metrics_from_network(net, sample_n = sample_n, taxa_n = taxa_n)
    metric_rows[[length(metric_rows) + 1]] <- cbind(
      data.frame(StudyID = dr$StudyID, Domain = dr$Domain, Group = group, GroupRole = role, stringsAsFactors = FALSE),
      metrics
    )

    if (nrow(net$edges) > 0) {
      stem <- safe_file(paste(dr$StudyID, dr$Domain, group, sep = "__"))
      write_csv_utf8(net$edges, file.path(edge_dir, paste0(stem, "_edges.csv")))
      write_csv_utf8(net$nodes, file.path(node_dir, paste0(stem, "_nodes.csv")))
    }
  }
}

eligible <- if (length(elig_rows) == 0) data.frame() else do.call(rbind, elig_rows)
metrics <- if (length(metric_rows) == 0) data.frame() else do.call(rbind, metric_rows)

write_csv_utf8(eligible, file.path(out_dir, "eligible_groups_over5.csv"))
write_csv_utf8(metrics, file.path(out_dir, "network_metrics_over5.csv"))
write_process_note(file.path(out_dir, "network_process_README.txt"))

message("[network-groups] Group rows: ", nrow(eligible))
message("[network-groups] Eligible groups: ", sum(eligible$EligibleOver5, na.rm = TRUE))
message("[network-groups] Metric rows: ", nrow(metrics))
message("[network-groups] Output: ", out_dir)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/run_menap_style_network_analysis.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "scripts/run_menap_style_network_analysis.R"
root_dir <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
setwd(root_dir)

source("R/utils.R")
source("R/00_config.R")

cfg <- load_config()
network_dir <- file.path(cfg$output_dir, "04_network_groups")
figure_dir <- file.path(network_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

paired_path <- file.path(network_dir, "paired_control_treatment_network_delta_metrics.csv")
if (!file.exists(paired_path)) {
  stop("Missing paired network table: ", paired_path, call. = FALSE)
}

min_mean_ra <- as.numeric(Sys.getenv("MENAP_MIN_MEAN_RA", unset = "0.0005"))
cor_cutoff <- as.numeric(Sys.getenv("MENAP_LOCAL_COR_CUTOFF", unset = "0.8"))
nperm <- as.integer(Sys.getenv("MENAP_COHESION_NPERM", unset = "50"))
random_remove_fraction <- as.numeric(Sys.getenv("MENAP_RANDOM_REMOVE_FRACTION", unset = "0.5"))
target_remove_n <- as.integer(Sys.getenv("MENAP_TARGET_REMOVE_N", unset = "4"))

safe_file <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

read_domain_pair_data <- function(study_id, domain) {
  domain_row <- data.frame(
    StudyID = study_id,
    Domain = domain,
    DomainPath = normalize_path(file.path(cfg$base_data_dir, study_id, domain)),
    ResultPath = normalize_path(file.path(cfg$base_data_dir, study_id, domain, "result")),
    stringsAsFactors = FALSE
  )
  otu <- read_otutab_rare(domain_row)
  md <- read_domain_metadata(domain_row)
  if (is.null(otu) || is.null(md)) return(NULL)
  sample_by_taxa <- otu_samples_by_rows(otu)
  common <- intersect(rownames(sample_by_taxa), md$SampleID)
  if (length(common) == 0) return(NULL)
  list(
    sample_by_taxa = sample_by_taxa[common, , drop = FALSE],
    metadata = md[match(common, md$SampleID), , drop = FALSE],
    taxonomy_path = file.path(domain_row$ResultPath, "taxonomy.txt")
  )
}

subset_group_otu <- function(pair_data, group) {
  sample_ids <- pair_data$metadata$SampleID[pair_data$metadata$Group == group]
  sample_ids <- intersect(sample_ids, rownames(pair_data$sample_by_taxa))
  if (length(sample_ids) == 0) return(NULL)
  pair_data$sample_by_taxa[sample_ids, , drop = FALSE]
}

filter_menap_taxa <- function(sample_by_taxa) {
  rel <- relative_abundance(sample_by_taxa)
  mean_ra <- colMeans(rel, na.rm = TRUE)
  keep <- names(mean_ra)[is.finite(mean_ra) & mean_ra > min_mean_ra]
  if (length(keep) < 3) return(NULL)
  sample_by_taxa[, keep, drop = FALSE]
}

make_pearson_network <- function(sample_by_taxa) {
  x <- filter_menap_taxa(sample_by_taxa)
  if (is.null(x) || nrow(x) < 3 || ncol(x) < 3) {
    return(list(status = "TooFewSamplesOrTaxa", filtered = x, adj = NULL, cor = NULL))
  }
  rel <- relative_abundance(x)
  transformed <- log1p(as.matrix(rel) * 10000)
  cmat <- suppressWarnings(cor(transformed, method = "pearson", use = "pairwise.complete.obs"))
  cmat[!is.finite(cmat)] <- 0
  diag(cmat) <- 0
  adj <- cmat
  adj[abs(adj) < cor_cutoff] <- 0
  diag(adj) <- 0
  keep <- rowSums(abs(adj)) > 0
  adj <- adj[keep, keep, drop = FALSE]
  if (nrow(adj) < 2 || sum(adj[upper.tri(adj)] != 0) == 0) {
    return(list(status = "NoEdges", filtered = x, adj = adj, cor = cmat))
  }
  list(status = "OK", filtered = x, adj = adj, cor = cmat)
}

natural_connectivity <- function(adj) {
  adj <- as.matrix(adj)
  if (nrow(adj) < 2) return(NA_real_)
  eig <- tryCatch(eigen(adj, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(eig)) return(NA_real_)
  log(mean(exp(eig)))
}

normalized_natural_connectivity <- function(adj) {
  raw <- natural_connectivity(adj)
  n <- nrow(adj)
  if (!is.finite(raw) || n < 2) return(NA_real_)
  denom <- n - log(n)
  if (!is.finite(denom) || denom == 0) raw else raw / denom
}

module_roles <- function(adj_binary) {
  g <- graph_from_adjacency_matrix(adj_binary, mode = "undirected", diag = FALSE)
  if (vcount(g) < 3 || ecount(g) == 0) {
    return(data.frame(Taxon = V(g)$name, Zi = 0, Pi = 0, Role = "Peripheral", Module = 1))
  }
  cl <- cluster_fast_greedy(g)
  module <- membership(cl)
  deg_total <- degree(g)
  modules <- sort(unique(module))
  zi <- numeric(vcount(g))
  pi <- numeric(vcount(g))
  for (m in modules) {
    idx <- which(module == m)
    subg <- induced_subgraph(g, vids = idx)
    k_in <- degree(subg)
    sig <- sd(k_in)
    zi[idx] <- if (is.finite(sig) && sig > 0) (k_in - mean(k_in)) / sig else 0
  }
  for (i in seq_len(vcount(g))) {
    ki <- deg_total[i]
    if (ki == 0) {
      pi[i] <- 0
    } else {
      neigh <- neighbors(g, i)
      counts <- vapply(modules, function(m) sum(module[neigh] == m), numeric(1))
      pi[i] <- 1 - sum((counts / ki)^2)
    }
  }
  role <- ifelse(zi > 2.5 & pi < 0.62, "Module hub",
                 ifelse(zi > 2.5 & pi > 0.62, "Network hub",
                        ifelse(zi < 2.5 & pi > 0.62, "Connector", "Peripheral")))
  data.frame(Taxon = V(g)$name, Zi = zi, Pi = pi, Role = role, Module = module, stringsAsFactors = FALSE)
}

random_remain_once <- function(adj_signed, remove_fraction) {
  n <- nrow(adj_signed)
  remove_n <- min(n, round(n * remove_fraction))
  id_rm <- if (remove_n > 0) sample(seq_len(n), remove_n) else integer(0)
  net <- adj_signed
  if (length(id_rm) > 0) {
    net[id_rm, ] <- 0
    net[, id_rm] <- 0
  }
  mean_interaction <- colMeans(net)
  failed <- which(mean_interaction <= 0)
  (n - length(failed)) / n
}

random_robustness <- function(adj_signed, remove_fraction = 0.5, nperm = 100) {
  if (is.null(adj_signed) || nrow(adj_signed) < 2) return(NA_real_)
  vals <- replicate(nperm, random_remain_once(adj_signed, remove_fraction))
  mean(vals, na.rm = TRUE)
}

targeted_remain <- function(adj_signed, roles, remove_n = 4) {
  if (is.null(adj_signed) || nrow(adj_signed) < 2) {
    return(list(value = NA_real_, removed = 0, target = "None"))
  }
  module_hubs <- roles$Taxon[roles$Role == "Module hub"]
  network_hubs <- roles$Taxon[roles$Role == "Network hub"]
  connectors <- roles$Taxon[roles$Role == "Connector"]
  targets <- unique(c(module_hubs, network_hubs, connectors))
  target_type <- "Keystone Zi-Pi"
  if (length(targets) == 0) {
    g <- graph_from_adjacency_matrix((abs(adj_signed) > 0) * 1, mode = "undirected", diag = FALSE)
    ord <- order(betweenness(g), degree(g), decreasing = TRUE)
    targets <- V(g)$name[ord]
    target_type <- "Betweenness fallback"
  }
  targets <- head(targets, remove_n)
  keep <- setdiff(rownames(adj_signed), targets)
  if (length(keep) < 2) {
    return(list(value = 0, removed = length(targets), target = target_type))
  }
  net <- adj_signed[keep, keep, drop = FALSE]
  mean_interaction <- colMeans(net)
  failed <- which(mean_interaction <= 0)
  list(value = (nrow(adj_signed) - length(targets) - length(failed)) / nrow(adj_signed),
       removed = length(targets), target = target_type)
}

betweenness_nc_slope <- function(adj_binary) {
  if (is.null(adj_binary) || nrow(adj_binary) < 5) return(NA_real_)
  g <- graph_from_adjacency_matrix(adj_binary, mode = "undirected", diag = FALSE)
  ord <- order(betweenness(g), degree(g), decreasing = TRUE)
  max_remove <- max(1, floor(vcount(g) * 0.8))
  steps <- unique(round(seq(0, max_remove, length.out = min(16, max_remove + 1))))
  nc <- vapply(steps, function(r) {
    if (r == 0) return(normalized_natural_connectivity(adj_binary))
    keep <- setdiff(seq_len(vcount(g)), ord[seq_len(r)])
    if (length(keep) < 2) return(NA_real_)
    normalized_natural_connectivity(adj_binary[keep, keep, drop = FALSE])
  }, numeric(1))
  ok <- is.finite(nc)
  if (sum(ok) < 3) return(NA_real_)
  coef(lm(nc[ok] ~ steps[ok]))[[2]]
}

cohesion_metrics <- function(sample_by_taxa, adj_taxa, nperm = 50) {
  x <- sample_by_taxa[, adj_taxa, drop = FALSE]
  rel <- as.matrix(relative_abundance(x))
  transformed <- log1p(rel * 10000)
  true_cor <- suppressWarnings(cor(transformed, method = "pearson", use = "pairwise.complete.obs"))
  true_cor[!is.finite(true_cor)] <- 0
  diag(true_cor) <- 0
  expected <- matrix(0, ncol(transformed), ncol(transformed), dimnames = dimnames(true_cor))

  set.seed(20260531)
  for (j in seq_len(ncol(transformed))) {
    perm_cor <- matrix(NA_real_, nrow = ncol(transformed), ncol = nperm)
    for (p in seq_len(nperm)) {
      y <- sample(transformed[, j])
      vals <- suppressWarnings(cor(transformed, y, method = "pearson", use = "pairwise.complete.obs"))
      vals[!is.finite(vals)] <- 0
      perm_cor[, p] <- vals
    }
    expected[, j] <- apply(perm_cor, 1, median, na.rm = TRUE)
  }
  obs_exp <- true_cor - expected
  diag(obs_exp) <- 0
  pos_conn <- apply(obs_exp, 2, function(v) {
    v <- v[v > 0 & is.finite(v)]
    if (length(v) == 0) 0 else mean(v)
  })
  neg_conn <- apply(obs_exp, 2, function(v) {
    v <- v[v < 0 & is.finite(v)]
    if (length(v) == 0) 0 else mean(v)
  })
  pos_cohesion <- as.vector(rel %*% pos_conn)
  neg_cohesion <- as.vector(rel %*% neg_conn)
  pos_mean <- mean(pos_cohesion, na.rm = TRUE)
  neg_mean <- mean(neg_cohesion, na.rm = TRUE)
  ratio <- if (is.finite(pos_mean) && pos_mean > 0) abs(neg_mean) / pos_mean else NA_real_
  list(PositiveCohesion = pos_mean, NegativeCohesion = neg_mean, AbsNegPosCohesionRatio = ratio)
}

network_metrics <- function(sample_by_taxa, study_id, domain, group, group_role, pair_id) {
  net <- make_pearson_network(sample_by_taxa)
  if (is.null(net$adj) || net$status != "OK") {
    return(list(
      metrics = data.frame(
        PairID = pair_id, StudyID = study_id, Domain = domain, Group = group, GroupRole = group_role,
        Samples = nrow(sample_by_taxa), TaxaAfterRAFilter = if (is.null(net$filtered)) 0 else ncol(net$filtered),
        Nodes = 0, Edges = 0, PositiveEdges = 0, NegativeEdges = 0, Connectance = NA_real_,
        MeanDegree = NA_real_, Modularity = NA_real_, NaturalConnectivity = NA_real_,
        BetweennessNCSlope = NA_real_, Random50Remain = NA_real_, Targeted4Remain = NA_real_,
        TargetedRemovedNodes = 0, TargetType = NA_character_, PositiveCohesion = NA_real_,
        NegativeCohesion = NA_real_, AbsNegPosCohesionRatio = NA_real_, ModuleHubs = 0,
        Connectors = 0, NetworkHubs = 0, Status = net$status
      ),
      nodes = data.frame(), edges = data.frame()
    ))
  }

  adj <- net$adj
  adj_binary <- (abs(adj) > 0) * 1
  upper <- adj[upper.tri(adj)]
  g <- graph_from_adjacency_matrix(adj_binary, mode = "undirected", diag = FALSE)
  cl <- tryCatch(cluster_fast_greedy(g), error = function(e) NULL)
  modularity <- if (is.null(cl)) NA_real_ else modularity(cl)
  deg <- degree(g)
  roles <- module_roles(adj_binary)
  target <- targeted_remain(adj, roles, target_remove_n)
  cohesion <- cohesion_metrics(net$filtered, rownames(adj), nperm = nperm)

  taxa <- rownames(adj)
  idx <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  edges <- data.frame(
    PairID = pair_id, StudyID = study_id, Domain = domain, Group = group, GroupRole = group_role,
    Source = taxa[idx[, 1]], Target = taxa[idx[, 2]], Rho = adj[idx],
    Sign = ifelse(adj[idx] > 0, "Positive", "Negative"), stringsAsFactors = FALSE
  )
  nodes <- roles %>%
    mutate(
      PairID = pair_id, StudyID = study_id, Domain = domain, Group = group, GroupRole = group_role,
      Degree = as.numeric(deg[Taxon])
    )

  metrics <- data.frame(
    PairID = pair_id, StudyID = study_id, Domain = domain, Group = group, GroupRole = group_role,
    Samples = nrow(sample_by_taxa), TaxaAfterRAFilter = ncol(net$filtered),
    Nodes = nrow(adj), Edges = sum(upper != 0), PositiveEdges = sum(upper > 0), NegativeEdges = sum(upper < 0),
    Connectance = sum(upper != 0) / (nrow(adj) * (nrow(adj) - 1) / 2),
    MeanDegree = mean(deg),
    Modularity = modularity,
    NaturalConnectivity = normalized_natural_connectivity(adj_binary),
    BetweennessNCSlope = betweenness_nc_slope(adj_binary),
    Random50Remain = random_robustness(adj, random_remove_fraction, nperm = 100),
    Targeted4Remain = target$value,
    TargetedRemovedNodes = target$removed,
    TargetType = target$target,
    PositiveCohesion = cohesion$PositiveCohesion,
    NegativeCohesion = cohesion$NegativeCohion %||% cohesion$NegativeCohesion,
    AbsNegPosCohesionRatio = cohesion$AbsNegPosCohesionRatio,
    ModuleHubs = sum(roles$Role == "Module hub"),
    Connectors = sum(roles$Role == "Connector"),
    NetworkHubs = sum(roles$Role == "Network hub"),
    Status = "OK",
    stringsAsFactors = FALSE
  )

  list(metrics = metrics, nodes = nodes, edges = edges)
}

paired <- read.csv(paired_path, check.names = FALSE)
paired <- paired[paired$ControlSamples >= 6 & paired$TreatmentSamples >= 6, , drop = FALSE]
paired$PairID <- paste(paired$StudyID, paired$Domain, paired$ControlGroup, paired$TreatmentGroup, sep = "__")

metric_rows <- list()
node_rows <- list()
edge_rows <- list()
skipped <- list()

for (i in seq_len(nrow(paired))) {
  row <- paired[i, ]
  pair_data <- read_domain_pair_data(row$StudyID, row$Domain)
  if (is.null(pair_data)) {
    skipped[[length(skipped) + 1]] <- data.frame(PairID = row$PairID, Reason = "MissingOTUOrMetadata")
    next
  }
  ctrl_otu <- subset_group_otu(pair_data, row$ControlGroup)
  trt_otu <- subset_group_otu(pair_data, row$TreatmentGroup)
  if (is.null(ctrl_otu) || is.null(trt_otu)) {
    skipped[[length(skipped) + 1]] <- data.frame(PairID = row$PairID, Reason = "MissingGroupSamples")
    next
  }

  ctrl <- network_metrics(ctrl_otu, row$StudyID, row$Domain, row$ControlGroup, "Control", row$PairID)
  trt <- network_metrics(trt_otu, row$StudyID, row$Domain, row$TreatmentGroup, "Treatment", row$PairID)

  metric_rows[[length(metric_rows) + 1]] <- ctrl$metrics
  metric_rows[[length(metric_rows) + 1]] <- trt$metrics
  if (nrow(ctrl$nodes) > 0) node_rows[[length(node_rows) + 1]] <- ctrl$nodes
  if (nrow(trt$nodes) > 0) node_rows[[length(node_rows) + 1]] <- trt$nodes
  if (nrow(ctrl$edges) > 0) edge_rows[[length(edge_rows) + 1]] <- ctrl$edges
  if (nrow(trt$edges) > 0) edge_rows[[length(edge_rows) + 1]] <- trt$edges
}

metrics <- if (length(metric_rows) == 0) data.frame() else bind_rows(metric_rows)
nodes <- if (length(node_rows) == 0) data.frame() else bind_rows(node_rows)
edges <- if (length(edge_rows) == 0) data.frame() else bind_rows(edge_rows)
skipped_df <- if (length(skipped) == 0) data.frame() else bind_rows(skipped)

write_csv_utf8(metrics, file.path(network_dir, "menap_style_local_network_metrics.csv"))
write_csv_utf8(nodes, file.path(network_dir, "menap_style_local_zipi_nodes.csv"))
write_csv_utf8(edges, file.path(network_dir, "menap_style_local_edges.csv"))
write_csv_utf8(skipped_df, file.path(network_dir, "menap_style_local_skipped_pairs.csv"))

metric_names <- c(
  "Nodes", "Edges", "PositiveEdges", "NegativeEdges", "Connectance", "MeanDegree", "Modularity",
  "NaturalConnectivity", "BetweennessNCSlope", "Random50Remain", "Targeted4Remain",
  "PositiveCohesion", "NegativeCohesion", "AbsNegPosCohesionRatio"
)

delta_rows <- list()
if (nrow(metrics) > 0) {
  for (pid in unique(metrics$PairID)) {
    x <- metrics[metrics$PairID == pid & metrics$Status == "OK", , drop = FALSE]
    c_row <- x[x$GroupRole == "Control", , drop = FALSE]
    t_row <- x[x$GroupRole == "Treatment", , drop = FALSE]
    if (nrow(c_row) != 1 || nrow(t_row) != 1) next
    out <- data.frame(
      PairID = pid, StudyID = t_row$StudyID, Domain = t_row$Domain,
      ControlGroup = c_row$Group, TreatmentGroup = t_row$Group,
      stringsAsFactors = FALSE
    )
    for (m in metric_names) {
      out[[paste0("Control_", m)]] <- suppressWarnings(as.numeric(c_row[[m]]))
      out[[paste0("Treatment_", m)]] <- suppressWarnings(as.numeric(t_row[[m]]))
      out[[paste0("Delta_", m)]] <- out[[paste0("Treatment_", m)]] - out[[paste0("Control_", m)]]
    }
    delta_rows[[length(delta_rows) + 1]] <- out
  }
}
delta <- if (length(delta_rows) == 0) data.frame() else bind_rows(delta_rows)
write_csv_utf8(delta, file.path(network_dir, "menap_style_local_pair_deltas.csv"))

test_rows <- list()
if (nrow(delta) > 0) {
  for (domain in sort(unique(delta$Domain))) {
    d <- delta[delta$Domain == domain, , drop = FALSE]
    for (m in metric_names) {
      c_col <- paste0("Control_", m)
      t_col <- paste0("Treatment_", m)
      ok <- is.finite(d[[c_col]]) & is.finite(d[[t_col]])
      p <- NA_real_
      stat <- NA_real_
      if (sum(ok) >= 3) {
        wt <- tryCatch(wilcox.test(d[[t_col]][ok], d[[c_col]][ok], paired = TRUE, exact = FALSE), error = function(e) NULL)
        if (!is.null(wt)) {
          p <- wt$p.value
          stat <- unname(wt$statistic)
        }
      }
      test_rows[[length(test_rows) + 1]] <- data.frame(
        Domain = domain, Metric = m, Npairs = sum(ok), WilcoxonV = stat, Pvalue = p,
        MeanControl = if (sum(ok) > 0) mean(d[[c_col]][ok]) else NA_real_,
        MeanTreatment = if (sum(ok) > 0) mean(d[[t_col]][ok]) else NA_real_,
        MeanDelta = if (sum(ok) > 0) mean(d[[paste0("Delta_", m)]][ok]) else NA_real_
      )
    }
  }
}
tests <- if (length(test_rows) == 0) data.frame() else bind_rows(test_rows)
if (nrow(tests) > 0) {
  tests$FDR_BH <- ave(tests$Pvalue, tests$Domain, FUN = function(p) p.adjust(p, method = "BH"))
}
write_csv_utf8(tests, file.path(network_dir, "menap_style_local_wilcoxon_tests.csv"))

method_comparison <- data.frame(
  MethodElement = c(
    "ASV filtering", "Correlation", "Threshold", "Network construction", "Module detection",
    "Zi-Pi node role", "Random robustness", "Targeted robustness", "Natural connectivity", "Cohesion",
    "Statistical test"
  ),
  NEE_Method = c(
    "ASVs with relative abundance >0.05%",
    "Pearson correlation",
    "RMT threshold in MENAP",
    "MENAP web pipeline",
    "Greedy modularity optimization",
    "Zi >2.5 and Pi cutoffs at 0.62",
    "Remaining species after 50% random node removal",
    "Remaining species after removing four module hubs",
    "Average eigenvalue from network spectrum; sequential betweenness removal",
    "Observed-minus-expected positive/negative cohesion",
    "Two-tailed Wilcoxon rank-sum/rank tests + BH FDR"
  ),
  LocalImplementation = c(
    paste0("Mean relative abundance > ", min_mean_ra * 100, "% within each group"),
    "Pearson on log1p(relative abundance * 10000)",
    paste0("Fixed |r| >= ", cor_cutoff, "; this is not a true RMT threshold"),
    "Local R/igraph implementation",
    "igraph fast-greedy modularity",
    "Same Zi/Pi cutoffs; module hubs, network hubs and connectors exported",
    paste0("Mean remaining proportion after ", random_remove_fraction * 100, "% random node removal, 100 simulations"),
    paste0("Remove up to ", target_remove_n, " Zi-Pi keystone nodes; if none exist, betweenness fallback is flagged"),
    "Normalized natural connectivity and betweenness attack slope",
    paste0("Permutation-based observed-minus-expected cohesion, nperm=", nperm),
    "Paired Wilcoxon signed-rank tests within Domain + BH FDR"
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(method_comparison, file.path(network_dir, "menap_style_method_comparison.csv"))

plot_metrics <- c("Random50Remain", "Targeted4Remain", "NaturalConnectivity", "AbsNegPosCohesionRatio")
long_rows <- list()
for (m in plot_metrics) {
  c_col <- paste0("Control_", m)
  t_col <- paste0("Treatment_", m)
  if (!(c_col %in% names(delta)) || !(t_col %in% names(delta))) next
  long_rows[[length(long_rows) + 1]] <- data.frame(
    PairID = delta$PairID, StudyID = delta$StudyID, Domain = delta$Domain,
    Metric = m, GroupRole = "Control", Value = delta[[c_col]], stringsAsFactors = FALSE
  )
  long_rows[[length(long_rows) + 1]] <- data.frame(
    PairID = delta$PairID, StudyID = delta$StudyID, Domain = delta$Domain,
    Metric = m, GroupRole = "Treatment", Value = delta[[t_col]], stringsAsFactors = FALSE
  )
}
plot_df <- if (length(long_rows) == 0) data.frame() else bind_rows(long_rows)
metric_labels <- c(
  Random50Remain = "Random removal\nrobustness",
  Targeted4Remain = "Targeted removal\nrobustness",
  NaturalConnectivity = "Natural\nconnectivity",
  AbsNegPosCohesionRatio = "|negative| / positive\ncohesion"
)
if (nrow(plot_df) > 0) {
  plot_df$MetricLabel <- metric_labels[plot_df$Metric]
  plot_df$GroupRole <- factor(plot_df$GroupRole, levels = c("Control", "Treatment"))
  cols <- c(Control = "#F2B43F", Treatment = "#4E9AD1")
  p <- ggplot(plot_df, aes(GroupRole, Value, fill = GroupRole)) +
    geom_line(aes(group = PairID), color = "grey65", linewidth = 0.22, alpha = 0.55) +
    geom_point(shape = 21, size = 1.8, color = "black", stroke = 0.18) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.55, fatten = 0, color = "black", linewidth = 0.25) +
    facet_grid(Domain ~ MetricLabel, scales = "free_y") +
    scale_fill_manual(values = cols) +
    labs(
      title = "MENAP-style local paired network stability metrics",
      subtitle = "Local approximation: Pearson, mean RA >0.05%, fixed correlation cutoff; paired CK-vs-treatment networks",
      x = NULL,
      y = "Metric value"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold")
    )
  png_path <- file.path(figure_dir, "fig_menap_style_local_paired_metrics.png")
  pdf_path <- file.path(figure_dir, "fig_menap_style_local_paired_metrics.pdf")
  png(png_path, width = 8.4, height = 4.6, units = "in", res = 450)
  print(p)
  dev.off()
  pdf(pdf_path, width = 8.4, height = 4.6, family = "sans", useDingbats = FALSE)
  print(p)
  dev.off()
}

readme <- c(
  "MENAP-style local network analysis",
  "",
  "This script was written after comparing the user's NEE methods paragraph with the current local network workflow.",
  "",
  "Important: this is a local approximation, not an actual MENAP/RMT run.",
  "For exact wording identical to the NEE method, the correlation matrices should be uploaded to MENAP and RMT thresholds/network properties should be exported from MENAP.",
  "",
  "Main choices:",
  paste0("- ASV filter: mean relative abundance > ", min_mean_ra * 100, "% within each group."),
  "- Correlation: Pearson on log1p(relative abundance * 10000).",
  paste0("- Local threshold: |r| >= ", cor_cutoff, " (RMT replacement, not true RMT)."),
  paste0("- Cohesion permutations: ", nperm, "."),
  paste0("- Random robustness: remaining taxa after ", random_remove_fraction * 100, "% random node removal."),
  paste0("- Targeted robustness: up to ", target_remove_n, " Zi-Pi keystone nodes; betweenness fallback is flagged when Zi-Pi keystone nodes are absent."),
  "",
  "Outputs:",
  file.path(network_dir, "menap_style_local_network_metrics.csv"),
  file.path(network_dir, "menap_style_local_pair_deltas.csv"),
  file.path(network_dir, "menap_style_local_wilcoxon_tests.csv"),
  file.path(network_dir, "menap_style_local_zipi_nodes.csv"),
  file.path(network_dir, "menap_style_local_edges.csv"),
  file.path(network_dir, "menap_style_method_comparison.csv"),
  file.path(figure_dir, "fig_menap_style_local_paired_metrics.png")
)
writeLines(readme, file.path(network_dir, "menap_style_local_README.txt"), useBytes = TRUE)

message("MENAP-style local metrics rows: ", nrow(metrics))
message("OK networks: ", sum(metrics$Status == "OK", na.rm = TRUE))
message("Pair delta rows: ", nrow(delta))
message("Wilcoxon rows: ", nrow(tests))
message("Output dir: ", network_dir)


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/build_paired_network_delta.R
# ================================================================

options(stringsAsFactors = FALSE)

pipeline_dir <- Sys.getenv(
  "PIPELINE_DIR",
  unset = "E:/BaiduSyncdisk/N_deposition1/Data/Global_Nitrogen_Pipeline"
)

network_dir <- file.path(pipeline_dir, "output", "04_network_groups")
metrics_path <- file.path(network_dir, "network_metrics_over5.csv")

if (!file.exists(metrics_path)) {
  stop("Missing network metrics file: ", metrics_path, call. = FALSE)
}

metrics <- read.csv(metrics_path, check.names = FALSE)

required <- c("StudyID", "Domain", "Group", "GroupRole", "Samples", "Status")
missing_required <- setdiff(required, names(metrics))
if (length(missing_required) > 0) {
  stop(
    "network_metrics_over5.csv is missing required columns: ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

metrics <- metrics[metrics$Status == "OK" & metrics$Samples >= 6, , drop = FALSE]
controls <- metrics[metrics$GroupRole == "Control", , drop = FALSE]
treatments <- metrics[metrics$GroupRole == "Treatment", , drop = FALSE]

key <- function(x) paste(x$StudyID, x$Domain, sep = "\r")
control_keys <- key(controls)
treatment_keys <- key(treatments)
paired_keys <- intersect(unique(control_keys), unique(treatment_keys))

metric_cols <- c(
  "Taxa", "Nodes", "Edges", "PositiveEdges", "NegativeEdges",
  "Connectance", "MeanDegree", "Modularity", "NaturalConnectivity",
  "RobustnessSlope"
)
metric_cols <- intersect(metric_cols, names(metrics))

ratio_metric_cols <- intersect(
  c("Taxa", "Nodes", "Edges", "PositiveEdges", "NegativeEdges",
    "Connectance", "MeanDegree", "NaturalConnectivity"),
  metric_cols
)

safe_num <- function(x) suppressWarnings(as.numeric(x))

make_pair_row <- function(ctrl, trt) {
  out <- data.frame(
    StudyID = trt$StudyID,
    Domain = trt$Domain,
    ControlGroup = ctrl$Group,
    TreatmentGroup = trt$Group,
    ControlSamples = safe_num(ctrl$Samples),
    TreatmentSamples = safe_num(trt$Samples),
    stringsAsFactors = FALSE
  )

  for (metric in metric_cols) {
    cval <- safe_num(ctrl[[metric]])
    tval <- safe_num(trt[[metric]])
    out[[paste0("Control_", metric)]] <- cval
    out[[paste0("Treatment_", metric)]] <- tval
    out[[paste0("Delta_", metric)]] <- tval - cval

    if (metric %in% ratio_metric_cols) {
      out[[paste0("LogRatio_", metric)]] <- ifelse(
        is.finite(cval) & is.finite(tval) & cval > 0 & tval > 0,
        log(tval / cval),
        NA_real_
      )
    }
  }

  if ("NaturalConnectivity" %in% metric_cols) {
    d <- out$Delta_NaturalConnectivity
    out$NaturalConnectivityChange <- ifelse(
      is.na(d), NA_character_,
      ifelse(d > 0, "HigherThanControl", ifelse(d < 0, "LowerThanControl", "NoChange"))
    )
  }

  if ("RobustnessSlope" %in% metric_cols) {
    d <- out$Delta_RobustnessSlope
    out$RobustnessSlopeChange <- ifelse(
      is.na(d), NA_character_,
      ifelse(d > 0, "LessNegativeThanControl", ifelse(d < 0, "MoreNegativeThanControl", "NoChange"))
    )
  }

  out
}

pair_rows <- list()
for (one_key in paired_keys) {
  ctrl_rows <- controls[control_keys == one_key, , drop = FALSE]
  trt_rows <- treatments[treatment_keys == one_key, , drop = FALSE]

  if (nrow(ctrl_rows) == 0 || nrow(trt_rows) == 0) {
    next
  }

  # Current metadata has one CK-like control per StudyID/Domain. If a future
  # dataset has multiple controls, keep all pairwise combinations visible.
  for (i in seq_len(nrow(ctrl_rows))) {
    for (j in seq_len(nrow(trt_rows))) {
      pair_rows[[length(pair_rows) + 1]] <- make_pair_row(ctrl_rows[i, , drop = FALSE], trt_rows[j, , drop = FALSE])
    }
  }
}

if (length(pair_rows) == 0) {
  paired_delta <- data.frame()
} else {
  paired_delta <- do.call(rbind, pair_rows)
  paired_delta <- paired_delta[order(paired_delta$Domain, paired_delta$StudyID, paired_delta$TreatmentGroup), , drop = FALSE]
}

delta_path <- file.path(network_dir, "paired_control_treatment_network_delta_metrics.csv")
write.csv(paired_delta, delta_path, row.names = FALSE, fileEncoding = "UTF-8")

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

summary_rows <- list()
if (nrow(paired_delta) > 0) {
  for (domain in sort(unique(paired_delta$Domain))) {
    x <- paired_delta[paired_delta$Domain == domain, , drop = FALSE]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Domain = domain,
      PairedComparisons = nrow(x),
      StudyDomainPairs = length(unique(paste(x$StudyID, x$Domain, sep = "\r"))),
      MeanDelta_Edges = if ("Delta_Edges" %in% names(x)) mean_or_na(x$Delta_Edges) else NA_real_,
      MeanDelta_Connectance = if ("Delta_Connectance" %in% names(x)) mean_or_na(x$Delta_Connectance) else NA_real_,
      MeanDelta_MeanDegree = if ("Delta_MeanDegree" %in% names(x)) mean_or_na(x$Delta_MeanDegree) else NA_real_,
      MeanDelta_NaturalConnectivity = if ("Delta_NaturalConnectivity" %in% names(x)) mean_or_na(x$Delta_NaturalConnectivity) else NA_real_,
      MeanDelta_RobustnessSlope = if ("Delta_RobustnessSlope" %in% names(x)) mean_or_na(x$Delta_RobustnessSlope) else NA_real_,
      LowerNaturalConnectivity = if ("NaturalConnectivityChange" %in% names(x)) sum(x$NaturalConnectivityChange == "LowerThanControl", na.rm = TRUE) else NA_integer_,
      HigherNaturalConnectivity = if ("NaturalConnectivityChange" %in% names(x)) sum(x$NaturalConnectivityChange == "HigherThanControl", na.rm = TRUE) else NA_integer_,
      MoreNegativeRobustnessSlope = if ("RobustnessSlopeChange" %in% names(x)) sum(x$RobustnessSlopeChange == "MoreNegativeThanControl", na.rm = TRUE) else NA_integer_,
      LessNegativeRobustnessSlope = if ("RobustnessSlopeChange" %in% names(x)) sum(x$RobustnessSlopeChange == "LessNegativeThanControl", na.rm = TRUE) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}

summary_delta <- if (length(summary_rows) == 0) data.frame() else do.call(rbind, summary_rows)
summary_path <- file.path(network_dir, "paired_control_treatment_network_delta_summary_by_domain.csv")
write.csv(summary_delta, summary_path, row.names = FALSE, fileEncoding = "UTF-8")

notes_path <- file.path(network_dir, "paired_control_treatment_network_delta_README.txt")
notes <- c(
  "Paired CK-vs-treatment network delta tables",
  "",
  "Input:",
  paste0("- ", metrics_path),
  "",
  "Eligibility:",
  "- Network Status == OK.",
  "- Samples >= 6 in both the control group and the treatment group.",
  "- Control and treatment are paired only within the same StudyID and Domain.",
  "",
  "Interpretation:",
  "- Each row in paired_control_treatment_network_delta_metrics.csv is one treatment network compared with its CK/control network.",
  "- Delta_* columns are Treatment - Control.",
  "- LogRatio_* columns are log(Treatment / Control) and are NA when either value is missing or <= 0.",
  "- Positive Delta_NaturalConnectivity means the treatment network has higher natural connectivity than CK.",
  "- For RobustnessSlope, larger/less negative values are generally interpreted as slower decline under node removal.",
  "",
  "Outputs:",
  paste0("- ", delta_path),
  paste0("- ", summary_path)
)
writeLines(notes, notes_path, useBytes = TRUE)

message("Wrote: ", delta_path)
message("Wrote: ", summary_path)
message("Wrote: ", notes_path)
message("Paired comparisons: ", nrow(paired_delta))
if (nrow(summary_delta) > 0) {
  message("Domains: ", paste(summary_delta$Domain, collapse = ", "))
}


# ================================================================
# Source: Data/Global_Nitrogen_Pipeline/scripts/plot_paired_network_delta.R
# ================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

pipeline_dir <- Sys.getenv(
  "PIPELINE_DIR",
  unset = "E:/BaiduSyncdisk/N_deposition1/Data/Global_Nitrogen_Pipeline"
)

network_dir <- file.path(pipeline_dir, "output", "04_network_groups")
figure_dir <- file.path(network_dir, "figures")
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

delta_path <- file.path(network_dir, "paired_control_treatment_network_delta_metrics.csv")
if (!file.exists(delta_path)) {
  stop("Missing paired delta table: ", delta_path, call. = FALSE)
}

dat <- read.csv(delta_path, check.names = FALSE)
if (nrow(dat) == 0) {
  stop("Paired delta table is empty: ", delta_path, call. = FALSE)
}

num_cols <- c(
  "ControlSamples", "TreatmentSamples", "Delta_Edges",
  "Delta_Connectance", "Delta_MeanDegree",
  "Delta_NaturalConnectivity", "Delta_RobustnessSlope"
)
for (col in intersect(num_cols, names(dat))) {
  dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
}

dat <- dat %>%
  mutate(
    PairLabel = paste0(StudyID, " / ", Domain, " / ", TreatmentGroup),
    ComparisonID = paste0(StudyID, " | ", TreatmentGroup),
    Domain = factor(Domain, levels = c("Bacterial", "Fungi"))
  )

domain_cols <- c(Bacterial = "#0072B2", Fungi = "#D55E00")
direction_cols <- c(
  "Lower than CK" = "#0072B2",
  "Higher than CK" = "#D55E00",
  "No change" = "grey55"
)

base_theme <- theme_classic(base_size = 9) +
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
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))

  png(filename = png_path, width = width, height = height, units = "in", res = 450)
  print(plot)
  dev.off()

  cairo_pdf(filename = pdf_path, width = width, height = height, family = "sans")
  print(plot)
  dev.off()

  c(png = png_path, pdf = pdf_path)
}

# Figure 1: every eligible CK-vs-treatment comparison, sorted by natural connectivity change.
fig1_dat <- dat %>%
  arrange(Delta_NaturalConnectivity) %>%
  mutate(PairLabel = factor(PairLabel, levels = unique(PairLabel)))

p1 <- ggplot(fig1_dat, aes(x = Delta_NaturalConnectivity, y = PairLabel, color = Domain)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey35") +
  geom_segment(aes(x = 0, xend = Delta_NaturalConnectivity, yend = PairLabel),
               linewidth = 0.35, alpha = 0.65) +
  geom_point(size = 1.8) +
  scale_color_manual(values = domain_cols, drop = FALSE) +
  labs(
    title = "Network natural connectivity response to nitrogen addition",
    subtitle = "Each row is one treatment network compared with its CK network; values are treatment minus CK",
    x = "Delta natural connectivity",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 5.2),
    legend.position = "top"
  )

paths1 <- save_plot(p1, "fig_network_delta_natural_connectivity_lollipop", width = 7.2, height = 9.2)

# Figure 2: distribution of major network stability indicators by domain.
long_dat <- bind_rows(
  dat %>% transmute(Domain, StudyID, TreatmentGroup, Metric = "Edges", Delta = Delta_Edges),
  dat %>% transmute(Domain, StudyID, TreatmentGroup, Metric = "Natural connectivity", Delta = Delta_NaturalConnectivity),
  dat %>% transmute(Domain, StudyID, TreatmentGroup, Metric = "Robustness slope", Delta = Delta_RobustnessSlope)
) %>%
  mutate(
    Metric = factor(Metric, levels = c("Edges", "Natural connectivity", "Robustness slope"))
  )

p2 <- ggplot(long_dat, aes(x = Domain, y = Delta, fill = Domain, color = Domain)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey35") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 20260531),
             size = 1.5, alpha = 0.75) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = domain_cols, drop = FALSE) +
  scale_color_manual(values = domain_cols, drop = FALSE) +
  labs(
    title = "Paired CK-vs-treatment network metric changes",
    subtitle = "Delta values are treatment minus CK; each point is one eligible treatment group",
    x = NULL,
    y = "Delta metric"
  ) +
  base_theme +
  theme(legend.position = "none")

paths2 <- save_plot(p2, "fig_network_delta_metric_distributions", width = 7.2, height = 3.0)

# Figure 3: direction counts for natural connectivity and robustness slope.
count_dat <- bind_rows(
  dat %>%
    mutate(
      Metric = "Natural connectivity",
      Direction = case_when(
        Delta_NaturalConnectivity < 0 ~ "Lower than CK",
        Delta_NaturalConnectivity > 0 ~ "Higher than CK",
        TRUE ~ "No change"
      )
    ) %>%
    count(Domain, Metric, Direction, name = "Comparisons"),
  dat %>%
    mutate(
      Metric = "Robustness slope",
      Direction = case_when(
        Delta_RobustnessSlope < 0 ~ "Lower than CK",
        Delta_RobustnessSlope > 0 ~ "Higher than CK",
        TRUE ~ "No change"
      )
    ) %>%
    count(Domain, Metric, Direction, name = "Comparisons")
) %>%
  mutate(
    Metric = factor(Metric, levels = c("Natural connectivity", "Robustness slope")),
    Direction = factor(Direction, levels = c("Lower than CK", "Higher than CK", "No change"))
  )

p3 <- ggplot(count_dat, aes(x = Domain, y = Comparisons, fill = Direction)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  facet_wrap(~ Metric, nrow = 1) +
  scale_fill_manual(values = direction_cols, drop = FALSE) +
  labs(
    title = "Direction of network stability changes",
    subtitle = "Counts are eligible CK-vs-treatment comparisons",
    x = NULL,
    y = "Number of comparisons"
  ) +
  base_theme

paths3 <- save_plot(p3, "fig_network_delta_direction_counts", width = 6.4, height = 3.0)

figure_index <- data.frame(
  Figure = c(
    "Natural connectivity lollipop",
    "Metric distributions",
    "Direction counts"
  ),
  PNG = c(paths1[["png"]], paths2[["png"]], paths3[["png"]]),
  PDF = c(paths1[["pdf"]], paths2[["pdf"]], paths3[["pdf"]]),
  stringsAsFactors = FALSE
)

index_path <- file.path(figure_dir, "network_delta_figure_index.csv")
write.csv(figure_index, index_path, row.names = FALSE, fileEncoding = "UTF-8")

message("Wrote figure index: ", index_path)
for (i in seq_len(nrow(figure_index))) {
  message("PNG: ", figure_index$PNG[i])
  message("PDF: ", figure_index$PDF[i])
}


# ================================================================
# Source: Data_analysis/Remote_working_package/07_Network/plot_network_paired_metrics_pretty.R
# ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gghalves)

net_dir <- "E:/BaiduSyncdisk/N_deposition1/Data_analysis/outputs/network_analysis"
fig_dir <- file.path(net_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

wide_path <- file.path(net_dir, "network_metrics_paired_wide.csv")
stats_path <- file.path(net_dir, "network_paired_statistics.csv")

core_metrics <- data.frame(
  Metric = c(
    "natural_connectivity",
    "random_robustness",
    "targeted_robustness",
    "cohesion_ratio_abs_neg_pos"
  ),
  MetricTitle = c(
    "natural connectivity",
    "random robustness",
    "targeted robustness",
    "negative/positive cohesion"
  ),
  stringsAsFactors = FALSE
)

plot_colors <- c("CK" = "#A8A8A8", "N addition" = "#F2B84B")
plot_point_edge <- c("CK" = "#666666", "N addition" = "#B88016")

p_text <- function(p) {
  ifelse(is.na(p), "P = NA", ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p)))
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))
  }
  m <- mean(x)
  se <- if (length(x) > 1) sd(x) / sqrt(length(x)) else 0
  data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}

build_plot_data <- function(wide, stats, domains = c("Bacterial", "Fungi")) {
  set.seed(20260612)

  value_rows <- list()
  summary_rows <- list()
  label_rows <- list()
  panel_i <- 1

  for (domain_i in domains) {
    dat_domain <- wide %>% filter(Domain == domain_i)
    stat_domain <- stats %>% filter(Domain == domain_i)
    domain_lower <- ifelse(domain_i == "Bacterial", "bacterial", "fungal")

    for (metric_i in core_metrics$Metric) {
      metric_title <- core_metrics$MetricTitle[core_metrics$Metric == metric_i]
      panel_letter <- letters[panel_i]
      panel_label <- paste0(panel_letter, "\nPaired response of\n", domain_lower, " ", metric_title)
      ck_col <- paste0(metric_i, "_CK")
      n_col <- paste0(metric_i, "_N")

      vals <- dat_domain %>%
        select(PairID, all_of(ck_col), all_of(n_col)) %>%
        rename(CK = all_of(ck_col), `N addition` = all_of(n_col)) %>%
        mutate(across(c(CK, `N addition`), as.numeric)) %>%
        filter(!is.na(CK), !is.na(`N addition`))

      vals$CK_X <- runif(nrow(vals), 1.07, 1.19)
      vals$N_X <- runif(nrow(vals), 1.81, 1.93)

      long_vals <- vals %>%
        pivot_longer(cols = c(CK, `N addition`), names_to = "TreatmentClass", values_to = "Value") %>%
        mutate(
          Domain = domain_i,
          Metric = metric_i,
          Panel = panel_label,
          X = ifelse(TreatmentClass == "CK", CK_X, N_X),
          TreatmentClass = factor(TreatmentClass, levels = c("CK", "N addition"))
        )

      value_rows[[length(value_rows) + 1]] <- long_vals

      for (treatment_i in c("CK", "N addition")) {
        x <- vals[[treatment_i]]
        ci <- ci95(x)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Domain = domain_i,
          Metric = metric_i,
          Panel = panel_label,
          TreatmentClass = treatment_i,
          X = ifelse(treatment_i == "CK", 1, 2),
          y = ci$y,
          ymin = ci$ymin,
          ymax = ci$ymax
        )
      }

      stat_row <- stat_domain %>% filter(Metric == metric_i)
      if (nrow(stat_row) > 0) {
        label <- paste0(
          p_text(stat_row$Wilcoxon_P[1]), "\n",
          "FDR = ", sprintf("%.3f", stat_row$Wilcoxon_FDR[1]), "\n",
          "n = ", stat_row$Pairs[1]
        )
      } else {
        label <- paste0("n = ", nrow(vals))
      }

      label_rows[[length(label_rows) + 1]] <- data.frame(
        Domain = domain_i,
        Metric = metric_i,
        Panel = panel_label,
        X = 0.78,
        Y = Inf,
        Label = label
      )

      panel_i <- panel_i + 1
    }
  }

  list(
    values = bind_rows(value_rows),
    summary = bind_rows(summary_rows),
    labels = bind_rows(label_rows)
  )
}

draw_network_plot <- function(wide, stats, domains, out_prefix, width, height) {
  plot_data <- build_plot_data(wide, stats, domains)

  ncol_use <- ifelse(length(domains) == 1, 4, 4)

  violin_data <- plot_data$values %>%
    mutate(
      TreatmentX = ifelse(TreatmentClass == "CK", 1, 2),
      PointX = X
    )

  make_one_panel <- function(panel_name) {
    panel_values <- violin_data %>% filter(Panel == panel_name)
    panel_label <- plot_data$labels %>% filter(Panel == panel_name)

    ggplot() +
      geom_half_violin(
        data = panel_values %>% filter(TreatmentClass == "CK"),
        aes(x = TreatmentX, y = Value, fill = TreatmentClass),
        side = "l",
        color = NA,
        alpha = 0.35,
        trim = FALSE,
        width = 0.72
      ) +
      geom_half_violin(
        data = panel_values %>% filter(TreatmentClass == "N addition"),
        aes(x = TreatmentX, y = Value, fill = TreatmentClass),
        side = "r",
        color = NA,
        alpha = 0.35,
        trim = FALSE,
        width = 0.72
      ) +
      geom_line(
        data = panel_values,
        aes(x = PointX, y = Value, group = PairID),
        color = "#C7B8B1",
        linewidth = 0.36,
        alpha = 0.42
      ) +
      geom_point(
        data = panel_values,
        aes(x = PointX, y = Value, fill = TreatmentClass, color = TreatmentClass),
        shape = 21,
        size = 2.0,
        stroke = 0.55,
        alpha = 0.88
      ) +
      geom_half_boxplot(
        data = panel_values %>% filter(TreatmentClass == "CK"),
        aes(x = TreatmentX, y = Value, fill = TreatmentClass, group = TreatmentClass),
        side = "l",
        errorbar.draw = FALSE,
        width = 0.20,
        outlier.shape = NA,
        alpha = 0.78,
        linewidth = 0.62
      ) +
      geom_half_boxplot(
        data = panel_values %>% filter(TreatmentClass == "N addition"),
        aes(x = TreatmentX, y = Value, fill = TreatmentClass, group = TreatmentClass),
        side = "r",
        errorbar.draw = FALSE,
        width = 0.20,
        outlier.shape = NA,
        alpha = 0.78,
        linewidth = 0.62
      ) +
      geom_text(
        data = panel_label,
        aes(x = X, y = Y, label = Label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1.12,
        size = 2.8,
        family = "Times New Roman"
      ) +
      scale_color_manual(values = plot_point_edge, breaks = c("CK", "N addition")) +
      scale_fill_manual(values = plot_colors, breaks = c("CK", "N addition")) +
      scale_x_continuous(
        breaks = c(1, 2),
        labels = c("CK", "N addition"),
        limits = c(0.62, 2.38),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(
        title = panel_name,
        x = NULL,
        y = "Network metric value",
        color = NULL,
        fill = NULL
      ) +
      theme_classic(base_family = "Times New Roman", base_size = 10) +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 11, lineheight = 0.95, hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.55),
        axis.ticks = element_line(color = "black", linewidth = 0.45),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.4),
        panel.grid.major.x = element_line(color = "#EFEFEF", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.margin = margin(8, 8, 8, 8)
      ) +
      guides(
        color = "none",
        fill = guide_legend(override.aes = list(size = 4, alpha = 1))
      )
  }

  panel_order <- unique(plot_data$labels$Panel)
  panel_plots <- lapply(panel_order, make_one_panel)
  p <- wrap_plots(panel_plots, ncol = ncol_use, guides = "collect") &
    theme(legend.position = "top")

  png_path <- file.path(fig_dir, paste0(out_prefix, ".png"))
  pdf_path <- file.path(fig_dir, paste0(out_prefix, ".pdf"))
  ggsave(png_path, p, width = width, height = height, dpi = 300)
  ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf)
  message("Saved: ", png_path)
  message("Saved: ", pdf_path)
}

main <- function() {
  wide <- read.csv(wide_path, check.names = FALSE, stringsAsFactors = FALSE)
  stats <- read.csv(stats_path, check.names = FALSE, stringsAsFactors = FALSE)

  draw_network_plot(
    wide,
    stats,
    domains = c("Bacterial", "Fungi"),
    out_prefix = "network_paired_metrics_pretty_Bacterial_Fungi",
    width = 13.2,
    height = 7.0
  )

  draw_network_plot(
    wide,
    stats,
    domains = "Bacterial",
    out_prefix = "network_paired_metrics_pretty_Bacterial",
    width = 12.6,
    height = 3.8
  )

  draw_network_plot(
    wide,
    stats,
    domains = "Fungi",
    out_prefix = "network_paired_metrics_pretty_Fungi",
    width = 12.6,
    height = 3.8
  )
}

main()


# ================================================================
# Source: Data_analysis/Remote_working_package/07_Network/plot_representative_networks_phylum_style.R
# ================================================================

rm(list = ls())
setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

required_packages <- c("dplyr", "igraph")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(igraph)

network_dir <- file.path("outputs", "network_analysis")
edge_dir <- file.path(network_dir, "edges")
node_dir <- file.path(network_dir, "nodes")
fig_dir <- file.path(network_dir, "figures")
gephi_dir <- file.path(network_dir, "gephi_phylum_style")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gephi_dir, recursive = TRUE, showWarnings = FALSE)

selected_path <- file.path(network_dir, "representative_networks_pretty_selected_pairs.csv")
selected_pairs <- read.csv(selected_path, stringsAsFactors = FALSE, check.names = FALSE)

phylum_palette <- c(
  Proteobacteria = "#B58CE5",
  Actinobacteria = "#59C36A",
  Acidobacteria = "#7B4B1D",
  Firmicutes = "#5BAED6",
  Bacteroidetes = "#DE6BAA",
  Verrucomicrobia = "#E99627",
  Planctomycetes = "#9DB51A",
  Chloroflexi = "#56B4B8",
  Gemmatimonadetes = "#F06B6B",
  Myxococcota = "#6A82C2",
  Ascomycota = "#B58CE5",
  Basidiomycota = "#DE6BAA",
  Mortierellomycota = "#59C36A",
  Chytridiomycota = "#5BAED6",
  Rozellomycota = "#E99627",
  Mucoromycota = "#7B4B1D",
  Glomeromycota = "#9DB51A",
  Kickxellomycota = "#56B4B8",
  Zoopagomycota = "#6A82C2",
  Other = "#BDBDBD",
  Unclassified = "#D9D9D9"
)

treatment_colors <- c("CK" = "#A8A8A8", "N" = "#F2B84B")

normalize_phylum <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unclassified"
  x <- gsub("^p__", "", x)

  dplyr::case_when(
    x %in% c("Pseudomonadota", "Proteobacteria") ~ "Proteobacteria",
    x %in% c("Actinomycetota", "Actinobacteriota", "Actinobacteria") ~ "Actinobacteria",
    x %in% c("Bacillota", "Firmicutes") ~ "Firmicutes",
    x %in% c("Bacteroidota", "Bacteroidetes") ~ "Bacteroidetes",
    x %in% c("Acidobacteriota", "Acidobacteria") ~ "Acidobacteria",
    x %in% c("Verrucomicrobiota", "Verrucomicrobia") ~ "Verrucomicrobia",
    x %in% c("Planctomycetota", "Planctomycetes") ~ "Planctomycetes",
    x %in% c("Chloroflexota", "Chloroflexi") ~ "Chloroflexi",
    x %in% c("Gemmatimonadota", "Gemmatimonadetes") ~ "Gemmatimonadetes",
    x %in% c("Myxococcota") ~ "Myxococcota",
    x %in% names(phylum_palette) ~ x,
    TRUE ~ "Other"
  )
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

build_panel_info <- function() {
  bacterial <- selected_pairs %>% filter(Domain == "Bacterial") %>% slice(1)
  fungi <- selected_pairs %>% filter(Domain == "Fungi") %>% slice(1)

  bind_rows(bacterial, fungi) %>%
    tidyr::expand_grid(NetworkGroup = c("CK", "N")) %>%
    mutate(
      TreatmentLabel = ifelse(NetworkGroup == "CK", "CK", "N addition"),
      PanelTitle = paste(Domain, TreatmentLabel),
      FileGroup = NetworkGroup
    )
}

make_network_graph <- function(pair_id, network_group) {
  edge_path <- file.path(edge_dir, paste0(pair_id, "_", network_group, "_edges.csv"))
  node_path <- file.path(node_dir, paste0(pair_id, "_", network_group, "_nodes_roles.csv"))

  edges <- safe_read_csv(edge_path) %>%
    mutate(
      source = as.character(source),
      target = as.character(target),
      r = suppressWarnings(as.numeric(r)),
      sign = ifelse(!is.na(r) & r < 0, "negative", "positive")
    )

  nodes <- safe_read_csv(node_path) %>%
    mutate(
      name = as.character(OTUID),
      degree = suppressWarnings(as.numeric(degree)),
      PhylumClean = normalize_phylum(Phylum)
    ) %>%
    arrange(desc(degree)) %>%
    distinct(name, .keep_all = TRUE) %>%
    select(name, everything())

  graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  graph <- simplify(
    graph,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(
      r = "first",
      p = "first",
      sign = "first",
      "ignore"
    )
  )

  degree_values <- suppressWarnings(as.numeric(V(graph)$degree))
  if (all(is.na(degree_values))) degree_values <- degree(graph)
  degree_values[is.na(degree_values)] <- 0

  V(graph)$size <- 1.8 + 6.8 * sqrt(degree_values / max(degree_values, 1, na.rm = TRUE))
  V(graph)$color <- phylum_palette[V(graph)$PhylumClean]
  V(graph)$color[is.na(V(graph)$color)] <- phylum_palette["Other"]
  V(graph)$frame.color <- grDevices::adjustcolor("white", alpha.f = 0.55)
  V(graph)$label <- NA

  edge_df <- as_data_frame(graph, what = "edges")
  from_phy <- V(graph)$PhylumClean[match(edge_df$from, V(graph)$name)]
  to_phy <- V(graph)$PhylumClean[match(edge_df$to, V(graph)$name)]

  edge_col <- rep("#CFCFCF", nrow(edge_df))
  same_phy <- !is.na(from_phy) & from_phy == to_phy & from_phy %in% names(phylum_palette)
  edge_col[same_phy] <- phylum_palette[from_phy[same_phy]]
  if ("sign" %in% names(edge_df)) {
    edge_alpha <- ifelse(!is.na(edge_df$sign) & edge_df$sign == "negative", 0.16, 0.24)
  } else {
    edge_alpha <- rep(0.22, nrow(edge_df))
  }
  E(graph)$color <- vapply(
    seq_along(edge_col),
    function(i) grDevices::adjustcolor(edge_col[i], alpha.f = edge_alpha[i]),
    character(1)
  )
  E(graph)$width <- 0.35
  E(graph)$curved <- 0.12

  graph
}

plot_one <- function(graph, title_text, seed = 20260614, show_legend = TRUE) {
  set.seed(seed)
  layout <- layout_with_fr(graph, niter = 999, grid = "nogrid")

  par(mar = c(0.15, 0.15, 2.35, ifelse(show_legend, 2.05, 0.15)))
  plot(
    graph,
    layout = layout,
    vertex.label = NA,
    vertex.size = V(graph)$size,
    vertex.color = V(graph)$color,
    vertex.frame.color = V(graph)$frame.color,
    edge.color = E(graph)$color,
    edge.width = E(graph)$width,
    edge.curved = E(graph)$curved,
    asp = 0
  )

  title(title_text, line = 0.55, cex.main = 1.18, font.main = 2)

  if (show_legend) {
    phy <- sort(table(V(graph)$PhylumClean), decreasing = TRUE)
    phy <- names(phy)[seq_len(min(8, length(phy)))]
    legend(
      "right",
      legend = phy,
      pch = 21,
      pt.bg = phylum_palette[phy],
      pt.cex = 1.25,
      col = "white",
      bty = "n",
      cex = 0.92,
      xpd = NA,
      title = "Phylum"
    )
    legend(
      "bottomright",
      legend = c("Low degree", "High degree"),
      pch = 21,
      pt.bg = "grey70",
      col = "white",
      pt.cex = c(0.85, 1.75),
      bty = "n",
      cex = 0.92,
      xpd = NA,
      title = "Degree"
    )
  }
}

export_gephi_tables <- function(graph, panel_row, out_stub) {
  edge_df <- as_data_frame(graph, what = "edges")
  node_df <- as_data_frame(graph, what = "vertices")

  nodes_out <- node_df %>%
    transmute(
      Id = name,
      Label = name,
      Domain = panel_row$Domain,
      Treatment = ifelse(panel_row$FileGroup == "CK", "CK", "N addition"),
      PairID = panel_row$PairID,
      Phylum = PhylumClean,
      Color = color,
      Size = round(size, 4),
      Degree = suppressWarnings(as.numeric(degree)),
      module = if ("module" %in% names(node_df)) module else NA,
      Zi = if ("Zi" %in% names(node_df)) Zi else NA,
      Pi = if ("Pi" %in% names(node_df)) Pi else NA,
      role = if ("role" %in% names(node_df)) role else NA
    )

  edge_from_phy <- nodes_out$Phylum[match(edge_df$from, nodes_out$Id)]
  edge_to_phy <- nodes_out$Phylum[match(edge_df$to, nodes_out$Id)]
  edge_color <- rep("#D0D0D0", nrow(edge_df))
  same_phy <- !is.na(edge_from_phy) & edge_from_phy == edge_to_phy &
    edge_from_phy %in% names(phylum_palette)
  edge_color[same_phy] <- phylum_palette[edge_from_phy[same_phy]]

  edges_out <- edge_df %>%
    transmute(
      Source = from,
      Target = to,
      Type = "Undirected",
      Weight = if ("r" %in% names(edge_df)) pmax(abs(suppressWarnings(as.numeric(r))), 0.001, na.rm = TRUE) else 1,
      r = if ("r" %in% names(edge_df)) suppressWarnings(as.numeric(r)) else NA_real_,
      p = if ("p" %in% names(edge_df)) suppressWarnings(as.numeric(p)) else NA_real_,
      sign = if ("sign" %in% names(edge_df)) sign else ifelse(!is.na(r) & r < 0, "negative", "positive"),
      Color = edge_color
    )

  write.csv(
    nodes_out,
    file.path(gephi_dir, paste0(out_stub, "_nodes.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    edges_out,
    file.path(gephi_dir, paste0(out_stub, "_edges.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

plot_all <- function(panel_info) {
  graphs <- lapply(seq_len(nrow(panel_info)), function(i) {
    make_network_graph(panel_info$PairID[i], panel_info$FileGroup[i])
  })

  titles <- vapply(seq_along(graphs), function(i) {
    paste0(
      panel_info$PanelTitle[i], "\n",
      "n = ", vcount(graphs[[i]])
    )
  }, character(1))

  for (i in seq_along(graphs)) {
    out_stub <- paste0(
      "representative_network_phylum_style_",
      panel_info$Domain[i], "_",
      ifelse(panel_info$FileGroup[i] == "CK", "CK", "N")
    )
    export_gephi_tables(graphs[[i]], panel_info[i, ], out_stub)
    pdf_path <- file.path(fig_dir, paste0(out_stub, ".pdf"))
    png_path <- file.path(fig_dir, paste0(out_stub, ".png"))

    grDevices::pdf(pdf_path, width = 5.9, height = 5.55, family = "sans")
    plot_one(graphs[[i]], titles[i], seed = 20260614 + i)
    dev.off()

    grDevices::png(png_path, width = 1770, height = 1665, res = 300, type = "cairo")
    plot_one(graphs[[i]], titles[i], seed = 20260614 + i)
    dev.off()

    message("Wrote: ", normalizePath(pdf_path, winslash = "/"))
    message("Wrote: ", normalizePath(png_path, winslash = "/"))
  }

  combined_pdf <- file.path(fig_dir, "representative_networks_phylum_style_Bacterial_Fungi.pdf")
  combined_png <- file.path(fig_dir, "representative_networks_phylum_style_Bacterial_Fungi.png")

  combined_phyla <- unique(unlist(lapply(graphs, function(g) V(g)$PhylumClean)))
  combined_phyla <- combined_phyla[combined_phyla %in% names(phylum_palette)]
  combined_counts <- sort(table(unlist(lapply(graphs, function(g) V(g)$PhylumClean))), decreasing = TRUE)
  combined_phyla <- names(combined_counts)[names(combined_counts) %in% combined_phyla]
  combined_phyla <- combined_phyla[seq_len(min(12, length(combined_phyla)))]

  draw_combined_legend <- function() {
    plot.new()
    par(mar = c(0, 0, 0, 0))
    legend(
      "topleft",
      legend = combined_phyla,
      pch = 21,
      pt.bg = phylum_palette[combined_phyla],
      pt.cex = 1.35,
      col = "white",
      bty = "n",
      cex = 1.1,
      title = "Phylum"
    )
    legend(
      "bottomleft",
      legend = c("Low degree", "High degree"),
      pch = 21,
      pt.bg = "grey70",
      col = "white",
      pt.cex = c(0.9, 1.8),
      bty = "n",
      cex = 1.05,
      title = "Degree"
    )
  }

  grDevices::pdf(combined_pdf, width = 10.8, height = 7.0, family = "sans")
  layout(matrix(c(1, 2, 5,
                  3, 4, 5), nrow = 2, byrow = TRUE), widths = c(1, 1, 0.30))
  par(oma = c(0.1, 0.1, 0.1, 0.1))
  for (i in seq_along(graphs)) {
    plot_one(
      graphs[[i]], titles[i],
      seed = 20260614 + i,
      show_legend = FALSE
    )
  }
  draw_combined_legend()
  dev.off()

  grDevices::png(combined_png, width = 3240, height = 2100, res = 300, type = "cairo")
  layout(matrix(c(1, 2, 5,
                  3, 4, 5), nrow = 2, byrow = TRUE), widths = c(1, 1, 0.30))
  par(oma = c(0.1, 0.1, 0.1, 0.1))
  for (i in seq_along(graphs)) {
    plot_one(
      graphs[[i]], titles[i],
      seed = 20260614 + i,
      show_legend = FALSE
    )
  }
  draw_combined_legend()
  dev.off()

  message("Wrote: ", normalizePath(combined_pdf, winslash = "/"))
  message("Wrote: ", normalizePath(combined_png, winslash = "/"))
}

panel_info <- build_panel_info()
plot_all(panel_info)

