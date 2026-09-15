#!/usr/bin/env Rscript

# ============================================================
# Celltaminate.R script
#
# This script can be run using command-line arguments.
#
# Required arguments:
#   --input_files        Comma-separated Kraken/KrakenUniq report files
#                        OR use --input_list for one file per line
#   --output_dir         Output directory
#
# Optional arguments:
#   --metadata_tsv               Metadata table with columns:
#                                sample, group, sample_type, is_control
#   --input_list                 Text file with one report path per line
#   --tax_level                  G or S [default: S]
#   --host_species               Comma-separated: human,mouse,drosophila [default: human]
#   --host_taxids_manual         Additional host taxids, comma-separated
#   --include_protozoa           TRUE/FALSE [default: TRUE]
#   --collapse_species           TRUE/FALSE [default: TRUE]
#   --collapse_min_genus_reads   Integer [default: 30]
#   --collapse_top_frac          Numeric [default: 0.85]
#   --fp_falsepos_cutoff         Numeric [default: 75]
#   --fp_true_cutoff             Numeric [default: 1]
#   --ref_bg_tsv                 Reference background TSV
#   --kitome_blacklist_tsv       Kitome/background blacklist TSV
#   --clinical_panel_tsv         Clinically important pathogens TSV
#   --quick_remove               Comma-separated: false,unc,true [default: false,unc]
#   --protected_taxa_tsv         Optional TSV with columns sample and taxon
#   --keep_descendants           TRUE/FALSE [default: FALSE]
#   --show_fp_breakdown          TRUE/FALSE [default: FALSE]
#   --save_cleaned_reanalysis    TRUE/FALSE [default: TRUE]
#
# Example:
# Rscript celltaminate.R \
#   --input_files sample1.txt,sample2.txt \
#   --metadata_tsv metadata.tsv \
#   --output_dir results \
#   --ref_bg_tsv refined_cell.lines.tsv \
#   --kitome_blacklist_tsv kitome_and_background_blacklist.tsv \
#   --clinical_panel_tsv clinically_important_pathogens.tsv
# ============================================================

# ------------------- LIBRARIES -------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stats)
  library(ggplot2)
  library(ggrepel)
  library(fmsb)
  library(writexl)
  library(tibble)
  library(tidyr)
})

# ------------------- OPTIONAL PACKAGES -------------------
HAS_ZIP <- requireNamespace("zip", quietly = TRUE)

# ------------------- UTILITIES -------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

safe_log <- function(x, eps = 1e-9) {
  suppressWarnings(log(pmax(as.numeric(x), 0) + eps))
}

safe_log10 <- function(x, eps = 1e-9) {
  suppressWarnings(log10(pmax(as.numeric(x), 0) + eps))
}

scale01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

coalesce0 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}

cap <- function(x, lo = 0, hi = 1) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  pmin(pmax(x, lo), hi)
}

trim_left <- function(x) trimws(x, which = "left")
safe_sample_id <- function(x) gsub("[^A-Za-z0-9]", "_", x)

read_optional_lines <- function(path) {
  if (is.null(path) || !file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

first_existing_path <- function(paths) {
  if (is.null(paths) || length(paths) == 0) return(NULL)
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0) return(NULL)
  for (p in unique(paths)) {
    if (file.exists(p)) return(p)
  }
  NULL
}

ui_default <- function(defaults, name, fallback) {
  if (is.null(defaults) || !is.list(defaults) || is.null(defaults[[name]])) return(fallback)
  defaults[[name]]
}

call_bucket_label <- function(call) {
  dplyr::case_when(
    call == "Likely true" ~ "Likely true",
    call == "Likely false positive / background" ~ "Likely false positive/background",
    TRUE ~ "Uncertain"
  )
}

CALL_BUCKET_LEVELS <- c(
  "Likely true",
  "Uncertain",
  "Likely false positive/background"
)

CALL_PALETTE <- c(
  "Likely false positive/background" = "#8B0000",
  "Uncertain" = "#FF9800",
  "Likely true" = "#2E7D32"
)

make_prevalence_abundance_df <- function(gs_long, tax_level = "S", total_samples = NULL) {
  df <- gs_long %>%
    filter(rank == tax_level, call != "Non-microbial / Host")

  if (nrow(df) == 0) {
    return(tibble())
  }

  if ("detected" %in% colnames(df)) {
    df$display_detected <- df$detected %in% TRUE
  } else {
    df$display_detected <- is.finite(df$reads_clade) & df$reads_clade > 0
  }

  n_samples <- suppressWarnings(as.integer(total_samples %||% dplyr::n_distinct(gs_long$sample)))
  if (!is.finite(n_samples) || n_samples <= 0) {
    return(tibble())
  }

  prev_df <- df %>%
    group_by(name_clean) %>%
    summarise(
      prevalence = n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]) / n_samples,
      mean_rpmm_detected = if (any(display_detected & is.finite(rpmm), na.rm = TRUE)) {
        mean(rpmm[display_detected & is.finite(rpmm)], na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_log_rpmm_detected = if (any(display_detected & is.finite(log_rpmm), na.rm = TRUE)) {
        mean(log_rpmm[display_detected & is.finite(log_rpmm)], na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  bucket_df <- df %>%
    filter(is.finite(reads_clade), reads_clade > 0, is.finite(rpmm)) %>%
    mutate(call_bucket = factor(call_bucket_label(call), levels = CALL_BUCKET_LEVELS)) %>%
    count(name_clean, call_bucket, name = "n_bucket") %>%
    arrange(name_clean, desc(n_bucket), call_bucket) %>%
    group_by(name_clean) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(call_bucket = as.character(call_bucket)) %>%
    select(name_clean, call_bucket)

  prev_df %>%
    left_join(bucket_df, by = "name_clean") %>%
    mutate(
      call_bucket = dplyr::coalesce(call_bucket, "Uncertain"),
      prevalence = cap(prevalence, 0, 1)
    )
}

get_bioai_taxa_info <- function(gs, tax_level = "S") {
  if (is.null(gs) || nrow(gs) == 0) {
    return(list(
      candidates = character(0),
      true_taxa = character(0),
      initial_selected = character(0),
      initial_choices = character(0)
    ))
  }

  base <- gs %>%
    filter(rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm), nzchar(name_clean)) %>%
    mutate(
      call_priority = dplyr::case_when(
        call == "Likely true" ~ 1,
        call == "Uncertain" ~ 2,
        call == "Likely false positive / background" ~ 3,
        TRUE ~ 4
      )
    )

  candidates <- base %>%
    arrange(call_priority, desc(rpmm)) %>%
    pull(name_clean) %>%
    unique()

  true_taxa <- base %>%
    filter(call == "Likely true") %>%
    arrange(desc(rpmm)) %>%
    pull(name_clean) %>%
    unique()

  initial_selected <- if (length(true_taxa) > 0) true_taxa else head(candidates, 3)
  initial_choices <- unique(c(initial_selected, head(candidates, 100)))

  list(
    candidates = as.character(candidates),
    true_taxa = as.character(true_taxa),
    initial_selected = as.character(initial_selected),
    initial_choices = as.character(initial_choices)
  )
}

# ------------------- KRAKEN REPORT PARSING -------------------
detect_report_format <- function(df) {
  # Returns "krakenuniq_8col", "kraken2_6col", or NA
  if (ncol(df) >= 8) return("krakenuniq_8col")
  if (ncol(df) == 6) return("kraken2_6col")
  NA_character_
}

read_kraken_report <- function(file_path, sample_name = NULL) {
  raw <- tryCatch(
    read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", quote = ""),
    error = function(e) NULL
  )

  if (is.null(raw) || nrow(raw) == 0) {
    return(list(error = paste0("Failed to read file: ", basename(file_path))))
  }

  fmt <- detect_report_format(raw)

  if (is.na(fmt)) {
    return(list(error = paste0(
      "Invalid Kraken report format for: ", basename(file_path),
      ". Expected 6-column Kraken2 report or 8-column KrakenUniq report."
    )))
  }

  if (fmt == "krakenuniq_8col") {
    raw8 <- raw[, 1:8]
    colnames(raw8) <- c("pct", "reads_clade", "reads_direct", "uniq_kmers", "dup_kmers", "rank", "taxid", "name_raw")
    df <- raw8 %>%
      mutate(
        pct = as.numeric(pct),
        reads_clade = as.numeric(reads_clade),
        reads_direct = as.numeric(reads_direct),
        uniq_kmers = as.numeric(uniq_kmers),
        dup_kmers = as.numeric(dup_kmers),
        rank = as.character(rank),
        taxid = suppressWarnings(as.integer(taxid)),
        name_raw = as.character(name_raw)
      )
  } else {
    raw6 <- raw
    colnames(raw6) <- c("pct", "reads_clade", "reads_direct", "rank", "taxid", "name_raw")
    df <- raw6 %>%
      mutate(
        pct = as.numeric(pct),
        reads_clade = as.numeric(reads_clade),
        reads_direct = as.numeric(reads_direct),
        uniq_kmers = NA_real_,
        dup_kmers = NA_real_,
        rank = as.character(rank),
        taxid = suppressWarnings(as.integer(taxid)),
        name_raw = as.character(name_raw)
      )
  }

  df <- df %>%
    mutate(
      sample = sample_name %||% basename(file_path),
      name_clean = trim_left(name_raw)
    )

  lead_spaces <- nchar(df$name_raw) - nchar(df$name_clean)
  lead_spaces[!is.finite(lead_spaces)] <- 0
  df$depth <- floor(pmax(lead_spaces, 0) / 2)
  df$kmer_per_read <- ifelse(is.na(df$uniq_kmers) | df$reads_clade <= 0, NA_real_, df$uniq_kmers / df$reads_clade)

  list(error = NULL, raw = raw, df = df, report_format = fmt)
}

annotate_lineage <- function(df) {
  n <- nrow(df)
  lineage <- character(n)
  superkingdom <- character(n)
  genus_ancestor <- character(n)
  species_ancestor <- character(n)

  stack_names <- character(0)
  stack_ranks <- character(0)

  for (i in seq_len(n)) {
    d <- df$depth[i]
    if (!is.finite(d) || d < 0) d <- 0

    if (length(stack_names) < d) {
      stack_names <- c(stack_names, rep(NA_character_, d - length(stack_names)))
      stack_ranks <- c(stack_ranks, rep(NA_character_, d - length(stack_ranks)))
    }

    if (length(stack_names) > d) {
      stack_names <- stack_names[1:d]
      stack_ranks <- stack_ranks[1:d]
    }

    stack_names <- c(stack_names, df$name_clean[i])
    stack_ranks <- c(stack_ranks, df$rank[i])

    lineage[i] <- paste(stack_names, collapse = ";")

    if ("Bacteria" %in% stack_names) superkingdom[i] <- "Bacteria" else
      if ("Archaea" %in% stack_names) superkingdom[i] <- "Archaea" else
        if ("Viruses" %in% stack_names) superkingdom[i] <- "Viruses" else
          if ("Fungi" %in% stack_names) superkingdom[i] <- "Fungi" else
            if ("Eukaryota" %in% stack_names) superkingdom[i] <- "Eukaryota" else
              superkingdom[i] <- NA_character_

    g_idx <- which(stack_ranks == "G")
    s_idx <- which(stack_ranks == "S")

    genus_ancestor[i] <- if (length(g_idx) > 0) stack_names[max(g_idx)] else NA_character_
    species_ancestor[i] <- if (length(s_idx) > 0) stack_names[max(s_idx)] else NA_character_
  }

  df$lineage_path <- lineage
  df$superkingdom <- superkingdom
  df$genus_ancestor <- genus_ancestor
  df$species_ancestor <- species_ancestor
  df
}

compute_parent_idx_from_depth <- function(depth) {
  depth <- suppressWarnings(as.integer(depth))
  parent_idx <- rep(NA_integer_, length(depth))
  stack <- integer(0)

  for (i in seq_along(depth)) {
    d <- depth[i]
    if (!is.finite(d) || d < 0) d <- 0L

    if (length(stack) > d) {
      stack <- stack[seq_len(d)]
    }

    parent_idx[i] <- if (d <= 0 || length(stack) < d) NA_integer_ else stack[d]

    if (length(stack) < d + 1) {
      stack <- c(stack, rep(NA_integer_, d + 1 - length(stack)))
    }

    stack[d + 1] <- i
  }

  parent_idx
}

annotate_top_strain_rows <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  df$Top_strain <- NA_character_
  df$Top_strain_rank <- NA_character_
  df$Top_strain_reads <- NA_integer_

  if (!all(c("name_clean", "rank", "reads_clade", "depth") %in% colnames(df))) {
    return(df)
  }

  rank <- trimws(as.character(df$rank))
  name_clean <- trimws(as.character(df$name_clean))
  reads_clade <- suppressWarnings(as.numeric(df$reads_clade))
  depth <- suppressWarnings(as.integer(df$depth))
  parent_idx <- compute_parent_idx_from_depth(depth)

  species_idx <- which(
    rank == "S" &
      !is.na(name_clean) &
      nzchar(name_clean) &
      is.finite(reads_clade) &
      reads_clade > 0
  )

  strain_idx <- which(
    grepl("^S[0-9]+$", rank) &
      !is.na(name_clean) &
      nzchar(name_clean) &
      is.finite(reads_clade) &
      reads_clade > 0
  )

  if (length(species_idx) == 0 || length(strain_idx) == 0) {
    return(df)
  }

  best_idx_by_species <- rep(NA_integer_, nrow(df))

  for (i in strain_idx) {
    p <- parent_idx[i]
    species_parent <- NA_integer_

    while (is.finite(p)) {
      if (!is.na(rank[p]) && rank[p] == "S") {
        species_parent <- p
        break
      }
      p <- parent_idx[p]
    }

    if (!is.finite(species_parent)) {
      next
    }

    current_best <- best_idx_by_species[species_parent]

    if (
      is.na(current_best) ||
      reads_clade[i] > reads_clade[current_best] ||
      (reads_clade[i] == reads_clade[current_best] && name_clean[i] < name_clean[current_best])
    ) {
      best_idx_by_species[species_parent] <- i
    }
  }

  assigned_species <- which(is.finite(best_idx_by_species))

  if (length(assigned_species) > 0) {
    df$Top_strain[assigned_species] <- name_clean[best_idx_by_species[assigned_species]]
    df$Top_strain_rank[assigned_species] <- rank[best_idx_by_species[assigned_species]]
    df$Top_strain_reads[assigned_species] <- as.integer(round(reads_clade[best_idx_by_species[assigned_species]]))
  }

  df
}

infer_host_like <- function(df, host_taxids = c(9606L, 9605L, 10090L, 7227L), host_name_patterns = NULL) {
  if (is.null(host_name_patterns)) {
    host_name_patterns <- c(
      "\\bHomo\\b",
      "\\bHomo sapiens\\b",
      "\\bPrimates\\b",
      "\\bHominidae\\b",
      "\\bChordata\\b",
      "\\bMammalia\\b",
      "\\bMetazoa\\b",
      "\\bMus\\b",
      "\\bMus musculus\\b",
      "\\bRodentia\\b",
      "\\bDrosophila\\b",
      "\\bDrosophilidae\\b",
      "\\bInsecta\\b"
    )
  }

  host_taxid_flag <- (!is.na(df$taxid) & df$taxid %in% as.integer(host_taxids))

  if (length(host_name_patterns) > 0) {
    host_name_flag <- Reduce(`|`, lapply(host_name_patterns, function(p) grepl(p, df$name_clean, ignore.case = TRUE)))
  } else {
    host_name_flag <- rep(FALSE, nrow(df))
  }

  lineage <- df$lineage_path %||% rep("", nrow(df))
  lineage[is.na(lineage)] <- ""

  if (length(host_name_patterns) > 0) {
    host_lineage_flag <- Reduce(`|`, lapply(host_name_patterns, function(p) grepl(p, lineage, ignore.case = TRUE)))
  } else {
    host_lineage_flag <- rep(FALSE, nrow(df))
  }

  df$is_host <- host_taxid_flag | host_name_flag | host_lineage_flag
  df
}

infer_plant_like <- function(df) {
  lineage <- df$lineage_path %||% rep("", nrow(df))
  lineage[is.na(lineage)] <- ""
  df$is_plant <- grepl("\\bViridiplantae\\b|\\bEmbryophyta\\b|\\bStreptophyta\\b", lineage, ignore.case = TRUE)
  df
}

infer_microbial <- function(df, include_protozoa = TRUE) {
  sk <- df$superkingdom
  base_microbe <- sk %in% c("Bacteria", "Archaea", "Viruses", "Fungi")
  proto_microbe <- (include_protozoa & sk == "Eukaryota" & !df$is_host & !df$is_plant)
  df$is_microbial <- base_microbe | proto_microbe
  df
}

compute_reads_summary <- function(df) {
  root_reads <- df$reads_clade[df$rank %in% c("R") | df$taxid == 1L | df$name_clean == "root"]
  unclassified_reads <- df$reads_clade[df$rank %in% c("U") | df$taxid == 0L | df$name_clean == "unclassified"]

  root_reads <- root_reads[!is.na(root_reads)]
  unclassified_reads <- unclassified_reads[!is.na(unclassified_reads)]

  classified <- ifelse(length(root_reads) > 0, root_reads[1], NA_real_)
  unclassified <- ifelse(length(unclassified_reads) > 0, unclassified_reads[1], 0)
  total <- NA_real_

  if (!is.na(classified)) {
    total <- classified + unclassified
  } else {
    total <- sum(df$reads_direct, na.rm = TRUE)
  }

  get_reads_name <- function(nm) {
    vals <- df$reads_clade[df$name_clean == nm]
    if (length(vals) == 0) return(0)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(0)
    vals[1]
  }

  bacteria <- get_reads_name("Bacteria")
  archaea <- get_reads_name("Archaea")
  viruses <- get_reads_name("Viruses")
  eukaryota <- get_reads_name("Eukaryota")
  metazoa <- get_reads_name("Metazoa")
  viridiplantae <- get_reads_name("Viridiplantae")

  euk_microbe <- max(eukaryota - metazoa - viridiplantae, 0)
  microbial_reads <- bacteria + archaea + viruses + euk_microbe
  host_reads <- metazoa

  list(
    total_reads = as.numeric(total),
    classified_reads = as.numeric(classified %||% NA_real_),
    unclassified_reads = as.numeric(unclassified),
    microbial_reads = as.numeric(microbial_reads),
    host_reads = as.numeric(host_reads)
  )
}

calc_shannon <- function(p) {
  p <- as.numeric(p)
  p <- p[p > 0 & is.finite(p)]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

# ------------------- REFERENCE-FREE COHORT ANALYSIS -------------------
default_kitome <- function() {
  c(
    "Ralstonia",
    "Cutibacterium acnes",
    "Propionibacterium acnes",
    "Bradyrhizobium",
    "Sphingomonas",
    "Methylobacterium",
    "Acinetobacter",
    "Pseudomonas",
    "Bacillus",
    "Staphylococcus",
    "Corynebacterium",
    "Escherichia coli",
    "Enterococcus",
    "Micrococcus",
    "Delftia",
    "Burkholderia",
    "Comamonas",
    "Achromobacter",
    "Cupriavidus",
    "Serratia",
    "Klebsiella"
  )
}

KITOME_BLACKLIST_PATH <- "kitome_and_background_blacklist.tsv"

empty_kitome_blacklist_df <- function() {
  tibble(
    organism_name = character(0),
    evidence_class = character(0),
    doi_hyperlink = character(0),
    note = character(0),
    organism_name_lc = character(0)
  )
}

load_kitome_blacklist <- function(path = KITOME_BLACKLIST_PATH) {
  if (is.null(path) || !file.exists(path)) {
    return(empty_kitome_blacklist_df())
  }

  df <- tryCatch(
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "	", quote = ""),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0) {
    return(empty_kitome_blacklist_df())
  }

  nms <- colnames(df)
  if (!"organism_name" %in% nms) {
    if ("name" %in% nms) df$organism_name <- df$name
    if ("organism" %in% nms) df$organism_name <- df$organism
  }
  if (!"evidence_class" %in% nms) {
    df$evidence_class <- NA_character_
  }
  if (!"doi_hyperlink" %in% nms) {
    if ("doi" %in% nms) df$doi_hyperlink <- df$doi else df$doi_hyperlink <- NA_character_
  }
  if (!"note" %in% nms) {
    df$note <- NA_character_
  }

  if (!"organism_name" %in% colnames(df)) {
    return(empty_kitome_blacklist_df())
  }

  df %>%
    transmute(
      organism_name = trimws(as.character(organism_name)),
      evidence_class = trimws(as.character(evidence_class)),
      doi_hyperlink = trimws(as.character(doi_hyperlink)),
      note = trimws(as.character(note)),
      organism_name_lc = tolower(trimws(as.character(organism_name)))
    ) %>%
    filter(nzchar(organism_name)) %>%
    distinct(organism_name_lc, .keep_all = TRUE)
}

KITOME_BG_BLACKLIST <- empty_kitome_blacklist_df()

get_kitome_blacklist_names <- function(blacklist_df = KITOME_BG_BLACKLIST) {
  base_names <- default_kitome()
  file_names <- if (!is.null(blacklist_df) && nrow(blacklist_df) > 0) blacklist_df$organism_name else character(0)
  unique(c(base_names, file_names))
}

CLINICALLY_IMPORTANT_PATHOGENS_PATH <- "clinically_important_pathogens.tsv"

empty_clinically_important_pathogens_df <- function() {
  tibble(
    category = character(0),
    Name = character(0),
    clinical_category = character(0),
    Remarks = character(0),
    doi = character(0),
    name_lc = character(0)
  )
}

load_clinically_important_pathogens <- function(path = CLINICALLY_IMPORTANT_PATHOGENS_PATH) {
  if (is.null(path) || !file.exists(path)) {
    return(empty_clinically_important_pathogens_df())
  }

  df <- tryCatch(
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", quote = ""),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0) {
    return(empty_clinically_important_pathogens_df())
  }

  nms <- colnames(df)

  if (!"Name" %in% nms) {
    if ("name" %in% nms) df$Name <- df$name
    if ("organism_name" %in% nms) df$Name <- df$organism_name
    if ("organism" %in% nms) df$Name <- df$organism
  }
  if (!"clinical_category" %in% nms) {
    df$clinical_category <- NA_character_
  }
  if (!"Remarks" %in% nms) {
    if ("remark" %in% nms) df$Remarks <- df$remark else
    if ("remarks" %in% nms) df$Remarks <- df$remarks else
    if ("note" %in% nms) df$Remarks <- df$note else
      df$Remarks <- NA_character_
  }
  if (!"doi" %in% nms) {
    if ("doi_hyperlink" %in% nms) df$doi <- df$doi_hyperlink else df$doi <- NA_character_
  }
  if (!"category" %in% nms) {
    df$category <- NA_character_
  }

  if (!"Name" %in% colnames(df)) {
    return(empty_clinically_important_pathogens_df())
  }

  df %>%
    transmute(
      category = trimws(as.character(category)),
      Name = trimws(as.character(Name)),
      clinical_category = trimws(as.character(clinical_category)),
      Remarks = trimws(as.character(Remarks)),
      doi = trimws(as.character(doi)),
      name_lc = tolower(trimws(as.character(Name)))
    ) %>%
    filter(nzchar(Name)) %>%
    distinct(name_lc, .keep_all = TRUE)
}

CLINICALLY_IMPORTANT_PATHOGENS <- empty_clinically_important_pathogens_df()

get_clinically_important_pathogen_names <- function(pathogen_df = CLINICALLY_IMPORTANT_PATHOGENS) {
  if (is.null(pathogen_df) || nrow(pathogen_df) == 0) return(character(0))
  unique(pathogen_df$Name)
}

collapse_species_within_genus <- function(gs_long, min_genus_reads = 20, top_frac_keep = 0.80) {
  sp <- gs_long %>% filter(rank == "S")
  if (nrow(sp) == 0) return(gs_long)

  sp_sum <- sp %>%
    group_by(sample, genus) %>%
    summarize(
      total_sp_reads = sum(reads_clade, na.rm = TRUE),
      top_reads = max(reads_clade, na.rm = TRUE),
      top_name = name_clean[which.max(reads_clade)[1]],
      top_frac = ifelse(total_sp_reads > 0, top_reads / total_sp_reads, NA_real_),
      .groups = "drop"
    )

  sp <- sp %>%
    left_join(sp_sum, by = c("sample", "genus"))

  keep_sp <- (sp$total_sp_reads >= min_genus_reads) &
    (sp$name_clean == sp$top_name) &
    (sp$top_frac >= top_frac_keep)

  drop_sp <- (sp$total_sp_reads >= min_genus_reads) & (!keep_sp)
  sp$keep_flag <- ifelse(is.na(drop_sp), TRUE, !drop_sp)

  sp_filtered <- sp %>%
    filter(keep_flag) %>%
    select(-total_sp_reads, -top_reads, -top_name, -top_frac, -keep_flag)

  other <- gs_long %>% filter(rank != "S")
  bind_rows(other, sp_filtered)
}

# ------------------- BACKGROUND NEGATIVE CONTROL REFERENCE DATASET  -------------------
REF_CELL_LINES_PATH <- "refined_cell.lines.tsv"

load_ref_cell_lines <- function(path = REF_CELL_LINES_PATH, eps = 1e-9) {
  if (!file.exists(path)) {
    return(list(available = FALSE, path = path, ref = NULL, stats = NULL, comp_genus = NULL))
  }

  ref <- tryCatch(read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)

  if (is.null(ref) || nrow(ref) == 0) {
    return(list(available = FALSE, path = path, ref = NULL, stats = NULL, comp_genus = NULL))
  }

  if (!"name" %in% colnames(ref)) {
    if ("taxon" %in% colnames(ref)) ref$name <- ref$taxon
    if ("Name" %in% colnames(ref)) ref$name <- ref$Name
  }

  if (!"rank" %in% colnames(ref)) {
    if ("Rank" %in% colnames(ref)) ref$rank <- ref$Rank
  }

  if (!"rpmm" %in% colnames(ref)) {
    if ("RPMM" %in% colnames(ref)) ref$rpmm <- ref$RPMM
    if ("value" %in% colnames(ref)) ref$rpmm <- ref$value
  }

  if (!all(c("rank", "name", "rpmm") %in% colnames(ref))) {
    return(list(available = FALSE, path = path, ref = NULL, stats = NULL, comp_genus = NULL))
  }

  ref <- ref %>%
    mutate(
      rank = as.character(rank),
      name_clean = trimws(as.character(name)),
      rpmm = suppressWarnings(as.numeric(rpmm)),
      log_rpmm = safe_log(rpmm, eps = eps)
    ) %>%
    filter(rank %in% c("G", "S"), nzchar(name_clean))

  if (!"sample" %in% colnames(ref)) {
    ref$sample <- NA_character_
  }
  ref$sample <- as.character(ref$sample)

  n_ref_samples <- length(unique(ref$sample[!is.na(ref$sample) & nzchar(ref$sample)]))
  if (!is.finite(n_ref_samples) || n_ref_samples <= 0) n_ref_samples <- NA_real_

  if (nrow(ref) == 0) {
    return(list(available = FALSE, path = path, ref = NULL, stats = NULL, comp_genus = NULL))
  }

  stats <- ref %>%
    group_by(rank, name_clean) %>%
    summarise(
      ref_n = sum(is.finite(rpmm)),
      ref_prev = if_else(
        is.finite(n_ref_samples),
        n_distinct(sample[rpmm > 0]) / n_ref_samples,
        mean(rpmm > 0, na.rm = TRUE)
      ),
      ref_q50_log = if_else(ref_n >= 1, as.numeric(stats::quantile(log_rpmm, 0.50, na.rm = TRUE)), NA_real_),
      ref_q95_log = if_else(ref_n >= 1, as.numeric(stats::quantile(log_rpmm, 0.95, na.rm = TRUE)), NA_real_),
      ref_q99_log = if_else(ref_n >= 1, as.numeric(stats::quantile(log_rpmm, 0.99, na.rm = TRUE)), NA_real_),
      ref_med_rpmm = if_else(ref_n >= 1, as.numeric(stats::median(rpmm, na.rm = TRUE)), NA_real_),
      ref_q95_rpmm = if_else(ref_n >= 1, as.numeric(stats::quantile(rpmm, 0.95, na.rm = TRUE)), NA_real_),
      .groups = "drop"
    )

  comp_genus <- ref %>%
    filter(rank == "G", is.finite(rpmm), rpmm > 0) %>%
    group_by(name_clean) %>%
    summarise(ref_med = median(rpmm, na.rm = TRUE), .groups = "drop") %>%
    mutate(ref_med = coalesce0(ref_med)) %>%
    arrange(desc(ref_med))

  tot <- sum(comp_genus$ref_med, na.rm = TRUE)
  comp_genus <- comp_genus %>%
    mutate(ref_p = ifelse(is.finite(tot) & tot > 0, ref_med / tot, NA_real_))

  list(available = TRUE, path = path, ref = ref, stats = stats, comp_genus = comp_genus)
}

REF_BG <- list(available = FALSE, path = REF_CELL_LINES_PATH, ref = NULL, stats = NULL, comp_genus = NULL)

js_divergence <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  p[!is.finite(p)] <- 0
  q[!is.finite(q)] <- 0
  p <- p / sum(p + eps)
  q <- q / sum(q + eps)
  m <- 0.5 * (p + q)

  kl <- function(a, b) {
    a <- a + eps
    b <- b + eps
    sum(a * log2(a / b))
  }

  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

final_score_model <- function() {
  feature_names <- c(
    "reads_support",
    "rpmm_support",
    "uniq_kmer_support",
    "kmer_per_read_support",
    "ambiguity_unresolved_fraction",
    "ambiguity_species_nondominance",
    "reference_prevalence",
    "reference_median_abundance",
    "reference_relative_abundance",
    "reference_at_or_below_q99",
    "reference_at_or_below_q95",
    "reference_at_or_below_q50",
    "decontaminated_abundance",
    "decontamination_ratio",
    "kitome_only",
    "kitome_clinical_overlap",
    "clinical_membership",
    "low_biomass_context",
    "reference_like_sample_composition",
    "clinical_x_decon_support",
    "kitome_only_x_decon_support",
    "kitome_overlap_x_decon_support",
    "reference_prevalence_x_low_enrichment",
    "unresolved_ambiguity_x_kmer_support",
    "species_nondominance_x_kmer_support"
  )
  
  list(
    feature_names = feature_names,
    
    intercept = 1.4925028628700236,
    
    score_offset = -4.9116688158923365,
    
    impute_values = setNames(
      c(
        -4.037345601140845,
        -5.560973858316463,
        -6.791598579642281,
        -2.868097637369349,
        0.13617358617707387,
        0.3418085172884586,
        0.5147502612680851,
        3.9850698122986192,
        -1.5728884233552725,
        0.8623956126137355,
        0.8000747849931447,
        0.16708213885080395,
        -4.451715457237202,
        -0.5715219157760044,
        0.1524074074074074,
        0.057901234567901236,
        -0.1365432098765432,
        0.0,
        -0.8760988289857169,
        -0.6162082363065745,
        -0.6411394789276108,
        -0.21839854637573825,
        0.09053060438834876,
        -0.8806494090358415,
        -2.1563274408371043
      ),
      feature_names
    ),
    
    means = setNames(
      c(
        -4.037345601140845,
        -5.560973858316463,
        -6.791598579642281,
        -2.868097637369349,
        0.13617358617707384,
        0.3418085172884586,
        0.5147502612680849,
        3.9850698122986192,
        -1.5728884233552725,
        0.8623956126137355,
        0.8000747849931447,
        0.16708213885080392,
        -4.451715457237202,
        -0.5715219157760044,
        0.1524074074074074,
        0.057901234567901236,
        -0.1365432098765432,
        0.0,
        -0.8760988289857169,
        -0.6162082363065745,
        -0.6411394789276108,
        -0.21839854637573825,
        0.09053060438834876,
        -0.8806494090358415,
        -2.1563274408371043
      ),
      feature_names
    ),
    
    scales = setNames(
      c(
        1.5128497297793833,
        2.150800601166756,
        1.5801497288378767,
        0.7665211438955583,
        0.16611865990321514,
        0.283347787544235,
        0.2586009543782396,
        1.4482795250848732,
        2.159788578677823,
        0.3428430077184419,
        0.3980383930748458,
        0.37127182030409095,
        2.974864934519208,
        0.3576872461912685,
        0.3594153441003318,
        0.2335565918646145,
        0.3433644735745873,
        1.0,
        0.037336150006150516,
        2.1033078430172503,
        1.9433149640003147,
        1.3129930478842984,
        0.37318236981981967,
        1.1019121244936676,
        1.7770162314781006
      ),
      feature_names
    ),
    
    standardized_weights = setNames(
      c(
        0.25588616843532014,
        0.015857076698091055,
        0.23738110321065942,
        0.0,
        0.0,
        0.13855067199962626,
        0.011924245703707377,
        0.1542911530541554,
        0.1201126650864708,
        0.14623256321543165,
        0.11672991539965896,
        0.04072473708446774,
        0.0684245954292958,
        0.11686135228472111,
        0.0,
        0.027464346253092846,
        0.207776863730765,
        0.0,
        0.0,
        0.2429811424456025,
        0.1631094360713755,
        0.0,
        0.07293897105982951,
        0.13037149291281877,
        0.0
      ),
      feature_names
    )
  )
}

compute_sample_priors <- function(gs_long, meta_df, params, ref_comp_genus = NULL) {
  df0 <- gs_long %>%
    distinct(sample, total_reads, microbial_reads) %>%
    mutate(
      log_microbial_reads = safe_log10(microbial_reads + 1),
      biomass_scaled = scale01(log_microbial_reads)
    )
  
  js_df <- df0 %>%
    mutate(
      js_div = NA_real_
    )
  
  if (!is.null(ref_comp_genus) && nrow(ref_comp_genus) > 0) {
    ref_vec <- ref_comp_genus %>%
      filter(is.finite(ref_p), ref_p > 0) %>%
      select(
        genus = name_clean,
        ref_p
      )
    
    samp_genus <- gs_long %>%
      filter(
        rank == "G",
        is_microbial,
        !is_host,
        !is_plant,
        is.finite(rpmm),
        rpmm > 0
      ) %>%
      group_by(
        sample,
        genus = name_clean
      ) %>%
      summarise(
        rpmm = sum(rpmm, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (nrow(samp_genus) > 0) {
      js_calc <- samp_genus %>%
        group_by(sample) %>%
        group_modify(
          ~{
            d <- .x
            
            p <- d$rpmm
            p[!is.finite(p)] <- 0
            names(p) <- d$genus
            
            if (sum(p) <= 0) {
              return(
                tibble(
                  js_div = NA_real_
                )
              )
            }
            
            p <- p / sum(p)
            
            all_g <- union(
              names(p),
              ref_vec$genus
            )
            
            p2 <- rep(
              0,
              length(all_g)
            )
            names(p2) <- all_g
            
            q2 <- rep(
              0,
              length(all_g)
            )
            names(q2) <- all_g
            
            p2[names(p)] <- p
            q2[ref_vec$genus] <- ref_vec$ref_p
            
            qsum <- sum(q2)
            
            if (!is.finite(qsum) || qsum <= 0) {
              return(
                tibble(
                  js_div = NA_real_
                )
              )
            }
            
            q2 <- q2 / qsum
            
            tibble(
              js_div = js_divergence(
                p2,
                q2
              )
            )
          }
        ) %>%
        ungroup()
      
      js_df <- df0 %>%
        left_join(
          js_calc,
          by = "sample"
        )
    }
  }
  
  js_df %>%
    mutate(
      js_div = if_else(
        is.finite(js_div),
        js_div,
        NA_real_
      )
    ) %>%
    select(
      sample,
      biomass_scaled,
      js_div
    )
}

compute_ambiguity_index <- function(gs_long) {
  base <- gs_long %>%
    filter(
      rank %in% c("G", "S"),
      is_microbial,
      !is_host,
      !is_plant
    )
  
  if (nrow(base) == 0) {
    return(
      tibble(
        sample = character(0),
        rank = character(0),
        name_clean = character(0),
        ambiguity_unresolved_fraction = numeric(0),
        ambiguity_species_nondominance = numeric(0)
      )
    )
  }
  
  genus_rows <- base %>%
    filter(rank == "G") %>%
    transmute(
      sample,
      genus = name_clean,
      genus_clade = suppressWarnings(
        as.numeric(reads_clade)
      ),
      genus_direct = suppressWarnings(
        as.numeric(reads_direct)
      )
    )
  
  sp_sum <- base %>%
    filter(rank == "S") %>%
    mutate(
      species_direct_nonnegative = pmax(
        coalesce0(reads_direct),
        0
      )
    ) %>%
    group_by(
      sample,
      genus
    ) %>%
    summarise(
      sp_total_direct = sum(
        species_direct_nonnegative,
        na.rm = TRUE
      ),
      top_reads = max(
        species_direct_nonnegative,
        na.rm = TRUE
      ),
      top_frac = if_else(
        sp_total_direct > 0,
        top_reads / sp_total_direct,
        NA_real_
      ),
      .groups = "drop"
    )
  
  sp_out <- base %>%
    filter(rank == "S") %>%
    left_join(
      sp_sum,
      by = c(
        "sample",
        "genus"
      )
    ) %>%
    left_join(
      genus_rows,
      by = c(
        "sample",
        "genus"
      )
    ) %>%
    mutate(
      ambiguity_unresolved_fraction = if_else(
        is.finite(genus_clade) &
          genus_clade > 0 &
          is.finite(genus_direct),
        cap(
          genus_direct / genus_clade,
          0,
          1
        ),
        NA_real_
      ),
      
      ambiguity_species_nondominance = if_else(
        is.finite(top_frac),
        cap(
          1 - top_frac,
          0,
          1
        ),
        NA_real_
      )
    ) %>%
    select(
      sample,
      rank,
      name_clean,
      ambiguity_unresolved_fraction,
      ambiguity_species_nondominance
    )
  
  g_out <- genus_rows %>%
    mutate(
      rank = "G",
      name_clean = genus,
      
      ambiguity_unresolved_fraction = if_else(
        is.finite(genus_clade) &
          genus_clade > 0 &
          is.finite(genus_direct),
        cap(
          genus_direct / genus_clade,
          0,
          1
        ),
        NA_real_
      ),
      
      ambiguity_species_nondominance = NA_real_
    ) %>%
    select(
      sample,
      rank,
      name_clean,
      ambiguity_unresolved_fraction,
      ambiguity_species_nondominance
    )
  
  bind_rows(
    g_out,
    sp_out
  )
}

# ------------------- Celltaminate Algorithm -----------------------------------------------------------
compute_taxon_features <- function(gs_long, meta_df, params, user_contam = character(0), clinical_panel = character(0), kitome = default_kitome()) {
  user_contam <- unique(trimws(user_contam))
  user_contam <- user_contam[nzchar(user_contam)]

  clinical_panel <- unique(trimws(clinical_panel))
  clinical_panel <- clinical_panel[nzchar(clinical_panel)]

  gs_long <- gs_long %>%
    mutate(
      log_microbial_reads = safe_log10(microbial_reads + 1),
      log_rpmm = safe_log(rpmm, eps = params$eps),
      nonmicrobe_frac = if_else(
        is.finite(total_reads) & total_reads > 0,
        pmax(total_reads - microbial_reads, 0) / total_reads,
        NA_real_
      )
    )

  has_kmers <- any(!is.na(gs_long$uniq_kmers))
  det_kmer_ok <- if (has_kmers) {
    (
      (is.na(gs_long$uniq_kmers) & is.na(gs_long$kmer_per_read)) |
      (!is.na(gs_long$kmer_per_read) & gs_long$kmer_per_read >= params$min_kmer_per_read_prevalence) |
      (!is.na(gs_long$uniq_kmers) & gs_long$uniq_kmers >= params$min_uniq_kmers_prevalence)
    )
  } else {
    TRUE
  }

  gs_long <- gs_long %>%
    mutate(
      detected = is_microbial &
        !is_host &
        !is_plant &
        reads_clade >= params$min_reads_prevalence &
        rpmm >= params$min_rpmm_prevalence &
        det_kmer_ok
    )

  n_samples <- meta_df %>% distinct(sample) %>% nrow()

  control_samples <- character(0)
  if ("is_control" %in% colnames(meta_df)) {
    control_flag <- tolower(trimws(as.character(meta_df$is_control))) %in% c("true", "t", "1", "yes", "y")
    control_samples <- unique(as.character(meta_df$sample[!is.na(control_flag) & control_flag]))
    control_samples <- control_samples[nzchar(control_samples)]
  }
  n_controls <- length(unique(control_samples))
  
  use_cohort_prevalence <- is.finite(n_samples) && n_samples >= 3
  use_control_prevalence <- is.finite(n_controls) && n_controls >= 3
  
  base_filt <- gs_long$is_microbial & !gs_long$is_host & !gs_long$is_plant & gs_long$rank %in% c("G", "S")

  tax_base <- gs_long %>%
    filter(base_filt) %>%
    group_by(rank, name_clean) %>%
    summarise(
      n_detected = n_distinct(sample[(detected %in% TRUE)]),
      n_present = n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]),
      mean_rpmm_detected = if_else(any(detected %in% TRUE), mean(rpmm[detected %in% TRUE], na.rm = TRUE), NA_real_),
      median_rpmm_detected = if_else(any(detected %in% TRUE), median(rpmm[detected %in% TRUE], na.rm = TRUE), NA_real_),
      mean_log_rpmm_detected = if_else(any(detected %in% TRUE), mean(log_rpmm[detected %in% TRUE], na.rm = TRUE), NA_real_),
      sd_log_rpmm_detected = if_else(sum(detected %in% TRUE, na.rm = TRUE) >= 2, sd(log_rpmm[detected %in% TRUE], na.rm = TRUE), NA_real_),
      median_rpmm_present = if_else(any(is.finite(reads_clade) & reads_clade > 0), median(rpmm[is.finite(reads_clade) & reads_clade > 0], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      prevalence_any = ifelse(use_cohort_prevalence, n_present / n_samples, NA_real_),
      prevalence = ifelse(use_cohort_prevalence, n_detected / n_samples, NA_real_),
      cv_detected = if_else(
        is.finite(mean_rpmm_detected) & mean_rpmm_detected > 0,
        sqrt(pmax(exp(sd_log_rpmm_detected^2) - 1, 0)),
        NA_real_
      ),
      cv_cap = 3,
      ubiquity = ifelse(
        use_cohort_prevalence,
        cap(
          pmax(coalesce0(prevalence), coalesce0(prevalence_any)) *
            (1 - pmin(coalesce0(cv_detected) / cv_cap, 1)),
          0,
          1
        ),
        NA_real_
      ),
      cohort_bg_rpmm = ifelse(
        use_cohort_prevalence,
        0.25 * coalesce0(prevalence_any) * coalesce0(median_rpmm_present),
        NA_real_
      )
    )

  if (n_controls > 0) {
    control_tax <- gs_long %>%
      filter(base_filt, sample %in% control_samples) %>%
      group_by(rank, name_clean) %>%
      summarise(
        control_prevalence = if_else(
          use_control_prevalence,
          n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]) / n_controls,
          NA_real_
        ),
        mean_rpmm_control_detected = if_else(
          any(is.finite(reads_clade) & reads_clade > 0),
          mean(rpmm[is.finite(reads_clade) & reads_clade > 0], na.rm = TRUE),
          NA_real_
        ),
        median_rpmm_control_detected = if_else(
          any(is.finite(reads_clade) & reads_clade > 0),
          median(rpmm[is.finite(reads_clade) & reads_clade > 0], na.rm = TRUE),
          NA_real_
        ),
        .groups = "drop"
      )
    
    tax_base <- tax_base %>%
      left_join(control_tax, by = c("rank", "name_clean")) %>%
      mutate(
        control_prevalence = if_else(use_control_prevalence, coalesce(control_prevalence, 0), NA_real_),
        mean_rpmm_control_detected = case_when(
          use_control_prevalence & control_prevalence == 0 ~ 0,
          is.finite(mean_rpmm_control_detected) ~ mean_rpmm_control_detected,
          TRUE ~ NA_real_
        ),
        median_rpmm_control_detected = case_when(
          use_control_prevalence & control_prevalence == 0 ~ 0,
          is.finite(median_rpmm_control_detected) ~ median_rpmm_control_detected,
          TRUE ~ NA_real_
        )
      )
  } else {
    tax_base <- tax_base %>%
      mutate(
        control_prevalence = NA_real_,
        mean_rpmm_control_detected = NA_real_,
        median_rpmm_control_detected = NA_real_
      )
  }

  tax_cor <- gs_long %>%
    filter(base_filt) %>%
    group_by(rank, name_clean) %>%
    group_modify(~{
      d <- .x

      if (n_distinct(d$sample) < 3) {
        return(tibble(cor_with_biomass = NA_real_, cor_with_nonmicrobe_frac = NA_real_))
      }

      cor_biomass <- suppressWarnings(cor(d$log_microbial_reads, d$log_rpmm, use = "pairwise.complete.obs"))
      cor_nonmicrobe <- suppressWarnings(cor(d$nonmicrobe_frac, d$log_rpmm, use = "pairwise.complete.obs"))

      tibble(cor_with_biomass = cor_biomass, cor_with_nonmicrobe_frac = cor_nonmicrobe)
    }) %>%
    ungroup()

  tax_base <- tax_base %>%
    left_join(tax_cor, by = c("rank", "name_clean"))

  if (isTRUE(REF_BG$available) && !is.null(REF_BG$stats) && nrow(REF_BG$stats) > 0) {
    tax_base <- tax_base %>%
      left_join(REF_BG$stats, by = c("rank", "name_clean"))
  } else {
    tax_base <- tax_base %>%
      mutate(
        ref_n = NA_real_,
        ref_prev = NA_real_,
        ref_q50_log = NA_real_,
        ref_q95_log = NA_real_,
        ref_q99_log = NA_real_,
        ref_med_rpmm = NA_real_,
        ref_q95_rpmm = NA_real_
      )
  }

  tax_base <- tax_base %>%
    mutate(
      in_user_contam = name_clean %in% user_contam,
      in_kitome = name_clean %in% kitome,
      in_clinical_panel = name_clean %in% clinical_panel
    )

  if (has_kmers) {
    kmer_med <- gs_long %>%
      filter(base_filt) %>%
      group_by(rank, name_clean) %>%
      summarise(
        kmer_med = median(kmer_per_read, na.rm = TRUE),
        uniq_med = median(uniq_kmers, na.rm = TRUE),
        .groups = "drop"
      )

    tax_base <- tax_base %>%
      left_join(kmer_med, by = c("rank", "name_clean"))
  } else {
    tax_base$kmer_med <- NA_real_
    tax_base$uniq_med <- NA_real_
  }

  tax_base <- tax_base %>%
    mutate(bg_rpmm = pmax(coalesce0(median_rpmm_control_detected), coalesce0(ref_med_rpmm), coalesce0(cohort_bg_rpmm)))

  sample_priors <- compute_sample_priors(gs_long, meta_df, params, ref_comp_genus = REF_BG$comp_genus)

  list(
    tax_features = tax_base,
    gs_long = gs_long,
    sample_priors = sample_priors
  )
}


log1p_nonnegative_score <- function(x) {
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  out <- rep(
    NA_real_,
    length(x)
  )
  
  ok <- is.finite(x)
  
  out[ok] <- log1p(
    pmax(
      x[ok],
      0
    )
  )
  
  out
}

compute_final_score_features <- function(out) {
  log_reads <- log1p_nonnegative_score(
    out$reads_clade
  )
  
  log_rpmm_score <- log1p_nonnegative_score(
    out$rpmm
  )
  
  log_uniq <- log1p_nonnegative_score(
    out$uniq_kmers
  )
  
  log_kpr <- log1p_nonnegative_score(
    out$kmer_per_read
  )
  
  log_decon <- log1p_nonnegative_score(
    out$decon_rpmm
  )
  
  log_ref_median <- log1p_nonnegative_score(
    out$ref_med_rpmm
  )
  
  ref_available <-
    is.finite(out$log_rpmm) &
    is.finite(out$ref_q50_log) &
    is.finite(out$ref_q95_log) &
    is.finite(out$ref_q99_log)
  
  reference_at_or_below_q99 <- ifelse(
    ref_available,
    as.numeric(
      out$log_rpmm <=
        out$ref_q99_log
    ),
    NA_real_
  )
  
  reference_at_or_below_q95 <- ifelse(
    ref_available,
    as.numeric(
      out$log_rpmm <=
        out$ref_q95_log
    ),
    NA_real_
  )
  
  reference_at_or_below_q50 <- ifelse(
    ref_available,
    as.numeric(
      out$log_rpmm <=
        out$ref_q50_log
    ),
    NA_real_
  )
  
  reference_prevalence <- ifelse(
    is.finite(out$ref_prev),
    cap(
      out$ref_prev,
      0,
      1
    ),
    NA_real_
  )
  
  reference_relative_abundance <- ifelse(
    is.finite(out$rpmm) &
      is.finite(out$ref_med_rpmm),
    log_ref_median -
      log_rpmm_score,
    NA_real_
  )
  
  ambiguity_unresolved_fraction <- ifelse(
    is.finite(
      out$ambiguity_unresolved_fraction
    ),
    cap(
      out$ambiguity_unresolved_fraction,
      0,
      1
    ),
    NA_real_
  )
  
  ambiguity_species_nondominance <- ifelse(
    is.finite(
      out$ambiguity_species_nondominance
    ),
    cap(
      out$ambiguity_species_nondominance,
      0,
      1
    ),
    NA_real_
  )
  
  decontamination_ratio_raw <- ifelse(
    is.finite(out$decon_ratio),
    cap(
      out$decon_ratio,
      0,
      1
    ),
    NA_real_
  )
  
  biomass_scaled <- ifelse(
    is.finite(out$biomass_scaled),
    cap(
      out$biomass_scaled,
      0,
      1
    ),
    NA_real_
  )
  
  js_div <- ifelse(
    is.finite(out$js_div),
    cap(
      out$js_div,
      0,
      1
    ),
    NA_real_
  )
  
  in_clinical <-
    !is.na(out$in_clinical_panel) &
    out$in_clinical_panel
  
  in_kitome <-
    !is.na(out$in_kitome) &
    out$in_kitome
  
  clinical <- as.numeric(
    in_clinical
  )
  
  kitome_only <- as.numeric(
    in_kitome &
      !in_clinical
  )
  
  kitome_clinical_overlap <- as.numeric(
    in_kitome &
      in_clinical
  )
  
  low_enrichment <- ifelse(
    is.finite(
      reference_relative_abundance
    ),
    pmax(
      reference_relative_abundance,
      0
    ),
    NA_real_
  )
  
  tibble(
    reads_support =
      -log_reads,
    
    rpmm_support =
      -log_rpmm_score,
    
    uniq_kmer_support =
      -log_uniq,
    
    kmer_per_read_support =
      -log_kpr,
    
    ambiguity_unresolved_fraction =
      ambiguity_unresolved_fraction,
    
    ambiguity_species_nondominance =
      ambiguity_species_nondominance,
    
    reference_prevalence =
      reference_prevalence,
    
    reference_median_abundance =
      log_ref_median,
    
    reference_relative_abundance =
      reference_relative_abundance,
    
    reference_at_or_below_q99 =
      reference_at_or_below_q99,
    
    reference_at_or_below_q95 =
      reference_at_or_below_q95,
    
    reference_at_or_below_q50 =
      reference_at_or_below_q50,
    
    decontaminated_abundance =
      -log_decon,
    
    decontamination_ratio =
      -decontamination_ratio_raw,
    
    kitome_only =
      kitome_only,
    
    kitome_clinical_overlap =
      kitome_clinical_overlap,
    
    clinical_membership =
      -clinical,
    
    low_biomass_context =
      -biomass_scaled,
    
    reference_like_sample_composition =
      -js_div,
    
    clinical_x_decon_support =
      -clinical *
      log_decon,
    
    kitome_only_x_decon_support =
      -kitome_only *
      log_decon,
    
    kitome_overlap_x_decon_support =
      -kitome_clinical_overlap *
      log_decon,
    
    reference_prevalence_x_low_enrichment =
      reference_prevalence *
      low_enrichment,
    
    unresolved_ambiguity_x_kmer_support =
      -ambiguity_unresolved_fraction *
      log_uniq,
    
    species_nondominance_x_kmer_support =
      -ambiguity_species_nondominance *
      log_uniq
  )
}

compute_final_score <- function(out) {
  model <- final_score_model()
  
  features <- compute_final_score_features(
    out
  )
  
  x <- as.matrix(
    features[
      ,
      model$feature_names,
      drop = FALSE
    ]
  )
  
  storage.mode(x) <- "double"
  
  for (
    j in seq_along(
      model$feature_names
    )
  ) {
    bad <- !is.finite(
      x[, j]
    )
    
    if (any(bad)) {
      x[bad, j] <-
        model$impute_values[j]
    }
  }
  
  z <- sweep(
    x,
    2,
    model$means,
    FUN = "-"
  )
  
  z <- sweep(
    z,
    2,
    model$scales,
    FUN = "/"
  )
  
  contributions <- sweep(
    z,
    2,
    model$standardized_weights,
    FUN = "*"
  )
  
  eta <-
    model$intercept +
    rowSums(
      contributions
    )
  
  mapped_logit <-
    eta +
    model$score_offset
  
  contam_prob <- plogis(
    mapped_logit
  )
  
  fp_score <-
    100 *
    contam_prob
  
  colnames(
    contributions
  ) <- paste0(
    "score_contrib_",
    model$feature_names
  )
  
  list(
    features = features,
    
    standardized_features = z,
    
    contributions = as.data.frame(
      contributions,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    eta = eta,
    
    mapped_logit = mapped_logit,
    
    contam_prob = contam_prob,
    
    fp_score = fp_score
  )
}

apply_calls <- function(gs_long, tax_features, meta_df, params, sample_priors) {
  meta2 <- meta_df %>%
    mutate(group = if_else(is.na(group) | !nzchar(group), "Group 1", group))

  group_sizes <- meta2 %>%
    count(group, name = "n_group")

  tf_keep <- tax_features %>%
    select(
      rank,
      name_clean,
      prevalence,
      prevalence_any,
      control_prevalence,
      ubiquity,
      cor_with_biomass,
      cor_with_nonmicrobe_frac,
      in_user_contam,
      in_kitome,
      in_clinical_panel,
      ref_n,
      ref_prev,
      ref_q50_log,
      ref_q95_log,
      ref_q99_log,
      ref_med_rpmm,
      median_rpmm_control_detected,
      cohort_bg_rpmm,
      bg_rpmm
    )

  out <- gs_long %>%
    left_join(tf_keep, by = c("rank", "name_clean")) %>%
    left_join(group_sizes, by = "group") %>%
    left_join(sample_priors %>% select(sample, biomass_scaled, js_div), by = "sample")

  out <- out %>%
    mutate(
      ref_n = suppressWarnings(as.numeric(ref_n)),
      ref_prev = suppressWarnings(as.numeric(ref_prev)),
      ref_q50_log = suppressWarnings(as.numeric(ref_q50_log)),
      ref_q95_log = suppressWarnings(as.numeric(ref_q95_log)),
      ref_q99_log = suppressWarnings(as.numeric(ref_q99_log)),
      ref_med_rpmm = suppressWarnings(as.numeric(ref_med_rpmm)),
      median_rpmm_control_detected = suppressWarnings(as.numeric(median_rpmm_control_detected)),
      cohort_bg_rpmm = suppressWarnings(as.numeric(cohort_bg_rpmm)),
      bg_rpmm = suppressWarnings(as.numeric(bg_rpmm))
    )

  out <- out %>%
    mutate(
      reference_median_rpmm = if_else(is.finite(ref_med_rpmm), ref_med_rpmm, NA_real_),
      control_median_rpmm = if_else(is.finite(median_rpmm_control_detected), median_rpmm_control_detected, NA_real_),
      log2FC_vs_reference_median = log2((coalesce0(rpmm) + params$fc_pseudocount) / (coalesce0(reference_median_rpmm) + params$fc_pseudocount)),
      log2FC_vs_control_median = log2((coalesce0(rpmm) + params$fc_pseudocount) / (coalesce0(control_median_rpmm) + params$fc_pseudocount)),
      enriched_vs_controls = if_else(is.finite(control_median_rpmm), log2FC_vs_control_median > 0, NA),
      qc_pass = is_microbial & !is_host & !is_plant
    )

  amb <- compute_ambiguity_index(out)

  out <- out %>%
    left_join(amb, by = c("sample", "rank", "name_clean"))

  out <- out %>%
    mutate(
      bg_rpmm = coalesce0(bg_rpmm),
      decon_source_rpmm = coalesce0(rpmm),
      decon_rpmm = pmax(decon_source_rpmm - bg_rpmm, 0),
      decon_ratio = if_else(is.finite(rpmm) & rpmm > 0, decon_rpmm / rpmm, NA_real_),
      log_decon_rpmm = safe_log(decon_rpmm, eps = params$eps)
    )

  score_result <- compute_final_score(
    out
  )
  
  out <- bind_cols(
    out,
    as_tibble(
      score_result$contributions
    )
  ) %>%
    mutate(
      final_model_logit =
        score_result$eta,
      
      score_mapped_logit =
        score_result$mapped_logit,
      
      contam_prob =
        score_result$contam_prob,
      
      fp_score =
        score_result$fp_score,
      
      fp_falsepos_flag =
        fp_score >=
        params$fp_falsepos_cutoff,
      
      fp_true_flag =
        fp_score <=
        params$fp_true_cutoff
    )

  out <- out %>%
    mutate(
      call = case_when(
        !is_microbial | is_host | is_plant ~ "Non-microbial / Host",
        fp_true_flag ~ "Likely true",
        fp_falsepos_flag ~ "Likely false positive / background",
        TRUE ~ "Uncertain"
      ),
      call_reason = case_when(
        !is_microbial | is_host | is_plant ~
          "Excluded as non-microbial or host",
        
        call == "Likely true" ~ paste0(
          "Score ",
          round(fp_score, 1),
          " <= likely true cutoff ",
          round(params$fp_true_cutoff, 1),
          "; calibrated evidence supports a likely true microbial signal"
        ),
        
        call == "Likely false positive / background" ~ paste0(
          "Score ",
          round(fp_score, 1),
          " >= false positive cutoff ",
          round(params$fp_falsepos_cutoff, 1),
          "; calibrated evidence supports background or contamination"
        ),
        
        TRUE ~ paste0(
          "Score ",
          round(fp_score, 1),
          " is between ",
          round(params$fp_true_cutoff, 1),
          " and ",
          round(params$fp_falsepos_cutoff, 1),
          "; evidence is mixed between true-signal and background patterns"
        )
      )
    )

  out
}

# ------------------- DECONTAMINATION (KRAKEN REPORT CLEANING) -------------------
get_kraken_tree_info <- function(raw_df) {
  name_col <- raw_df[[ncol(raw_df)]]
  name_raw <- as.character(name_col)
  name_clean <- trim_left(name_raw)

  lead_spaces <- nchar(name_raw) - nchar(name_clean)
  lead_spaces[!is.finite(lead_spaces)] <- 0
  depth <- floor(pmax(lead_spaces, 0) / 2)

  parent_idx <- rep(NA_integer_, length(depth))
  stack <- integer(0)

  for (i in seq_along(depth)) {
    d <- depth[i]

    if (length(stack) > d) {
      stack <- stack[seq_len(d)]
    }

    parent_idx[i] <- if (d <= 0 || length(stack) < d) NA_integer_ else stack[d]

    if (length(stack) < d + 1) {
      stack <- c(stack, rep(NA_integer_, d + 1 - length(stack)))
    }

    stack[d + 1] <- i
  }

  list(
    name_raw = name_raw,
    name_clean = name_clean,
    depth = depth,
    parent_idx = parent_idx
  )
}

compute_remove_taxa_for_sample <- function(gs, quick, protected_taxa = character(0)) {
  remove_taxa <- character(0)

  if ("false" %in% quick) {
    remove_taxa <- c(remove_taxa, gs$name_clean[gs$call == "Likely false positive / background"])
  }
  if ("unc" %in% quick) {
    remove_taxa <- c(remove_taxa, gs$name_clean[gs$call == "Uncertain"])
  }
  if ("true" %in% quick) {
    remove_taxa <- c(remove_taxa, gs$name_clean[gs$call == "Likely true"])
  }

  remove_taxa <- unique(remove_taxa)
  remove_taxa <- remove_taxa[!is.na(remove_taxa) & nzchar(remove_taxa)]
  setdiff(remove_taxa, protected_taxa)
}

compute_keep_taxa_for_sample <- function(gs, remove_taxa, protected_taxa = character(0)) {
  keep_taxa <- gs %>%
    filter(rank %in% c("G", "S"), call == "Likely true") %>%
    pull(name_clean) %>%
    unique()

  keep_taxa <- setdiff(keep_taxa, remove_taxa)
  keep_taxa <- unique(c(keep_taxa, protected_taxa))
  keep_taxa <- keep_taxa[!is.na(keep_taxa) & nzchar(keep_taxa)]
  keep_taxa
}

retain_kraken_raw <- function(raw_df, keep_taxa, keep_descendants = FALSE, keep_root = TRUE, keep_unclassified = TRUE) {
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    return(raw_df)
  }

  info <- get_kraken_tree_info(raw_df)
  name_clean <- info$name_clean
  depth <- info$depth
  parent_idx <- info$parent_idx

  keep_taxa <- unique(trimws(as.character(keep_taxa)))
  keep_taxa <- keep_taxa[nzchar(keep_taxa)]

  keep_idx <- rep(FALSE, nrow(raw_df))
  target_idx <- which(name_clean %in% keep_taxa)

  for (i in target_idx) {
    keep_idx[i] <- TRUE

    p <- parent_idx[i]
    while (is.finite(p)) {
      keep_idx[p] <- TRUE
      p <- parent_idx[p]
    }

    if (isTRUE(keep_descendants)) {
      d <- depth[i]
      j <- i + 1
      while (j <= nrow(raw_df) && depth[j] > d) {
        keep_idx[j] <- TRUE
        j <- j + 1
      }
    }
  }

  rank_col <- if (ncol(raw_df) >= 8) 6 else if (ncol(raw_df) >= 6) 4 else NA_integer_
  taxid_col <- if (ncol(raw_df) >= 8) 7 else if (ncol(raw_df) >= 6) 5 else NA_integer_

  ranks <- if (is.finite(rank_col)) as.character(raw_df[[rank_col]]) else rep("", nrow(raw_df))
  taxids <- if (is.finite(taxid_col)) suppressWarnings(as.integer(raw_df[[taxid_col]])) else rep(NA_integer_, nrow(raw_df))

  special_idx <- which(
    (keep_root & (ranks %in% c("R") | taxids == 1L | name_clean == "root")) |
      (keep_unclassified & (ranks %in% c("U") | taxids == 0L | name_clean == "unclassified"))
  )
  keep_idx[special_idx] <- TRUE

  reads_direct <- suppressWarnings(as.numeric(raw_df[[3]]))
  reads_direct[!is.finite(reads_direct)] <- 0
  reads_direct[!keep_idx] <- 0

  reads_clade <- reads_direct
  for (i in rev(seq_len(nrow(raw_df)))) {
    if (!keep_idx[i]) next
    p <- parent_idx[i]
    if (is.finite(p) && keep_idx[p]) {
      reads_clade[p] <- reads_clade[p] + reads_clade[i]
    }
  }

  kept_rows <- which(keep_idx)
  if (length(kept_rows) == 0) {
    return(raw_df[0, , drop = FALSE])
  }

  special_keep <- kept_rows %in% special_idx
  nonzero_keep <- reads_clade[kept_rows] > 0 | reads_direct[kept_rows] > 0 | special_keep
  kept_rows <- kept_rows[nonzero_keep]
  if (length(kept_rows) == 0) {
    return(raw_df[0, , drop = FALSE])
  }

  root_rows <- kept_rows[ranks[kept_rows] %in% c("R") | taxids[kept_rows] == 1L | name_clean[kept_rows] == "root"]
  unclassified_rows <- kept_rows[ranks[kept_rows] %in% c("U") | taxids[kept_rows] == 0L | name_clean[kept_rows] == "unclassified"]

  total_reads <- 0
  if (length(root_rows) > 0) {
    total_reads <- total_reads + reads_clade[root_rows[1]]
  }
  if (length(unclassified_rows) > 0) {
    total_reads <- total_reads + reads_clade[unclassified_rows[1]]
  }
  if (!is.finite(total_reads) || total_reads <= 0) {
    total_reads <- sum(reads_direct[kept_rows], na.rm = TRUE)
  }
  if (!is.finite(total_reads) || total_reads <= 0) {
    total_reads <- 1
  }

  cleaned <- raw_df[kept_rows, , drop = FALSE]
  cleaned[[1]] <- round(100 * reads_clade[kept_rows] / total_reads, 6)
  cleaned[[2]] <- as.numeric(reads_clade[kept_rows])
  cleaned[[3]] <- as.numeric(reads_direct[kept_rows])

  cleaned
}

build_cleaned_report_for_sample <- function(raw_df, gs, quick, protected_taxa = character(0), keep_descendants = FALSE) {
  remove_taxa <- compute_remove_taxa_for_sample(gs, quick = quick, protected_taxa = protected_taxa)
  keep_taxa <- compute_keep_taxa_for_sample(gs, remove_taxa = remove_taxa, protected_taxa = protected_taxa)
  retain_kraken_raw(raw_df, keep_taxa = keep_taxa, keep_descendants = keep_descendants)
}

zip_files_safely <- function(files, zipfile) {
  files <- files[file.exists(files)]
  if (length(files) == 0) stop("No files to zip.")

  if (HAS_ZIP) {
    zip::zipr(zipfile, files = files)
    return(zipfile)
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)

  tmpdir <- tempfile("celltaminate_zip_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  file.copy(files, tmpdir, overwrite = TRUE)
  setwd(tmpdir)
  utils::zip(zipfile = zipfile, files = basename(files))
  zipfile
}



spider_plot_layout <- function(taxa) {
  taxa <- as.character(taxa)
  taxa <- taxa[!is.na(taxa)]
  max_chars <- if (length(taxa) > 0) max(nchar(taxa), na.rm = TRUE) else 10
  label_cex <- max(0.45, min(0.72, 20 / max(20, max_chars)))
  label_radius <- min(1.28, 1.12 + (max_chars / 120))
  list(
    mar = c(2.2, 2.2, 4.0, 2.2),
    label_cex = label_cex,
    label_radius = label_radius
  )
}

draw_spider_labels <- function(taxa, radius = 1.16, cex = 0.6) {
  n <- length(taxa)
  if (n == 0) return(invisible(NULL))

  angles_deg <- seq(90, 90 - 360 + 360 / n, length.out = n)
  angles_rad <- angles_deg * pi / 180

  x <- radius * cos(angles_rad)
  y <- radius * sin(angles_rad)

  srt <- ifelse(angles_deg < -90, angles_deg + 180, ifelse(angles_deg > 90, angles_deg - 180, angles_deg))
  adj_x <- ifelse(cos(angles_rad) > 0.15, 0, ifelse(cos(angles_rad) < -0.15, 1, 0.5))

  for (i in seq_along(taxa)) {
    text(
      x = x[i],
      y = y[i],
      labels = taxa[i],
      srt = srt[i],
      cex = cex,
      adj = c(adj_x[i], 0.5),
      xpd = NA
    )
  }

  invisible(NULL)
}


# ------------------- COMMAND LINE HELPERS -------------------
empty_plot <- function(msg) {
  ggplot() +
    theme_void() +
    xlim(0, 1) +
    ylim(0, 1) +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 5)
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  x <- tolower(trimws(as.character(x)[1]))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

split_csv <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  vals <- unlist(strsplit(as.character(x)[1], ",", fixed = TRUE))
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

HELP_TEXT <- paste(
  "Usage:",
  "  Rscript celltaminate_cli.R --input_files sample1.txt,sample2.txt --output_dir results [options]",
  "",
  "Required:",
  "  --output_dir PATH",
  "  --input_files PATH1,PATH2,...   or   --input_list file_with_paths.txt",
  "",
  "Optional:",
  "  --metadata_tsv PATH",
  "  --tax_level G|S",
  "  --host_species human,mouse,drosophila",
  "  --host_taxids_manual 9606,10090",
  "  --include_protozoa TRUE|FALSE",
  "  --collapse_species TRUE|FALSE",
  "  --collapse_min_genus_reads INTEGER",
  "  --collapse_top_frac NUMERIC",
  "  --fp_falsepos_cutoff NUMERIC",
  "  --fp_true_cutoff NUMERIC",
  "  --ref_bg_tsv PATH",
  "  --kitome_blacklist_tsv PATH",
  "  --clinical_panel_tsv PATH",
  "  --quick_remove false,unc,true",
  "  --protected_taxa_tsv PATH",
  "  --keep_descendants TRUE|FALSE",
  "  --show_fp_breakdown TRUE|FALSE",
  "  --save_cleaned_reanalysis TRUE|FALSE",
  sep = "\n"
)

parse_cli_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)

  if (length(x) == 0 || any(x %in% c("-h", "--help"))) {
    cat(HELP_TEXT, "\n")
    quit(save = "no", status = 0)
  }

  if (length(x) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }

  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    val <- x[i + 1]
    if (!startsWith(key, "--")) {
      stop(paste0("Invalid argument name: ", key), call. = FALSE)
    }
    key <- sub("^--", "", key)
    out[[key]] <- val
    i <- i + 2
  }
  out
}

require_arg <- function(args, key) {
  val <- args[[key]] %||% ""
  if (!nzchar(val)) {
    stop(paste0("Missing required argument --", key), call. = FALSE)
  }
  val
}

make_unique_names <- function(nms) {
  if (length(nms) == 0) return(nms)
  dup <- duplicated(nms)
  if (!any(dup)) return(nms)

  out <- nms
  counts <- list()

  for (i in seq_along(out)) {
    nm <- out[i]
    if (is.null(counts[[nm]])) counts[[nm]] <- 0
    counts[[nm]] <- counts[[nm]] + 1
    if (counts[[nm]] > 1) {
      stem <- tools::file_path_sans_ext(nm)
      ext <- tools::file_ext(nm)
      if (nzchar(ext)) {
        out[i] <- paste0(stem, "_", counts[[nm]], ".", ext)
      } else {
        out[i] <- paste0(stem, "_", counts[[nm]])
      }
    }
  }

  out
}

collect_input_reports <- function(input_files = NULL, input_list = NULL) {
  files <- character(0)

  if (!is.null(input_list) && nzchar(input_list)) {
    if (!file.exists(input_list)) {
      stop(paste0("Input list file not found: ", input_list), call. = FALSE)
    }
    files <- c(files, read_optional_lines(input_list))
  }

  if (!is.null(input_files) && nzchar(input_files)) {
    files <- c(files, split_csv(input_files))
  }

  files <- unique(trimws(files))
  files <- files[nzchar(files)]

  if (length(files) == 0) {
    stop("No input report files were provided.", call. = FALSE)
  }

  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop(
      paste0(
        "The following input files do not exist:\n",
        paste(missing, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  sample_names <- make_unique_names(basename(files))
  out <- as.list(files)
  names(out) <- sample_names
  out
}

build_host_taxids <- function(host_species = c("human"), host_taxids_manual = "") {
  ids <- integer(0)

  if ("human" %in% host_species) ids <- c(ids, 9606L, 9605L)
  if ("mouse" %in% host_species) ids <- c(ids, 10090L)
  if ("drosophila" %in% host_species) ids <- c(ids, 7227L)

  manual <- gsub("[^0-9,]", "", host_taxids_manual %||% "")
  if (nzchar(manual)) {
    parts <- unlist(strsplit(manual, ",", fixed = TRUE))
    parts <- parts[nzchar(parts)]
    ids <- c(ids, suppressWarnings(as.integer(parts)))
  }

  ids <- ids[is.finite(ids)]
  unique(ids)
}

build_host_name_patterns <- function(host_species = c("human")) {
  pats <- character(0)

  if ("human" %in% host_species) {
    pats <- c(
      pats,
      "\\bHomo\\b",
      "\\bHomo sapiens\\b",
      "\\bPrimates\\b",
      "\\bHominidae\\b",
      "\\bChordata\\b",
      "\\bMammalia\\b",
      "\\bMetazoa\\b"
    )
  }

  if ("mouse" %in% host_species) {
    pats <- c(
      pats,
      "\\bMus\\b",
      "\\bMus musculus\\b",
      "\\bRodentia\\b"
    )
  }

  if ("drosophila" %in% host_species) {
    pats <- c(
      pats,
      "\\bDrosophila\\b",
      "\\bDrosophilidae\\b",
      "\\bInsecta\\b"
    )
  }

  unique(pats)
}

default_metadata <- function(sample_names) {
  tibble(
    sample = as.character(sample_names),
    group = "Group 1",
    sample_type = "Clinical / sterile",
    is_control = FALSE
  )
}

load_metadata <- function(metadata_tsv, sample_names) {
  if (is.null(metadata_tsv) || !nzchar(metadata_tsv)) {
    return(default_metadata(sample_names))
  }

  if (!file.exists(metadata_tsv)) {
    stop(paste0("Metadata file not found: ", metadata_tsv), call. = FALSE)
  }

  meta <- tryCatch(
    read.delim(metadata_tsv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", quote = ""),
    error = function(e) NULL
  )

  if (is.null(meta) || nrow(meta) == 0) {
    stop("Metadata table could not be read or is empty.", call. = FALSE)
  }

  if (!"sample" %in% colnames(meta)) {
    stop("Metadata table must contain a 'sample' column.", call. = FALSE)
  }

  meta$sample <- as.character(meta$sample)

  if (!"group" %in% colnames(meta)) meta$group <- "Group 1"
  if (!"sample_type" %in% colnames(meta)) meta$sample_type <- "Clinical / sterile"
  if (!"is_control" %in% colnames(meta)) meta$is_control <- FALSE

  meta <- meta %>%
    transmute(
      sample = as.character(sample),
      group = as.character(group),
      sample_type = as.character(sample_type),
      is_control = tolower(trimws(as.character(is_control))) %in% c("true", "t", "1", "yes", "y")
    ) %>%
    distinct(sample, .keep_all = TRUE)

  base_meta <- default_metadata(sample_names)
  meta <- base_meta %>%
    left_join(meta, by = "sample", suffix = c(".default", "")) %>%
    transmute(
      sample = sample,
      group = dplyr::coalesce(group, group.default),
      sample_type = dplyr::coalesce(sample_type, sample_type.default),
      is_control = dplyr::coalesce(is_control, is_control.default)
    )

  meta
}

load_protected_taxa_map <- function(path) {
  if (is.null(path) || !nzchar(path)) return(list())
  if (!file.exists(path)) {
    stop(paste0("Protected taxa file not found: ", path), call. = FALSE)
  }

  df <- tryCatch(
    read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", quote = ""),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0) return(list())
  if (!"sample" %in% colnames(df)) stop("Protected taxa TSV must contain a 'sample' column.", call. = FALSE)

  taxon_col <- NULL
  for (nm in c("taxon", "name_clean", "name", "organism_name")) {
    if (nm %in% colnames(df)) {
      taxon_col <- nm
      break
    }
  }

  if (is.null(taxon_col)) {
    stop("Protected taxa TSV must contain one of: taxon, name_clean, name, organism_name.", call. = FALSE)
  }

  df <- df %>%
    transmute(
      sample = as.character(sample),
      taxon = trimws(as.character(.data[[taxon_col]]))
    ) %>%
    filter(nzchar(sample), nzchar(taxon))

  split(df$taxon, df$sample)
}

build_params <- function(args) {
  list(
    analysis_mode = "Cell-line",
    
    tax_level = toupper(
      args$tax_level %||% "S"
    ),
    
    eps = 1e-9,
    
    min_reads_prevalence = 5,
    
    min_rpmm_prevalence = 0.2,
    
    min_uniq_kmers_prevalence = 96,
    
    min_kmer_per_read_prevalence = 0.4,
    
    fp_falsepos_cutoff = as.numeric(
      args$fp_falsepos_cutoff %||% 75
    ),
    
    fp_true_cutoff = as.numeric(
      args$fp_true_cutoff %||% 1
    ),
    
    fc_pseudocount = 0.1
  )
}

init_reference_tables <- function(args) {
  REF_BG <<- if (!is.null(args$ref_bg_tsv) && nzchar(args$ref_bg_tsv)) {
    load_ref_cell_lines(path = args$ref_bg_tsv, eps = 1e-9)
  } else {
    list(available = FALSE, path = NULL, ref = NULL, stats = NULL, comp_genus = NULL)
  }

  KITOME_BG_BLACKLIST <<- if (!is.null(args$kitome_blacklist_tsv) && nzchar(args$kitome_blacklist_tsv)) {
    load_kitome_blacklist(path = args$kitome_blacklist_tsv)
  } else {
    empty_kitome_blacklist_df()
  }

  CLINICALLY_IMPORTANT_PATHOGENS <<- if (!is.null(args$clinical_panel_tsv) && nzchar(args$clinical_panel_tsv)) {
    load_clinically_important_pathogens(path = args$clinical_panel_tsv)
  } else {
    empty_clinically_important_pathogens_df()
  }
}

parse_reports <- function(file_map, host_taxids, host_name_patterns, include_protozoa = TRUE) {
  sample_list <- list()
  errors <- character(0)

  for (nm in names(file_map)) {
    fpath <- file_map[[nm]]
    parsed <- read_kraken_report(fpath, sample_name = nm)

    if (!is.null(parsed$error)) {
      errors <- c(errors, parsed$error)
      next
    }

    df <- parsed$df
    df <- annotate_lineage(df)
    df <- infer_host_like(df, host_taxids = host_taxids, host_name_patterns = host_name_patterns)
    df <- infer_plant_like(df)
    df <- infer_microbial(df, include_protozoa = include_protozoa)
    df <- annotate_top_strain_rows(df)
    summ <- compute_reads_summary(df)

    sample_list[[nm]] <- list(
      sample = nm,
      path = fpath,
      report_format = parsed$report_format,
      raw = parsed$raw,
      df = df,
      summary = summ
    )
  }

  list(sample_list = sample_list, errors = errors)
}

analyze_sample_list <- function(sample_list, meta, prm, collapse_species_flag = TRUE,
                                collapse_min_genus_reads = 30, collapse_top_frac = 0.85,
                                clinical_orgs = character(0), kitome = character(0)) {
  gs_list <- list()
  sample_summary <- list()

  for (nm in names(sample_list)) {
    obj <- sample_list[[nm]]
    df <- obj$df
    summ <- obj$summary

    gs <- df %>%
      filter(rank %in% c("G", "S")) %>%
      mutate(
        sample = nm,
        total_reads = summ$total_reads,
        microbial_reads = summ$microbial_reads,
        host_reads = summ$host_reads,
        rpm = ifelse(is.finite(total_reads) & total_reads > 0, reads_clade / total_reads * 1e6, NA_real_),
        rpmm = ifelse(is_microbial & is.finite(microbial_reads) & microbial_reads > 0, reads_clade / microbial_reads * 1e6, NA_real_),
        log_rpmm = safe_log(rpmm, eps = prm$eps),
        genus = ifelse(rank == "G", name_clean, genus_ancestor)
      )

    gs_list[[nm]] <- gs

    n_genus_raw <- gs %>%
      filter(rank == "G", is_microbial, !is_host, !is_plant, reads_clade > 0) %>%
      pull(name_clean) %>%
      unique() %>%
      length()

    n_species_raw <- gs %>%
      filter(rank == "S", is_microbial, !is_host, !is_plant, reads_clade > 0) %>%
      pull(name_clean) %>%
      unique() %>%
      length()

    sp_rpmm <- gs %>%
      filter(rank == "S", is_microbial, !is_host, !is_plant, !is.na(rpmm), rpmm > 0)

    shannon <- calc_shannon((sp_rpmm$rpmm %||% numeric(0)) / 1e6)

    sample_summary[[nm]] <- tibble(
      sample = nm,
      report_format = obj$report_format,
      total_reads = summ$total_reads,
      microbial_reads = summ$microbial_reads,
      n_genus_raw = n_genus_raw,
      n_species_raw = n_species_raw,
      shannon = shannon
    )
  }

  gs_long <- bind_rows(gs_list)
  cohort_summary <- bind_rows(sample_summary)

  meta$sample <- as.character(meta$sample)
  meta$group <- as.character(meta$group)
  meta$sample_type <- as.character(meta$sample_type)
  if (!("is_control" %in% colnames(meta))) meta$is_control <- FALSE
  meta$is_control <- tolower(trimws(as.character(meta$is_control))) %in% c("true", "t", "1", "yes", "y")

  gs_long <- gs_long %>%
    left_join(meta, by = "sample")

  if (isTRUE(collapse_species_flag)) {
    gs_long <- collapse_species_within_genus(
      gs_long,
      min_genus_reads = collapse_min_genus_reads,
      top_frac_keep = collapse_top_frac
    )
  }

  engine <- compute_taxon_features(
    gs_long,
    meta,
    prm,
    user_contam = character(0),
    clinical_panel = clinical_orgs,
    kitome = kitome
  )

  tax_features <- engine$tax_features
  gs_long <- engine$gs_long
  sample_priors <- engine$sample_priors

  gs_long <- apply_calls(gs_long, tax_features, meta, prm, sample_priors)

  n_cohort_samples <- dplyr::n_distinct(meta$sample)
  prevalence_display <- gs_long %>%
    filter(rank %in% c("G", "S"), call != "Non-microbial / Host") %>%
    group_by(rank, name_clean) %>%
    summarise(
      cohort_prevalence_display = ifelse(
        n_cohort_samples > 0,
        n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]) / n_cohort_samples,
        NA_real_
      ),
      .groups = "drop"
    )
  gs_long <- gs_long %>%
    left_join(prevalence_display, by = c("rank", "name_clean"))

  tax_level_now <- prm$tax_level %||% "S"

  qc_counts2 <- gs_long %>%
    filter(call != "Non-microbial / Host") %>%
    group_by(sample) %>%
    summarise(
      n_genus_qc = n_distinct(name_clean[rank == "G" & qc_pass]),
      n_species_qc = n_distinct(name_clean[rank == "S" & qc_pass]),
      shannon = {
        sp <- rpmm[rank == "S" & qc_pass]
        sp <- sp[is.finite(sp) & sp > 0]
        if (length(sp) == 0) NA_real_ else calc_shannon(sp / 1e6)
      },
      n_true = sum(rank == tax_level_now & call == "Likely true", na.rm = TRUE),
      n_uncertain = sum(rank == tax_level_now & call == "Uncertain", na.rm = TRUE),
      n_contaminant = sum(rank == tax_level_now & call == "Likely false positive / background", na.rm = TRUE),
      .groups = "drop"
    )

  cohort_summary <- cohort_summary %>%
    select(-any_of(c("n_genus_qc", "n_species_qc", "shannon", "n_true", "n_uncertain", "n_contaminant"))) %>%
    left_join(qc_counts2, by = "sample")

  list(
    cohort_summary = cohort_summary,
    meta = meta,
    gs_long = gs_long,
    tax_features = tax_features,
    sample_priors = sample_priors
  )
}

run_pipeline <- function(file_map, meta, args, host_taxids, host_name_patterns) {
  parse_out <- parse_reports(
    file_map = file_map,
    host_taxids = host_taxids,
    host_name_patterns = host_name_patterns,
    include_protozoa = parse_bool(args$include_protozoa, TRUE)
  )

  sample_list <- parse_out$sample_list
  errors <- parse_out$errors

  if (length(sample_list) == 0) {
    stop(paste(c("No valid report files were parsed.", errors), collapse = "\n"), call. = FALSE)
  }

  meta <- load_metadata(args$metadata_tsv %||% "", names(sample_list))

  prm <- build_params(args)

  result <- analyze_sample_list(
    sample_list = sample_list,
    meta = meta,
    prm = prm,
    collapse_species_flag = parse_bool(args$collapse_species, TRUE),
    collapse_min_genus_reads = as.numeric(args$collapse_min_genus_reads %||% 30),
    collapse_top_frac = as.numeric(args$collapse_top_frac %||% 0.85),
    clinical_orgs = get_clinically_important_pathogen_names(CLINICALLY_IMPORTANT_PATHOGENS),
    kitome = get_kitome_blacklist_names(KITOME_BG_BLACKLIST)
  )

  result$sample_list <- sample_list
  result$parse_errors <- errors
  result$params <- prm
  result
}

# ------------------- OUTPUT HELPERS -------------------
build_display_taxa_table <- function(gs, show_fp_breakdown = FALSE) {
  if (nrow(gs) == 0) return(tibble())

  if (!"prevalence_any" %in% names(gs)) {
    gs$prevalence_any <- gs$prevalence
  }

  tax_level <- if (all(gs$rank == "S", na.rm = TRUE)) "S" else "G"

  if (tax_level == "S") {
    if (!all(c("Top_strain", "Top_strain_rank", "Top_strain_reads") %in% names(gs))) {
      gs <- gs %>%
        mutate(
          Top_strain = NA_character_,
          Top_strain_rank = NA_character_,
          Top_strain_reads = NA_integer_
        )
    }
  } else {
    gs <- gs %>%
      mutate(
        Top_strain = NA_character_,
        Top_strain_rank = NA_character_,
        Top_strain_reads = NA_integer_
      )
  }

  gs_view <- gs %>%
    mutate(
      Name = as.character(name_clean),
      Score_level = ifelse(tax_level == "S", "Species", "Genus"),
      Top_strain = as.character(Top_strain),
      Top_strain_rank = as.character(Top_strain_rank),
      Top_strain_reads = as.integer(Top_strain_reads),
      Reads = as.integer(round(reads_clade)),
      RPM = round(rpm, 1),
      RPMM = round(rpmm, 3),
      Decontaminated_RPMM = round(decon_rpmm, 3),
      Reference_median_RPMM = ifelse(
        is.finite(suppressWarnings(as.numeric(reference_median_rpmm))),
        round(suppressWarnings(as.numeric(reference_median_rpmm)), 3),
        NA_real_
      ),
      log2FC_vs_reference_median = round(log2FC_vs_reference_median, 2),
      Celltaminate_Score = round(fp_score, 1),
      Cohort_prevalence = round(cohort_prevalence_display, 3),
      Reference_prevalence = ifelse(
        is.finite(suppressWarnings(as.numeric(ref_prev))),
        round(suppressWarnings(as.numeric(ref_prev)), 3),
        NA_real_
      ),
      In_kitome = in_kitome,
      In_clinical_panel = in_clinical_panel,
      Call = as.character(call),
      Call_reason = ifelse(
        tax_level == "S" & !is.na(Top_strain) & nzchar(Top_strain),
        paste0(call_reason, " Display note: the top strain shown here inherits the species-level Celltaminate score."),
        call_reason
      )
    )

  base_cols <- c(
    "Name",
    "Score_level",
    if (tax_level == "S") c("Top_strain", "Top_strain_rank", "Top_strain_reads") else NULL,
    "Reads",
    "RPM",
    "RPMM",
    "Decontaminated_RPMM",
    "Reference_median_RPMM",
    "log2FC_vs_reference_median",
    "Celltaminate_Score",
    "Cohort_prevalence",
    "Reference_prevalence",
    "In_kitome",
    "In_clinical_panel",
    "Call",
    "Call_reason"
  )

  if (isTRUE(show_fp_breakdown)) {
    base_cols <- c(
      base_cols,
      "score_contrib_reads_support",
      "score_contrib_rpmm_support",
      "score_contrib_uniq_kmer_support",
      "score_contrib_ambiguity_species_nondominance",
      "score_contrib_reference_prevalence",
      "score_contrib_reference_median_abundance",
      "score_contrib_reference_relative_abundance",
      "score_contrib_reference_at_or_below_q99",
      "score_contrib_reference_at_or_below_q95",
      "score_contrib_reference_at_or_below_q50",
      "score_contrib_decontaminated_abundance",
      "score_contrib_decontamination_ratio",
      "score_contrib_kitome_clinical_overlap",
      "score_contrib_clinical_membership",
      "score_contrib_clinical_x_decon_support",
      "score_contrib_kitome_only_x_decon_support",
      "score_contrib_reference_prevalence_x_low_enrichment",
      "score_contrib_unresolved_ambiguity_x_kmer_support"
    )
  }

  gs_view %>% select(any_of(base_cols))
}

save_tables <- function(result, out_dir, show_fp_breakdown = FALSE) {
  tables_dir <- file.path(out_dir, "tables")
  per_sample_dir <- file.path(tables_dir, "per_sample")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)

  write.table(result$cohort_summary,
              file = file.path(tables_dir, "cohort_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result$meta,
              file = file.path(tables_dir, "sample_metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result$gs_long,
              file = file.path(tables_dir, "all_taxa_calls.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result$tax_features,
              file = file.path(tables_dir, "tax_features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(result$sample_priors,
              file = file.path(tables_dir, "sample_priors.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  prev_df <- make_prevalence_abundance_df(
    result$gs_long,
    tax_level = result$params$tax_level,
    total_samples = n_distinct(result$meta$sample)
  )
  if (ncol(prev_df) == 0) {
    prev_df <- tibble(
      name_clean = character(0),
      prevalence = numeric(0),
      mean_rpmm_detected = numeric(0),
      mean_log_rpmm_detected = numeric(0),
      call_bucket = character(0)
    )
  }
  write.table(prev_df,
              file = file.path(tables_dir, "prevalence_abundance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  for (samp in unique(result$gs_long$sample)) {
    sid <- safe_sample_id(samp)
    gs <- result$gs_long %>%
      filter(sample == samp, rank == result$params$tax_level, call != "Non-microbial / Host")
    display_df <- build_display_taxa_table(gs, show_fp_breakdown = show_fp_breakdown)
    write.table(display_df,
                file = file.path(per_sample_dir, paste0(sid, "_taxa_table.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  workbook_path <- file.path(out_dir, "celltaminate_summary.xlsx")
  writexl::write_xlsx(
    list(
      cohort_summary = result$cohort_summary,
      sample_metadata = result$meta,
      tax_features = result$tax_features,
      sample_priors = result$sample_priors,
      prevalence_abundance = prev_df
    ),
    path = workbook_path
  )
}

save_cohort_pca_plot <- function(result, file_path) {
  tax_level <- result$params$tax_level %||% "S"

  df <- result$gs_long %>%
    filter(rank == tax_level, call != "Non-microbial / Host", qc_pass) %>%
    group_by(sample, name_clean) %>%
    summarize(val = max(log_rpmm, na.rm = TRUE), .groups = "drop")

  p <- if (nrow(df) == 0 || n_distinct(df$sample) < 2) {
    empty_plot("No QC-pass taxa to build PCA.")
  } else {
    wide <- tidyr::pivot_wider(df, names_from = name_clean, values_from = val, values_fill = 0)
    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$sample

    if (nrow(m) < 2 || ncol(m) < 1) {
      empty_plot("No QC-pass taxa to build PCA.")
    } else {
      keep_cols <- apply(m, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
      if (!any(keep_cols)) {
        empty_plot("No variable taxa to build PCA.")
      } else {
        m <- m[, keep_cols, drop = FALSE]
        pca <- stats::prcomp(m, center = TRUE, scale. = TRUE)
        pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
        pcs$sample <- rownames(pcs)
        pcs <- left_join(pcs, result$meta, by = "sample")

        ggplot(pcs, aes(x = PC1, y = PC2, label = sample, color = group)) +
          geom_point(size = 3, alpha = 0.9) +
          ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
          theme_minimal() +
          labs(
            title = paste0("Cohort PCA (", ifelse(tax_level == "G", "Genus", "Species"), " level)"),
            color = "Group"
          )
      }
    }
  }

  ggsave(file_path, p, width = 8, height = 6, dpi = 300)
}

save_prevalence_plot <- function(result, file_path, title_suffix = "") {
  tax_level <- result$params$tax_level %||% "S"
  tf <- make_prevalence_abundance_df(
    result$gs_long,
    tax_level = tax_level,
    total_samples = n_distinct(result$meta$sample)
  )

  p <- if (nrow(tf) == 0) {
    empty_plot("No cohort taxa available.")
  } else {
    top_labels <- tf %>%
      filter(is.finite(prevalence), is.finite(mean_log_rpmm_detected)) %>%
      arrange(desc(mean_rpmm_detected), desc(prevalence)) %>%
      head(20)

    ggplot(tf, aes(x = prevalence, y = mean_log_rpmm_detected, color = call_bucket)) +
      geom_point(
        position = position_jitter(width = 0.015, height = 0.02),
        size = 2.2,
        alpha = 0.8
      ) +
      ggrepel::geom_text_repel(
        data = top_labels,
        aes(label = name_clean),
        size = 3,
        max.overlaps = 40
      ) +
      scale_color_manual(values = CALL_PALETTE, breaks = CALL_BUCKET_LEVELS) +
      theme_minimal() +
      labs(
        title = paste0("Prevalence vs abundance", title_suffix),
        x = "Prevalence across samples",
        y = "Mean log(RPMM) where present",
        color = NULL
      )
  }

  ggsave(file_path, p, width = 8, height = 6, dpi = 300)
}

save_heatmap_plot <- function(result, file_path, title_prefix = "") {
  tax_level <- result$params$tax_level %||% "S"

  gs <- result$gs_long %>%
    filter(rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm))

  focus <- gs %>% filter(call %in% c("Likely true"))
  if (nrow(focus) == 0) {
    focus <- gs %>% filter(call != "Likely false positive / background")
  }

  p <- if (nrow(focus) == 0) {
    empty_plot("No taxa available for heatmap.")
  } else {
    top_taxa <- focus %>%
      group_by(name_clean) %>%
      summarize(mean_rpmm = mean(rpmm, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_rpmm)) %>%
      head(30)

    heat <- gs %>%
      filter(name_clean %in% top_taxa$name_clean) %>%
      select(sample, name_clean, rpmm) %>%
      mutate(log_rpmm_tile = log10(rpmm + 1))

    heat <- left_join(heat, result$meta %>% select(sample, group), by = "sample")

    samp_order <- result$meta %>%
      arrange(group, sample) %>%
      pull(sample)

    heat$sample <- factor(heat$sample, levels = samp_order)
    heat$name_clean <- factor(heat$name_clean, levels = rev(top_taxa$name_clean))

    ggplot(heat, aes(x = sample, y = name_clean, fill = log_rpmm_tile)) +
      geom_tile() +
      scale_fill_gradient(low = "grey95", high = "navy") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste0(title_prefix, "Heatmap of top taxa (", ifelse(tax_level == "G", "Genus", "Species"), ")"),
        x = NULL,
        y = NULL,
        fill = "log10(RPMM+1)"
      )
  }

  ggsave(file_path, p, width = 10, height = 8, dpi = 300)
}

save_fp_by_sample_plot <- function(result, file_path, title_prefix = "") {
  tax_level <- result$params$tax_level %||% "S"

  df <- result$gs_long %>%
    filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score), qc_pass) %>%
    group_by(sample) %>%
    summarize(
      median_fp = median(fp_score, na.rm = TRUE),
      mean_fp = mean(fp_score, na.rm = TRUE),
      n_taxa = dplyr::n(),
      .groups = "drop"
    ) %>%
    left_join(result$meta, by = "sample")

  p <- if (nrow(df) == 0) {
    empty_plot("No Celltaminate scores available.")
  } else {
    samp_order <- result$meta %>%
      arrange(group, sample) %>%
      pull(sample)
    df$sample <- factor(df$sample, levels = samp_order)

    ggplot(df, aes(x = sample, y = median_fp, color = group)) +
      geom_point(size = 3, alpha = 0.9) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = paste0(title_prefix, "Median Celltaminate score per sample (QC-pass taxa)"),
        x = NULL,
        y = "Median Celltaminate score (0-100)",
        color = "Group"
      )
  }

  ggsave(file_path, p, width = 9, height = 5.5, dpi = 300)
}

save_fp_by_call_plot <- function(result, file_path, title_suffix = "") {
  tax_level <- result$params$tax_level %||% "S"

  df <- result$gs_long %>%
    filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score)) %>%
    mutate(call_bucket = call_bucket_label(call))

  p <- if (nrow(df) == 0) {
    empty_plot("No Celltaminate scores available.")
  } else {
    ggplot(df, aes(x = call_bucket, y = fp_score)) +
      geom_boxplot(outlier.size = 0.6) +
      theme_minimal() +
      labs(
        title = paste0("Celltaminate score distribution by call category", title_suffix),
        x = NULL,
        y = "Celltaminate score (0-100)"
      )
  }

  ggsave(file_path, p, width = 8, height = 5.5, dpi = 300)
}

save_top_taxa_plot <- function(result, sample_name, file_path) {
  tax_level <- result$params$tax_level %||% "S"
  sel_calls <- c("Likely true", "Uncertain")

  gs <- result$gs_long %>%
    filter(sample == sample_name, rank == tax_level, call %in% sel_calls, is.finite(rpmm))

  p <- if (nrow(gs) == 0) {
    empty_plot("No taxa in selected categories.")
  } else {
    gs <- gs %>%
      mutate(
        call_bucket = call_bucket_label(call),
        y = log(rpmm + 1)
      )

    topn <- gs %>%
      mutate(
        call_priority = dplyr::case_when(
          call == "Likely true" ~ 1L,
          call == "Uncertain" ~ 2L,
          call == "Likely false positive / background" ~ 3L,
          TRUE ~ 4L
        )
      ) %>%
      arrange(call_priority, desc(y)) %>%
      head(12) %>%
      mutate(name_clean = factor(name_clean, levels = rev(name_clean)))

    ggplot(topn, aes(x = name_clean, y = y, fill = call_bucket)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      scale_fill_manual(values = CALL_PALETTE, breaks = CALL_BUCKET_LEVELS) +
      labs(x = NULL, y = "log(RPMM+1)", fill = NULL, title = paste0("Top taxa: ", sample_name))
  }

  ggsave(file_path, p, width = 7.5, height = 5.5, dpi = 300)
}

save_fp_scatter_plot <- function(result, sample_name, file_path) {
  tax_level <- result$params$tax_level %||% "S"

  gs <- result$gs_long %>%
    filter(sample == sample_name, rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm), is.finite(fp_score)) %>%
    mutate(
      call_bucket = call_bucket_label(call),
      log10_rpmm = log10(rpmm + 1)
    )

  p <- if (nrow(gs) == 0) {
    empty_plot("No taxa available.")
  } else {
    ggplot(gs, aes(x = log10_rpmm, y = fp_score, color = call_bucket)) +
      geom_point(alpha = 0.7) +
      theme_minimal() +
      scale_color_manual(values = CALL_PALETTE, breaks = CALL_BUCKET_LEVELS) +
      labs(
        title = paste0("Celltaminate score vs abundance: ", sample_name),
        x = "log10(RPMM+1)",
        y = "Celltaminate score (0-100)",
        color = NULL
      )
  }

  ggsave(file_path, p, width = 7.5, height = 5.5, dpi = 300)
}

save_radar_plot <- function(result, sample_name, file_path, n_show = 5) {
  tax_level <- result$params$tax_level %||% "S"
  sel_calls <- c("Likely true", "Uncertain")
  n_show <- max(3, min(10, as.integer(n_show)))

  png(file_path, width = 2200, height = 1800, res = 250)
  on.exit(dev.off(), add = TRUE)

  gs_samp <- result$gs_long %>%
    filter(sample == sample_name, rank == tax_level, call %in% sel_calls, is.finite(log_rpmm))

  if (nrow(gs_samp) == 0) {
    plot.new()
    text(0.5, 0.5, "No taxa in selected categories.")
    return(invisible(NULL))
  }

  top <- gs_samp %>%
    arrange(desc(rpmm)) %>%
    head(n_show)

  taxa <- top$name_clean

  if (length(taxa) < 3) {
    plot.new()
    text(0.5, 0.5, "Not enough taxa for spider plot.")
    return(invisible(NULL))
  }

  tf <- result$tax_features %>%
    filter(rank == tax_level, name_clean %in% taxa)

  ref_med <- tf$ref_q50_log[match(taxa, tf$name_clean)]
  ref_med[!is.finite(ref_med)] <- 0

  samp_vals <- top$log_rpmm
  samp_vals[!is.finite(samp_vals)] <- 0

  all_vals <- c(ref_med, samp_vals)
  mx <- max(all_vals, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) mx <- 1

  ref_scaled <- ref_med / mx
  samp_scaled <- samp_vals / mx

  df_radar <- as.data.frame(rbind(
    rep(1, length(taxa)),
    rep(0, length(taxa)),
    ref_scaled,
    samp_scaled
  ))
  colnames(df_radar) <- taxa
  rownames(df_radar) <- c("max", "min", "Reference median", "User sample")

  layout_cfg <- spider_plot_layout(taxa)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mar = layout_cfg$mar, xpd = NA)
  fmsb::radarchart(
    df_radar,
    axistype = 1,
    pcol = c("black", "red"),
    plty = c(1, 1),
    plwd = c(2, 2),
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    vlcex = 0,
    title = paste0("Reference median vs sample (scaled log-RPMM): ", sample_name)
  )

  draw_spider_labels(
    taxa = taxa,
    radius = layout_cfg$label_radius,
    cex = layout_cfg$label_cex
  )

  legend(
    x = "topright",
    legend = c("Reference median", "User sample"),
    col = c("black", "red"),
    lty = 1,
    bty = "n",
    cex = 0.9
  )
}

save_all_plots <- function(result, out_dir, prefix = NULL) {
  plots_dir <- file.path(out_dir, "plots")
  per_sample_dir <- file.path(plots_dir, "per_sample")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)

  pref <- if (is.null(prefix) || !nzchar(prefix)) "" else paste0(prefix, "_")
  title_prefix <- if (is.null(prefix) || !nzchar(prefix)) "" else paste0(prefix, " ")
  title_suffix <- if (is.null(prefix) || !nzchar(prefix)) "" else paste0(" (", prefix, ")")

  save_cohort_pca_plot(result, file.path(plots_dir, paste0(pref, "cohort_pca.png")))
  save_prevalence_plot(result, file.path(plots_dir, paste0(pref, "cohort_prevalence_abundance.png")), title_suffix = title_suffix)
  save_heatmap_plot(result, file.path(plots_dir, paste0(pref, "cohort_heatmap.png")), title_prefix = title_prefix)
  save_fp_by_sample_plot(result, file.path(plots_dir, paste0(pref, "cohort_fp_by_sample.png")), title_prefix = title_prefix)
  save_fp_by_call_plot(result, file.path(plots_dir, paste0(pref, "cohort_fp_by_call.png")), title_suffix = title_suffix)

  for (samp in unique(result$gs_long$sample)) {
    sid <- safe_sample_id(samp)
    save_top_taxa_plot(result, samp, file.path(per_sample_dir, paste0(pref, sid, "_top_taxa.png")))
    save_fp_scatter_plot(result, samp, file.path(per_sample_dir, paste0(pref, sid, "_fp_scatter.png")))
    save_radar_plot(result, samp, file.path(per_sample_dir, paste0(pref, sid, "_radar.png")))
  }
}

write_cleaned_reports <- function(result, out_dir, quick = c("false", "unc"),
                                  protected_taxa_map = list(), keep_descendants = FALSE) {
  cleaned_dir <- file.path(out_dir, "cleaned_reports")
  dir.create(cleaned_dir, recursive = TRUE, showWarnings = FALSE)

  out_files <- list()

  for (samp in names(result$sample_list)) {
    sid <- safe_sample_id(samp)
    gs_samp <- result$gs_long %>% filter(sample == samp)
    protected_taxa <- protected_taxa_map[[samp]] %||% character(0)
    raw_df <- result$sample_list[[samp]]$raw

    cleaned <- build_cleaned_report_for_sample(
      raw_df = raw_df,
      gs = gs_samp,
      quick = quick,
      protected_taxa = protected_taxa,
      keep_descendants = keep_descendants
    )

    out_path <- file.path(cleaned_dir, paste0("cleaned_", sid, ".txt"))
    write.table(cleaned, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    out_files[[samp]] <- out_path
  }

  out_files
}

save_run_summary <- function(result, out_dir, args, label = "original") {
  path <- file.path(out_dir, paste0(label, "_run_summary.txt"))
  lines <- c(
    paste0("label\t", label),
    paste0("date\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("tax_level\t", result$params$tax_level),
    paste0("n_samples\t", n_distinct(result$meta$sample)),
    paste0("n_rows_gs_long\t", nrow(result$gs_long)),
    paste0("parse_errors\t", if (length(result$parse_errors) > 0) paste(result$parse_errors, collapse = " | ") else "none"),
    "",
    "[arguments]"
  )

  if (length(args) > 0) {
    arg_lines <- vapply(names(args), function(k) paste0(k, "\t", args[[k]]), character(1))
    lines <- c(lines, arg_lines)
  }

  writeLines(lines, con = path)
}

# ------------------- MAIN -------------------
main <- function() {
  args <- parse_cli_args()

  input_files <- args$input_files %||% ""
  input_list <- args$input_list %||% ""
  output_dir <- require_arg(args, "output_dir")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  init_reference_tables(args)

  file_map <- collect_input_reports(input_files = input_files, input_list = input_list)

  host_species <- split_csv(args$host_species %||% "human")
  if (length(host_species) == 0) host_species <- c("human")

  host_taxids <- build_host_taxids(
    host_species = host_species,
    host_taxids_manual = args$host_taxids_manual %||% ""
  )
  host_name_patterns <- build_host_name_patterns(host_species = host_species)

  meta <- load_metadata(args$metadata_tsv %||% "", names(file_map))

  result <- run_pipeline(
    file_map = file_map,
    meta = meta,
    args = args,
    host_taxids = host_taxids,
    host_name_patterns = host_name_patterns
  )

  save_tables(
    result = result,
    out_dir = output_dir,
    show_fp_breakdown = parse_bool(args$show_fp_breakdown, FALSE)
  )
  save_all_plots(result = result, out_dir = output_dir, prefix = "original")
  save_run_summary(result = result, out_dir = output_dir, args = args, label = "original")

  protected_taxa_map <- load_protected_taxa_map(args$protected_taxa_tsv %||% "")
  quick_remove <- split_csv(args$quick_remove %||% "false,unc")
  keep_descendants <- parse_bool(args$keep_descendants, FALSE)

  cleaned_files <- write_cleaned_reports(
    result = result,
    out_dir = output_dir,
    quick = quick_remove,
    protected_taxa_map = protected_taxa_map,
    keep_descendants = keep_descendants
  )

  if (parse_bool(args$save_cleaned_reanalysis, TRUE) && length(cleaned_files) > 0) {
    cleaned_result <- run_pipeline(
      file_map = cleaned_files,
      meta = result$meta,
      args = args,
      host_taxids = host_taxids,
      host_name_patterns = host_name_patterns
    )

    cleaned_out_dir <- file.path(output_dir, "cleaned_reanalysis")
    dir.create(cleaned_out_dir, recursive = TRUE, showWarnings = FALSE)

    save_tables(
      result = cleaned_result,
      out_dir = cleaned_out_dir,
      show_fp_breakdown = parse_bool(args$show_fp_breakdown, FALSE)
    )
    save_all_plots(result = cleaned_result, out_dir = cleaned_out_dir, prefix = "cleaned")
    save_run_summary(result = cleaned_result, out_dir = cleaned_out_dir, args = args, label = "cleaned")
  }

  message("Done. Results written to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
