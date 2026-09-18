
# ------------------- LIBRARIES -------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(dplyr)
  library(stats)
  library(ggplot2)
  library(ggrepel)
  library(fmsb)
  library(writexl)
  library(DT)
  library(tibble)
})

HAS_HTTR <- requireNamespace("httr", quietly = TRUE)
HAS_JSONLITE <- requireNamespace("jsonlite", quietly = TRUE)
HAS_ZIP <- requireNamespace("zip", quietly = TRUE)
HAS_XML2 <- requireNamespace("xml2", quietly = TRUE)

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
    call == "Prioritized" ~ "Prioritized",
    TRUE ~ "Not prioritized"
  )
}

CALL_BUCKET_LEVELS <- c(
  "Prioritized",
  "Not prioritized"
)

CALL_PALETTE <- c(
  "Prioritized" = "#2E7D32",
  "Not prioritized" = "#8B0000"
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
      call_bucket = dplyr::coalesce(call_bucket, "Not prioritized"),
      prevalence = cap(prevalence, 0, 1)
    )
}

get_bioai_taxa_info <- function(gs, tax_level = "S") {
  if (is.null(gs) || nrow(gs) == 0) {
    return(list(
      candidates = character(0),
      prioritized_taxa = character(0),
      initial_selected = character(0),
      initial_choices = character(0)
    ))
  }

  base <- gs %>%
    filter(rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm), nzchar(name_clean)) %>%
    mutate(
      call_priority = dplyr::case_when(
        call == "Prioritized" ~ 1,
        TRUE ~ 2
      )
    )

  candidates <- base %>%
    arrange(call_priority, desc(rpmm)) %>%
    pull(name_clean) %>%
    unique()

  prioritized_taxa <- base %>%
    filter(call == "Prioritized") %>%
    arrange(desc(rpmm)) %>%
    pull(name_clean) %>%
    unique()

  initial_selected <- if (length(prioritized_taxa) > 0) prioritized_taxa else head(candidates, 3)
  initial_choices <- unique(c(initial_selected, head(candidates, 100)))

  list(
    candidates = as.character(candidates),
    prioritized_taxa = as.character(prioritized_taxa),
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
  lines <- tryCatch(
    readLines(file_path, warn = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(lines) || length(lines) == 0) {
    return(list(error = paste0("Failed to read file: ", basename(file_path))))
  }
  
  lines <- gsub("\r$", "", lines)
  lines <- lines[nzchar(lines)]
  
  if (length(lines) == 0) {
    return(list(error = paste0("Failed to read file: ", basename(file_path))))
  }
  
  parts <- strsplit(lines, "\t", fixed = TRUE)
  n_fields <- vapply(parts, length, integer(1))
  max_fields <- max(n_fields)
  
  fmt <- if (max_fields >= 8) {
    "krakenuniq_8col"
  } else if (max_fields == 6) {
    "kraken2_6col"
  } else {
    NA_character_
  }
  
  if (is.na(fmt)) {
    return(list(error = paste0(
      "Invalid Kraken report format for: ", basename(file_path),
      ". Expected 6-column Kraken2 report or 8-column Kraken2/KrakenUniq report."
    )))
  }
  
  use_n <- if (fmt == "krakenuniq_8col") 8 else 6
  
  raw_mat <- do.call(
    rbind,
    lapply(parts, function(x) {
      if (length(x) < use_n) {
        x <- c(x, rep("", use_n - length(x)))
      }
      x[seq_len(use_n)]
    })
  )
  
  raw <- as.data.frame(raw_mat, stringsAsFactors = FALSE, check.names = FALSE)
  raw[] <- lapply(raw, as.character)
  
  if (fmt == "krakenuniq_8col") {
    colnames(raw) <- c("pct", "reads_clade", "reads_direct", "uniq_kmers", "dup_kmers", "rank", "taxid", "name_raw")
    
    df <- raw %>%
      mutate(
        pct = suppressWarnings(as.numeric(pct)),
        reads_clade = suppressWarnings(as.numeric(reads_clade)),
        reads_direct = suppressWarnings(as.numeric(reads_direct)),
        uniq_kmers = suppressWarnings(as.numeric(uniq_kmers)),
        dup_kmers = suppressWarnings(as.numeric(dup_kmers)),
        rank = trimws(as.character(rank)),
        taxid = suppressWarnings(as.integer(taxid)),
        name_raw = as.character(name_raw)
      )
  } else {
    colnames(raw) <- c("pct", "reads_clade", "reads_direct", "rank", "taxid", "name_raw")
    
    df <- raw %>%
      mutate(
        pct = suppressWarnings(as.numeric(pct)),
        reads_clade = suppressWarnings(as.numeric(reads_clade)),
        reads_direct = suppressWarnings(as.numeric(reads_direct)),
        uniq_kmers = NA_real_,
        dup_kmers = NA_real_,
        rank = trimws(as.character(rank)),
        taxid = suppressWarnings(as.integer(taxid)),
        name_raw = as.character(name_raw)
      )
  }
  
  name_left_trim <- trimws(df$name_raw, which = "left")
  lead_spaces <- nchar(df$name_raw) - nchar(name_left_trim)
  lead_spaces[!is.finite(lead_spaces)] <- 0
  
  df <- df %>%
    mutate(
      sample = sample_name %||% basename(file_path),
      name_clean = trimws(name_raw)
    )
  
  df$depth <- floor(pmax(lead_spaces, 0) / 2)
  df$kmer_per_read <- ifelse(
    is.na(df$uniq_kmers) | df$reads_clade <= 0,
    NA_real_,
    df$uniq_kmers / df$reads_clade
  )
  
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

  use_cohort_prevalence <- is.finite(n_samples) && n_samples >= 3

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
    mutate(bg_rpmm = pmax(coalesce0(ref_med_rpmm), coalesce0(cohort_bg_rpmm)))

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

taxa_table_columns <- function(tax_level = "S", show_fp_breakdown = FALSE) {
  cols <- c(
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
    cols <- c(
      cols,
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
  
  cols
}

prepare_taxa_table_data <- function(gs, tax_level = "S", show_fp_breakdown = FALSE, show_kitome = TRUE, show_clinical_panel = TRUE) {
  if (!"prevalence_any" %in% names(gs)) {
    gs$prevalence_any <- gs$prevalence
  }

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
      In_kitome = if (isTRUE(show_kitome)) in_kitome else NA,
      In_clinical_panel = if (isTRUE(show_clinical_panel)) in_clinical_panel else NA,
      Call = as.character(call),
      Call_reason = ifelse(
        tax_level == "S" & !is.na(Top_strain) & nzchar(Top_strain),
        paste0(call_reason, " Display note: the top strain shown here inherits the species-level Celltaminate score."),
        call_reason
      )
    ) %>%
    select(any_of(taxa_table_columns(tax_level = tax_level, show_fp_breakdown = show_fp_breakdown))) %>%
    mutate(across(everything(), ~ if (is.factor(.x)) as.character(.x) else .x))

  gs_view <- as.data.frame(gs_view, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(gs_view) <- NULL
  gs_view
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
      cohort_bg_rpmm = suppressWarnings(as.numeric(cohort_bg_rpmm)),
      bg_rpmm = suppressWarnings(as.numeric(bg_rpmm))
    )

  out <- out %>%
    mutate(
      reference_median_rpmm = if_else(is.finite(ref_med_rpmm), ref_med_rpmm, NA_real_),
      log2FC_vs_reference_median = log2((coalesce0(rpmm) + params$fc_pseudocount) / (coalesce0(reference_median_rpmm) + params$fc_pseudocount)),
      qc_pass = is_microbial & !is_host & !is_plant
    )

  amb <- compute_ambiguity_index(out)
  
  out <- out %>%
    left_join(
      amb,
      by = c(
        "sample",
        "rank",
        "name_clean"
      )
    )
  
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
      
      fp_true_flag =
        fp_score <=
        params$fp_true_cutoff
    )

  out <- out %>%
    mutate(
      call = case_when(
        !is_microbial | is_host | is_plant ~ "Non-microbial / Host",
        fp_true_flag ~ "Prioritized",
        TRUE ~ "Not prioritized"
      ),
      call_reason = case_when(
        !is_microbial | is_host | is_plant ~ "Excluded as non-microbial or host",
        call == "Prioritized" ~ paste0(
          "Score ", round(fp_score, 1),
          " <= prioritization cutoff ", round(params$fp_true_cutoff, 1)
        ),
        TRUE ~ paste0(
          "Score ", round(fp_score, 1),
          " > prioritization cutoff ", round(params$fp_true_cutoff, 1)
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

compute_remove_taxa_for_sample <- function(gs, protected_taxa = character(0)) {
  remove_taxa <- unique(gs$name_clean[gs$call == "Not prioritized"])
  remove_taxa <- remove_taxa[!is.na(remove_taxa) & nzchar(remove_taxa)]
  setdiff(remove_taxa, protected_taxa)
}

compute_keep_taxa_for_sample <- function(gs, remove_taxa, protected_taxa = character(0)) {
  keep_taxa <- gs %>%
    filter(rank %in% c("G", "S"), call == "Prioritized") %>%
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

build_cleaned_report_for_sample <- function(raw_df, gs, protected_taxa = character(0), keep_descendants = FALSE) {
  remove_taxa <- compute_remove_taxa_for_sample(gs, protected_taxa = protected_taxa)
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

# ------------------- BIOAI -----------------------------
bioai_lite_explain <- function(sample_name, sample_type, top_df) {
  if (is.null(top_df) || nrow(top_df) == 0) {
    return(HTML("<b>BioAI:</b> No high-confidence microbial taxa detected after Celltaminate filters."))
  }

  top_df <- top_df %>% head(8)

  bullets <- paste0(
    "<li><b>", htmltools::htmlEscape(top_df$name_clean), "</b> (", top_df$rank, ") — ",
    "reads: ", format(top_df$reads_clade, scientific = FALSE, big.mark = ","), "; ",
    "rpmm: ", format(round(top_df$rpmm, 3), nsmall = 3), "; ",
    "call: <i>", htmltools::htmlEscape(top_df$call), "</i>",
    ifelse(!is.na(top_df$log2FC_vs_reference_median), paste0("; log2FC vs reference median: ", round(top_df$log2FC_vs_reference_median, 2)), ""),
    "</li>"
  )

  HTML(paste0(
    "<b>BioAI:</b><br/>",
    "Specimen context: <b>", htmltools::htmlEscape(sample_type %||% "Unknown"), "</b>.<br/>",
    "Prioritized taxa (after Celltaminate scoring):",
    "<ul>", paste(bullets, collapse = ""), "</ul>",
    "<b>Interpretation tips:</b><br/>",
    "• Any microbial signal should be interpreted alongside other orthogonal validations.<br/>",
    "• The Celltaminate score integrates quantitative sequence evidence, reference background, decontamination, and curated contextual evidence to prioritize microbial taxa.<br/>",
    "• Consider orthogonal confirmation before clinical decisions."
  ))
}

text_is_present <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  y <- suppressWarnings(as.character(x)[1])
  if (length(y) == 0 || is.na(y)) return(FALSE)
  nzchar(trimws(y))
}

lookup_kitome_blacklist_entry <- function(names_to_check, blacklist_df = KITOME_BG_BLACKLIST) {
  if (is.null(blacklist_df) || nrow(blacklist_df) == 0) return(NULL)

  names_to_check <- as.character(names_to_check)
  names_to_check <- trimws(names_to_check)
  names_to_check <- names_to_check[!is.na(names_to_check) & nzchar(names_to_check)]
  if (length(names_to_check) == 0) return(NULL)

  hit <- blacklist_df %>%
    filter(organism_name_lc %in% tolower(names_to_check)) %>%
    slice(1)

  if (nrow(hit) == 0) return(NULL)
  hit
}

lookup_clinically_important_pathogen_entry <- function(names_to_check, pathogen_df = CLINICALLY_IMPORTANT_PATHOGENS) {
  if (is.null(pathogen_df) || nrow(pathogen_df) == 0) return(NULL)

  names_to_check <- as.character(names_to_check)
  names_to_check <- trimws(names_to_check)
  names_to_check <- names_to_check[!is.na(names_to_check) & nzchar(names_to_check)]
  if (length(names_to_check) == 0) return(NULL)

  hit <- pathogen_df %>%
    filter(name_lc %in% tolower(names_to_check)) %>%
    slice(1)

  if (nrow(hit) == 0) return(NULL)
  hit
}

make_celltaminate_remark_html <- function(names_to_check, blacklist_df = KITOME_BG_BLACKLIST, pathogen_df = CLINICALLY_IMPORTANT_PATHOGENS) {
  out <- character(0)

  kit_hit <- lookup_kitome_blacklist_entry(names_to_check, blacklist_df = blacklist_df)
  if (!is.null(kit_hit) && nrow(kit_hit) > 0) {
    note <- as.character(kit_hit$note[1])
    doi <- as.character(kit_hit$doi_hyperlink[1])

    if (!text_is_present(note)) note <- "Listed in kitome/background blacklist"

    note_html <- htmltools::htmlEscape(note)

    if (text_is_present(doi)) {
      doi_html <- htmltools::htmlEscape(doi)
      note_html <- paste0(
        "<a href='", doi_html, "' target='_blank' rel='noopener noreferrer'>",
        note_html,
        "</a>"
      )
    }

    out <- c(out, paste0("<b>Celltaminate Remark:</b> ", note_html, "<br/>"))
  }

  clin_hit <- lookup_clinically_important_pathogen_entry(names_to_check, pathogen_df = pathogen_df)
  if (!is.null(clin_hit) && nrow(clin_hit) > 0) {
    clinical_category <- as.character(clin_hit$clinical_category[1])
    doi <- as.character(clin_hit$doi[1])

    if (!text_is_present(clinical_category)) clinical_category <- "clinically_important_pathogen"

    clinical_html <- htmltools::htmlEscape(clinical_category)

    if (text_is_present(doi)) {
      doi_url <- if (grepl("^https?://", doi, ignore.case = TRUE)) doi else paste0("https://doi.org/", doi)
      doi_html <- htmltools::htmlEscape(doi_url)
      clinical_html <- paste0(
        "<a href='", doi_html, "' target='_blank' rel='noopener noreferrer'>",
        clinical_html,
        "</a>"
      )
    }

    out <- c(out, paste0("<b>Celltaminate Remark:</b> ", clinical_html, "<br/>"))
  }

  paste0(out, collapse = "")
}

ncbi_taxonomy_lookup <- function(query, timeout_sec = 10) {
  if (!HAS_HTTR || !HAS_XML2) return(NULL)

  esearch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  r1 <- tryCatch(
    httr::GET(esearch_url, query = list(db = "taxonomy", term = query, retmode = "xml"), httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  if (is.null(r1) || httr::status_code(r1) != 200) return(NULL)

  x1 <- tryCatch(xml2::read_xml(httr::content(r1, as = "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(x1)) return(NULL)

  ids <- xml2::xml_text(xml2::xml_find_all(x1, ".//Id"))
  if (length(ids) == 0) return(NULL)

  taxid <- ids[1]

  efetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  r2 <- tryCatch(
    httr::GET(efetch_url, query = list(db = "taxonomy", id = taxid, retmode = "xml"), httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  if (is.null(r2) || httr::status_code(r2) != 200) return(NULL)

  x2 <- tryCatch(xml2::read_xml(httr::content(r2, as = "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(x2)) return(NULL)

  sci <- xml2::xml_text(xml2::xml_find_first(x2, ".//ScientificName"))
  rank <- xml2::xml_text(xml2::xml_find_first(x2, ".//Rank"))
  lin_nodes <- xml2::xml_find_all(x2, ".//LineageEx/Taxon/ScientificName")
  lineage <- xml2::xml_text(lin_nodes)
  lineage <- lineage[!is.na(lineage) & nzchar(lineage)]

  list(taxid = taxid, scientific_name = sci, rank = rank, lineage = lineage)
}

wikipedia_summary <- function(title, timeout_sec = 10) {
  if (!HAS_HTTR || !HAS_JSONLITE) return(NULL)

  base <- "https://en.wikipedia.org/api/rest_v1/page/summary/"
  url <- paste0(base, URLencode(title, reserved = TRUE))
  r <- tryCatch(httr::GET(url, httr::timeout(timeout_sec)), error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200) return(NULL)

  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  js <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (is.null(js)) return(NULL)

  extract <- js$extract %||% ""
  extract <- as.character(extract)[1]
  if (!text_is_present(extract)) return(NULL)
  extract
}

wikidata_summary <- function(query, timeout_sec = 10) {
  if (!HAS_HTTR || !HAS_JSONLITE) return(NULL)

  s1 <- tryCatch(
    httr::GET(
      "https://www.wikidata.org/w/api.php",
      query = list(action = "wbsearchentities", search = query, language = "en", format = "json", limit = 1),
      httr::timeout(timeout_sec)
    ),
    error = function(e) NULL
  )
  if (is.null(s1) || httr::status_code(s1) != 200) return(NULL)

  js1 <- tryCatch(jsonlite::fromJSON(httr::content(s1, as = "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(js1) || is.null(js1$search) || length(js1$search) == 0 || is.null(js1$search$id) || length(js1$search$id) == 0) return(NULL)

  qid <- suppressWarnings(as.character(js1$search$id[1]))
  if (!text_is_present(qid)) return(NULL)

  s2 <- tryCatch(
    httr::GET(
      "https://www.wikidata.org/w/api.php",
      query = list(action = "wbgetentities", ids = qid, props = "descriptions", languages = "en", format = "json"),
      httr::timeout(timeout_sec)
    ),
    error = function(e) NULL
  )
  if (is.null(s2) || httr::status_code(s2) != 200) return(NULL)

  js2 <- tryCatch(jsonlite::fromJSON(httr::content(s2, as = "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(js2) || is.null(js2$entities) || is.null(js2$entities[[qid]]) || is.null(js2$entities[[qid]]$descriptions) || is.null(js2$entities[[qid]]$descriptions$en) || is.null(js2$entities[[qid]]$descriptions$en$value)) return(NULL)

  out <- js2$entities[[qid]]$descriptions$en$value
  out <- suppressWarnings(as.character(out)[1])
  if (!text_is_present(out)) return(NULL)
  out
}

bioai_web_explain <- function(sample_name, sample_type, top_df) {
  if (is.null(top_df) || nrow(top_df) == 0) {
    return(HTML("<b>BioAI-web:</b> No taxa selected."))
  }

  out_sections <- c()

  for (i in seq_len(nrow(top_df))) {
    nm <- as.character(top_df$name_clean[i])
    rk <- as.character(top_df$rank[i])
    reads <- top_df$reads_clade[i]
    rpmm <- top_df$rpmm[i]
    fp_score <- suppressWarnings(as.numeric(top_df$fp_score[i]))
    call <- as.character(top_df$call[i])
    call_reason <- as.character(top_df$call_reason[i])
    if (is.na(call_reason) || !nzchar(call_reason)) call_reason <- ""

    tx <- ncbi_taxonomy_lookup(nm)
    sci <- if (!is.null(tx) && text_is_present(tx$scientific_name)) tx$scientific_name[1] else nm
    lineage_txt <- if (!is.null(tx) && length(tx$lineage) > 0) paste(tx$lineage, collapse = " > ") else ""

    wk <- wikipedia_summary(sci)
    if (!text_is_present(wk)) wk <- wikipedia_summary(nm)

    wd <- NULL
    if (!text_is_present(wk)) wd <- wikidata_summary(sci)
    if (!text_is_present(wd)) wd <- wikidata_summary(nm)

    info_txt <- wk
    if (!text_is_present(info_txt)) info_txt <- wd
    if (!text_is_present(info_txt)) info_txt <- "(No Wikipedia or Wikidata summary found.)"

    section <- paste0(
      "<h4>", htmltools::htmlEscape(sci), " <span style='font-size:12px;color:#666;'>(", htmltools::htmlEscape(rk), ")</span></h4>",
      "<div style='color:#444;'>",
      "<b>Sample:</b> ", htmltools::htmlEscape(sample_name), " &nbsp; | &nbsp; ",
      "<b>Specimen context:</b> ", htmltools::htmlEscape(sample_type %||% "Unknown"), "<br/>",
      "<b>Evidence:</b> reads=", format(reads, scientific = FALSE, big.mark = ","), ", rpmm=", round(as.numeric(rpmm), 3), "<br/>",
      "<b>Celltaminate score:</b> ", ifelse(is.finite(fp_score), format(round(fp_score, 1), nsmall = 1), "NA"), "<br/>",
      "<b>Celltaminate call:</b> ", htmltools::htmlEscape(call), "<br/>",
      if (text_is_present(call_reason)) paste0("<b>Call reason:</b> ", htmltools::htmlEscape(call_reason), "<br/>") else "",
      make_celltaminate_remark_html(c(nm, sci)),
      if (text_is_present(lineage_txt)) paste0("<b>NCBI lineage:</b> ", htmltools::htmlEscape(lineage_txt), "<br/>") else "",
      "</div>",
      "<div style='margin-top:6px; white-space:pre-wrap;'>", htmltools::htmlEscape(info_txt), "</div>",
      "<hr/>"
    )

    out_sections <- c(out_sections, section)
  }

  HTML(paste0(
    "<b>BioAI :</b><br/>",
    "<div class='small-note'>This tool uses summaries from NCBI Taxonomy and Wikipedia/WikiData. It is informational and should not be solely used for clinical decision.</div><br/>",
    paste(out_sections, collapse = "\n")
  ))
}

safe_bioai_web_explain <- function(sample_name, sample_type, top_df) {
  tryCatch(
    bioai_web_explain(sample_name, sample_type, top_df),
    error = function(e) {
      HTML(paste0(
        "<b>BioAI :</b><br/>",
        "<div class='small-note'>BioAI interpretation could not be generated for the selected taxa.</div><br/>",
        "<div style='color:#444;'>",
        htmltools::htmlEscape(conditionMessage(e)),
        "</div>"
      ))
    }
  )
}

can_run_pca <- function(wide_df) {
  if (is.null(wide_df) || nrow(wide_df) < 2 || ncol(wide_df) < 3) return(FALSE)
  m <- as.matrix(wide_df[, -1, drop = FALSE])
  if (nrow(m) < 2 || ncol(m) < 2) return(FALSE)
  any(apply(m, 2, function(x) stats::var(x, na.rm = TRUE) > 0))
}

spider_plot_layout <- function(taxa) {
  taxa <- as.character(taxa)
  taxa <- taxa[!is.na(taxa)]
  max_chars <- if (length(taxa) > 0) max(nchar(taxa), na.rm = TRUE) else 10
  label_cex <- max(0.5, min(0.78, 20 / max(20, max_chars)))
  label_radius <- max(1.16, min(1.34, 1.12 + 0.010 * max_chars))
  list(
    mai = c(0.7, 0.7, 1.0, 0.7),
    label_cex = label_cex,
    label_radius = label_radius
  )
}

# ------------------- UI HELPERS -------------------
analysis_settings_ui <- function(include_analyze_button = TRUE, include_file_inputs = TRUE, defaults = NULL) {
  tagList(
    if (include_file_inputs) {
      tagList(
        tags$hr(),
        h4("Inputs"),
        fileInput(
          "files",
          "Upload Kraken report files (.txt/.tsv or .zip of many)",
          multiple = TRUE,
          accept = c(".txt", ".tsv", ".report", ".zip")
        ),
        if (include_analyze_button) {
          tagList(
            actionButton("analyze_button", "Run analysis", icon = icon("play"), class = "btn-primary run-analysis-button"),
            progressBar(id = "analysis_progress", value = 0, display_pct = TRUE),
            tags$small(class = "text-muted", "The progress bar updates as uploaded reports are parsed and scored.")
          )
        } else {
          NULL
        }
      )
    } else {
      NULL
    },
    
    tags$details(
      class = "advanced-parameters-box",
      tags$summary(tags$strong("Show advanced parameters")),
      
      selectInput(
        "tax_level",
        "Taxonomy level to show in plots",
        choices = c("Genus" = "G", "Species" = "S"),
        selected = ui_default(defaults, "tax_level", "S")
      ),
      tags$small(class = "text-muted", "A single taxonomy level is applied across all tables and plots."),
      br(),
      
      checkboxGroupInput(
        "host_species",
        "Host species to exclude",
        choices = c("human", "mouse", "drosophila"),
        selected = ui_default(defaults, "host_species", c("human"))
      ),
      textInput("host_taxids_manual", "Additional host taxids (comma-separated, optional)", value = ui_default(defaults, "host_taxids_manual", "")),
      checkboxInput("include_protozoa", "Include non-host eukaryotic microbes (e.g. protozoa)", value = ui_default(defaults, "include_protozoa", TRUE)),
      
      tags$hr(),
      h4("Merge Species Parameters"),
      checkboxInput("collapse_species", "Merge noisy species within a genus", value = ui_default(defaults, "collapse_species", TRUE)),
      sliderInput("collapse_min_genus_reads", "Min genus reads for species collapse", min = 0, max = 5000, value = ui_default(defaults, "collapse_min_genus_reads", 30), step = 1),
      sliderInput("collapse_top_frac", "Dominance fraction to keep top species", min = 0.5, max = 1.0, value = ui_default(defaults, "collapse_top_frac", 0.85), step = 0.05),
      
      tags$hr(),
      tags$hr(),
      
      h4(
        "Celltaminate score"
      ),
      sliderInput(
        "fp_true_cutoff",
        "Prioritization cutoff (≤)",
        min = 0,
        max = 100,
        value = ui_default(
          defaults,
          "fp_true_cutoff",
          1
        ),
        step = 1
      ),
      
      tags$details(
        tags$summary(
          tags$strong(
            "Advanced prevalence inputs"
          )
        ),
        
        h5(
          "Prevalence and ubiquity inputs"
        ),
        
        sliderInput(
          "min_reads_prevalence",
          "Min reads for prevalence",
          min = 0,
          max = 500,
          value = ui_default(
            defaults,
            "min_reads_prevalence",
            5
          ),
          step = 1
        ),
        
        sliderInput(
          "min_rpmm_prevalence",
          "Min RPMM for prevalence",
          min = 0,
          max = 100,
          value = ui_default(
            defaults,
            "min_rpmm_prevalence",
            0.2
          ),
          step = 0.1
        ),
        
        sliderInput(
          "min_uniq_kmers_prevalence",
          "Min unique k-mers for prevalence",
          min = 0,
          max = 5000,
          value = ui_default(
            defaults,
            "min_uniq_kmers_prevalence",
            96
          ),
          step = 1
        ),
        
        sliderInput(
          "min_kmer_per_read_prevalence",
          "Min unique k-mers/read for prevalence",
          min = 0,
          max = 50,
          value = ui_default(
            defaults,
            "min_kmer_per_read_prevalence",
            0.4
          ),
          step = 0.1
        )
      ),
      checkboxInput(
        "show_fp_breakdown",
        "Show score breakdown columns in taxa tables",
        value = ui_default(
          defaults,
          "show_fp_breakdown",
          FALSE
        )
      )
    )
  )
}

# ------------------- UI -------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Arial, sans-serif; }
      .well { border-radius: 8px; }
      .celltaminate-title { font-weight: 800; letter-spacing: 0.2px; }
      .celltaminate-subtitle { font-size: 18px; font-weight: 600; color: #555; margin-bottom: 8px; }
      .subtle { color: #444; }
      .call-pill { padding: 2px 8px; border-radius: 12px; font-size: 12px; display: inline-block; }
      .pill-true { background: #e6f4ea; }
      .pill-cont { background: #fdecea; }
      .pill-low { background: #fff4e5; }
      .pill-unc { background: #eef2ff; }
      .pill-host { background: #f3f4f6; }
      .small-note { font-size: 12px; color: #666; }
      .intro-card { background: #f8fafc; border: 1px solid #e5e7eb; border-radius: 10px; padding: 14px; min-height: 130px; margin-bottom: 12px; }
      .intro-card h4 { margin-top: 0; margin-bottom: 8px; font-size: 18px; }
      .intro-card p { margin-bottom: 0; color: #444; }
      .preview-card { background: #ffffff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 12px; margin-bottom: 18px; }
      .preview-title { text-align:center; font-weight:700; margin-top:10px; margin-bottom:6px; }
      .preview-text { text-align:center; color:#555; font-size:14px; line-height:1.4; min-height: 38px; }

      .run-analysis-button {
        width: 100%;
        margin-top: 0px;
        margin-bottom: 10px;
        font-weight: 700;
      }

      .advanced-parameters-box {
        margin-top: 14px;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        background: #f8fafc;
        padding: 10px 12px 12px 12px;
      }

      .advanced-parameters-box > summary {
        cursor: pointer;
        color: #2c3e50;
        margin-bottom: 8px;
      }
    "))
  ),
  uiOutput("main_page")
)

server <- function(input, output, session) {
  app_dir_candidates <- unique(c(
    tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE), error = function(e) NA_character_),
    tryCatch({
      ofiles <- unlist(lapply(sys.frames(), function(x) {
        if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_
      }), use.names = FALSE)
      ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
      if (length(ofiles) > 0) dirname(normalizePath(ofiles[1], winslash = "/", mustWork = FALSE)) else NA_character_
    }, error = function(e) NA_character_),
    tryCatch({
      if (requireNamespace("rstudioapi", quietly = TRUE)) {
        p <- rstudioapi::getSourceEditorContext()$path
        if (!is.null(p) && nzchar(p)) dirname(normalizePath(p, winslash = "/", mustWork = FALSE)) else NA_character_
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)
  ))
  app_dir_candidates <- unique(app_dir_candidates[!is.na(app_dir_candidates) & nzchar(app_dir_candidates)])
  app_dir <- first_existing_path(app_dir_candidates) %||% normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  asset_src <- function(fname) {
    as.character(fname)[1]
  }

  resolve_app_file <- function(fname) {
    first_existing_path(c(
      fname,
      file.path("www", fname),
      file.path("data", "panels", fname),
      file.path("data", "reference", fname),
      unlist(lapply(app_dir_candidates, function(d) c(
        file.path(d, fname),
        file.path(d, "www", fname),
        file.path(d, "data", "panels", fname),
        file.path(d, "data", "reference", fname),
        file.path(d, "..", "data", "panels", fname),
        file.path(d, "..", "data", "reference", fname)
      )), use.names = FALSE),
      file.path(app_dir, fname),
      file.path(app_dir, "www", fname),
      file.path(app_dir, "..", "data", "panels", fname),
      file.path(app_dir, "..", "data", "reference", fname)
    ))
  }

  ensure_reference_assets_loaded <- local({
    loaded <- FALSE
    function(force = FALSE) {
      if (loaded && !isTRUE(force)) return(invisible(TRUE))

      ref_path <- resolve_app_file("refined_cell.lines.tsv")
      if (!is.null(ref_path)) {
        REF_BG <<- load_ref_cell_lines(path = ref_path, eps = 1e-9)
      } else {
        REF_BG <<- load_ref_cell_lines(path = "refined_cell.lines.tsv", eps = 1e-9)
      }

      kitome_blacklist_path <- resolve_app_file("kitome_and_background_blacklist.tsv")
      if (!is.null(kitome_blacklist_path)) {
        KITOME_BG_BLACKLIST <<- load_kitome_blacklist(path = kitome_blacklist_path)
      } else {
        KITOME_BG_BLACKLIST <<- empty_kitome_blacklist_df()
      }

      clinically_important_pathogens_path <- resolve_app_file("clinically_important_pathogens.tsv")
      if (!is.null(clinically_important_pathogens_path)) {
        CLINICALLY_IMPORTANT_PATHOGENS <<- load_clinically_important_pathogens(path = clinically_important_pathogens_path)
      } else {
        CLINICALLY_IMPORTANT_PATHOGENS <<- empty_clinically_important_pathogens_df()
      }

      loaded <<- TRUE
      invisible(TRUE)
    }
  })

  capture_ui_defaults <- function() {
    list(
      tax_level = input$tax_level %||% "S",
      
      host_species =
        input$host_species %||%
        c("human"),
      
      host_taxids_manual =
        input$host_taxids_manual %||%
        "",
      
      include_protozoa =
        isTRUE(
          input$include_protozoa
        ),
      
      collapse_species =
        isTRUE(
          input$collapse_species
        ),
      
      collapse_min_genus_reads =
        as.numeric(
          input$collapse_min_genus_reads %||%
            30
        ),
      
      collapse_top_frac =
        as.numeric(
          input$collapse_top_frac %||%
            0.85
        ),
      
      fp_true_cutoff =
        as.numeric(
          input$fp_true_cutoff %||%
            1
        ),
      
      min_reads_prevalence =
        as.numeric(
          input$min_reads_prevalence %||%
            5
        ),
      
      min_rpmm_prevalence =
        as.numeric(
          input$min_rpmm_prevalence %||%
            0.2
        ),
      
      min_uniq_kmers_prevalence =
        as.numeric(
          input$min_uniq_kmers_prevalence %||%
            96
        ),
      
      min_kmer_per_read_prevalence =
        as.numeric(
          input$min_kmer_per_read_prevalence %||%
            0.4
        ),
      
      show_fp_breakdown =
        isTRUE(
          input$show_fp_breakdown
        )
    )
  }

  analysis_defaults <- reactiveVal(NULL)
  landing_page <- reactiveVal(TRUE)

  default_organisms <- reactive({
    get_kitome_blacklist_names(KITOME_BG_BLACKLIST)
  })

  clinical_panel_orgs <- reactive({
    get_clinically_important_pathogen_names(CLINICALLY_IMPORTANT_PATHOGENS)
  })

  has_kitome_display <- reactive({
    length(default_organisms() %||% character(0)) > 0
  })

  has_clinical_panel_display <- reactive({
    length(clinical_panel_orgs() %||% character(0)) > 0
  })

  parsed_samples <- reactiveVal(NULL)
  meta_df <- reactiveVal(NULL)
  last_error <- reactiveVal(NULL)
  cleaned_parsed_samples <- reactiveVal(NULL)
  cleaned_last_error <- reactiveVal(NULL)

  host_taxids <- reactive({
    ids <- integer(0)
    if ("human" %in% (input$host_species %||% character(0))) ids <- c(ids, 9606L, 9605L)
    if ("mouse" %in% (input$host_species %||% character(0))) ids <- c(ids, 10090L)
    if ("drosophila" %in% (input$host_species %||% character(0))) ids <- c(ids, 7227L)

    manual <- input$host_taxids_manual %||% ""
    manual <- gsub("[^0-9,]", "", manual)

    if (nzchar(manual)) {
      parts <- unlist(strsplit(manual, ","))
      parts <- parts[nzchar(parts)]
      ids <- c(ids, suppressWarnings(as.integer(parts)))
    }

    ids <- ids[is.finite(ids)]
    unique(ids)
  })

  host_name_patterns <- reactive({
    pats <- character(0)
    if ("human" %in% (input$host_species %||% character(0))) {
      pats <- c(pats, "\\bHomo\\b", "\\bHomo sapiens\\b", "\\bPrimates\\b", "\\bHominidae\\b", "\\bChordata\\b", "\\bMammalia\\b", "\\bMetazoa\\b")
    }
    if ("mouse" %in% (input$host_species %||% character(0))) {
      pats <- c(pats, "\\bMus\\b", "\\bMus musculus\\b", "\\bRodentia\\b")
    }
    if ("drosophila" %in% (input$host_species %||% character(0))) {
      pats <- c(pats, "\\bDrosophila\\b", "\\bDrosophilidae\\b", "\\bInsecta\\b")
    }
    unique(pats)
  })

  params <- reactive({
    list(
      analysis_mode = "Cell-line",
      
      tax_level =
        input$tax_level %||%
        "S",
      
      eps = 1e-9,
      
      min_reads_prevalence =
        as.numeric(
          input$min_reads_prevalence %||%
            5
        ),
      
      min_rpmm_prevalence =
        as.numeric(
          input$min_rpmm_prevalence %||%
            0.2
        ),
      
      min_uniq_kmers_prevalence =
        as.numeric(
          input$min_uniq_kmers_prevalence %||%
            96
        ),
      
      min_kmer_per_read_prevalence =
        as.numeric(
          input$min_kmer_per_read_prevalence %||%
            0.4
        ),
      
      fp_true_cutoff =
        as.numeric(
          input$fp_true_cutoff %||%
            1
        ),
      
      fc_pseudocount = 0.1
    )
  })

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
        out[i] <- paste0(tools::file_path_sans_ext(nm), "_", counts[[nm]], ".", tools::file_ext(nm))
      }
    }

    out
  }

  collect_uploaded_reports <- function(file_df) {
    paths <- character(0)
    names_vec <- character(0)

    if (is.null(file_df) || nrow(file_df) == 0) return(list())

    for (i in seq_len(nrow(file_df))) {
      fname <- file_df$name[i]
      fpath <- file_df$datapath[i]
      if (is.na(fname) || is.na(fpath)) next

      if (grepl("\\.zip$", fname, ignore.case = TRUE)) {
        exdir <- tempfile("celltaminate_zip_")
        dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
        utils::unzip(fpath, exdir = exdir)
        extracted <- list.files(exdir, recursive = TRUE, full.names = TRUE)
        extracted <- extracted[grepl("\\.(txt|tsv|report)$", extracted, ignore.case = TRUE)]
        if (length(extracted) == 0) next
        paths <- c(paths, extracted)
        names_vec <- c(names_vec, basename(extracted))
      } else {
        paths <- c(paths, fpath)
        names_vec <- c(names_vec, fname)
      }
    }

    names_vec <- make_unique_names(names_vec)
    out <- as.list(paths)
    names(out) <- names_vec
    out
  }

  observeEvent(input$analyze_button, {
    analysis_defaults(isolate(capture_ui_defaults()))
    updateProgressBar(session, "analysis_progress", value = 0)
    req(input$files)

    cleaned_parsed_samples(NULL)
    cleaned_last_error(NULL)
    last_error(NULL)
    ensure_reference_assets_loaded()

    files <- collect_uploaded_reports(input$files)

    if (length(files) == 0) {
      last_error("No valid report files found (txt/tsv or zip).")
      return()
    }

    sample_list <- list()
    errors <- c()
    n_files <- length(files)
    idx <- 0

    for (nm in names(files)) {
      idx <- idx + 1
      fpath <- files[[nm]]
      prog <- round((idx - 1) / n_files * 100)
      updateProgressBar(session, "analysis_progress", value = prog)

      parsed <- read_kraken_report(fpath, sample_name = nm)

      if (!is.null(parsed$error)) {
        errors <- c(errors, parsed$error)
        next
      }

      df <- parsed$df
      df <- annotate_lineage(df)
      df <- infer_host_like(df, host_taxids = host_taxids(), host_name_patterns = host_name_patterns())
      df <- infer_plant_like(df)
      df <- infer_microbial(df, include_protozoa = isTRUE(input$include_protozoa))
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

    if (length(sample_list) == 0) {
      last_error(paste(errors, collapse = "\n"))
      return()
    }

    meta0 <- tibble(
      sample = names(sample_list),
      group = "Group 1",
      sample_type = "Clinical / sterile"
    )

    parsed_samples(sample_list)
    meta_df(as.data.frame(meta0))
    last_error(if (length(errors) > 0) paste(errors, collapse = "\n") else NULL)

    updateProgressBar(session, "analysis_progress", value = 100)
    landing_page(FALSE)
  })

  cohort_analysis <- reactive({
    req(parsed_samples())
    ensure_reference_assets_loaded()
    meta <- meta_df()
    req(meta)

    sample_list <- parsed_samples()
    prm <- params()

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

    gs_long <- gs_long %>%
      left_join(meta, by = "sample")

    if (isTRUE(input$collapse_species)) {
      gs_long <- collapse_species_within_genus(
        gs_long,
        min_genus_reads = as.numeric(input$collapse_min_genus_reads %||% 30),
        top_frac_keep = as.numeric(input$collapse_top_frac %||% 0.85)
      )
    }

    user_contam <- character(0)
    clinical_orgs <- clinical_panel_orgs() %||% character(0)
    kitome <- default_organisms() %||% default_kitome()

    engine <- compute_taxon_features(
      gs_long,
      meta,
      prm,
      user_contam = user_contam,
      clinical_panel = clinical_orgs,
      kitome = kitome
    )

    tax_features <- engine$tax_features
    gs_long <- engine$gs_long
    sample_priors <- engine$sample_priors

    gs_long <- apply_calls(gs_long, tax_features, meta, prm, sample_priors)
    for (dbg_nm in unique(gs_long$sample)) {
    }

    n_cohort_samples <- dplyr::n_distinct(meta$sample)
    prevalence_display <- gs_long %>%
      filter(rank %in% c("G", "S"), call != "Non-microbial / Host") %>%
      group_by(rank, name_clean) %>%
      summarise(
        cohort_prevalence_display = ifelse(n_cohort_samples > 0, n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]) / n_cohort_samples, NA_real_),
        .groups = "drop"
      )
    gs_long <- gs_long %>% left_join(prevalence_display, by = c("rank", "name_clean"))

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
        n_prioritized = sum(rank == tax_level_now & call == "Prioritized", na.rm = TRUE),
        n_not_prioritized = sum(rank == tax_level_now & call == "Not prioritized", na.rm = TRUE),
        .groups = "drop"
      )

    cohort_summary <- cohort_summary %>%
      select(-any_of(c("n_genus_qc", "n_species_qc", "shannon", "n_prioritized", "n_not_prioritized"))) %>%
      left_join(qc_counts2, by = "sample")

    list(
      cohort_summary = cohort_summary,
      meta = meta,
      gs_long = gs_long,
      tax_features = tax_features
    )
  })

  output$main_page <- renderUI({
    if (landing_page()) {
      fluidPage(
        div(
          style = "max-width: 1500px; margin: 0 auto; padding: 10px 12px 0 12px;",
          div(
            style = "background: linear-gradient(135deg, #f8fbff 0%, #eef4ff 100%); border: 1px solid #dbe7ff; border-radius: 18px; padding: 22px 24px; margin-top: 16px; margin-bottom: 18px; box-shadow: 0 4px 16px rgba(37, 99, 235, 0.08);",
            fluidRow(
              column(
                12,
                div(
                  style = "text-align:center;",
                  h1(
                    "Celltaminate",
                    class = "celltaminate-title",
                    style = "margin-top: 0; margin-bottom: 8px; font-size: 40px; color: #1e3a8a;"
                  ),
                  div(
                    style = "font-size: 18px; font-weight: 600; color: #475569; margin-bottom: 6px;",
                    "An app for rapid detection of true microorganisms"
                  ),
                  div(
                    style = "font-size: 15px; color: #64748b; max-width: 900px; margin: 0 auto;",
                    "Celltaminate calls true microbial cells from kraken report generated from both illumina sequence and oxford nanopore"
                  )
                )
              )
            )
          ),
          
          fluidRow(
            column(
              4,
              wellPanel(
                style = "background: #ffffff; border: 1px solid #e5e7eb; border-radius: 18px; padding: 18px; box-shadow: 0 4px 14px rgba(0,0,0,0.05);",
                analysis_settings_ui(include_analyze_button = TRUE, defaults = analysis_defaults())
              )
            ),
            column(
              8,
              wellPanel(
                style = "background: #ffffff; border: 1px solid #e5e7eb; border-radius: 18px; padding: 18px; box-shadow: 0 4px 14px rgba(0,0,0,0.05);",
                div(
                  style = "background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%); border: 1px solid #e5e7eb; border-radius: 16px; padding: 16px; margin-bottom: 16px;",
                  fluidRow(
                    column(
                      12,
                      div(
                        style = "text-align:center;",
                        img(
                          src = asset_src("app.jpg"),
                          style = "max-width:720px; width:100%; height:auto; border:1px solid #dbe7eb; border-radius:12px; box-shadow: 0 6px 18px rgba(0,0,0,0.08);"
                        )
                      )
                    )
                  )
                ),
              
              div(
                class = "intro-card",
                style = "background: linear-gradient(135deg, #f8fbff 0%, #eef4ff 100%); border: 1px solid #dbe7ff; border-left: 6px solid #2563eb; border-radius: 14px; padding: 18px 20px; box-shadow: 0 2px 10px rgba(37, 99, 235, 0.08);",
                h4(style = "margin-top:0; margin-bottom:4px; color:#1e3a8a;", "Before you start"),
                tags$p("To run Celltaminate, follow the steps below from left to right.")
              ),
              
              fluidRow(
                column(
                  4,
                  div(
                    class = "intro-card",
                    style = "border-top: 5px solid #2563eb; border-radius: 14px; min-height: 190px; box-shadow: 0 2px 10px rgba(0,0,0,0.05);",
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#1d4ed8; background:#dbeafe; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "Step 1"
                    ),
                    h4("Upload reports"),
                    tags$p(
                      "Upload one or more Kraken reports. You can generate a Kraken report by running Kraken2 locally or on the ",
                      tags$a(
                        "Galaxy server",
                        href = "https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fkraken2%2Fkraken2%2F2.17.1%2Bgalaxy0&version=latest",
                        target = "_blank"
                      ),
                      "."
                    )
                  )
                ),
                column(
                  4,
                  div(
                    class = "intro-card",
                    style = "border-top: 5px solid #2563eb; border-radius: 14px; min-height: 190px; box-shadow: 0 2px 10px rgba(0,0,0,0.05);",
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#1d4ed8; background:#dbeafe; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "Step 2"
                    ),
                    h4("Choose the view"),
                    tags$p("Choose whether to view results at the genus level or the species level.")
                  )
                ),
                column(
                  4,
                  div(
                    class = "intro-card",
                    style = "border-top: 5px solid #2563eb; border-radius: 14px; min-height: 190px; box-shadow: 0 2px 10px rgba(0,0,0,0.05);",
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#1d4ed8; background:#dbeafe; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "Step 3"
                    ),
                    h4("Run the analysis"),
                    tags$p("Click the Run analysis button at the bottom left of the page.")
                  )
                )
              ),
              
              fluidRow(
                column(
                  6,
                  div(
                    class = "intro-card",
                    style = "border-top: 5px solid #2563eb; border-radius: 14px; min-height: 170px; box-shadow: 0 2px 10px rgba(0,0,0,0.05);",
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#1d4ed8; background:#dbeafe; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "Step 4"
                    ),
                    h4("Review the result summary"),
                    tags$p("After the analysis is complete, the app will show a summary overview of all uploaded Kraken reports.")
                  )
                ),
                column(
                  6,
                  div(
                    class = "intro-card",
                    style = "border-top: 5px solid #2563eb; border-radius: 14px; min-height: 170px; box-shadow: 0 2px 10px rgba(0,0,0,0.05);",
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#1d4ed8; background:#dbeafe; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "Step 5"
                    ),
                    h4("Review each sample"),
                    tags$p("Open each sample tab to review the result for each Kraken report."),
                    tags$p("A lower Celltaminate score means the organism behaves more like a true signal in this dataset. A higher Celltaminate score means it behaves more like background.")
                  )
                )
              ),
              
              div(
                style = "background:#f8fafc; border:1px solid #dbeafe; border-left:5px solid #2563eb; border-radius:12px; padding:14px 16px; margin-top:10px; margin-bottom:14px;",
                
                HTML(
                  "<b>Conceptual equation</b><br/>",
                  "Celltaminate score = 100 × plogis(fitted intercept + standardized weighted evidence + score-mapping offset)"
                )
              ),
              
              tags$p(
                HTML(
                  "Celltaminate evaluates quantitative and biological evidence for each organism using a calibrated model. ",
                  "Each fitted feature is transformed and standardized using parameters learned from the development benchmark, multiplied by its fitted coefficient, and combined with the fitted intercept."
                )
              ),
              
              tags$p(
                HTML(
                  "The resulting value is mapped to the 0–100 Celltaminate score. ",
                  "Lower scores indicate stronger support for a prioritized organism, while higher scores indicate stronger background or contamination evidence."
                )
                ),
              
              fluidRow(
                column(
                  6,
                  div(
                    class = "intro-card",
                    style = "background:#f0fdf4; border-left:6px solid #16a34a; border-radius:14px; padding:18px; box-shadow:0 2px 10px rgba(0,0,0,0.06); min-height: 520px;",
                    h4(style = "margin-top:0; color:#166534;", "Organism A"),
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#166534; background:#dcfce7; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "More strongly supported for prioritization because"
                    ),
                    tags$ul(
                      tags$li(
                        "read support is stronger"
                      ),
                      
                      tags$li(
                        "RPMM or abundance is higher"
                      ),
                      
                      tags$li(
                        "unique k-mer support is stronger"
                      ),
                      
                      tags$li(
                        "species-level taxonomic support is less ambiguous"
                      ),
                      
                      tags$li(
                        "the organism is less consistent with abundance patterns in the sterile reference dataset"
                      ),
                      
                      tags$li(
                        "more microbial abundance remains after background subtraction"
                      ),
                      
                      tags$li(
                        "the decontaminated fraction is higher"
                      ),
                      
                      tags$li(
                        "the organism is represented in the clinically important pathogen panel when supported by quantitative sequence evidence"
                      )
                    ),
                    tags$p("These findings support prioritization of this organism over background-like taxa."),
                    div(
                      style = "margin-top:10px; padding:10px 12px; border-radius:10px; font-weight:600; background:#dcfce7; color:#166534;",
                      "Interpretation: the final Celltaminate score is more likely to fall in the prioritized range."
                    )
                  )
                ),
                column(
                  6,
                  div(
                    class = "intro-card",
                    style = "background:#fef2f2; border-left:6px solid #dc2626; border-radius:14px; padding:18px; box-shadow:0 2px 10px rgba(0,0,0,0.06); min-height: 520px;",
                    h4(style = "margin-top:0; color:#991b1b;", "Organism B"),
                    tags$span(
                      style = "display:inline-block; font-size:12px; font-weight:700; color:#991b1b; background:#fee2e2; border-radius:999px; padding:4px 10px; margin-bottom:8px;",
                      "More consistent with a not-prioritized background-like signal because"
                    ),
                    tags$ul(
                      tags$li(
                        "read support is weaker"
                      ),
                      
                      tags$li(
                        "RPMM or abundance is lower"
                      ),
                      
                      tags$li(
                        "unique k-mer support is weaker"
                      ),
                      
                      tags$li(
                        "species-level taxonomic support is more ambiguous"
                      ),
                      
                      tags$li(
                        "the observed abundance is more consistent with the sterile reference background"
                      ),
                      
                      tags$li(
                        "less microbial abundance remains after background subtraction"
                      ),
                      
                      tags$li(
                        "the decontaminated fraction is lower"
                      ),
                      
                      tags$li(
                        "kitome membership can contribute together with weak decontaminated support"
                      )
                    ),
                    tags$p("These findings are more consistent with a background-like signal and therefore do not support prioritization."),
                    div(
                      style = "margin-top:10px; padding:10px 12px; border-radius:10px; font-weight:600; background:#fee2e2; color:#991b1b;",
                      "Interpretation: the final Celltaminate score is more likely to fall above the prioritization cutoff."
                    )
                  )
                )
              )
            )
          )
        ),
        br(),
        fluidRow(
          column(
            12,
            h2("About Celltaminate"),
            tags$p(
              class = "subtle",
              "Celltaminate integrates quantitative sequence evidence, sterile-reference background, decontamination, and curated contextual evidence to summarize each sample in a way that is easier to interpret than a long Kraken report."
            ),
            div(
              class = "intro-card",
              h4("Sample Results"),
              p("These example panels show how the app summarizes each sample. Together they help explain how many taxa were found, how abundant they are, how they compare across samples, and how strongly the score supports a prioritized versus not-prioritized microbial signal.")
            ),
            fluidRow(
              column(
                6,
                div(
                  class = "preview-card",
                  img(src = asset_src("Summary.png"), style = "width:100%; border-radius:6px;"),
                  div(class = "preview-title", "Sample Summary Analysis"),
                  div(class = "preview-text", "A quick overview of reads, microbial burden, diversity, and how many organisms fall into each score-based category.")
                )
              ),
              column(
                6,
                div(
                  class = "preview-card",
                  img(src = asset_src("barplot.jpg"), style = "width:100%; border-radius:6px;"),
                  div(class = "preview-title", "Sample Barplot"),
                  div(class = "preview-text", "Highlights the top organisms in a sample so users can quickly see which taxa dominate the microbial signal.")
                )
              )
            ),
            fluidRow(
              column(
                6,
                div(
                  class = "preview-card",
                  img(src = asset_src("boxplot.png"), style = "width:100%; border-radius:6px;"),
                  div(class = "preview-title", "Sample Boxplot"),
                  div(class = "preview-text", "Shows how score distributions or abundance patterns compare across call categories and across the cohort.")
                )
              ),
              column(
                6,
                div(
                  class = "preview-card",
                  img(src = asset_src("spiderplot.png"), style = "width:100%; border-radius:6px;"),
                  div(class = "preview-title", "Sample Spider Plot"),
                  div(class = "preview-text", "Compares a sample against reference or cohort context across several organisms at the same time.")
                )
              )
            )
          )
        )
      )
    )
    } else {
      fluidPage(
        fluidRow(
          column(
            12,
            div(
              style = "margin:20px; text-align:right;",
              actionButton("return_home_btn", "Return to Home", icon = icon("home"))
            )
          )
        ),
        uiOutput("results_tabs")
      )
    }
  })

  output$results_tabs <- renderUI({
    req(parsed_samples())

    sample_names <- names(parsed_samples())
    tabs <- list()

    tabs <- c(tabs, list(
      tabPanel(
        title = "Cohort Overview",
        fluidRow(
          column(
            4,
            wellPanel(
              h3("Settings"),
              tags$p(class = "small-note", "The selected taxonomy level is applied across all cohort and sample outputs."),
              analysis_settings_ui(include_analyze_button = FALSE, include_file_inputs = FALSE, defaults = analysis_defaults()),
              tags$hr(),
              h4("Groups"),
              tags$p(class = "small-note", "You can edit sample groups in the metadata table or use auto-clustering."),
              actionButton("auto_cluster_btn", "Auto-cluster samples", icon = icon("project-diagram")),
              numericInput("n_clusters", "Number of clusters", value = 2, min = 2, max = 10, step = 1),
              tags$hr(),
              h4("Downloads"),
              downloadButton("download_cohort_summary", "Download cohort summary (.xlsx)"),
              br(), br(),
              downloadButton("download_all_cleaned_zip", "Download all cleaned reports (.zip)"),
              tags$p(class = "small-note", "Cleaned reports keep prioritized taxa by default and preserve any additional taxa you explicitly keep in each sample tab.")
            )
          ),
          column(
            8,
            wellPanel(
              tabsetPanel(
                tabPanel(
                  "Original cohort",
                  h2("Cohort Summary"),
                  verbatimTextOutput("cohort_errors"),
                  br(),
                  dataTableOutput("cohort_summary_table"),
                  br(),
                  h3("Cohort PCA"),
                  plotOutput("cohort_pca_plot", height = "420px"),
                  br(),
                  h3("Prevalence vs abundance"),
                  plotOutput("cohort_prev_plot", height = "420px"),
                  br(),
                  h3("Cohort heatmap (top prioritized taxa)"),
                  plotOutput("cohort_heatmap_plot", height = "520px"),
                  br(),
                  h3("Celltaminate score summaries"),
                  plotOutput("cohort_fp_by_sample_plot", height = "360px"),
                  br(),
                  plotOutput("cohort_fp_by_call_plot", height = "360px"),
                  br(),
                  h3("Sample metadata (editable)"),
                  dataTableOutput("meta_table")
                )
              )
            )
          )
        )
      )
    ))

    for (samp in sample_names) {
      sid <- safe_sample_id(samp)

      tabs <- c(tabs, list(
        tabPanel(
          title = samp,
          h3(paste("Results for:", samp)),
          htmlOutput(paste0("sample_summary_", sid)),
          br(),
          fluidRow(
            column(
              8,
              h4("Taxa table"),
              selectInput(
                inputId = paste0("taxa_view_", sid),
                label = "Show",
                choices = c(
                  "All" = "all",
                  "Prioritized" = "prioritized",
                  "Not prioritized" = "not_prioritized"
                ),
                selected = "all"
              ),
              dataTableOutput(paste0("taxa_table_", sid))
            ),
            column(
              4,
              h4("Top taxa"),
              checkboxGroupInput(
                inputId = paste0("top_taxa_calls_", sid),
                label = "Show categories",
                choices = c(
                  "Prioritized" = "Prioritized",
                  "Not prioritized" = "Not prioritized"
                ),
                selected = c("Prioritized")
              ),
              plotOutput(paste0("top_taxa_plot_", sid), height = "320px"),
              br(),
              h4("Celltaminate score vs abundance"),
              plotOutput(paste0("fp_scatter_plot_", sid), height = "320px"),
              br(),
              h4("Spider Plot"),
              sliderInput(
                inputId = paste0("spider_n_", sid),
                label = "Max taxa to show",
                min = 3,
                max = 10,
                value = 5,
                step = 1
              ),
              checkboxGroupInput(
                inputId = paste0("spider_calls_", sid),
                label = "Show categories",
                choices = c(
                  "Prioritized" = "Prioritized",
                  "Not prioritized" = "Not prioritized"
                ),
                selected = c("Prioritized")
              ),
              plotOutput(paste0("radar_plot_", sid), height = "340px")
            )
          ),
          tags$hr(),
          h3("Decontamination Feature"),
          helpText("Not-prioritized taxa are removed from cleaned reports by default. Use the search box below to keep specific taxa if needed."),
          uiOutput(paste0("decontam_ui_", sid)),
          tags$hr(),
          h3("BioAI interpretation"),
          uiOutput(paste0("bioai_ui_", sid)),
          uiOutput(paste0("bioai_out_", sid))
        )
      ))
    }

    do.call(tabsetPanel, tabs)
  })

  output$cohort_errors <- renderText({
    last_error() %||% ""
  })

  output$cohort_summary_table <- DT::renderDataTable({
    ca <- cohort_analysis()
    req(ca)

    df <- ca$cohort_summary %>%
      mutate(
        total_reads = as.integer(round(total_reads)),
        microbial_reads = as.integer(round(microbial_reads)),
        n_genus_raw = as.integer(round(n_genus_raw)),
        n_species_raw = as.integer(round(n_species_raw)),
        n_genus_qc = as.integer(round(n_genus_qc)),
        n_species_qc = as.integer(round(n_species_qc)),
        n_prioritized = as.integer(round(n_prioritized)),
        n_not_prioritized = as.integer(round(n_not_prioritized)),
        n_not_prioritized = as.integer(round(n_not_prioritized)),
        shannon = round(shannon, 3)
      )

    DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$download_cohort_summary <- downloadHandler(
    filename = function() paste0("celltaminate_cohort_summary_", Sys.Date(), ".xlsx"),
    content = function(file) {
      ca <- cohort_analysis()
      req(ca)

      writexl::write_xlsx(
        list(
          cohort_summary = ca$cohort_summary,
          sample_metadata = ca$meta %>% select(sample, group, sample_type)
        ),
        path = file
      )
    }
  )

  output$meta_table <- DT::renderDataTable({
    df <- meta_df()
    req(df)

    df_view <- df %>%
      select(sample, group, sample_type)

    DT::datatable(
      df_view,
      rownames = FALSE,
      editable = list(target = "cell", disable = list(columns = c(0))),
      options = list(pageLength = 20)
    )
  })

  observeEvent(input$meta_table_cell_edit, {
    info <- input$meta_table_cell_edit
    df <- meta_df()
    req(df)

    i <- info$row
    j <- info$col
    v <- info$value

    if (j == 1) {
      df$group[i] <- as.character(v)
    } else if (j == 2) {
      df$sample_type[i] <- as.character(v)
    }

    meta_df(df)
  })

  observeEvent(input$auto_cluster_btn, {
    ca <- cohort_analysis()
    df <- ca$gs_long
    meta <- meta_df()
    req(meta)

    tax_level <- params()$tax_level %||% "S"

    mat_df <- df %>%
      filter(rank == tax_level, call != "Non-microbial / Host", qc_pass) %>%
      group_by(sample, name_clean) %>%
      summarize(val = max(log_rpmm, na.rm = TRUE), .groups = "drop")

    if (nrow(mat_df) == 0) return()

    wide <- tidyr::pivot_wider(mat_df, names_from = name_clean, values_from = val, values_fill = 0)
    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$sample

    k <- as.integer(input$n_clusters %||% 2)
    k <- max(2, min(k, 10))

    set.seed(1)
    km <- stats::kmeans(m, centers = k)
    clusters <- paste0("Cluster ", km$cluster)

    meta2 <- meta
    meta2$group <- clusters[match(meta2$sample, names(km$cluster))]
    meta_df(meta2)
  })

  output$cohort_pca_plot <- renderPlot({
    ca <- cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", qc_pass) %>%
      group_by(sample, name_clean) %>%
      summarize(val = max(log_rpmm, na.rm = TRUE), .groups = "drop")

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No QC-pass taxa to build PCA.")
      return()
    }

    wide <- tidyr::pivot_wider(df, names_from = name_clean, values_from = val, values_fill = 0)

    if (!can_run_pca(wide)) {
      plot.new()
      text(0.5, 0.5, "PCA needs at least 2 samples and at least 2 variable taxa.")
      return()
    }

    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$sample

    pca <- stats::prcomp(m, center = TRUE, scale. = TRUE)
    pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    pcs$sample <- rownames(pcs)
    pcs <- left_join(pcs, ca$meta, by = "sample")

    ggplot(pcs, aes(x = PC1, y = PC2, label = sample, color = group)) +
      geom_point(size = 3, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
      theme_minimal() +
      labs(
        title = paste0("Cohort PCA (", ifelse(tax_level == "G", "Genus", "Species"), " level)"),
        color = "Group"
      )
  })

  output$cohort_prev_plot <- renderPlot({
    ca <- cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    tf <- make_prevalence_abundance_df(ca$gs_long, tax_level = tax_level, total_samples = n_distinct(ca$meta$sample))

    if (nrow(tf) == 0) {
      plot.new()
      text(0.5, 0.5, "No cohort taxa available.")
      return()
    }

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
        title = "Prevalence vs abundance",
        x = "Prevalence across samples",
        y = "Mean log(RPMM) where present",
        color = NULL
      )
  })

  output$cohort_heatmap_plot <- renderPlot({
    ca <- cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    calls_prioritized <- c("Prioritized")

    gs <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm))

    focus <- gs %>%
      filter(call %in% calls_prioritized)

    if (nrow(focus) == 0) {
      focus <- gs %>% filter(call != "Not prioritized")
    }

    if (nrow(focus) == 0) {
      plot.new()
      text(0.5, 0.5, "No taxa available for heatmap.")
      return()
    }

    top_taxa <- focus %>%
      group_by(name_clean) %>%
      summarize(mean_rpmm = mean(rpmm, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_rpmm)) %>%
      head(30)

    heat <- gs %>%
      filter(name_clean %in% top_taxa$name_clean) %>%
      select(sample, name_clean, rpmm) %>%
      mutate(log_rpmm_tile = log10(rpmm + 1))

    heat <- left_join(heat, ca$meta %>% select(sample, group), by = "sample")

    samp_order <- ca$meta %>%
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
        title = paste0("Heatmap of top taxa (", ifelse(tax_level == "G", "Genus", "Species"), ")"),
        x = NULL,
        y = NULL,
        fill = "log10(RPMM+1)"
      )
  })

  output$cohort_fp_by_sample_plot <- renderPlot({
    ca <- cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score), qc_pass) %>%
      group_by(sample) %>%
      summarize(
        median_fp = median(fp_score, na.rm = TRUE),
        mean_fp = mean(fp_score, na.rm = TRUE),
        n_taxa = dplyr::n(),
        .groups = "drop"
      ) %>%
      left_join(ca$meta, by = "sample")

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No Celltaminate scores available.")
      return()
    }

    samp_order <- ca$meta %>%
      arrange(group, sample) %>%
      pull(sample)

    df$sample <- factor(df$sample, levels = samp_order)

    ggplot(df, aes(x = sample, y = median_fp, color = group)) +
      geom_point(size = 3, alpha = 0.9) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Median Celltaminate score per sample (QC-pass taxa)",
        x = NULL,
        y = "Median Celltaminate score (0-100)",
        color = "Group"
      )
  })

  output$cohort_fp_by_call_plot <- renderPlot({
    ca <- cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    calls_prioritized <- c("Prioritized")

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score)) %>%
      mutate(call_bucket = call_bucket_label(call))

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No Celltaminate scores available.")
      return()
    }

    ggplot(df, aes(x = call_bucket, y = fp_score)) +
      geom_boxplot(outlier.size = 0.6) +
      theme_minimal() +
      labs(
        title = "Celltaminate score distribution by call category",
        x = NULL,
        y = "Celltaminate score (0-100)"
      )
  })

  observeEvent(input$analyze_cleaned_button, {
    req(parsed_samples())
    ca <- cohort_analysis()
    req(ca)

    cleaned_last_error(NULL)
    sample_list <- parsed_samples()

    if (is.null(sample_list) || length(sample_list) == 0) {
      cleaned_last_error("No samples loaded.")
      return()
    }

    withProgress(message = "Analyzing decontaminated cohort...", value = 0, {
      tmpdir <- tempfile("celltaminate_cleaned_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

      cleaned_list <- list()
      errs <- c()
      n <- length(sample_list)

      for (idx in seq_along(sample_list)) {
        samp <- names(sample_list)[idx]
        incProgress(1 / max(n, 1), detail = samp)

        sid <- safe_sample_id(samp)
        sel_id <- paste0("decontam_select_", sid)
        keep_desc_id <- paste0("decontam_keep_desc_", sid)

        gs <- ca$gs_long %>% filter(sample == samp)

        protected_taxa <- input[[sel_id]] %||% character(0)
        keep_desc <- isTRUE(input[[keep_desc_id]] %||% FALSE)

        raw_df <- sample_list[[samp]]$raw
        cleaned_raw <- build_cleaned_report_for_sample(
          raw_df = raw_df,
          gs = gs,
          protected_taxa = protected_taxa,
          keep_descendants = keep_desc
        )

        out_path <- file.path(tmpdir, paste0("cleaned_", sid, ".txt"))
        write.table(cleaned_raw, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

        parsed <- tryCatch(
          read_kraken_report(out_path, sample_name = samp),
          error = function(e) {
            errs <<- c(errs, paste0(samp, ": ", conditionMessage(e)))
            return(NULL)
          }
        )

        if (is.null(parsed)) next

        df <- parsed$df
        df <- annotate_lineage(df)
        df <- infer_host_like(df, host_taxids = host_taxids(), host_name_patterns = host_name_patterns())
        df <- infer_plant_like(df)
        df <- infer_microbial(df, include_protozoa = isTRUE(input$include_protozoa))
        df <- annotate_top_strain_rows(df)
        summ <- compute_reads_summary(df)

        cleaned_list[[samp]] <- list(
          path = out_path,
          raw = parsed$raw,
          df = df,
          report_format = parsed$report_format,
          summary = summ
        )
      }

      cleaned_parsed_samples(cleaned_list)
      if (length(errs) > 0) cleaned_last_error(paste(head(errs, 10), collapse = "\n"))
    })
  })

  observeEvent(input$clear_cleaned_button, {
    cleaned_parsed_samples(NULL)
    cleaned_last_error(NULL)
  })

  cleaned_cohort_analysis <- reactive({
    req(cleaned_parsed_samples())
    req(meta_df())
    ensure_reference_assets_loaded()

    meta <- meta_df()
    sample_list <- cleaned_parsed_samples()
    prm <- params()

    user_contam <- character(0)
    clinical_orgs <- clinical_panel_orgs() %||% character(0)
    kitome <- default_organisms() %||% default_kitome()

    gs_list <- list()
    cohort_rows <- list()

    for (nm in names(sample_list)) {
      df <- sample_list[[nm]]$df
      summ <- sample_list[[nm]]$summary

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

      cohort_rows[[nm]] <- tibble(
        sample = nm,
        report_format = sample_list[[nm]]$report_format,
        total_reads = summ$total_reads,
        microbial_reads = summ$microbial_reads,
        host_reads = summ$host_reads,
        n_genus_raw = n_distinct(df$name_clean[df$rank == "G" & df$reads_clade > 0]),
        n_species_raw = n_distinct(df$name_clean[df$rank == "S" & df$reads_clade > 0])
      )
    }

    gs_long <- bind_rows(gs_list)

    meta <- meta %>%
      mutate(
        group = as.character(group)
      )

    gs_long <- left_join(gs_long, meta, by = "sample")

    if (isTRUE(input$collapse_species)) {
      gs_long <- collapse_species_within_genus(
        gs_long,
        min_genus_reads = as.numeric(input$collapse_min_genus_reads %||% 30),
        top_frac_keep = as.numeric(input$collapse_top_frac %||% 0.85)
      )
    }

    engine2 <- compute_taxon_features(
      gs_long,
      meta_df = meta,
      params = prm,
      user_contam = user_contam,
      clinical_panel = clinical_orgs,
      kitome = kitome
    )

    gs_long2 <- engine2$gs_long

    gs_called <- apply_calls(
      gs_long2,
      tax_features = engine2$tax_features,
      meta_df = meta,
      params = prm,
      sample_priors = engine2$sample_priors
    )

    n_cohort_samples <- dplyr::n_distinct(meta$sample)
    prevalence_display <- gs_called %>%
      filter(rank %in% c("G", "S"), call != "Non-microbial / Host") %>%
      group_by(rank, name_clean) %>%
      summarise(
        cohort_prevalence_display = ifelse(n_cohort_samples > 0, n_distinct(sample[is.finite(reads_clade) & reads_clade > 0]) / n_cohort_samples, NA_real_),
        .groups = "drop"
      )
    gs_called <- gs_called %>% left_join(prevalence_display, by = c("rank", "name_clean"))

    tax_features <- engine2$tax_features
    tax_level_now <- prm$tax_level %||% "S"

    qc_counts <- gs_called %>%
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
        n_prioritized = sum(rank == tax_level_now & call == "Prioritized", na.rm = TRUE),
        n_not_prioritized = sum(rank == tax_level_now & call == "Not prioritized", na.rm = TRUE),
        .groups = "drop"
      )

    cohort_summary <- bind_rows(cohort_rows) %>%
      left_join(qc_counts, by = "sample")

    list(
      gs_long = gs_called,
      cohort_summary = cohort_summary,
      meta = meta,
      tax_features = tax_features
    )
  })

  output$cleaned_cohort_status <- renderUI({
    if (is.null(cleaned_parsed_samples())) {
      return(tags$div(class = "small-note", "No decontaminated cohort results yet. Click 'Analyze decontaminated cohort' after selecting taxa to remove in each sample tab."))
    }

    err <- cleaned_last_error() %||% ""

    tagList(
      h2("Decontaminated Cohort Summary"),
      if (nzchar(err)) tags$pre(style = "color:#b91c1c;", err) else NULL,
      tags$p(class = "small-note", "These outputs are computed after applying your per-sample decontamination selections, then re-running the Celltaminate algorithm.")
    )
  })

  output$cleaned_summary_table <- DT::renderDataTable({
    ca <- cleaned_cohort_analysis()
    req(ca)

    df <- ca$cohort_summary %>%
      mutate(
        total_reads = as.integer(round(total_reads)),
        microbial_reads = as.integer(round(microbial_reads)),
        host_reads = as.integer(round(host_reads)),
        n_genus_raw = as.integer(round(n_genus_raw)),
        n_species_raw = as.integer(round(n_species_raw)),
        n_genus_qc = as.integer(round(n_genus_qc)),
        n_species_qc = as.integer(round(n_species_qc)),
        n_prioritized = as.integer(round(n_prioritized)),
        n_not_prioritized = as.integer(round(n_not_prioritized)),
        n_not_prioritized = as.integer(round(n_not_prioritized)),
        shannon = round(shannon, 3)
      )

    DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$cleaned_pca_plot <- renderPlot({
    ca <- cleaned_cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", qc_pass) %>%
      group_by(sample, name_clean) %>%
      summarize(val = max(log_rpmm, na.rm = TRUE), .groups = "drop")

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No QC-pass taxa to build PCA.")
      return()
    }

    wide <- tidyr::pivot_wider(df, names_from = name_clean, values_from = val, values_fill = 0)

    if (!can_run_pca(wide)) {
      plot.new()
      text(0.5, 0.5, "PCA needs at least 2 samples and at least 2 variable taxa.")
      return()
    }

    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$sample

    pca <- stats::prcomp(m, center = TRUE, scale. = TRUE)
    pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    pcs$sample <- rownames(pcs)
    pcs <- left_join(pcs, ca$meta, by = "sample")

    ggplot(pcs, aes(x = PC1, y = PC2, label = sample, color = group)) +
      geom_point(size = 3, alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
      theme_minimal() +
      labs(
        title = paste0("Decontaminated cohort PCA (", ifelse(tax_level == "G", "Genus", "Species"), " level)"),
        color = "Group"
      )
  })

  output$cleaned_prev_plot <- renderPlot({
    ca <- cleaned_cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    tf <- make_prevalence_abundance_df(ca$gs_long, tax_level = tax_level, total_samples = n_distinct(ca$meta$sample))

    if (nrow(tf) == 0) {
      plot.new()
      text(0.5, 0.5, "No cohort taxa available.")
      return()
    }

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
        title = "Prevalence vs abundance (decontaminated)",
        x = "Prevalence across samples",
        y = "Mean log(RPMM) where present",
        color = NULL
      )
  })

  output$cleaned_heatmap_plot <- renderPlot({
    ca <- cleaned_cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    calls_prioritized <- c("Prioritized")

    gs <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm))

    focus <- gs %>% filter(call %in% calls_prioritized)

    if (nrow(focus) == 0) {
      focus <- gs %>% filter(call != "Not prioritized")
    }

    if (nrow(focus) == 0) {
      plot.new()
      text(0.5, 0.5, "No taxa available for heatmap.")
      return()
    }

    top_taxa <- focus %>%
      group_by(name_clean) %>%
      summarize(mean_rpmm = mean(rpmm, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_rpmm)) %>%
      head(30)

    heat <- gs %>%
      filter(name_clean %in% top_taxa$name_clean) %>%
      select(sample, name_clean, rpmm) %>%
      mutate(log_rpmm_tile = log10(rpmm + 1))

    heat <- left_join(heat, ca$meta %>% select(sample, group), by = "sample")

    samp_order <- ca$meta %>%
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
        title = paste0("Decontaminated heatmap (", ifelse(tax_level == "G", "Genus", "Species"), ")"),
        x = NULL,
        y = NULL,
        fill = "log10(RPMM+1)"
      )
  })

  output$cleaned_fp_by_sample_plot <- renderPlot({
    ca <- cleaned_cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score), qc_pass) %>%
      group_by(sample) %>%
      summarize(
        median_fp = median(fp_score, na.rm = TRUE),
        mean_fp = mean(fp_score, na.rm = TRUE),
        n_taxa = dplyr::n(),
        .groups = "drop"
      ) %>%
      left_join(ca$meta, by = "sample")

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No Celltaminate scores available.")
      return()
    }

    samp_order <- ca$meta %>%
      arrange(group, sample) %>%
      pull(sample)

    df$sample <- factor(df$sample, levels = samp_order)

    ggplot(df, aes(x = sample, y = median_fp, color = group)) +
      geom_point(size = 3, alpha = 0.9) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Median Celltaminate score per sample (QC-pass taxa; decontaminated)",
        x = NULL,
        y = "Median Celltaminate score (0-100)",
        color = "Group"
      )
  })

  output$cleaned_fp_by_call_plot <- renderPlot({
    ca <- cleaned_cohort_analysis()
    req(ca)

    tax_level <- params()$tax_level %||% "S"
    calls_prioritized <- c("Prioritized")

    df <- ca$gs_long %>%
      filter(rank == tax_level, call != "Non-microbial / Host", is.finite(fp_score)) %>%
      mutate(call_bucket = call_bucket_label(call))

    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No Celltaminate scores available.")
      return()
    }

    ggplot(df, aes(x = call_bucket, y = fp_score)) +
      geom_boxplot(outlier.size = 0.6) +
      theme_minimal() +
      labs(
        title = "Celltaminate score distribution by call category (decontaminated)",
        x = NULL,
        y = "Celltaminate score (0-100)"
      )
  })

  output$download_all_cleaned_zip <- downloadHandler(
    filename = function() paste0("celltaminate_cleaned_reports_", Sys.Date(), ".zip"),
    content = function(file) {
      req(parsed_samples())
      ca <- cohort_analysis()
      req(ca)

      tmpdir <- tempfile("celltaminate_cleaned_zip_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

      gs <- ca$gs_long
      sample_list <- parsed_samples()
      out_files <- c()

      for (samp in names(sample_list)) {
        sid <- safe_sample_id(samp)
        sel_id <- paste0("decontam_select_", sid)
        keep_desc_id <- paste0("decontam_keep_desc_", sid)

        gs_samp <- gs %>% filter(sample == samp)
        protected_taxa <- input[[sel_id]] %||% character(0)
        keep_desc <- isTRUE(input[[keep_desc_id]] %||% FALSE)

        raw_df <- sample_list[[samp]]$raw
        cleaned <- build_cleaned_report_for_sample(
          raw_df = raw_df,
          gs = gs_samp,
          protected_taxa = protected_taxa,
          keep_descendants = keep_desc
        )

        out_path <- file.path(tmpdir, paste0("cleaned_", samp, ".txt"))
        write.table(cleaned, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        out_files <- c(out_files, out_path)
      }

      zip_files_safely(out_files, zipfile = file)
    }
  )

  observeEvent(cohort_analysis(), {
    ca <- cohort_analysis()
    sample_list <- parsed_samples()
    req(sample_list)

    for (samp in names(sample_list)) {
      local({
        sample_name <- samp
        sid <- safe_sample_id(sample_name)
        bioai_text <- reactiveVal(NULL)

        output[[paste0("sample_summary_", sid)]] <- renderUI({
          ca <- cohort_analysis()
          cs <- ca$cohort_summary %>% filter(sample == sample_name)
          meta <- ca$meta %>% filter(sample == sample_name)

          if (nrow(cs) == 0) return(NULL)

          n_prioritized <- cs$n_prioritized %||% 0
          n_cont <- cs$n_not_prioritized %||% 0
          n_unc <- cs$n_not_prioritized %||% 0

          pills <- tagList(
            span(class = "call-pill pill-true", paste0("Prioritized: ", n_prioritized)),
            span(style = "margin-left:8px;", class = "call-pill pill-cont", paste0("Not prioritized: ", n_cont))
          )

          tagList(
            HTML(paste0(
              "<div class='subtle'>",
              "<b>Microbial reads:</b> ", format(as.integer(round(cs$microbial_reads)), scientific = FALSE, big.mark = ","), " &nbsp; | &nbsp; ",
              "<b>Total reads:</b> ", format(as.integer(round(cs$total_reads)), scientific = FALSE, big.mark = ","), "<br/>",
              "<b>Unique taxa (raw):</b> Genus=", as.integer(round(cs$n_genus_raw)), ", Species=", as.integer(round(cs$n_species_raw)), " &nbsp; | &nbsp; ",
              "<b>Unique taxa (scored):</b> Genus=", as.integer(round(cs$n_genus_qc)), ", Species=", as.integer(round(cs$n_species_qc)), "<br/>",
              "<b>Group:</b> ", htmltools::htmlEscape(meta$group[1] %||% "Group 1"), " &nbsp; | &nbsp; ",
              "<b>Shannon (species RPMM):</b> ", round(cs$shannon, 3),
              "</div>"
            )),
            pills
          )
        })

        output[[paste0("taxa_table_", sid)]] <- DT::renderDataTable({
          ca <- cohort_analysis()
          req(ca)
          req(parsed_samples())
          
          tax_level <- params()$tax_level %||% "S"
          
          gs <- ca$gs_long %>%
            filter(sample == sample_name, rank == tax_level, call != "Non-microbial / Host")
          
          view_mode <- input[[paste0("taxa_view_", sid)]] %||% "all"
          
          gs <- gs %>% mutate(call = as.character(call))
          
          if (view_mode == "prioritized") {
            gs <- gs %>% filter(call == "Prioritized")
          } else if (view_mode == "not_prioritized") {
            gs <- gs %>% filter(call == "Not prioritized")
          }
          
          gs_view <- prepare_taxa_table_data(
            gs = gs,
            tax_level = tax_level,
            show_fp_breakdown = isTRUE(input$show_fp_breakdown),
            show_kitome = isTRUE(has_kitome_display()),
            show_clinical_panel = isTRUE(has_clinical_panel_display())
          )

          DT::datatable(
            gs_view,
            rownames = FALSE,
            escape = TRUE,
            selection = "none",
            options = list(
              pageLength = 25,
              scrollX = TRUE,
              autoWidth = TRUE,
              deferRender = FALSE
            )
          )
        }, server = FALSE)

        output[[paste0("top_taxa_plot_", sid)]] <- renderPlot({
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"

          sel_calls <- input[[paste0("top_taxa_calls_", sid)]] %||% character(0)

          if (length(sel_calls) == 0) {
            plot.new()
            text(0.5, 0.5, "No categories selected.")
            return()
          }

          gs <- ca$gs_long %>%
            filter(sample == sample_name, rank == tax_level, call %in% sel_calls, is.finite(rpmm))

          if (nrow(gs) == 0) {
            plot.new()
            text(0.5, 0.5, "No taxa in selected categories.")
            return()
          }

          gs <- gs %>%
            mutate(
              call_bucket = call_bucket_label(call),
              y = log(rpmm + 1)
            )

          topn <- gs %>%
            mutate(
              call_priority = dplyr::case_when(
                call == "Prioritized" ~ 1L,
                call == "Not prioritized" ~ 2L,
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
            labs(x = NULL, y = "log(RPMM+1)", fill = NULL, title = "Top taxa")
        })

        output[[paste0("fp_scatter_plot_", sid)]] <- renderPlot({
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"

          gs <- ca$gs_long %>%
            filter(sample == sample_name, rank == tax_level, call != "Non-microbial / Host", is.finite(rpmm), is.finite(fp_score)) %>%
            mutate(
              call_bucket = call_bucket_label(call),
              log10_rpmm = log10(rpmm + 1)
            )

          if (nrow(gs) == 0) {
            plot.new()
            text(0.5, 0.5, "No taxa available.")
            return()
          }

          ggplot(gs, aes(x = log10_rpmm, y = fp_score, color = call_bucket)) +
            geom_point(alpha = 0.7) +
            theme_minimal() +
            scale_color_manual(values = CALL_PALETTE, breaks = CALL_BUCKET_LEVELS) +
            labs(
              title = "Celltaminate score vs abundance",
              x = "log10(RPMM+1)",
              y = "Celltaminate score (0-100)",
              color = NULL
            )
        })

        output[[paste0("radar_plot_", sid)]] <- renderPlot({
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"

          sel_calls <- input[[paste0("spider_calls_", sid)]] %||% character(0)

          n_show <- as.integer(input[[paste0("spider_n_", sid)]] %||% 5)
          n_show <- max(3, min(10, n_show))

          gs_samp <- ca$gs_long %>%
            filter(sample == sample_name, rank == tax_level, call %in% sel_calls, is.finite(log_rpmm))

          if (nrow(gs_samp) == 0) {
            plot.new()
            text(0.5, 0.5, "No taxa in selected categories.")
            return()
          }

          top <- gs_samp %>%
            arrange(desc(rpmm)) %>%
            head(n_show)

          taxa <- top$name_clean

          if (length(taxa) < 3) {
            plot.new()
            text(0.5, 0.5, "Not enough taxa for spider plot.")
            return()
          }

          tf <- ca$tax_features %>%
            filter(rank == tax_level, name_clean %in% taxa)

          ref_med <- tf$ref_q50_log[match(taxa, tf$name_clean)]
          ref_med[!is.finite(ref_med)] <- 0

          ctrl_val <- rep(0, length(taxa))
          samp_vals <- top$log_rpmm
          samp_vals[!is.finite(samp_vals)] <- 0

          all_vals <- c(ref_med, ctrl_val, samp_vals)
          mx <- max(all_vals, na.rm = TRUE)
          if (!is.finite(mx) || mx <= 0) mx <- 1

          ref_scaled <- ref_med / mx
          samp_scaled <- samp_vals / mx

          rows <- list(
            rep(1, length(taxa)),
            rep(0, length(taxa)),
            ref_scaled,
            samp_scaled
          )
          row_names <- c("max", "min", "Reference median", "User sample")
          cols <- c("black", "red")

          df_radar <- as.data.frame(do.call(rbind, rows))
          colnames(df_radar) <- taxa
          rownames(df_radar) <- row_names

          layout_cfg <- spider_plot_layout(taxa)

          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar), add = TRUE)

          par(mai = layout_cfg$mai, xpd = NA)
          fmsb::radarchart(
            df_radar,
            axistype = 1,
            pcol = cols,
            plty = rep(1, length(cols)),
            plwd = rep(2, length(cols)),
            cglcol = "grey",
            cglty = 1,
            axislabcol = "grey",
            vlcex = 0,
            vlabels = rep("", length(taxa)),
            title = "Reference median vs sample (scaled log-RPMM)"
          )

          n_taxa <- length(taxa)
          theta <- pi / 2 - 2 * pi * ((seq_len(n_taxa) - 1) / n_taxa)
          label_r <- layout_cfg$label_radius
          x_lab <- label_r * cos(theta)
          y_lab <- label_r * sin(theta)
          angle_deg <- theta * 180 / pi

          label_rot <- ifelse(angle_deg < -90 | angle_deg > 90, angle_deg + 180, angle_deg)
          label_adj <- ifelse(angle_deg < -90 | angle_deg > 90, 1, 0)

          for (i in seq_len(n_taxa)) {
            graphics::text(
              x = x_lab[i],
              y = y_lab[i],
              labels = taxa[i],
              srt = label_rot[i],
              adj = c(label_adj[i], 0.5),
              cex = layout_cfg$label_cex
            )
          }

          legend(
            x = "topright",
            legend = row_names[-c(1, 2)],
            col = cols,
            lty = 1,
            bty = "n",
            cex = 0.8
          )
        })

        output[[paste0("decontam_ui_", sid)]] <- renderUI({
          ca <- cohort_analysis()
          gs <- ca$gs_long %>% filter(sample == sample_name)
          tagList(
            checkboxInput(paste0("decontam_keep_desc_", sid), "For kept genera include descendant species", value = FALSE),
            selectizeInput(
              inputId = paste0("decontam_select_", sid),
              label = "Select additional taxa to keep",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(placeholder = "Type to search...", maxOptions = 10000)
            ),
            tags$p(class = "small-note", paste0("Cleaned reports remove not-prioritized taxa by default. Prioritized taxa are kept.")),
            downloadButton(paste0("downloadDecontaminated_", sid), "Download decontaminated report")
          )
        })

        outputOptions(output, paste0("decontam_ui_", sid), suspendWhenHidden = FALSE)

        observeEvent(cohort_analysis(), {
          ca <- cohort_analysis()
          gs <- ca$gs_long %>% filter(sample == sample_name)

          choices <- gs %>%
            filter(is.na(call) | call != "Non-microbial / Host") %>%
            pull(name_clean) %>%
            unique()

          choices <- sort(choices)
          current <- input[[paste0("decontam_select_", sid)]]
          sel <- if (is.null(current)) character(0) else current
          choices <- unname(as.character(choices))

          updateSelectizeInput(
            session,
            inputId = paste0("decontam_select_", sid),
            choices = choices,
            selected = sel,
            server = TRUE
          )
        }, ignoreInit = FALSE)

        output[[paste0("downloadDecontaminated_", sid)]] <- downloadHandler(
          filename = function() paste0("decontaminated_", sample_name, "_", Sys.Date(), ".txt"),
          content = function(file) {
            req(parsed_samples())

            sample_list <- parsed_samples()
            raw_df <- sample_list[[sample_name]]$raw
            ca <- cohort_analysis()
            gs <- ca$gs_long %>% filter(sample == sample_name)

            protected_taxa <- input[[paste0("decontam_select_", sid)]] %||% character(0)
                keep_desc <- isTRUE(input[[paste0("decontam_keep_desc_", sid)]] %||% FALSE)

            cleaned <- build_cleaned_report_for_sample(
              raw_df = raw_df,
              gs = gs,
              protected_taxa = protected_taxa,
              keep_descendants = keep_desc
            )
            write.table(cleaned, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          }
        )

        output[[paste0("bioai_ui_", sid)]] <- renderUI({
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"
          gs <- ca$gs_long %>% filter(sample == sample_name)
          bioai_info <- get_bioai_taxa_info(gs, tax_level = tax_level)

          tagList(
            tags$p(class = "small-note", "BioAI uses NCBI Taxonomy and Wikipedia/WikiData summaries. It is informational and should not be solely used for clinical decision."),
            selectizeInput(
              inputId = paste0("bioai_taxa_", sid),
              label = "Taxa to interpret",
              choices = as.character(bioai_info$initial_choices),
              selected = as.character(bioai_info$initial_selected),
              multiple = TRUE,
              options = list(placeholder = "Type to search taxa...", maxOptions = 10000)
            ),
            actionButton(paste0("bioai_go_", sid), "Generate interpretation", icon = icon("magic"))
          )
        })

        outputOptions(output, paste0("bioai_ui_", sid), suspendWhenHidden = FALSE)

        observeEvent(cohort_analysis(), {
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"
          gs <- ca$gs_long %>% filter(sample == sample_name)
          bioai_info <- get_bioai_taxa_info(gs, tax_level = tax_level)
          current_sel <- input[[paste0("bioai_taxa_", sid)]] %||% character(0)
          selected_taxa <- if (length(current_sel) > 0) current_sel else bioai_info$initial_selected
          if (length(selected_taxa) == 0 && length(bioai_info$candidates) > 0) {
            selected_taxa <- head(bioai_info$candidates, 3)
          }

          session$onFlushed(function() {
            updateSelectizeInput(
              session,
              inputId = paste0("bioai_taxa_", sid),
              choices = as.character(bioai_info$candidates),
              selected = as.character(selected_taxa),
              server = TRUE
            )
          }, once = TRUE)
        }, ignoreInit = FALSE)

        observeEvent(input[[paste0("bioai_go_", sid)]], {
          ca <- cohort_analysis()
          tax_level <- params()$tax_level %||% "S"

          gs <- ca$gs_long %>% filter(sample == sample_name)
          md <- ca$meta %>% filter(sample == sample_name)
          taxa_sel <- input[[paste0("bioai_taxa_", sid)]] %||% character(0)

          if (length(taxa_sel) == 0) {
            bioai_text(HTML("<b>Please select one or more taxa.</b>"))
            return()
          }

          top_df <- gs %>%
            filter(rank == tax_level, name_clean %in% taxa_sel) %>%
            arrange(desc(rpmm)) %>%
            select(rank, name_clean, reads_clade, rpmm, fp_score, call, call_reason, log2FC_vs_reference_median) %>%
            distinct()

          sample_type <- as.character(md$sample_type[1])
          if (!text_is_present(sample_type)) sample_type <- "Clinical / sterile"
          bioai_text(safe_bioai_web_explain(sample_name, sample_type, top_df))
        })

        output[[paste0("bioai_out_", sid)]] <- renderUI({
          bioai_text()
        })

        outputOptions(output, paste0("bioai_out_", sid), suspendWhenHidden = FALSE)
      })
    }
  }, ignoreInit = FALSE)

  observeEvent(input$return_home_btn, {
    parsed_samples(NULL)
    meta_df(NULL)
    last_error(NULL)
    cleaned_parsed_samples(NULL)
    cleaned_last_error(NULL)
    updateProgressBar(session, "analysis_progress", value = 0)
    landing_page(TRUE)
  })
}

shinyApp(ui = ui, server = server)