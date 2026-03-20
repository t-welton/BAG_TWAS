# libraries

library(tidyverse)
library(biomaRt)
library(colorRamps)
library(colorspace)
library(ggmanh)
library(gtexr)
library(corrplot)
library(coloc)
library(data.table)
library(arrow)
library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)
library(grid)
library(WGCNA)

# setup paths

dir_root        <- file.path([[PLACEHOLDER]])
dir_results     <- file.path(dir_root,"results")
dir_twas        <- file.path(dir_results,"twas")
dir_gwas        <- file.path(dir_results,"harmonised_imputed_merged_gwas")

datasets <- c(
  "wen_davatzikos_combined",
  "wen_davatzikos_wm",
  "wen_davatzikos_gm",
  "wen_davatzikos_fc",
  "jawinski_markett_gm",
  "jawinski_markett_gwm",
  "jawinski_markett_wm",
  "kaufmann_westlye",
  "kim_lee",
  "leonardsen_wang",
  "pineiro_fjell",
  "yi_huang",
  "zhao_zhao"
)
model_type_expression_splicing <- c("eqtl","sqtl")
model_type_multi_individual <- c("smultixcan","spredixcan")
model_type_brain_all <- c("brain_tissues","all_tissues")
gtex_brain_tissues <- c(
  "Amygdala",
  "Anterior_cingulate_cortex_BA24",
  "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere",
  "Cerebellum",
  "Cortex",
  "Frontal_Cortex_BA9",
  "Hippocampus",
  "Hypothalamus",
  "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia",
  "Spinal_cord_cervical_c-1",
  "Substantia_nigra"
)
gtex_all_tissues <- c(
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum",
  "Adrenal_Gland",
  "Artery_Aorta",
  "Artery_Coronary",
  "Artery_Tibial",
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra",
  "Breast_Mammary_Tissue",
  "Cells_Cultured_fibroblasts",
  "Cells_EBV-transformed_lymphocytes",
  "Colon_Sigmoid",
  "Colon_Transverse",
  "Esophagus_Gastroesophageal_Junction",
  "Esophagus_Mucosa",
  "Esophagus_Muscularis",
  "Heart_Atrial_Appendage",
  "Heart_Left_Ventricle",
  "Kidney_Cortex",
  "Liver",
  "Lung",
  "Minor_Salivary_Gland",
  "Muscle_Skeletal",
  "Nerve_Tibial",
  "Ovary",
  "Pancreas",
  "Pituitary",
  "Prostate",
  "Skin_Not_Sun_Exposed_Suprapubic",
  "Skin_Sun_Exposed_Lower_leg",
  "Small_Intestine_Terminal_Ileum",
  "Spleen",
  "Stomach",
  "Testis",
  "Thyroid",
  "Uterus",
  "Vagina",
  "Whole_Blood"
)

# load twas results data
files <- list.files(dir_twas, pattern = "\\.csv$", full.names = TRUE)

n_files <- length(files)
cat("Found", n_files, "CSV files in", dir_twas, "\n")

df_twas_results <- list()
problems <- character(0)

for (i in seq_along(files)) {
  
  f <- files[i]
  fname <- basename(f)
  
  cat(sprintf("[%d/%d] Loading %s ... ", i, n_files, fname))
  
  stem <- sub("\\.csv$", "", fname)
  
  m <- regexec("^(eqtl|sqtl)_(smultixcan|spredixcan)_(.+)$", stem)
  parts <- regmatches(stem, m)[[1]]
  if (length(parts) == 0) {
    cat("SKIPPED (bad prefix)\n")
    problems <- c(problems, paste0("Bad prefix: ", fname))
    next
  }
  
  qtl   <- parts[2]
  model <- parts[3]
  rest  <- parts[4]
  
  dataset <- NA_character_
  scope   <- NA_character_
  tissue  <- NA_character_
  key4    <- NA_character_
  
  if (model == "smultixcan") {
    
    m2 <- regexec(paste0("^(", paste(datasets, collapse="|"), ")_(brain_tissues|all_tissues)$"), rest)
    p2 <- regmatches(rest, m2)[[1]]
    if (length(p2) == 0) {
      cat("SKIPPED (smultixcan parse fail)\n")
      problems <- c(problems, paste0("smultixcan parse failed: ", fname))
      next
    }
    
    dataset <- p2[2]
    scope   <- p2[3]
    key4    <- scope
    
  } else if (model == "spredixcan") {
    
    m2 <- regexec(paste0("^(", paste(datasets, collapse="|"), ")_mashr_(.+)$"), rest)
    p2 <- regmatches(rest, m2)[[1]]
    if (length(p2) == 0) {
      cat("SKIPPED (spredixcan parse fail)\n")
      problems <- c(problems, paste0("spredixcan parse failed: ", fname))
      next
    }
    
    dataset <- p2[2]
    tissue  <- p2[3]
    key4    <- tissue
    
    if (!(tissue %in% gtex_all_tissues)) {
      cat("SKIPPED (unknown tissue)\n")
      problems <- c(problems, paste0("Unknown tissue: ", tissue, " in ", fname))
      next
    }
    
  } else {
    cat("SKIPPED (unknown model)\n")
    problems <- c(problems, paste0("Unknown model: ", fname))
    next
  }
 
  sep_guess <- if (model == "smultixcan") "	" else ","
  
  df <- 
    data.table::fread(
      file = f,
      sep = sep_guess,
      data.table = FALSE,
      showProgress = FALSE,
      fill = TRUE,
      quote = "\"",
      check.names = FALSE
    )
  
  df$._dataset <- dataset
  df$._qtl     <- qtl
  df$._model   <- model
  df$._scope   <- scope
  df$._tissue  <- tissue
  df$._file    <- fname
  
  df_twas_results[[dataset]][[qtl]][[model]][[key4]] <- df
  
  cat("OK\n")
}

# usage example
# df_twas_results[["kim_lee"]][["eqtl"]][["smultixcan"]][["brain_tissues"]]

# annotate additional columns

n_slices <- 0

for (dataset in names(df_twas_results)) {
  for (qtl in names(df_twas_results[[dataset]])) {
    for (model in names(df_twas_results[[dataset]][[qtl]])) {
      for (key4 in names(df_twas_results[[dataset]][[qtl]][[model]])) {
        
        df <- df_twas_results[[dataset]][[qtl]][[model]][[key4]]
        if (is.null(df) || !is.data.frame(df)) next
        
        n_slices <- n_slices + 1
        cat(sprintf("[%d] %s | %s | %s | %s ... ", n_slices, dataset, qtl, model, key4))
        
        if (!("pvalue" %in% names(df))) {
          cat("SKIP (no pvalue column)
")
          next
        }
        df$pvalue <- as.numeric(df$pvalue)
        
        n0 <- nrow(df)
        df <- df %>% dplyr::filter(!is.na(pvalue))
        n1 <- nrow(df)
        
        if ("gene" %in% names(df)) {
          df$gene_noversion <- sub("\\.\\d+$", "", df$gene)
        } else {
          df$gene_noversion <- NA_character_
        }
        
        df$pvalue_fdr <- p.adjust(df$pvalue, method = "BH")
        df$fdr_significant <- df$pvalue_fdr < 0.05
        
        n_tests <- nrow(df)
        if (is.na(n_tests) || n_tests == 0) {
          df$bonf_significant <- FALSE
        } else {
          df$bonf_significant <- df$pvalue < (0.05 / n_tests)
        }
        
       # write back
        df_twas_results[[dataset]][[qtl]][[model]][[key4]] <- df
      }
    }
  }
}

# annotate positions eqtl

required_cols <- c("chromosome_name", "start_position", "end_position", "annotation_genome_version")
missing_genes_all <- character(0)

slice_counter <- 0
for (dataset in names(df_twas_results)) {
  if (!"eqtl" %in% names(df_twas_results[[dataset]])) next
  for (model in names(df_twas_results[[dataset]][["eqtl"]])) {
    for (key4 in names(df_twas_results[[dataset]][["eqtl"]][[model]])) {
      
      df <- df_twas_results[[dataset]][["eqtl"]][[model]][[key4]]
      if (is.null(df) || !is.data.frame(df)) next
      
      slice_counter <- slice_counter + 1
      cat(sprintf("[%d] %s | eqtl | %s | %s ... ", slice_counter, dataset, model, key4))
      
      if (!"gene_noversion" %in% names(df)) {
        cat("SKIP (no gene_noversion)
")
        next
      }
      
      for (cc in required_cols) {
        if (!cc %in% names(df)) df[[cc]] <- NA
      }
      
      df$chromosome_name <- as.character(df$chromosome_name)
      df$start_position <- suppressWarnings(as.integer(df$start_position))
      df$end_position <- suppressWarnings(as.integer(df$end_position))
      df$annotation_genome_version <- as.character(df$annotation_genome_version)
      
      miss_idx <- which(is.na(df$chromosome_name) & !is.na(df$gene_noversion) & df$gene_noversion != "")
      if (length(miss_idx) > 0) {
        missing_genes_all <- unique(c(missing_genes_all, df$gene_noversion[miss_idx]))
      }
      
      df_twas_results[[dataset]][["eqtl"]][[model]][[key4]] <- df
      cat(sprintf("OK (missing genes in slice: %d)
", length(unique(df$gene_noversion[miss_idx]))))
    }
  }
}

missing_genes_all <- unique(missing_genes_all)
cat("
Total unique eQTL genes needing position annotation: ", length(missing_genes_all), "
", sep = "")

# resolve positions from biomaRt
annotation_map <- data.frame(
  gene_noversion = character(0),
  chromosome_name = character(0),
  start_position = integer(0),
  end_position = integer(0),
  annotation_genome_version = character(0),
  stringsAsFactors = FALSE
)

if (length(missing_genes_all) > 0) {
  
  mart_versions <- listEnsemblArchives()$version
  mart_version_index <- 1
  
  repeat {
    
    if (length(missing_genes_all) == 0) break
    if (mart_version_index > (length(mart_versions) - 2)) break
    
    mart_version <- mart_versions[mart_version_index]
    cat("
Trying Ensembl mart version: ", mart_version,
        " | genes remaining: ", length(missing_genes_all), "
", sep = "")
    
    mart <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = mart_version),
      error = function(e) NULL
    )
    
    if (is.null(mart)) {
      cat("  -> failed to connect; skipping this version
")
      mart_version_index <- mart_version_index + 1
      next
    }
    
   # split missing genes by type
    is_ensg <- grepl("^ENSG", missing_genes_all)
    genes_ensg <- missing_genes_all[is_ensg]
    genes_sym <- missing_genes_all[!is_ensg]
    
   # query Ensembl IDs
    res1 <- NULL
    if (length(genes_ensg) > 0) {
      res1 <- tryCatch(
        getBM(
          attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
          filters = "ensembl_gene_id",
          values = genes_ensg,
          mart = mart
        ),
        error = function(e) NULL
      )
      if (!is.null(res1) && nrow(res1) > 0) {
        names(res1)[names(res1) == "ensembl_gene_id"] <- "gene_noversion"
      }
    }
    
    # query gene symbols
    res2 <- NULL
    if (length(genes_sym) > 0) {
      res2 <- tryCatch(
        getBM(
          attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
          filters = "external_gene_name",
          values = genes_sym,
          mart = mart
        ),
        error = function(e) NULL
      )
      if (!is.null(res2) && nrow(res2) > 0) {
        names(res2)[names(res2) == "external_gene_name"] <- "gene_noversion"
      }
    }
    
    # combine
    mart_result <- NULL
    if (!is.null(res1) && nrow(res1) > 0 && !is.null(res2) && nrow(res2) > 0) {
      mart_result <- dplyr::bind_rows(res1, res2)
    } else if (!is.null(res1) && nrow(res1) > 0) {
      mart_result <- res1
    } else if (!is.null(res2) && nrow(res2) > 0) {
      mart_result <- res2
    }
    
    if (!is.null(mart_result) && nrow(mart_result) > 0) {
      
      mart_result$annotation_genome_version <- mart_version
      mart_result <- mart_result %>%
        mutate(
          chromosome_name = as.character(chromosome_name),
          start_position = as.integer(start_position),
          end_position = as.integer(end_position),
          annotation_genome_version = as.character(annotation_genome_version)
        )
      
      mart_result <- mart_result %>%
        arrange(gene_noversion) %>%
        group_by(gene_noversion) %>%
        slice(1) %>%
        ungroup()
      
      new_genes <- setdiff(mart_result$gene_noversion, annotation_map$gene_noversion)
      if (length(new_genes) > 0) {
        annotation_map <- dplyr::bind_rows(
          annotation_map,
          mart_result %>% filter(gene_noversion %in% new_genes)
        )
      }
      missing_genes_all <- setdiff(missing_genes_all, mart_result$gene_noversion)
      cat("  -> mapped this round: ", nrow(mart_result),
          " | remaining: ", length(missing_genes_all), "
", sep = "")
    } else {
      cat("  -> no matches found in this mart version
")
    }
    
    mart_version_index <- mart_version_index + 1
  }
}

write_counter <- 0
for (dataset in names(df_twas_results)) {
  if (!"eqtl" %in% names(df_twas_results[[dataset]])) next
  for (model in names(df_twas_results[[dataset]][["eqtl"]])) {
    for (key4 in names(df_twas_results[[dataset]][["eqtl"]][[model]])) {
      
      df <- df_twas_results[[dataset]][["eqtl"]][[model]][[key4]]
      if (is.null(df) || !is.data.frame(df)) next
      if (!"gene_noversion" %in% names(df)) next
      
      write_counter <- write_counter + 1
      cat(sprintf("[%d] Writing positions: %s | eqtl | %s | %s ... ", write_counter, dataset, model, key4))
      
      # ensure columns exist
      for (cc in required_cols) {
        if (!cc %in% names(df)) df[[cc]] <- NA
      }
      
      idx <- match(df$gene_noversion, annotation_map$gene_noversion)
      can_fill <- !is.na(idx) & is.na(df$chromosome_name)
      
      df$chromosome_name[can_fill] <- annotation_map$chromosome_name[idx[can_fill]]
      df$start_position[can_fill] <- annotation_map$start_position[idx[can_fill]]
      df$end_position[can_fill] <- annotation_map$end_position[idx[can_fill]]
      df$annotation_genome_version[can_fill] <- annotation_map$annotation_genome_version[idx[can_fill]]
      
      # normalize types again
      df$chromosome_name <- as.character(df$chromosome_name)
      df$start_position <- suppressWarnings(as.integer(df$start_position))
      df$end_position <- suppressWarnings(as.integer(df$end_position))
      df$annotation_genome_version <- as.character(df$annotation_genome_version)
      
      df_twas_results[[dataset]][["eqtl"]][[model]][[key4]] <- df
      
      cat(sprintf("OK (filled rows: %d)
", sum(can_fill)))
    }
  }
}

# annotate positions and genes sqtl

slice_counter <- 0
for (dataset in names(df_twas_results)) {
  if (!"sqtl" %in% names(df_twas_results[[dataset]])) next
  for (model in names(df_twas_results[[dataset]][["sqtl"]])) {
    for (key4 in names(df_twas_results[[dataset]][["sqtl"]][[model]])) {
      
      df <- df_twas_results[[dataset]][["sqtl"]][[model]][[key4]]
      if (is.null(df) || !is.data.frame(df)) next
      
      slice_counter <- slice_counter + 1
      cat(sprintf("[%d] %s | sqtl | %s | %s ... ", slice_counter, dataset, model, key4))

      cand_cols <- c("intron_id", "phenotype_id", "gene")
      id_col <- cand_cols[cand_cols %in% names(df)][1]
      if (is.na(id_col)) {
        cat("SKIP (no intron_id/phenotype_id/gene column)
")
        df_twas_results[[dataset]][["sqtl"]][[model]][[key4]] <- df
        next
      }
      
      to_fill <- which(
        is.na(df$intron_chromosome_name) &
          !is.na(df[[id_col]]) &
          df[[id_col]] != ""
      )
      
      parsed_now <- 0L
      if (length(to_fill) > 0) {
        x <- as.character(df[[id_col]][to_fill])
        
        parts <- strsplit(x, "_", fixed = TRUE)
        last3 <- lapply(parts, function(v) if (length(v) >= 3) tail(v, 3) else c(NA, NA, NA))
        mat <- do.call(rbind, last3)
        
        chr <- gsub("^chr", "", mat[, 1])
        st  <- suppressWarnings(as.integer(mat[, 2]))
        en  <- suppressWarnings(as.integer(mat[, 3]))
        
        ok <- !is.na(chr) & chr != "" & !is.na(st) & !is.na(en)
        if (any(ok)) {
          df$intron_chromosome_name[to_fill[ok]] <- as.character(chr[ok])
          df$intron_start_position[to_fill[ok]] <- as.integer(st[ok])
          df$intron_end_position[to_fill[ok]] <- as.integer(en[ok])
          parsed_now <- sum(ok)
        }
      }
      
      df_twas_results[[dataset]][["sqtl"]][[model]][[key4]] <- df
      cat(sprintf("OK (id_col=%s; parsed %d/%d)
", id_col, parsed_now, length(to_fill)))
    }
  }
}

# count num of predicted assocs per twas
counts_list <- list()
row_i <- 0L

for (dataset in names(df_twas_results)) {
  for (qtl in names(df_twas_results[[dataset]])) {
    for (model in names(df_twas_results[[dataset]][[qtl]])) {
      for (key4 in names(df_twas_results[[dataset]][[qtl]][[model]])) {
        
        df <- df_twas_results[[dataset]][[qtl]][[model]][[key4]]
        if (is.null(df) || !is.data.frame(df)) next
        
        row_i <- row_i + 1L
        
        n_assoc <- nrow(df)
        n_fdr  <- if ("fdr_significant" %in% names(df)) sum(df$fdr_significant %in% TRUE, na.rm = TRUE) else NA_integer_
        n_bonf <- if ("bonf_significant" %in% names(df)) sum(df$bonf_significant %in% TRUE, na.rm = TRUE) else NA_integer_
        
        counts_list[[row_i]] <- data.frame(
          dataset = dataset,
          qtl = qtl,
          model = model,
          key4 = key4,
          n_assoc = n_assoc,
          n_fdr_significant = n_fdr,
          n_bonf_significant = n_bonf,
          stringsAsFactors = FALSE
        )
        
        if ("._file" %in% names(df)) {
          counts_list[[row_i]]$file <- as.character(df$._file[1])
        } else {
          counts_list[[row_i]]$file <- NA_character_
        }
        
        cat(sprintf("[%d] %s | %s | %s | %s ... n=%d\n", row_i, dataset, qtl, model, key4, n_assoc))
      }
    }
  }
}

counts_df <- dplyr::bind_rows(counts_list) %>%
  arrange(qtl, model, dataset, key4)


# PCA

dir_pca_base <- [[PLACEHOLDER]]
dir_pca_p    <- file.path(dir_pca_base, "pca_p")
dir_pca_z    <- file.path(dir_pca_base, "pca_z")
dir.create(dir_pca_p, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_pca_z, recursive = TRUE, showWarnings = FALSE)

dir_pca_leg_p <- file.path(dir_pca_p, "legends")
dir_pca_leg_z <- file.path(dir_pca_z, "legends")
dir.create(dir_pca_leg_p, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_pca_leg_z, recursive = TRUE, showWarnings = FALSE)

primary_dataset <- "wen_davatzikos_combined"
secondary_datasets <- c("wen_davatzikos_gm", "wen_davatzikos_wm", "wen_davatzikos_fc")
replication_datasets <- setdiff(datasets, c(primary_dataset, secondary_datasets))

# palette 
OKABE_ITO <- c(
  black  = "#000000",
  orange = "#E69F00",
  sky    = "#56B4E9",
  green  = "#009E73",
  yellow = "#F0E442",
  blue   = "#0072B2",
  verm   = "#D55E00",
  purple = "#CC79A7",
  grey   = "#7F7F7F"
)

# theme
theme_nature <- function(base_size = 9, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2, lineheight = 0.95),
      plot.subtitle = element_text(hjust = 0, size = base_size, lineheight = 1.1, margin = margin(b = 4)),
      axis.title = element_text(face = "plain", size = base_size + 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6)
    )
}

# helpers 
is_brain_gtex_label <- function(tissue_label) {
  if (is.na(tissue_label) || tissue_label == "") return(FALSE)
  if (grepl("^Brain_", tissue_label)) return(TRUE)
  if (tissue_label %in% gtex_brain_tissues) return(TRUE)
  FALSE
}

# extractors
extract_id_p <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  if (!("pvalue" %in% names(df))) return(NULL)
  
  id_col <- dplyr::case_when(
    "gene_noversion" %in% names(df) ~ "gene_noversion",
    "gene" %in% names(df) ~ "gene",
    "phenotype_id" %in% names(df) ~ "phenotype_id",
    "intron_id" %in% names(df) ~ "intron_id",
    TRUE ~ NA_character_
  )
  if (is.na(id_col)) return(NULL)
  
  out <- df %>%
    dplyr::transmute(
      id_raw = as.character(.data[[id_col]]),
      pvalue = suppressWarnings(as.numeric(pvalue))
    ) %>%
    dplyr::filter(!is.na(id_raw), id_raw != "", !is.na(pvalue), pvalue > 0)
  
  if (nrow(out) == 0) return(NULL)
  
  out <- out %>%
    dplyr::mutate(id = sub("[.][0-9]+$", "", id_raw)) %>%
    dplyr::select(id, pvalue)
  
  out %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(pvalue = min(pvalue, na.rm = TRUE), .groups = "drop")
}

extract_id_z <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  
  id_col <- dplyr::case_when(
    "gene_noversion" %in% names(df) ~ "gene_noversion",
    "gene" %in% names(df) ~ "gene",
    "phenotype_id" %in% names(df) ~ "phenotype_id",
    "intron_id" %in% names(df) ~ "intron_id",
    TRUE ~ NA_character_
  )
  if (is.na(id_col)) return(NULL)
  
  if ("zscore" %in% names(df)) {
    z_col <- "zscore"
  } else if ("z_mean" %in% names(df)) {
    z_col <- "z_mean"
  } else {
    return(NULL)
  }
  
  out <- df %>%
    dplyr::transmute(
      id_raw = as.character(.data[[id_col]]),
      zscore = suppressWarnings(as.numeric(.data[[z_col]]))
    ) %>%
    dplyr::filter(!is.na(id_raw), id_raw != "", is.finite(zscore))
  
  if (nrow(out) == 0) return(NULL)
  
  out <- out %>%
    dplyr::mutate(id = sub("[.][0-9]+$", "", id_raw)) %>%
    dplyr::select(id, zscore)
  
  out %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(zscore = zscore[which.max(abs(zscore))][1], .groups = "drop")
}

build_pca_matrix <- function(slice_list_named, mode = c("p", "z")) {
  mode <- match.arg(mode)
  stopifnot(is.list(slice_list_named), length(slice_list_named) >= 2)
  
  extractor <- if (mode == "p") extract_id_p else extract_id_z
  value_nm  <- if (mode == "p") "pvalue" else "zscore"
  
  cleaned <- purrr::imap(slice_list_named, function(df, nm) {
    x <- extractor(df)
    if (is.null(x)) return(NULL)
    x <- x %>% dplyr::rename(!!nm := .data[[value_nm]])
    x
  })
  cleaned <- cleaned[!vapply(cleaned, is.null, logical(1))]
  if (length(cleaned) < 2) stop("Too few valid slices after cleaning.")
  
  common <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "id"), cleaned)
  if (nrow(common) == 0) stop("No common ids across slices for PCA.")
  
  mat <- common %>% dplyr::select(-id) %>% as.matrix()
  if (mode == "p") {
    mat <- -log10(mat)
  }
  
  mat <- t(mat) # rows=points, cols=features
  list(mat = mat, n_features = ncol(mat))
}

run_pca_scores <- function(mat, meta_df = NULL) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)
  
  # keep finite columns only
  is_finite_col <- apply(mat, 2, function(v) all(is.finite(v)))
  mat <- mat[, is_finite_col, drop = FALSE]
  
  # remove zero-variance columns
  sds <- apply(mat, 2, stats::sd)
  keep <- !is.na(sds) & sds > 0
  n_dropped <- sum(!keep)
  if (n_dropped > 0) message("Dropped ", n_dropped, " constant/invalid feature columns before PCA.")
  mat <- mat[, keep, drop = FALSE]
  
  p <- prcomp(mat, scale. = TRUE)
  scores <- as.data.frame(p$x)
  scores$Point <- rownames(scores)
  
  ve <- (p$sdev^2) / sum(p$sdev^2)
  scores$PC1_var <- ve[1]
  scores$PC2_var <- ve[2]
  
  if (!is.null(meta_df)) scores <- dplyr::left_join(scores, meta_df, by = "Point")
  list(scores = scores, pca = p)
}

plot_pca <- function(scores_df, color_col, title, subtitle = NULL, colors = NULL, show_legend = FALSE) {
  stopifnot(color_col %in% names(scores_df))
  
  pc1_var <- round(unique(scores_df$PC1_var)[1] * 100, 1)
  pc2_var <- round(unique(scores_df$PC2_var)[1] * 100, 1)
  
  title_w <- stringr::str_wrap(title, width = 55)
  subtitle_w <- if (!is.null(subtitle)) stringr::str_wrap(subtitle, width = 70) else NULL
  
  g <- ggplot(scores_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = .data[[color_col]]), size = 2.3, alpha = 0.9) +
    labs(
      title = title_w,
      subtitle = subtitle_w,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)"),
      color = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    scale_y_continuous(expand = expansion(mult = 0.06)) +
    coord_cartesian(clip = "on") +
    theme_nature() +
    theme(legend.position = if (isTRUE(show_legend)) "right" else "none")
  
  if (!is.null(colors)) g <- g + scale_color_manual(values = colors)
  g
}

extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

save_legend_tiff <- function(legend_grob, filename, width_mm = 70, height_mm = 40) {
  if (is.null(legend_grob)) {
    message("No legend found for: ", filename)
    return(invisible(NULL))
  }
  grDevices::tiff(filename, units = "mm", res = 600, compression = "lzw", width = width_mm, height = height_mm)
  grid::grid.newpage()
  grid::grid.draw(legend_grob)
  grDevices::dev.off()
  invisible(filename)
}

get_spredixcan_slices <- function(dataset, qtl, tissues = NULL) {
  x <- df_twas_results[[dataset]][[qtl]][["spredixcan"]]
  if (is.null(x)) return(list())
  if (!is.null(tissues)) x <- x[intersect(names(x), tissues)]
  x
}

get_smultixcan_slices <- function(datasets_vec, qtl, scope, verbose = TRUE) {
  qtl <- trimws(qtl)
  scope <- trimws(scope)
  datasets_vec <- trimws(datasets_vec)
  
  out <- list()
  missing <- character(0)
  reasons <- character(0)
  
  for (ds in datasets_vec) {
    if (is.null(df_twas_results[[ds]])) {
      missing <- c(missing, ds)
      reasons <- c(reasons, "no dataset tier")
      next
    }
    if (is.null(df_twas_results[[ds]][[qtl]])) {
      missing <- c(missing, ds)
      reasons <- c(reasons, paste0("no qtl tier: ", qtl))
      next
    }
    if (is.null(df_twas_results[[ds]][[qtl]][["smultixcan"]])) {
      missing <- c(missing, ds)
      reasons <- c(reasons, "no smultixcan tier")
      next
    }
    if (is.null(df_twas_results[[ds]][[qtl]][["smultixcan"]][[scope]])) {
      missing <- c(missing, ds)
      reasons <- c(reasons, paste0("no scope tier: ", scope))
      next
    }
    
    x <- df_twas_results[[ds]][[qtl]][["smultixcan"]][[scope]]
    
    ok_df <- is.data.frame(x)
    ok_p  <- ok_df && ("pvalue" %in% names(x))
    ok_n  <- ok_df && nrow(x) > 0
    
    if (ok_df && ok_p && ok_n) {
      out[[ds]] <- x
    } else {
      missing <- c(missing, ds)
      if (!ok_df) {
        reasons <- c(reasons, paste0("not data.frame (class=", paste(class(x), collapse = ","), ")"))
      } else if (!ok_p) {
        reasons <- c(reasons, "missing pvalue column")
      } else if (!ok_n) {
        reasons <- c(reasons, "0 rows")
      } else {
        reasons <- c(reasons, "unknown")
      }
    }
  }
  out
}

# Brain vs Other palette
colors_brain_other <- c(
  Brain = unname(OKABE_ITO["verm"]),
  Other = unname(OKABE_ITO["grey"])
)

# dataset palette 
make_dataset_colors <- function(datasets_vec) {
  cols <- colorspace::qualitative_hcl(length(datasets_vec), palette = "Dark 3")
  setNames(cols, datasets_vec)
}

make_individual_tissue_pca <- function(dataset, qtl, tissues = gtex_all_tissues,
                                       brain_only = FALSE,
                                       mode = c("p", "z"),
                                       dir_out,
                                       file_stub = NULL,
                                       title_stub = NULL) {
  mode <- match.arg(mode)
  
  slices <- get_spredixcan_slices(dataset, qtl, tissues = tissues)
  if (length(slices) < 3) {
    return(NULL)
  }
  
  if (brain_only) {
    keep <- names(slices)[vapply(names(slices), is_brain_gtex_label, logical(1))]
    slices <- slices[keep]
    if (length(slices) < 3) {
      return(NULL)
    }
  }
  
  m <- build_pca_matrix(slices, mode = mode)
  
  meta <- data.frame(
    Point = rownames(m$mat),
    TissueGroup = ifelse(vapply(rownames(m$mat), is_brain_gtex_label, logical(1)), "Brain", "Other"),
    stringsAsFactors = FALSE
  )
  
  r <- run_pca_scores(m$mat, meta)
  
  title <- if (!is.null(title_stub)) title_stub else paste0("PCA: ", dataset)
  metric <- if (mode == "p") "-log10(P)" else "signed Z"
  subtitle <- paste0(qtl, " | spredixcan | ", metric, " | n_features=", m$n_features, " | n_points=", nrow(r$scores))
  
  g <- plot_pca(r$scores, "TissueGroup", title, subtitle, colors = colors_brain_other, show_legend = FALSE)
  
  if (!is.null(file_stub)) {
    fn <- file.path(dir_out, paste0(file_stub, "_", qtl, ".tiff"))
    ggsave(fn, g,
           device = "tiff",
           dpi = 600,
           compression = "lzw",
           width = 85,
           height = 70,
           units = "mm")
  }
  
  g
}

make_multidataset_pca <- function(datasets_vec, qtl, scope,
                                  mode = c("p", "z"),
                                  dir_out,
                                  file_stub = NULL,
                                  title_stub = NULL,
                                  dataset_colors = NULL) {
  mode <- match.arg(mode)
  
  slices <- get_smultixcan_slices(datasets_vec, qtl, scope, verbose = TRUE)
  if (length(slices) < 3) {
    return(NULL)
  }
  
  slices <- slices[intersect(datasets_vec, names(slices))]
  
  m <- build_pca_matrix(slices, mode = mode)
  
  meta <- data.frame(
    Point = rownames(m$mat),
    Dataset = rownames(m$mat),
    stringsAsFactors = FALSE
  )
  
  r <- run_pca_scores(m$mat, meta)
  
  title <- if (!is.null(title_stub)) title_stub else paste0("PCA: smultixcan ", scope)
  metric <- if (mode == "p") "-log10(P)" else "signed Z"
  subtitle <- paste0(qtl, " | smultixcan | ", metric, " | n_features=", m$n_features, " | n_points=", nrow(r$scores))
  
  if (is.null(dataset_colors)) dataset_colors <- make_dataset_colors(unique(r$scores$Dataset))
  
  g <- plot_pca(r$scores, "Dataset", title, subtitle, colors = dataset_colors, show_legend = FALSE)
  
  if (!is.null(file_stub)) {
    fn <- file.path(dir_out, paste0(file_stub, "_", qtl, ".tiff"))
    ggsave(fn, g,
           device = "tiff",
           dpi = 600,
           compression = "lzw",
           width = 85,
           height = 70,
           units = "mm")
  }
  
  g
}

# save legends 
legend_demo_ind <- data.frame(
  PC1 = c(0, 1), PC2 = c(0, 1),
  TissueGroup = c("Brain", "Other"),
  PC1_var = 0.5, PC2_var = 0.5
)

p_leg_ind <- plot_pca(legend_demo_ind, "TissueGroup",
                      title = "Legend", subtitle = NULL,
                      colors = colors_brain_other,
                      show_legend = TRUE) +
  theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

save_legend_tiff(extract_legend_grob(p_leg_ind), file.path(dir_pca_leg_p, "legend_individual_tissues_brain_vs_other.tiff"), width_mm = 70, height_mm = 35)
save_legend_tiff(extract_legend_grob(p_leg_ind), file.path(dir_pca_leg_z, "legend_individual_tissues_brain_vs_other.tiff"), width_mm = 70, height_mm = 35)

datasets_ps <- c(primary_dataset, secondary_datasets)
cols_ps <- make_dataset_colors(datasets_ps)
legend_demo_ps <- data.frame(
  PC1 = seq_along(datasets_ps), PC2 = seq_along(datasets_ps),
  Dataset = datasets_ps,
  PC1_var = 0.5, PC2_var = 0.5
)

p_leg_ps <- plot_pca(legend_demo_ps, "Dataset",
                     title = "Legend", subtitle = NULL,
                     colors = cols_ps,
                     show_legend = TRUE) +
  theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

save_legend_tiff(extract_legend_grob(p_leg_ps), file.path(dir_pca_leg_p, "legend_datasets_primary_secondary.tiff"), width_mm = 85, height_mm = 70)
save_legend_tiff(extract_legend_grob(p_leg_ps), file.path(dir_pca_leg_z, "legend_datasets_primary_secondary.tiff"), width_mm = 85, height_mm = 70)

cols_all <- make_dataset_colors(datasets)
legend_demo_all <- data.frame(
  PC1 = seq_along(datasets), PC2 = seq_along(datasets),
  Dataset = datasets,
  PC1_var = 0.5, PC2_var = 0.5
)

p_leg_all <- plot_pca(legend_demo_all, "Dataset",
                      title = "Legend", subtitle = NULL,
                      colors = cols_all,
                      show_legend = TRUE) +
  theme(plot.title = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

save_legend_tiff(extract_legend_grob(p_leg_all), file.path(dir_pca_leg_p, "legend_datasets_all.tiff"), width_mm = 85, height_mm = 120)
save_legend_tiff(extract_legend_grob(p_leg_all), file.path(dir_pca_leg_z, "legend_datasets_all.tiff"), width_mm = 85, height_mm = 120)


pca_plots_p <- list()

for (qtl in c("eqtl", "sqtl")) {
  
  cat("[PCA-P 01] Primary dataset – individual tissues (all) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("01_primary_individual_all_", qtl)]] <- make_individual_tissue_pca(
    dataset = primary_dataset,
    qtl = qtl,
    tissues = gtex_all_tissues,
    brain_only = FALSE,
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "01_primary_individual_all_tissues",
    title_stub = "Primary (wen_davatzikos_combined): individual tissues"
  )
  
  for (ds in secondary_datasets) {
    cat("[PCA-P 02–04] Secondary dataset – ", ds, " | individual tissues (all) | ", qtl, "
", sep = "")
    pca_plots_p[[paste0("02to04_", ds, "_all_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = FALSE,
      mode = "p",
      dir_out = dir_pca_p,
      file_stub = paste0("02_secondary_", ds, "_individual_all_tissues"),
      title_stub = paste0("Secondary (", ds, "): individual tissues")
    )
  }
  
  cat("[PCA-P 05] Primary dataset – individual tissues (brain only) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("05_primary_individual_brain_only_", qtl)]] <- make_individual_tissue_pca(
    dataset = primary_dataset,
    qtl = qtl,
    tissues = gtex_all_tissues,
    brain_only = TRUE,
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "05_primary_individual_brain_only",
    title_stub = "Primary: individual tissues (brain only)"
  )
  
  for (ds in secondary_datasets) {
    cat("[PCA-P 06] Secondary dataset – ", ds, " | individual tissues (brain only) | ", qtl, "
", sep = "")
    pca_plots_p[[paste0("06_", ds, "_brain_only_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = TRUE,
      mode = "p",
      dir_out = dir_pca_p,
      file_stub = paste0("06_secondary_", ds, "_individual_brain_only"),
      title_stub = paste0("Secondary (", ds, "): individual tissues (brain only)")
    )
  }
  
  cat("[PCA-P 07] Primary + secondary – multi-dataset (brain) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("07_primary_secondary_multidataset_brain_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "brain_tissues",
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "07_primary_secondary_multidataset_brain",
    title_stub = "Primary + secondary: smultixcan (brain tissues)",
    dataset_colors = cols_ps
  )
  
  cat("[PCA-P 08] Primary + secondary – multi-dataset (all) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("08_primary_secondary_multidataset_all_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "all_tissues",
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "08_primary_secondary_multidataset_all",
    title_stub = "Primary + secondary: smultixcan (all tissues)",
    dataset_colors = cols_ps
  )
  
  cat("[PCA-P 09] All datasets – multi-dataset (brain) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("09_all_multidataset_brain_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "brain_tissues",
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "09_all_datasets_multidataset_brain",
    title_stub = "All datasets: smultixcan (brain tissues)",
    dataset_colors = cols_all
  )
  
  cat("[PCA-P 10] All datasets – multi-dataset (all) | ", qtl, "
", sep = "")
  pca_plots_p[[paste0("10_all_multidataset_all_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "all_tissues",
    mode = "p",
    dir_out = dir_pca_p,
    file_stub = "10_all_datasets_multidataset_all",
    title_stub = "All datasets: smultixcan (all tissues)",
    dataset_colors = cols_all
  )
  
  for (ds in replication_datasets) {
    cat("[PCA-P 11+] Replication dataset – ", ds, " | individual tissues (all) | ", qtl, "
", sep = "")
    pca_plots_p[[paste0("11_replication_", ds, "_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = FALSE,
      mode = "p",
      dir_out = dir_pca_p,
      file_stub = paste0("11_replication_", ds, "_individual_all_tissues"),
      title_stub = paste0("Replication (", ds, "): individual tissues")
    )
  }
}

pca_plots_z <- list()

for (qtl in c("eqtl", "sqtl")) {

  cat("[PCA-Z 01] Primary dataset – individual tissues (all) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("01_primary_individual_all_", qtl)]] <- make_individual_tissue_pca(
    dataset = primary_dataset,
    qtl = qtl,
    tissues = gtex_all_tissues,
    brain_only = FALSE,
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "01_primary_individual_all_tissues",
    title_stub = "Primary (wen_davatzikos_combined): individual tissues"
  )
  
  for (ds in secondary_datasets) {
    cat("[PCA-Z 02–04] Secondary dataset – ", ds, " | individual tissues (all) | ", qtl, "
", sep = "")
    pca_plots_z[[paste0("02to04_", ds, "_all_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = FALSE,
      mode = "z",
      dir_out = dir_pca_z,
      file_stub = paste0("02_secondary_", ds, "_individual_all_tissues"),
      title_stub = paste0("Secondary (", ds, "): individual tissues")
    )
  }
  
  cat("[PCA-Z 05] Primary dataset – individual tissues (brain only) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("05_primary_individual_brain_only_", qtl)]] <- make_individual_tissue_pca(
    dataset = primary_dataset,
    qtl = qtl,
    tissues = gtex_all_tissues,
    brain_only = TRUE,
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "05_primary_individual_brain_only",
    title_stub = "Primary: individual tissues (brain only)"
  )
  
  for (ds in secondary_datasets) {
    cat("[PCA-Z 06] Secondary dataset – ", ds, " | individual tissues (brain only) | ", qtl, "
", sep = "")
    pca_plots_z[[paste0("06_", ds, "_brain_only_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = TRUE,
      mode = "z",
      dir_out = dir_pca_z,
      file_stub = paste0("06_secondary_", ds, "_individual_brain_only"),
      title_stub = paste0("Secondary (", ds, "): individual tissues (brain only)")
    )
  }
  
  cat("[PCA-Z 07] Primary + secondary – multi-dataset (brain) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("07_primary_secondary_multidataset_brain_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "brain_tissues",
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "07_primary_secondary_multidataset_brain",
    title_stub = "Primary + secondary: smultixcan (brain tissues)",
    dataset_colors = cols_ps
  )
  
  cat("[PCA-Z 08] Primary + secondary – multi-dataset (all) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("08_primary_secondary_multidataset_all_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "all_tissues",
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "08_primary_secondary_multidataset_all",
    title_stub = "Primary + secondary: smultixcan (all tissues)",
    dataset_colors = cols_ps
  )
  
  cat("[PCA-Z 09] All datasets – multi-dataset (brain) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("09_all_multidataset_brain_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "brain_tissues",
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "09_all_datasets_multidataset_brain",
    title_stub = "All datasets: smultixcan (brain tissues)",
    dataset_colors = cols_all
  )
  
  cat("[PCA-Z 10] All datasets – multi-dataset (all) | ", qtl, "
", sep = "")
  pca_plots_z[[paste0("10_all_multidataset_all_", qtl)]] <- make_multidataset_pca(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "all_tissues",
    mode = "z",
    dir_out = dir_pca_z,
    file_stub = "10_all_datasets_multidataset_all",
    title_stub = "All datasets: smultixcan (all tissues)",
    dataset_colors = cols_all
  )
  
  for (ds in replication_datasets) {
    cat("[PCA-Z 11+] Replication dataset – ", ds, " | individual tissues (all) | ", qtl, "
", sep = "")
    pca_plots_z[[paste0("11_replication_", ds, "_", qtl)]] <- make_individual_tissue_pca(
      dataset = ds,
      qtl = qtl,
      tissues = gtex_all_tissues,
      brain_only = FALSE,
      mode = "z",
      dir_out = dir_pca_z,
      file_stub = paste0("11_replication_", ds, "_individual_all_tissues"),
      title_stub = paste0("Replication (", ds, "): individual tissues")
    )
  }
}

# Correlation

# output folders 
dir_fig_base  <- [[PLACEHOLDER]]
dir_corr_base <- file.path(dir_fig_base, "correlation")
dir_corr_eqtl <- file.path(dir_corr_base, "eqtl")
dir_corr_sqtl <- file.path(dir_corr_base, "sqtl")
dir.create(dir_corr_eqtl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_corr_sqtl, recursive = TRUE, showWarnings = FALSE)

datasets_ps <- c(
  "wen_davatzikos_combined",   # primary
  "wen_davatzikos_gm",
  "wen_davatzikos_wm",
  "wen_davatzikos_fc"
)


#  color scale 
col_div_corr <- colorRampPalette(c(unname(OKABE_ITO["blue"]), "white", unname(OKABE_ITO["verm"])))
wrap_title <- function(x, width = 55) stringr::str_wrap(x, width = width)

# extract slices from df_twas_results 
get_spredixcan_slices <- function(dataset, qtl, tissues = NULL) {
  x <- df_twas_results[[dataset]][[qtl]][["spredixcan"]]
  if (is.null(x)) return(list())
  
  x
}

get_smultixcan_slices <- function(datasets_vec, qtl, scope, verbose = TRUE) {
  out <- list()
  missing <- character(0)
  
  for (ds in datasets_vec) {
    x <- df_twas_results[[ds]][[qtl]][["smultixcan"]][[scope]]
    if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
      missing <- c(missing, ds)
      next
    }
    out[[ds]] <- x
  }

  out[intersect(datasets_vec, names(out))]
}

build_z_wide_union <- function(slices_named, qtl, mode = c("individual", "multi")) {
  mode <- match.arg(mode)
  
  if (!is.list(slices_named) || length(slices_named) < 3) {
    message("  Too few points (<3); skip.")
    return(NULL)
  }
  
  id_col <- if (qtl == "eqtl") "gene_noversion" else "gene_name"  # always
  z_col  <- if (mode == "individual") "zscore" else "z_mean"
  
  message("    • Building wide Z matrix (Option A: union_all + pairwise.complete.obs)")
  message("      - qtl=", qtl, "; id_col=", id_col, "; z_col=", z_col)
  
  vecs <- purrr::imap(slices_named, function(df, nm) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      message("      - ", nm, ": empty; skip")
      return(NULL)
    }
    if (!(id_col %in% names(df))) {
      message("      - ", nm, ": missing ", id_col, "; skip")
      return(NULL)
    }
    if (!(z_col %in% names(df))) {
      message("      - ", nm, ": missing ", z_col, "; skip")
      return(NULL)
    }
    
    x <- df %>%
      dplyr::transmute(
        id_raw = as.character(.data[[id_col]]),
        z = suppressWarnings(as.numeric(.data[[z_col]]))
      ) %>%
      dplyr::filter(!is.na(id_raw), id_raw != "", is.finite(z)) %>%
      dplyr::mutate(id = sub("[.][0-9]+$", "", id_raw)) %>%
      dplyr::select(id, z)
    
    if (nrow(x) == 0) {
      message("      - ", nm, ": no finite Z after filtering; skip")
      return(NULL)
    }
    
    x <- x %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(z = z[which.max(abs(z))][1], .groups = "drop")
    
    v <- x$z
    names(v) <- x$id
    
    v
  })
  
  vecs <- vecs[!vapply(vecs, is.null, logical(1))]
  if (length(vecs) < 3) {
    return(NULL)
  }
  
  pts <- names(vecs)
  ids <- unique(unlist(lapply(vecs, names), use.names = FALSE))
  if (length(ids) < 2) {
    return(NULL)
  }
  
  mat <- do.call(rbind, lapply(pts, function(p) {
    v <- vecs[[p]]
    v[ids]  # NA for missing
  }))
  rownames(mat) <- pts
  
  list(mat = mat, n_features = length(ids), n_points = length(pts))
}

#  correlation + plotting 
calc_cor <- function(mat) {
  mat <- as.matrix(mat)

  keep <- apply(mat, 2, function(v) any(is.finite(v)))
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)

  stats::cor(t(mat), use = "pairwise.complete.obs", method = "pearson")
}

save_corrplot_tiff <- function(cmat, filename, title, subtitle = NULL,
                               tl_cex = 0.35, mar = c(0, 15, 2.5, 0)) {
  if (is.null(cmat) || !is.matrix(cmat) || nrow(cmat) < 2) return(NULL)
  
  message("    • Writing TIFF: ", filename)
  grDevices::tiff(filename, units = "mm", res = 600, compression = "lzw", width = 120, height = 120)
  op <- par(mar = mar)
  on.exit({par(op); grDevices::dev.off()}, add = TRUE)
  
  main_title <- wrap_title(title, 55)
  sub_title  <- if (!is.null(subtitle)) wrap_title(subtitle, 80) else NULL
  
  message("cmat dims:")
  cat("cmat dims:", paste(dim(cmat), collapse=" x "), "\n")
  
  message("corrplot")
  corrplot::corrplot(
    cmat,
    method = "color",
    type = "upper",
    order = "original",
    tl.col = "black",
    tl.cex = tl_cex,
    tl.srt = 45,
    tl.pos = "lt",
    addgrid.col = "white",
    col = col_div_corr(201),
    col.lim = c(-1, 1),
    diag = TRUE,
    title = if (is.null(sub_title)) main_title else paste0(main_title, "\n", sub_title)
  )
  
  invisible(filename)
}

#  wrappers 
run_within_dataset_tissue_corr <- function(dataset, qtl, tissues_vec, scope_label, out_dir) {
  
  all_names <- names(df_twas_results[[dataset]][[qtl]][["spredixcan"]])
  
  hits <- intersect(all_names, tissues_vec)
  if (length(hits) == 0 && identical(scope_label, "brain_only")) {
    hits <- intersect(all_names, paste0("Brain_", tissues_vec))
  }
  
  if (length(hits) < 3) {
    return(NULL)
  }
  
  slices <- df_twas_results[[dataset]][[qtl]][["spredixcan"]][hits]
  
  z <- build_z_wide_union(slices, qtl = qtl, mode = "individual")
  if (is.null(z)) return(NULL)
  
  cmat <- calc_cor(z$mat)
  if (is.null(cmat)) return(NULL)
  
  fn <- file.path(out_dir, paste0("within_", dataset, "_", scope_label, "_", qtl, ".tiff"))
  title <- paste0("Correlation across tissues: ", dataset)
  subtitle <- paste0(dataset, " | ", qtl, " | n_features=", z$n_features)
  
  save_corrplot_tiff(cmat, fn, title, subtitle)
}

run_multistudy_corr <- function(datasets_vec, qtl, scope, scope_label, out_dir, file_stub) {
  message("  → Multi-study study×study: ", file_stub, " | ", qtl, " | ", scope_label)
  slices <- get_smultixcan_slices(datasets_vec, qtl, scope, verbose = TRUE)
  if (length(slices) < 3) {
    return(NULL)
  }
  
  z <- build_z_wide_union(slices, qtl = qtl, mode = "multi")
  if (is.null(z)) return(NULL)
  
  cmat <- calc_cor(z$mat)
  if (is.null(cmat)) {
    return(NULL)
  }
  
  fn <- file.path(out_dir, paste0(file_stub, "_", scope_label, "_", qtl, ".tiff"))
  title <- "Correlation across studies"
  subtitle <- paste0("multi", " | ", qtl, " | n_features=", z$n_features)
  
  save_corrplot_tiff(cmat, fn, title, subtitle)
}

for (qtl in c("eqtl", "sqtl")) {
  out_dir <- if (qtl == "eqtl") dir_corr_eqtl else dir_corr_sqtl

  for (ds in datasets) {
    run_within_dataset_tissue_corr(ds, qtl, gtex_all_tissues, scope_label = "gtex49", out_dir = out_dir)
    
    run_within_dataset_tissue_corr(ds, qtl, gtex_brain_tissues, scope_label = "brain_only", out_dir = out_dir)
  }
  
  run_multistudy_corr(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "all_tissues",
    scope_label = "all_tissues",
    out_dir = out_dir,
    file_stub = "multistudy_primary_secondary"
  )
  
  run_multistudy_corr(
    datasets_vec = datasets_ps,
    qtl = qtl,
    scope = "brain_tissues",
    scope_label = "brain_tissues",
    out_dir = out_dir,
    file_stub = "multistudy_primary_secondary"
  )
  
  run_multistudy_corr(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "all_tissues",
    scope_label = "all_tissues",
    out_dir = out_dir,
    file_stub = "multistudy_all_datasets"
  )
  
  run_multistudy_corr(
    datasets_vec = datasets,
    qtl = qtl,
    scope = "brain_tissues",
    scope_label = "brain_tissues",
    out_dir = out_dir,
    file_stub = "multistudy_all_datasets"
  )
  
}

# significance

alpha_fdr_primary      <- 0.05
alpha_nominal_support  <- 0.05   

primary_dataset_fixed  <- "wen_davatzikos_combined"
primary_model_fixed    <- "smultixcan"
primary_scope_fixed    <- "brain_tissues"

id_col_fixed <- "gene_name"
secondary_datasets <- c("wen_davatzikos_fc", "wen_davatzikos_gm", "wen_davatzikos_wm")

# Replication datasets 
all_datasets <- names(df_twas_results)
replication_datasets <- setdiff(all_datasets, c(primary_dataset_fixed, secondary_datasets))

extract_id_pz <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    if (!(id_col_fixed %in% names(df))) return(NULL)
    if (!("pvalue" %in% names(df))) return(NULL)
    
    z_col <- if ("zscore" %in% names(df)) "zscore" else if ("z_mean" %in% names(df)) "z_mean" else NA_character_
    
    out <- df %>%
        dplyr::transmute(
            id = as.character(.data[[id_col_fixed]]),
            pvalue = suppressWarnings(as.numeric(pvalue)),
            zscore = if (!is.na(z_col)) suppressWarnings(as.numeric(.data[[z_col]])) else NA_real_
        ) %>%
        dplyr::filter(!is.na(id), id != "", !is.na(pvalue), pvalue > 0)
    
    if (nrow(out) == 0) return(NULL)
    
    out %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(
            pvalue = min(pvalue, na.rm = TRUE),
            zscore = zscore[which.min(pvalue)][1],
            .groups = "drop"
        )
}

primary_sig_ids    <- list(eqtl = character(0), sqtl = character(0))
primary_sig_tables <- list(eqtl = NULL,           sqtl = NULL)

for (qtl in c("eqtl", "sqtl")) {
    
    df_primary <- df_twas_results[[primary_dataset_fixed]][[qtl]][[primary_model_fixed]][[primary_scope_fixed]]
    
    if (is.null(df_primary)) {
        next
    }
    
    tab <- extract_id_pz(df_primary)
    if (is.null(tab) || nrow(tab) == 0) {
        next
    }
    
    tab <- tab %>%
        dplyr::mutate(
            p_fdr_primary     = p.adjust(pvalue, method = "BH"),
            primary_discovery = p_fdr_primary < alpha_fdr_primary
        )
    
    primary_sig_tables[[qtl]] <- tab
    primary_sig_ids[[qtl]]    <- tab$id[tab$primary_discovery %in% TRUE]
    
    message("Primary discovery (", qtl, "): ",
            sum(tab$primary_discovery %in% TRUE), " / ", nrow(tab),
            " IDs at BH q<", alpha_fdr_primary)
}

support_lookup <- list(eqtl = list(), sqtl = list())

for (qtl in c("eqtl", "sqtl")) {
    
    ids0 <- primary_sig_ids[[qtl]]
    if (length(ids0) == 0) next
    
    tab0 <- primary_sig_tables[[qtl]] %>% dplyr::select(id, z0 = zscore)
    
    for (ds in all_datasets) {
        
        df_sup <- df_twas_results[[ds]][[qtl]][[primary_model_fixed]][[primary_scope_fixed]]
        tab_sup <- extract_id_pz(df_sup)
        if (is.null(tab_sup)) next
        
        tab_sup <- tab_sup %>%
            dplyr::filter(id %in% ids0) %>%
            dplyr::left_join(tab0, by = "id") %>%
            dplyr::mutate(
                support_nominal   = pvalue < alpha_nominal_support,
                support_direction = ifelse(is.finite(zscore) & is.finite(z0), sign(zscore) == sign(z0), NA),
                # if direction is missing, allow nominal-only support
                support_both      = support_nominal & (is.na(support_direction) | support_direction)
            )
        
        support_lookup[[qtl]][[ds]] <- tab_sup
    }
}

slice_i <- 0L

for (dataset in names(df_twas_results)) {
    for (qtl in names(df_twas_results[[dataset]])) {
        for (model in names(df_twas_results[[dataset]][[qtl]])) {
            for (key4 in names(df_twas_results[[dataset]][[qtl]][[model]])) {
                
                df <- df_twas_results[[dataset]][[qtl]][[model]][[key4]]
                if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) next
                
                slice_i <- slice_i + 1L
                
                id_vec <- as.character(df[[id_col_fixed]])
                
                df$primary_discovery     <- id_vec %in% primary_sig_ids[[qtl]]
                df$secondary_support     <- FALSE
                df$replication_support   <- FALSE
                df$support_nominal_p     <- NA_real_
                df$support_direction     <- NA
                
                sup <- support_lookup[[qtl]][[dataset]]
                if (!is.null(sup)) {
                    sup2 <- sup %>% dplyr::select(id, pvalue, support_direction, support_both)
                    m <- match(id_vec, sup2$id)
                    
                    df$support_nominal_p[!is.na(m)] <- sup2$pvalue[m[!is.na(m)]]
                    df$support_direction[!is.na(m)] <- sup2$support_direction[m[!is.na(m)]]
                    
                    if (dataset %in% secondary_datasets) {
                        df$secondary_support[!is.na(m)] <- sup2$support_both[m[!is.na(m)]] %in% TRUE
                    }
                    if (dataset %in% replication_datasets) {
                        df$replication_support[!is.na(m)] <- sup2$support_both[m[!is.na(m)]] %in% TRUE
                    }
                }
                
                df_twas_results[[dataset]][[qtl]][[model]][[key4]] <- df
            }
        }
    }
}

# print sig primary / sec/ rep summary

cat("\n==== SUMMARY: primary discovery and cross-dataset support ====\n")

summary_list <- list()

for (qtl in c("eqtl", "sqtl")) {
    
    cat("\n--", toupper(qtl), "--\n")
    
    # primary discovery table
    df_primary <- df_twas_results[[primary_dataset_fixed]][[qtl]][[primary_model_fixed]][[primary_scope_fixed]]
    
    primary_ids <- unique(df_primary$gene_name[df_primary$primary_discovery %in% TRUE])
    n_primary <- length(primary_ids)
    
    cat("Primary discoveries:", n_primary, "\n")
    
    if (n_primary == 0) next
    
    sec_supported <- logical(n_primary)
    names(sec_supported) <- primary_ids
    
    for (ds in secondary_datasets) {
        df_sec <- df_twas_results[[ds]][[qtl]][[primary_model_fixed]][[primary_scope_fixed]]
        if (is.null(df_sec)) next
        
        m <- match(primary_ids, df_sec$gene_name)
        ok <- !is.na(m) & (df_sec$secondary_support[m] %in% TRUE)
        sec_supported[ok] <- TRUE
    }
    
    n_sec <- sum(sec_supported)
    cat("Secondary support (≥1 modality):", n_sec, "/", n_primary, "\n")
    
    rep_supported <- logical(n_primary)
    names(rep_supported) <- primary_ids
    
    for (ds in replication_datasets) {
        df_rep <- df_twas_results[[ds]][[qtl]][[primary_model_fixed]][[primary_scope_fixed]]
        if (is.null(df_rep)) next
        
        m <- match(primary_ids, df_rep$gene_name)
        ok <- !is.na(m) & (df_rep$replication_support[m] %in% TRUE)
        rep_supported[ok] <- TRUE
    }
    
    n_rep <- sum(rep_supported)
    cat("Replication support (≥1 study):", n_rep, "/", n_primary, "\n")
    
    # combined
    n_both <- sum(sec_supported & rep_supported)
    cat("Secondary + replication support:", n_both, "/", n_primary, "\n")
    
    summary_list[[qtl]] <- data.frame(
        qtl = qtl,
        primary = n_primary,
        secondary_any = n_sec,
        replication_any = n_rep,
        secondary_and_replication = n_both
    )
}

summary_df <- dplyr::bind_rows(summary_list)

cat("\n Tabular summary \n")
print(summary_df)

# two-sided manhattan plots

dir_manh_base <- [[PLACEHOLDER]]
dir_manh_eqtl <- file.path(dir_manh_base, "eqtl")
dir_manh_sqtl <- file.path(dir_manh_base, "sqtl")
dir.create(dir_manh_eqtl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_manh_sqtl, recursive = TRUE, showWarnings = FALSE)

primary_dataset <- "wen_davatzikos_combined"
secondary_datasets <- c("wen_davatzikos_fc", "wen_davatzikos_gm", "wen_davatzikos_wm")
datasets <- names(df_twas_results)
replication_datasets <- setdiff(datasets, c(primary_dataset, secondary_datasets))

theme_nature_manh <- function(base_size = 9, base_family = "") {
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
        theme(
            plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
            plot.subtitle = element_text(hjust = 0, size = base_size, margin = margin(b = 4)),
            axis.title = element_text(size = base_size + 1),
            axis.text.x = element_text(size = base_size - 1, colour = "black", margin = margin(t = 2)),
            axis.text.y = element_text(size = base_size, colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(6, 6, 6, 6),
            legend.position = "none"
        )
}

get_thresholds <- function(pvec) {
    pvec <- as.numeric(pvec)
    pvec <- pvec[is.finite(pvec) & !is.na(pvec) & pvec > 0]
    if (length(pvec) == 0) return(list(bonf_p = NA_real_, fdr_p = NA_real_))
    bonf_p <- 0.05 / length(pvec)
    padj <- p.adjust(pvec, method = "BH")
    fdr_p <- if (any(is.finite(padj) & padj < 0.05)) max(pvec[padj < 0.05], na.rm = TRUE) else NA_real_
    list(bonf_p = bonf_p, fdr_p = fdr_p)
}

prep_twosided_manhattan <- function(df, qtl) {
    if (!("pvalue" %in% names(df))) stop("Missing pvalue column.")
    if (!("z_mean" %in% names(df))) stop("Missing z_mean column (expected for smultixcan).")
    
    if (qtl == "eqtl") {
        CHR_COL <- "chromosome_name"
        POS_COL <- "start_position"
    } else if (qtl == "sqtl") {
        CHR_COL <- "intron_chromosome_name"
        POS_COL <- "intron_start_position"
    } else stop("qtl must be eqtl or sqtl")
    
    if (!(CHR_COL %in% names(df))) stop(paste0("Missing ", CHR_COL))
    if (!(POS_COL %in% names(df))) stop(paste0("Missing ", POS_COL))
    
    out <- df %>%
        dplyr::transmute(
            pvalue = suppressWarnings(as.numeric(pvalue)),
            z_mean = suppressWarnings(as.numeric(z_mean)),
            chr_raw = trimws(as.character(.data[[CHR_COL]])),
            pos = suppressWarnings(as.numeric(.data[[POS_COL]]))
        ) %>%
        dplyr::mutate(
            chr_raw = gsub("^chr", "", chr_raw, ignore.case = TRUE),
            chr_num = suppressWarnings(as.integer(chr_raw))
        ) %>%
        dplyr::filter(
            is.finite(pvalue), pvalue > 0,
            is.finite(z_mean),
            is.finite(chr_num), chr_num >= 1, chr_num <= 22,
            is.finite(pos), pos > 0
        ) %>%
        dplyr::mutate(
            signed_logp = sign(z_mean) * (-log10(pvalue)),
            direction = ifelse(z_mean >= 0, "pos", "neg")
        )
    
    if (nrow(out) == 0) return(NULL)
    
    chr_lens <- out %>%
        dplyr::group_by(chr_num) %>%
        dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(chr_num) %>%
        dplyr::mutate(
            chr_len = as.numeric(chr_len),
            cum_start = dplyr::lag(cumsum(as.numeric(chr_len)), default = 0)
        )
    
    axis_set <- chr_lens %>%
        dplyr::mutate(center = cum_start + chr_len / 2)
    
    out <- out %>%
        dplyr::left_join(chr_lens, by = "chr_num") %>%
        dplyr::mutate(
            bp_cum = as.numeric(pos) + as.numeric(cum_start)
        )
    
    list(df = out, axis_set = axis_set)
}

compute_global_ylim <- function(qtl) {
    yvals <- c()
    
    for (ds in datasets) {
        df <- df_twas_results[[ds]][[qtl]][["smultixcan"]][["brain_tissues"]]
        if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) next
        if (!("pvalue" %in% names(df)) || !("z_mean" %in% names(df))) next
        
        p <- suppressWarnings(as.numeric(df$pvalue))
        z <- suppressWarnings(as.numeric(df$z_mean))
        
        ok <- is.finite(p) & p > 0 & is.finite(z)
        if (!any(ok)) next
        
        yvals <- c(yvals, sign(z[ok]) * (-log10(p[ok])))
    }
    
    if (length(yvals) == 0) return(c(-5, 5))
    
    lim <- max(abs(yvals), na.rm = TRUE)
    lim <- ceiling(lim * 1.05)  # small headroom
    
    c(-lim, lim)
}

ylim_eqtl <- compute_global_ylim("eqtl")
ylim_sqtl <- compute_global_ylim("sqtl")

cat("Fixed y-limits:\n",
    "  eQTL:", paste(ylim_eqtl, collapse = " to "), "\n",
    "  sQTL:", paste(ylim_sqtl, collapse = " to "), "\n")

plot_twosided_manhattan <- function(df_in, axis_set, thr,
                                    title, subtitle = NULL,
                                    point_size = 0.45,
                                    ylim_fixed = NULL) {
    
    chr_cols <- rep(c("#6E6E6E", "#B0B0B0"), length.out = 22)
    names(chr_cols) <- as.character(1:22)
    
    df_in <- df_in %>%
        dplyr::mutate(
            bonf_sig = is.finite(thr$bonf_p) & (pvalue < thr$bonf_p),
            fdr_sig  = is.finite(thr$fdr_p)  & (pvalue < thr$fdr_p)
        )
    
    g <- ggplot(df_in, aes(x = bp_cum, y = signed_logp)) +
        geom_point(aes(colour = as.character(chr_num)), alpha = 0.65, size = point_size) +
        scale_colour_manual(values = chr_cols) +
        geom_point(
            data = df_in %>% dplyr::filter(fdr_sig),
            aes(fill = direction),
            shape = 21, colour = "black", stroke = 0.15,
            alpha = 0.90, size = point_size + 0.25
        ) +
        geom_point(
            data = df_in %>% dplyr::filter(bonf_sig),
            aes(fill = direction),
            shape = 21, colour = "black", stroke = 0.20,
            alpha = 0.95, size = point_size + 0.55
        ) +
        scale_fill_manual(values = c(pos = unname(OKABE_ITO["verm"]), neg = unname(OKABE_ITO["blue"]))) +
        guides(colour = "none", fill = "none") +
        scale_x_continuous(
            breaks = axis_set$center,
            labels = axis_set$chr_num,
            expand = expansion(mult = 0.01),
            guide  = guide_axis(check.overlap = FALSE)
        ) +
        scale_y_continuous(
            limits = ylim_fixed,
            labels = function(x) abs(x),
            expand = expansion(mult = 0)
        ) +
        labs(
            x = "Chromosome",
            y = expression(paste("Signed ", -log[10], italic("(P)"))),
            title = title,
            subtitle = subtitle
        ) +
        theme_nature_manh()
    
    if (is.finite(thr$fdr_p)) {
        fdr_lp <- -log10(thr$fdr_p)
        g <- g +
            geom_hline(yintercept =  fdr_lp, linetype = "dotted", linewidth = 0.35) +
            geom_hline(yintercept = -fdr_lp, linetype = "dotted", linewidth = 0.35)
    }
    
    if (is.finite(thr$bonf_p)) {
        bonf_lp <- -log10(thr$bonf_p)
        g <- g +
            geom_hline(yintercept =  bonf_lp, linetype = "dashed", linewidth = 0.45) +
            geom_hline(yintercept = -bonf_lp, linetype = "dashed", linewidth = 0.45)
    }
    
    g
}
run_all_twosided_manhattan_brain_only <- function() {
    
    for (qtl in c("eqtl", "sqtl")) {
        
        out_dir  <- if (qtl == "eqtl") dir_manh_eqtl else dir_manh_sqtl
        ylim_use <- if (qtl == "eqtl") ylim_eqtl     else ylim_sqtl
        
        for (ds in datasets) {
            
            df <- df_twas_results[[ds]][[qtl]][["smultixcan"]][["brain_tissues"]]
            
            if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
                message("Skip: ", ds, " | ", qtl, " (missing/empty smultixcan brain_tissues)")
                next
            }
            
            tier <- if (ds == primary_dataset) "Primary"
            else if (ds %in% secondary_datasets) "Secondary"
            else "Replication"
            
            thr  <- get_thresholds(df$pvalue)
            prep <- prep_twosided_manhattan(df, qtl)
            if (is.null(prep)) {
                message("Skip: ", ds, " | ", qtl, " (no valid autosomal points after filtering)")
                next
            }
            
            subtitle <- paste0(
                tier, " | smultixcan brain_tissues | n=", nrow(df),
                if (is.finite(thr$fdr_p)) paste0(" | FDR q<0.05: p<", signif(thr$fdr_p, 3)) else " | FDR q<0.05: none",
                if (is.finite(thr$bonf_p)) paste0(" | Bonf: p<", signif(thr$bonf_p, 3)) else " | Bonf: NA"
            )
            
            g <- plot_twosided_manhattan(
                df_in = prep$df,
                axis_set = prep$axis_set,
                thr = thr,
                title = paste0("BAG TWAS (", qtl, "): ", ds),
                subtitle = subtitle,
                point_size = 0.45,
                ylim_fixed = ylim_use
            )
            
            fn <- file.path(out_dir, paste0("twosided_manhattan_", qtl, "_", ds, "_smultixcan_brain_tissues.tiff"))
            
            ggsave(
                filename = fn,
                plot = g,
                device = "tiff",
                dpi = 600,
                compression = "lzw",
                width = 180, height = 95,
                units = "mm"
            )
            
        }
    }
    
}

run_all_twosided_manhattan_brain_only()

dir_manh_base_ps <- [[PLACEHOLDER]]
dir_manh_eqtl_ps <- file.path(dir_manh_base_ps, "eqtl")
dir_manh_sqtl_ps <- file.path(dir_manh_base_ps, "sqtl")
dir.create(dir_manh_eqtl_ps, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_manh_sqtl_ps, recursive = TRUE, showWarnings = FALSE)

# Dataset tiering
primary_dataset <- "wen_davatzikos_combined"
secondary_datasets <- c("wen_davatzikos_fc", "wen_davatzikos_gm", "wen_davatzikos_wm")
datasets <- names(df_twas_results)

compute_primary_ylim <- function(qtl, primary_ds, quant = 0.999, max_cap = 60) {
    
    df <- df_twas_results[[primary_ds]][[qtl]][["smultixcan"]][["brain_tissues"]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(c(-5, 5))
    if (!("pvalue" %in% names(df)) || !("z_mean" %in% names(df))) return(c(-5, 5))
    
    p <- suppressWarnings(as.numeric(df$pvalue))
    z <- suppressWarnings(as.numeric(df$z_mean))
    
    ok <- is.finite(p) & p > 0 & is.finite(z)
    if (!any(ok)) return(c(-5, 5))
    
    y <- sign(z[ok]) * (-log10(p[ok]))
    y <- y[is.finite(y)]
    if (length(y) == 0) return(c(-5, 5))
    
    lim <- stats::quantile(abs(y), probs = quant, na.rm = TRUE, names = FALSE)
    lim <- ceiling(lim * 1.10)   # headroom
    lim <- min(lim, max_cap)     # cap
    
    c(-lim, lim)
}

ylim_eqtl_primary <- compute_primary_ylim("eqtl", primary_dataset, quant = 0.999, max_cap = 60)
ylim_sqtl_primary <- compute_primary_ylim("sqtl", primary_dataset, quant = 0.999, max_cap = 60)

clamp_to_ylim <- function(df_in, ylim_fixed) {
    ylo <- ylim_fixed[1]
    yhi <- ylim_fixed[2]
    
    df_in %>%
        dplyr::mutate(
            signed_logp_raw = signed_logp,
            overflow = ifelse(signed_logp_raw < ylo | signed_logp_raw > yhi, TRUE, FALSE),
            signed_logp = pmax(pmin(signed_logp_raw, yhi), ylo)
        )
}

run_all_twosided_manhattan_primary_ylim_clamped <- function() {
    
    for (qtl in c("eqtl", "sqtl")) {
        
        out_dir  <- if (qtl == "eqtl") dir_manh_eqtl_ps else dir_manh_sqtl_ps
        ylim_use <- if (qtl == "eqtl") ylim_eqtl_primary else ylim_sqtl_primary
        
        for (ds in datasets) {
            
            df <- df_twas_results[[ds]][[qtl]][["smultixcan"]][["brain_tissues"]]
            if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
                message("Skip: ", ds, " | ", qtl, " (missing/empty smultixcan brain_tissues)")
                next
            }
            
            tier <- if (ds == primary_dataset) "Primary"
            else if (ds %in% secondary_datasets) "Secondary"
            else "Replication"
            
            thr  <- get_thresholds(df$pvalue)
            prep <- prep_twosided_manhattan(df, qtl)
            if (is.null(prep)) {
                message("Skip: ", ds, " | ", qtl, " (no valid autosomal points after filtering)")
                next
            }
            
            # Clamp to primary-derived y-limits
            prep$df <- clamp_to_ylim(prep$df, ylim_use)
            
            subtitle <- paste0(
                tier, " | smultixcan brain_tissues | n=", nrow(df),
                if (is.finite(thr$fdr_p)) paste0(" | FDR q<0.05: p<", signif(thr$fdr_p, 3)) else " | FDR q<0.05: none",
                if (is.finite(thr$bonf_p)) paste0(" | Bonf: p<", signif(thr$bonf_p, 3)) else " | Bonf: NA",
                " | y-lim from primary (clamped)"
            )
            
            g <- plot_twosided_manhattan(
                df_in = prep$df,
                axis_set = prep$axis_set,
                thr = thr,
                title = paste0("BAG TWAS (", qtl, "): ", ds),
                subtitle = subtitle,
                point_size = 0.45,
                ylim_fixed = ylim_use
            )
            
            fn <- file.path(out_dir, paste0("twosided_manhattan_", qtl, "_", ds, "_smultixcan_brain_tissues_primary_ylim_clamped.tiff"))
            
            ggsave(
                filename = fn,
                plot = g,
                device = "tiff",
                dpi = 600,
                compression = "lzw",
                width = 180, height = 95,
                units = "mm"
            )
            
            n_over <- sum(prep$df$overflow %in% TRUE, na.rm = TRUE)
            if (n_over > 0) message("  (clamped overflow points: ", n_over, ")")
            
            message("Wrote: ", fn)
        }
    }
    
}

run_all_twosided_manhattan_primary_ylim_clamped()




# tables: PRIMARY significant features

primary_dataset <- "wen_davatzikos_combined"
primary_model   <- "smultixcan"
primary_scope   <- "brain_tissues"
alpha_fdr_primary <- 0.05

dir_tables <- file.path([[PLACEHOLDER]])
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

make_primary_sig_table <- function(qtl) {
    
    df <- df_twas_results[[primary_dataset]][[qtl]][[primary_model]][[primary_scope]]
    stopifnot(is.data.frame(df), nrow(df) > 0)
    
    p <- suppressWarnings(as.numeric(df$pvalue))
    z <- suppressWarnings(as.numeric(df$z_mean))
    
    out <- df %>%
        dplyr::mutate(
            feature_id = as.character(gene_name),
            pvalue_num = p,
            z_mean_num = z,
            fdr_q = p.adjust(pvalue_num, method = "BH"),
            primary_sig = fdr_q < alpha_fdr_primary,
            direction = dplyr::if_else(is.finite(z_mean_num) & z_mean_num >= 0, "Positive", "Negative")
        ) %>%
        dplyr::filter(is.finite(pvalue_num), pvalue_num > 0, primary_sig %in% TRUE)
    
    if (qtl == "eqtl") {
        out <- out %>%
            dplyr::transmute(
                rank = dplyr::row_number(),
                feature_id = feature_id,
                chr = as.character(chromosome_name),
                pos = suppressWarnings(as.integer(start_position)),
                direction = direction,
                z_mean = z_mean_num,
                pvalue = pvalue_num,
                fdr_q = fdr_q
            )
    } else if (qtl == "sqtl") {
        out <- out %>%
            dplyr::transmute(
                rank = dplyr::row_number(),
                feature_id = feature_id,
                chr = as.character(intron_chromosome_name),
                pos = suppressWarnings(as.integer(intron_start_position)),
                direction = direction,
                z_mean = z_mean_num,
                pvalue = pvalue_num,
                fdr_q = fdr_q
            )
    }
    
    out %>% dplyr::arrange(pvalue) %>% dplyr::mutate(rank = dplyr::row_number())
}

# Build tables
tab_primary_eqtl <- make_primary_sig_table("eqtl")
tab_primary_sqtl <- make_primary_sig_table("sqtl")

# Write CSVs
fn_eqtl <- file.path(dir_tables, "primary_sig_eqtl_smultixcan_brain_tissues.csv")
fn_sqtl <- file.path(dir_tables, "primary_sig_sqtl_smultixcan_brain_tissues.csv")

readr::write_csv(tab_primary_eqtl, fn_eqtl)
readr::write_csv(tab_primary_sqtl, fn_sqtl)

# Venn

dir_venn <- [[PLACEHOLDER]]
dir.create(dir_venn, recursive = TRUE, showWarnings = FALSE)

# Datasets 
ds_primary <- "wen_davatzikos_combined"
ds_gm      <- "wen_davatzikos_gm"
ds_wm      <- "wen_davatzikos_wm"
ds_fc      <- "wen_davatzikos_fc"

model <- "smultixcan"
scope <- "brain_tissues"

alpha_fdr_primary <- 0.05   
alpha_nominal      <- 0.05 

get_sig_ids_primary_fdr <- function(dataset, qtl, alpha_fdr = 0.05) {
    df <- df_twas_results[[dataset]][[qtl]][[model]][[scope]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(character(0))
    if (!("gene_name" %in% names(df)) || !("pvalue" %in% names(df))) return(character(0))
    
    df %>%
        dplyr::transmute(
            id = as.character(gene_name),
            p  = suppressWarnings(as.numeric(pvalue))
        ) %>%
        dplyr::filter(!is.na(id), id != "", is.finite(p), p > 0) %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(p = min(p, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(q = p.adjust(p, method = "BH")) %>%
        dplyr::filter(is.finite(q), q < alpha_fdr) %>%
        dplyr::pull(id) %>%
        unique()
}

get_sig_ids_nominal <- function(dataset, qtl, alpha_p = 0.05) {
    df <- df_twas_results[[dataset]][[qtl]][[model]][[scope]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(character(0))
    if (!("gene_name" %in% names(df)) || !("pvalue" %in% names(df))) return(character(0))
    
    df %>%
        dplyr::transmute(
            id = as.character(gene_name),
            p  = suppressWarnings(as.numeric(pvalue))
        ) %>%
        dplyr::filter(!is.na(id), id != "", is.finite(p), p > 0) %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(p = min(p, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(p), p < alpha_p) %>%
        dplyr::pull(id) %>%
        unique()
}

#  plot
plot_venn_nature <- function(qtl) {
    
    primary_ids <- get_sig_ids_primary_fdr(ds_primary, qtl, alpha_fdr_primary)
    gm_ids      <- get_sig_ids_nominal(ds_gm, qtl, alpha_nominal)
    wm_ids      <- get_sig_ids_nominal(ds_wm, qtl, alpha_nominal)
    fc_ids      <- get_sig_ids_nominal(ds_fc, qtl, alpha_nominal)
    
    use_fc <- length(fc_ids) > 0
    
    subtitle_txt <- paste0(
        "Primary: BH-FDR q<", alpha_fdr_primary,
        "; GM/WM", if (use_fc) "/FC" else "",
        ": nominal p<", alpha_nominal
    )
    col_invis <- grDevices::adjustcolor("black", alpha.f = 0)
    
    n_primary <- length(primary_ids)
    n_gm      <- length(gm_ids)
    n_wm      <- length(wm_ids)
    n_fc      <- length(fc_ids)
    
    n_pg  <- length(intersect(primary_ids, gm_ids))
    n_pw  <- length(intersect(primary_ids, wm_ids))
    n_pf  <- length(intersect(primary_ids, fc_ids))
    
    n_pgw <- length(Reduce(intersect, list(primary_ids, gm_ids, wm_ids)))
    n_pgf <- length(Reduce(intersect, list(primary_ids, gm_ids, fc_ids)))
    n_pwf <- length(Reduce(intersect, list(primary_ids, wm_ids, fc_ids)))
    
    n_pgwf <- length(Reduce(intersect, list(primary_ids, gm_ids, wm_ids, fc_ids)))
    
    grid.newpage()
    
    if (use_fc) {
        n_gw  <- length(intersect(gm_ids, wm_ids))
        n_gf  <- length(intersect(gm_ids, fc_ids))
        n_wf  <- length(intersect(wm_ids, fc_ids))
        
        n_gwf <- length(Reduce(intersect, list(gm_ids, wm_ids, fc_ids)))
        draw.quad.venn(
            fontfamily = "Arial",
            family = "Arial",
            cat.fontfamily = "Arial",
            cat.family = "Arial",
            area1 = n_primary,
            area2 = n_gm,
            area3 = n_wm,
            area4 = n_fc,
            
            n12   = n_pg,
            n13   = n_pw,
            n14   = n_pf,
            
            n23   = n_gw,
            n24   = n_gf,
            n34   = n_wf,
            
            n123  = n_pgw,
            n124  = n_pgf,
            n134  = n_pwf,
            n234  = n_gwf,
            
            n1234 = n_pgwf,
            
            category = c("Primary", "GM", "WM", "FC"),
            
            label.col = c(
                col_invis, col_invis, col_invis, "black",
                "black", "black", col_invis,
                col_invis, "black", "black",
                "black", "black", col_invis,
                col_invis,
                "black"
            ),
            
            label.cex = c(
                1.2, 0, 0, 0,
                1.2, 1.2, 1.2,
                0, 0, 0,
                1.2, 1.2, 1.2,
                0,
                1.2
            ),
            
            fill  = c("#B2182B", "#BDBDBD", "#BDBDBD", "#BDBDBD"),
            alpha = c(0.6, 0.3, 0.3, 0.3),
            
            cat.cex = 1.2,
            margin = 0.1
        )
        
    } else {
        
        draw.triple.venn(
            area1 = n_primary,
            area2 = n_gm,
            area3 = n_wm,
            
            n12  = n_pg,
            n13  = n_pw,
            n123 = n_pgw,
            
            category = c("Primary", "GM", "WM"),
            
            label.col = c(
                "black",
                col_invis,
                col_invis,
                "black",
                "black",
                col_invis,
                "black"
            ),
            
            label.cex = c(1.2, 1.2, 1.2, 0, 0, 1.2),
            
            fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
            alpha = c(0.6, 0.3, 0.3),
            
            cat.cex = 1.2
        )
    }
    
    grid.text(
        subtitle_txt,
        y = 0.05,
        gp = gpar(fontsize = 9)
    )
}

save_venn_tiff <- function(qtl, filename, width_mm = 110, height_mm = 85, dpi = 600) {
    grDevices::tiff(filename, width = width_mm, height = height_mm, units = "mm",
                    res = dpi, compression = "lzw", type="cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_venn_nature(qtl)
}

venn_filename <- function(qtl) {
    fc_n <- length(get_sig_ids_nominal(ds_fc, qtl, alpha_nominal))
    if (fc_n > 0) {
        file.path(dir_venn, paste0("venn_primary_GM_WM_FC_", qtl, "_nature.tiff"))
    } else {
        file.path(dir_venn, paste0("venn_primary_GM_WM_", qtl, "_nature.tiff"))
    }
}

p_eqtl <- plot_venn_nature("eqtl")
fn_eqtl <- venn_filename("eqtl")
save_venn_tiff("eqtl", fn_eqtl)
cat("Wrote: ", fn_eqtl, "\n", sep = "")

p_sqtl <- plot_venn_nature("sqtl")
fn_sqtl <- venn_filename("sqtl")
save_venn_tiff("sqtl", fn_sqtl)
cat("Wrote: ", fn_sqtl, "\n", sep = "")

sets_eqtl <- list(
    Primary = get_sig_ids_primary_fdr(ds_primary, "eqtl", alpha_fdr_primary),
    GM      = get_sig_ids_nominal(ds_gm, "eqtl", alpha_nominal),
    WM      = get_sig_ids_nominal(ds_wm, "eqtl", alpha_nominal)
)
fc_eqtl <- get_sig_ids_nominal(ds_fc, "eqtl", alpha_nominal)
if (length(fc_eqtl) > 0) sets_eqtl$FC <- fc_eqtl

sets_sqtl <- list(
    Primary = get_sig_ids_primary_fdr(ds_primary, "sqtl", alpha_fdr_primary),
    GM      = get_sig_ids_nominal(ds_gm, "sqtl", alpha_nominal),
    WM      = get_sig_ids_nominal(ds_wm, "sqtl", alpha_nominal)
)
fc_sqtl <- get_sig_ids_nominal(ds_fc, "sqtl", alpha_nominal)
if (length(fc_sqtl) > 0) sets_sqtl$FC <- fc_sqtl

cat("\nSet sizes (eQTL):",
    paste(names(sets_eqtl), vapply(sets_eqtl, length, integer(1)), sep="=", collapse="  "),
    "\n", sep = " ")
cat("Set sizes (sQTL):",
    paste(names(sets_sqtl), vapply(sets_sqtl, length, integer(1)), sep="=", collapse="  "),
    "\n", sep = " ")



# WGCNA

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
datasets_all <- names(df_twas_results)

## output dirs
dir_wgcna_res <- file.path(dir_root, "results", "wgcna")
dir_wgcna_fig <- file.path(dir_root, "images",  "wgcna")
dir.create(dir_wgcna_res, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_wgcna_fig, recursive = TRUE, showWarnings = FALSE)
subdirs <- c("W1_eqtl","W1_sqtl","W2_eqtl","W2_sqtl")
for (sd in subdirs) {
    dir.create(file.path(dir_wgcna_res, sd), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dir_wgcna_fig, sd), recursive = TRUE, showWarnings = FALSE)
}

# conf
primary_dataset <- "wen_davatzikos_combined"
secondary_datasets <- c("wen_davatzikos_gm","wen_davatzikos_wm","wen_davatzikos_fc")

# 49 tissues
W1_MAX_BLOCK <- 5000
W1_MIN_MODULE <- 30
W1_DEEPSPLIT  <- 2

SOFT_POWERS <- c(1:10, seq(12, 30, 2))
NETWORK_TYPE <- "signed"
COR_TYPE <- "pearson"
MERGE_CUT_HEIGHT <- 0.25

.build_feature_by_point_intersection <- function(named_slices,
                                                 id_col = "gene_name",
                                                 z_col,
                                                 verbose = TRUE) {
    stopifnot(is.list(named_slices), length(named_slices) >= 2)
    
    vecs <- lapply(named_slices, function(df) {
        if (is.null(df) || !is.data.frame(df)) return(NULL)
        if (!(id_col %in% names(df))) return(NULL)
        if (!(z_col %in% names(df))) return(NULL)
        
        ids <- as.character(df[[id_col]])
        z   <- suppressWarnings(as.numeric(df[[z_col]]))
        
        ok <- !is.na(ids) & ids != "" & is.finite(z)
        ids <- ids[ok]; z <- z[ok]
        if (length(ids) == 0) return(NULL)
        
        tapply(z, ids, function(v) v[which.max(abs(v))][1])
    })
    
    vecs <- vecs[!vapply(vecs, is.null, logical(1))]
    if (length(vecs) < 2) stop("Too few valid points after extracting ", z_col)
    
    common_ids <- Reduce(intersect, lapply(vecs, names))
    common_ids <- unique(common_ids)
    
    mat <- sapply(vecs, function(v) v[common_ids])
    mat <- as.matrix(mat)
    rownames(mat) <- common_ids
    colnames(mat) <- names(vecs)
    
    if (sum(is.na(mat)) > 0) stop("Intersection matrix has NA values; check inputs.")
    mat
}

.to_wgcna_datExpr <- function(feature_by_point, verbose = TRUE) {
    datExpr <- t(feature_by_point)
    
    ok_cols <- apply(datExpr, 2, function(v) all(is.finite(v)))
    datExpr <- datExpr[, ok_cols, drop = FALSE]
    
    sds <- apply(datExpr, 2, stats::sd)
    keep <- !is.na(sds) & sds > 0
    datExpr <- datExpr[, keep, drop = FALSE]

    datExpr
}

build_wgcna_matrix_spredixcan <- function(df_twas_results, dataset, qtl = c("eqtl","sqtl"),
                                          tissues = NULL, verbose = TRUE) {
    qtl <- match.arg(qtl)
    
    x <- df_twas_results[[dataset]][[qtl]][["spredixcan"]]
    if (is.null(x)) stop()
    
    if (!is.null(tissues)) x <- x[intersect(names(x), tissues)]
    x <- x[vapply(x, is.data.frame, logical(1))]
    if (length(x) < 3) stop()
    
    mat <- .build_feature_by_point_intersection(
        named_slices = x,
        id_col = "gene_name",
        z_col  = "zscore",
        verbose = verbose
    )
    
    list(
        feature_by_point = mat,
        datExpr = .to_wgcna_datExpr(mat, verbose = verbose)
    )
}

build_wgcna_matrix_smultixcan <- function(df_twas_results, datasets_vec, qtl = c("eqtl","sqtl"),
                                          scope = c("brain_tissues","all_tissues"),
                                          verbose = TRUE) {
    qtl <- match.arg(qtl)
    scope <- match.arg(scope)
    
    slices <- list()
    missing <- character(0)
    
    for (ds in datasets_vec) {
        df <- df_twas_results[[ds]][[qtl]][["smultixcan"]][[scope]]
        if (is.data.frame(df) && nrow(df) > 0) slices[[ds]] <- df else missing <- c(missing, ds)
    }
    
    if (length(slices) < 3) {
        stop("Need >=3 datasets for WGCNA. Found: ", length(slices),
             if (length(missing)>0) paste0(" (missing: ", paste(missing, collapse=", "), ")") else "")
    }
    
    mat <- .build_feature_by_point_intersection(
        named_slices = slices,
        id_col = "gene_name",
        z_col  = "z_mean",
        verbose = verbose
    )
    
    list(
        feature_by_point = mat,
        datExpr = .to_wgcna_datExpr(mat, verbose = verbose)
    )
}

## helpers
make_traits_W1 <- function(datExpr) {
    data.frame(
        is_brain = as.integer(grepl("^Brain_", rownames(datExpr))),
        row.names = rownames(datExpr),
        stringsAsFactors = FALSE
    )
}

.save_tiff <- function(fn, w_mm = 180, h_mm = 120, res = 600) {
    tiff(fn, units="mm", width=w_mm, height=h_mm, res=res, compression="lzw")
}

run_wgcna_export <- function(datExpr, traits, tag, out_res_dir, out_fig_dir,
                             networkType = NETWORK_TYPE,
                             corType = COR_TYPE,
                             powers = SOFT_POWERS,
                             deepSplit = 2,
                             minModuleSize = 30,
                             mergeCutHeight = MERGE_CUT_HEIGHT,
                             maxBlockSize = 5000,
                             verbose = 3) {
    
    dir.create(out_res_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
    
    datExpr <- as.data.frame(datExpr)
    traits  <- as.data.frame(traits)
    
    # align traits to samples
    traits <- traits[rownames(datExpr), , drop = FALSE]
    
    # QC
    gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = verbose)
    if (!gsg$allOK) {
        datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
        traits  <- traits[rownames(datExpr), , drop = FALSE]
    }
    
    # Soft threshold
    sft <- WGCNA::pickSoftThreshold(
        datExpr,
        powerVector = powers,
        networkType = networkType,
        corFnc = if (corType == "pearson") "cor" else "bicor",
        verbose = verbose
    )
    
    fit <- sft$fitIndices
    idx <- which(fit$SFT.R.sq >= 0.80)
    softPower <- if (length(idx) > 0) fit$Power[min(idx)] else fit$Power[which.max(fit$SFT.R.sq)]
    
    #  Soft-threshold figure 
    fn_sft <- file.path(out_fig_dir, paste0(tag, "_soft_threshold.tiff"))
    tiff(fn_sft, units="mm", width=180, height=90, res=600, compression="lzw")
    par(mfrow=c(1,2), mar=c(4,4,2,1))
    plot(fit$Power, fit$SFT.R.sq, type="n",
         xlab="Soft threshold (power)", ylab="Scale-free fit (R^2)",
         main=paste0(tag, " scale-free fit"))
    text(fit$Power, fit$SFT.R.sq, labels=fit$Power, cex=0.7)
    abline(h=0.80, col="red", lty=2)
    
    plot(fit$Power, fit$mean.k., type="n",
         xlab="Soft threshold (power)", ylab="Mean connectivity",
         main=paste0(tag, " mean connectivity"))
    text(fit$Power, fit$mean.k., labels=fit$Power, cex=0.7)
    dev.off()
    
    #  Modules 
    net <- WGCNA::blockwiseModules(
        datExpr,
        power = softPower,
        networkType = networkType,
        TOMType = if (networkType == "signed") "signed" else "unsigned",
        corType = corType,
        maxBlockSize = maxBlockSize,
        minModuleSize = minModuleSize,
        deepSplit = deepSplit,
        mergeCutHeight = mergeCutHeight,
        numericLabels = TRUE,
        pamRespectsDendro = TRUE,
        saveTOMs = FALSE,
        verbose = verbose
    )
    
    moduleLabels <- net$colors
    moduleColors <- WGCNA::labels2colors(moduleLabels)
    MEs <- WGCNA::orderMEs(net$MEs)
    
    #  Dendrogram + colors 
    fn_dend <- file.path(out_fig_dir, paste0(tag, "_dendrogram_modules.tiff"))
    tiff(fn_dend, units="mm", width=180, height=140, res=600, compression="lzw")
    WGCNA::plotDendroAndColors(
        net$dendrograms[[1]],
        moduleColors[net$blockGenes[[1]]],
        groupLabels = "Modules",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05,
        main = paste0(tag, " dendrogram & modules")
    )
    dev.off()
    
    #  Module sizes 
    tab_sizes <- sort(table(moduleColors), decreasing = TRUE)
    df_sizes <- data.frame(
        module_color = names(tab_sizes),
        n_genes = as.integer(tab_sizes),
        stringsAsFactors = FALSE
    )
    write.csv(df_sizes, file.path(out_res_dir, paste0(tag, "_module_sizes.csv")), row.names = FALSE)
    
    #  Module–trait correlation & p 
    moduleTraitCor <- stats::cor(MEs, traits, use = "p")
    moduleTraitP   <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples=nrow(datExpr))
    
    #  Heatmap  
    fn_mth <- file.path(out_fig_dir, paste0(tag, "_module_trait_heatmap.tiff"))
    tiff(fn_mth, units="mm", width=180, height=140, res=600, compression="lzw")
    
    textMatrix <- paste0(
        signif(moduleTraitCor, 2), "\n(",
        format(moduleTraitP, digits=2, scientific=TRUE), ")"
    )
    dim(textMatrix) <- dim(moduleTraitCor)
    
    yLabs <- rownames(moduleTraitCor)
    
    WGCNA::labeledHeatmap(
        Matrix = moduleTraitCor,
        xLabels = colnames(traits),
        yLabels = yLabs,
        ySymbols = yLabs,
        colorLabels = FALSE,
        colors = WGCNA::blueWhiteRed(50),
        textMatrix = textMatrix,
        setStdMargins = FALSE,
        cex.text = 0.55,
        zlim = c(-1,1),
        main = paste0(tag, " module–trait relationships")
    )
    dev.off()
    
    #  Eigengene clustering 
    fn_meclust <- file.path(out_fig_dir, paste0(tag, "_eigengene_clustering.tiff"))
    tiff(fn_meclust, units="mm", width=180, height=120, res=600, compression="lzw")
    MEcor <- stats::cor(MEs, use = "p")
    MEtree <- stats::hclust(stats::as.dist(1 - MEcor), method="average")
    plot(MEtree, main=paste0(tag, " eigengene clustering"), xlab="", sub="")
    dev.off()
    
    #  Long module–trait table
    ME_names <- rownames(moduleTraitCor)
    trait_names <- colnames(moduleTraitCor)
    
    df_mt <- expand.grid(ME = ME_names, trait = trait_names, stringsAsFactors = FALSE)
    
    r <- match(df_mt$ME, ME_names)
    c <- match(df_mt$trait, trait_names)
    
    df_mt$cor <- moduleTraitCor[cbind(r, c)]
    df_mt$p   <- moduleTraitP[cbind(r, c)]
    
    df_mt$module_label <- as.integer(sub("^ME", "", df_mt$ME))
    df_mt$module_color <- WGCNA::labels2colors(df_mt$module_label)
    
    df_mt <- df_mt[order(df_mt$p), ]
    write.csv(df_mt, file.path(out_res_dir, paste0(tag, "_module_trait_long.csv")), row.names = FALSE)
    
    
    #  Best module per trait 
    df_best <- do.call(rbind, lapply(trait_names, function(tr) {
        sub <- df_mt[df_mt$trait == tr, ]
        sub[which.min(sub$p), c("trait","ME","module_label","module_color","cor","p")]
    }))
    rownames(df_best) <- NULL
    write.csv(df_best, file.path(out_res_dir, paste0(tag, "_best_module_per_trait.csv")), row.names = FALSE)
    
    #  Hub tables 
    kME <- stats::cor(datExpr, MEs, use="p")
    colnames(kME) <- paste0("kME_", colnames(MEs))
    
    get_hubs_for_ME <- function(MEname, top_n) {
        lbl <- as.integer(sub("^ME","", MEname))
        in_mod <- which(moduleLabels == lbl)
        if (length(in_mod) == 0) return(NULL)
        
        genes <- colnames(datExpr)[in_mod]
        kvec <- kME[in_mod, paste0("kME_", MEname)]
        o <- order(abs(kvec), decreasing = TRUE)
        
        head_n <- min(top_n, length(o))
        data.frame(
            module_ME = rep(MEname, head_n),
            module_label = rep(lbl, head_n),
            module_color = rep(WGCNA::labels2colors(lbl), head_n),
            gene = genes[o][seq_len(head_n)],
            kME = as.numeric(kvec[o][seq_len(head_n)]),
            stringsAsFactors = FALSE
        )
    }
    
    hubs20 <- do.call(rbind, lapply(colnames(MEs), get_hubs_for_ME, top_n = 20))
    hubs100 <- do.call(rbind, lapply(colnames(MEs), get_hubs_for_ME, top_n = 100))
    
    if (!is.null(hubs20))  write.csv(hubs20,  file.path(out_res_dir, paste0(tag, "_hubs_top20.csv")), row.names = FALSE)
    if (!is.null(hubs100)) write.csv(hubs100, file.path(out_res_dir, paste0(tag, "_hubs_top100_SUPP.csv")), row.names = FALSE)
    
    #  Genemodule membership 
    df_mem <- data.frame(
        gene = colnames(datExpr),
        module_label = moduleLabels,
        module_color = moduleColors,
        stringsAsFactors = FALSE
    )
    write.csv(df_mem, file.path(out_res_dir, paste0(tag, "_gene_module_membership_SUPP.csv")), row.names = FALSE)
    
    # Return object
    list(
        tag = tag,
        datExpr = datExpr,
        traits = traits,
        softPower = softPower,
        sft = sft,
        net = net,
        moduleLabels = moduleLabels,
        moduleColors = moduleColors,
        MEs = MEs,
        moduleTraitCor = moduleTraitCor,
        moduleTraitP = moduleTraitP,
        files = list(
            soft_threshold = fn_sft,
            dendrogram = fn_dend,
            module_trait_heatmap = fn_mth,
            eigengene_clustering = fn_meclust
        )
    )
}

## RUN

# W1 eQTL (primary, cross-tissue)
w1_eqtl <- build_wgcna_matrix_spredixcan(df_twas_results, dataset=primary_dataset, qtl="eqtl", verbose=TRUE)
traits_w1_eqtl <- make_traits_W1(w1_eqtl$datExpr)
res_W1_eqtl <- run_wgcna_export(
    datExpr = w1_eqtl$datExpr,
    traits  = traits_w1_eqtl,
    tag = "W1_eqtl",
    out_res_dir = file.path(dir_wgcna_res, "W1_eqtl"),
    out_fig_dir = file.path(dir_wgcna_fig, "W1_eqtl"),
    deepSplit = W1_DEEPSPLIT,
    minModuleSize = W1_MIN_MODULE,
    maxBlockSize = W1_MAX_BLOCK
)

# W1 sQTL
w1_sqtl <- build_wgcna_matrix_spredixcan(df_twas_results, dataset=primary_dataset, qtl="sqtl", verbose=TRUE)
traits_w1_sqtl <- make_traits_W1(w1_sqtl$datExpr)
res_W1_sqtl <- run_wgcna_export(
    datExpr = w1_sqtl$datExpr,
    traits  = traits_w1_sqtl,
    tag = "W1_sqtl",
    out_res_dir = file.path(dir_wgcna_res, "W1_sqtl"),
    out_fig_dir = file.path(dir_wgcna_fig, "W1_sqtl"),
    deepSplit = W1_DEEPSPLIT,
    minModuleSize = W1_MIN_MODULE,
    maxBlockSize = W1_MAX_BLOCK
)

# SUMMARY TABLE
wgcna_overview <- function(res) {
    tibble(
        network = res$tag,
        n_points = nrow(res$datExpr),
        n_features = ncol(res$datExpr),
        softPower = res$softPower,
        n_modules_incl_grey = length(unique(res$moduleColors)),
        n_modules_excl_grey = length(setdiff(unique(res$moduleColors), "grey")),
        largest_module_n = max(table(res$moduleColors))
    )
}

overview_all <- bind_rows(
    wgcna_overview(res_W1_eqtl),
    wgcna_overview(res_W1_sqtl),
    wgcna_overview(res_W2_eqtl),
    wgcna_overview(res_W2_sqtl)
)

write_csv(overview_all, file.path(dir_wgcna_res, "WGCNA_overview_all.csv"))
overview_all

best_W1_eqtl_isbrain <- read.csv(
    file.path(dir_wgcna_res, "W1_eqtl", "W1_eqtl_best_module_per_trait.csv")
) |>
    dplyr::filter(trait == "is_brain")
best_W1_sqtl_isbrain <- read.csv(
    file.path(dir_wgcna_res, "W1_sqtl", "W1_sqtl_best_module_per_trait.csv")
) |>
    dplyr::filter(trait == "is_brain")

best_W1_eqtl_isbrain
best_W1_sqtl_isbrain

# Extract top hubs for best modules
get_top_hubs_for_best <- function(res, trait_name, top_n=20) {
    best <- res$best_by_trait %>% filter(trait==trait_name)
    if (nrow(best)==0) return(tibble())
    MEname <- best$module_ME[1]
    res$hubs_top20 %>% filter(module_ME==MEname) %>% slice_head(n=top_n)
}

hubs_W1_eqtl_20 <- read.csv(file.path(dir_wgcna_res, "W1_eqtl", "W1_eqtl_hubs_top20.csv"))
top20_W1_eqtl <- hubs_W1_eqtl_20[hubs_W1_eqtl_20$module_ME == best_W1_eqtl_isbrain$ME, ]

hubs_W1_sqtl_20 <- read.csv(file.path(dir_wgcna_res, "W1_sqtl", "W1_sqtl_hubs_top20.csv"))
top20_W1_sqtl <- hubs_W1_sqtl_20[hubs_W1_sqtl_20$module_ME == best_W1_sqtl_isbrain$ME, ]

top20_W1_eqtl
top20_W1_sqtl

# replication

primary_dataset    <- "wen_davatzikos_combined"
secondary_datasets <- c("wen_davatzikos_gm", "wen_davatzikos_wm", "wen_davatzikos_fc")

all_datasets <- names(df_twas_results)
replication_datasets <- setdiff(all_datasets, c(primary_dataset, secondary_datasets))

cat("Primary:   ", primary_dataset, "\n")
cat("Secondary: ", paste(secondary_datasets, collapse=", "), "\n")
cat("Replication (n=", length(replication_datasets), "): ", paste(replication_datasets, collapse=", "), "\n", sep="")

pick_first <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NA_character_)
    hit[1]
}

get_slice <- function(dataset, qtl, model, key4) {
    x <- df_twas_results[[dataset]][[qtl]][[model]][[key4]]
    if (is.null(x) || !is.data.frame(x)) return(NULL)
    x
}
resolve_feature_col <- function(df, qtl) {
    if (qtl == "eqtl") {
        pick_first(df, c("gene_noversion","gene","gene_id","feature","id"))
    } else {
        pick_first(df, c("intron","intron_id","event","feature","id","gene_noversion","gene"))
    }
}

resolve_z_col <- function(df) {
    pick_first(df, c("zscore","z","z_mean","zscore_mean","Z","STAT"))
}

resolve_p_col <- function(df) {
    pick_first(df, c("pvalue","p","P","pval","PVAL"))
}

# 1) Feature-level replication 
build_primary_hits <- function(qtl) {
    dfp <- get_slice(primary_dataset, qtl, "smultixcan", "brain_tissues")
    stopifnot(!is.null(dfp))
    
    fcol <- resolve_feature_col(dfp, qtl)
    zcol <- resolve_z_col(dfp)
    pcol <- resolve_p_col(dfp)
    
    stopifnot(!is.na(fcol), !is.na(zcol), !is.na(pcol))
    stopifnot("pvalue_fdr" %in% names(dfp))
    
    dfp %>%
        mutate(
            feature = .data[[fcol]],
            z_primary = as.numeric(.data[[zcol]]),
            p_primary = as.numeric(.data[[pcol]]),
            q_primary = as.numeric(pvalue_fdr)
        ) %>%
        filter(!is.na(feature), !is.na(z_primary), !is.na(p_primary), !is.na(q_primary)) %>%
        filter(q_primary < 0.05) %>%
        distinct(feature, .keep_all = TRUE) %>%
        dplyr::select(feature, z_primary, p_primary, q_primary)
}

get_replication_stats_one_dataset <- function(dataset, qtl, primary_hits_df) {
    dfr <- get_slice(dataset, qtl, "smultixcan", "brain_tissues")
    if (is.null(dfr)) return(NULL)
    
    fcol <- resolve_feature_col(dfr, qtl)
    zcol <- resolve_z_col(dfr)
    pcol <- resolve_p_col(dfr)
    if (is.na(fcol) || is.na(zcol) || is.na(pcol)) return(NULL)
    
    dfr2 <- dfr %>%
        mutate(
            feature = .data[[fcol]],
            z_rep = as.numeric(.data[[zcol]]),
            p_rep = as.numeric(.data[[pcol]])
        ) %>%
        dplyr::select(feature, z_rep, p_rep) %>%
        filter(!is.na(feature), !is.na(z_rep), !is.na(p_rep)) %>%
        distinct(feature, .keep_all = TRUE) %>%
        mutate(dataset = dataset)
    
    # Join onto primary hits and compute replication criteria
    primary_hits_df %>%
        left_join(dfr2, by = "feature") %>%
        mutate(
            rep_nominal = !is.na(p_rep) & (p_rep < 0.05),
            rep_consistent = rep_nominal & !is.na(z_rep) & (sign(z_rep) == sign(z_primary))
        )
}

replicate_feature_level <- function(qtl) {
    primary_hits <- build_primary_hits(qtl)
    
    per_dataset <- lapply(replication_datasets, function(ds) {
        get_replication_stats_one_dataset(ds, qtl, primary_hits)
    })
    per_dataset <- bind_rows(per_dataset)
    
    out_long <- per_dataset
    
    out_summary <- per_dataset %>%
        group_by(feature, z_primary, p_primary, q_primary) %>%
        summarise(
            n_rep_available   = sum(!is.na(p_rep)),
            n_rep_nominal     = sum(rep_nominal, na.rm = TRUE),
            n_rep_consistent  = sum(rep_consistent, na.rm = TRUE),
            any_rep_nominal   = any(rep_nominal, na.rm = TRUE),
            any_rep_consistent= any(rep_consistent, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(n_rep_consistent), desc(n_rep_nominal), q_primary)
    
    list(primary_hits = primary_hits, long = out_long, summary = out_summary)
}

rep_eqtl <- replicate_feature_level("eqtl")
rep_sqtl <- replicate_feature_level("sqtl")

dir_rep <- file.path(dir_results, "replication")
dir.create(dir_rep, showWarnings = FALSE, recursive = TRUE)

write.csv(rep_eqtl$summary, file.path(dir_rep, "replication_feature_level_eqtl_summary.csv"), row.names = FALSE)
write.csv(rep_eqtl$long,    file.path(dir_rep, "replication_feature_level_eqtl_long.csv"),    row.names = FALSE)

write.csv(rep_sqtl$summary, file.path(dir_rep, "replication_feature_level_sqtl_summary.csv"), row.names = FALSE)
write.csv(rep_sqtl$long,    file.path(dir_rep, "replication_feature_level_sqtl_long.csv"),    row.names = FALSE)

cat("Wrote feature-level replication CSVs to: ", dir_rep, "\n", sep="")


# replication of wgcna modules

get_twas_slice <- function(dataset, qtl, model = "smultixcan", scope = "brain_tissues") {
    
    df <- df_twas_results[[dataset]][[qtl]][[model]][[scope]]
    if (is.null(df) || !is.data.frame(df)) {
        stop(sprintf("Missing TWAS slice: %s | %s | %s | %s", dataset, qtl, model, scope))
    }
    
    id_col <- "gene_name"
    
    z_col <- if (model == "smultixcan") "z_mean" else if (model == "spredixcan") "zscore" else NA_character_
    p_col <- "pvalue"
    if (is.na(z_col) || is.na(p_col)) stop("Could not find zscore/z_mean and pvalue/pval columns in TWAS slice.")
    
    out <- df %>%
        transmute(
            feature = .data[[id_col]],
            z = as.numeric(.data[[z_col]]),
            p = as.numeric(.data[[p_col]])
        ) %>%
        filter(!is.na(feature), !is.na(z), !is.na(p))
    
    out
}

read_wgcna_exports <- function(network_dir) {
    
    prefix <- basename(network_dir)
    
    files <- list(
        membership = file.path(network_dir, paste0(prefix, "_gene_module_membership_SUPP.csv")),
        best_trait = file.path(network_dir, paste0(prefix, "_best_module_per_trait.csv")),
        hubs20     = file.path(network_dir, paste0(prefix, "_hubs_top20.csv")),
        hubs100    = file.path(network_dir, paste0(prefix, "_hubs_top100_SUPP.csv"))
    )
    
    missing <- names(files)[!file.exists(unlist(files))]
    if (length(missing) > 0) {
        stop(sprintf("Missing WGCNA export(s) in %s: %s", network_dir, paste(missing, collapse=", ")))
    }
    
    membership <- readr::read_csv(files$membership, show_col_types = FALSE)
    best_trait <- readr::read_csv(files$best_trait, show_col_types = FALSE)
    hubs20     <- readr::read_csv(files$hubs20, show_col_types = FALSE)
    hubs100    <- readr::read_csv(files$hubs100, show_col_types = FALSE)
    
    membership <- membership %>%
        rename(feature = gene)
    
    hubs20 <- hubs20 %>%
        rename(feature = gene)
    
    hubs100 <- hubs100 %>%
        rename(feature = gene)
    
    list(
        membership = membership,
        best_trait = best_trait,
        hubs20 = hubs20,
        hubs100 = hubs100
    )
}


replicate_primary_brain_module <- function(qtl, network_dir, scope = "brain_tissues",
                                           use_top_hubs = 20, p_rep = 0.05) {
    
    w <- read_wgcna_exports(network_dir)
    
    bm <- w$best_trait %>% filter(trait == "is_brain")
    if (nrow(bm) != 1) stop("Expected exactly one row for trait == 'is_brain' in best_module_per_trait.csv")
    
    brain_color <- bm$module_color[[1]]
    brain_ME    <- bm$ME[[1]]
    
    module_features <- w$membership %>%
        filter(module_color == brain_color) %>%
        pull(feature) %>%
        unique()
    
    # choose hub list 
    hubs <- w$hubs100 %>%
        filter(module_color == brain_color) %>%
        arrange(desc(kME)) %>%
        slice_head(n = use_top_hubs) %>%
        pull(feature) %>%
        unique()
    
    tw_primary <- get_twas_slice(primary_dataset, qtl, model = "smultixcan", scope = scope) %>%
        filter(feature %in% module_features) %>%
        dplyr::select(feature, z_primary = z, p_primary = p)
    
    out <- lapply(replication_datasets, function(ds) {
        
        tw_rep <- get_twas_slice(ds, qtl, model = "smultixcan", scope = scope) %>%
            filter(feature %in% module_features) %>%
            dplyr::select(feature, z_rep = z, p_rep = p)
        
        m <- inner_join(tw_primary, tw_rep, by = "feature")
        
        if (nrow(m) < 20) {
            return(tibble(
                qtl = qtl,
                module_ME = brain_ME,
                module_color = brain_color,
                dataset = ds,
                n_overlap = nrow(m),
                cor_z_primary_vs_rep = NA_real_,
                p_cor = NA_real_,
                n_rep_nominal_concordant = NA_integer_,
                prop_rep_nominal_concordant = NA_real_,
                n_hubs_overlap = sum(hubs %in% m$feature),
                n_hubs_rep_nominal_concordant = NA_integer_
            ))
        }
        
        # correlation of Z vectors across module members (primary vs replication)
        ct <- suppressWarnings(cor.test(m$z_primary, m$z_rep, method = "pearson"))
        
        rep_flag <- (m$p_rep < p_rep) & (sign(m$z_rep) == sign(m$z_primary))
        
        # hub replication
        mh <- m %>% filter(feature %in% hubs)
        hub_flag <- (mh$p_rep < p_rep) & (sign(mh$z_rep) == sign(mh$z_primary))
        
        tibble(
            qtl = qtl,
            module_ME = brain_ME,
            module_color = brain_color,
            dataset = ds,
            n_overlap = nrow(m),
            cor_z_primary_vs_rep = unname(ct$estimate),
            p_cor = ct$p.value,
            n_rep_nominal_concordant = sum(rep_flag),
            prop_rep_nominal_concordant = mean(rep_flag),
            n_hubs_overlap = nrow(mh),
            n_hubs_rep_nominal_concordant = sum(hub_flag)
        )
    }) %>% bind_rows()
    
    out
}

# Run both eQTL + sQTL 

dir_wgcna <- file.path(dir_results, "wgcna")

res_rep_eqtl <- replicate_primary_brain_module(
    qtl = "eqtl",
    network_dir = file.path(dir_wgcna, "W1_eqtl"),
    scope = "brain_tissues",
    use_top_hubs = 20,
    p_rep = 0.05
)

res_rep_sqtl <- replicate_primary_brain_module(
    qtl = "sqtl",
    network_dir = file.path(dir_wgcna, "W1_sqtl"),
    scope = "brain_tissues",
    use_top_hubs = 20,
    p_rep = 0.05
)

res_rep_wgcna <- bind_rows(res_rep_eqtl, res_rep_sqtl)

dir_out <- file.path(dir_results, "wgcna_replication")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(res_rep_wgcna, file.path(dir_out, "W1_primary_brain_module_replication_eqtl_sqtl.csv"))

print(res_rep_wgcna)

# replication barplot

#  inputs 
eqtl_file <- [[PLACEHOLDER]]
sqtl_file <- [[PLACEHOLDER]]

rep_col <- "n_rep_consistent"

library(tidyverse)
library(colorspace)
library(patchwork)  

eqtl <- readr::read_csv(eqtl_file, show_col_types = FALSE) %>% mutate(qtl = "eQTL")
sqtl <- readr::read_csv(sqtl_file, show_col_types = FALSE) %>% mutate(qtl = "sQTL")

dat <- bind_rows(eqtl, sqtl)

stopifnot(rep_col %in% names(dat))

max_rep <- max(dat[[rep_col]], na.rm = TRUE)

plot_df <- dat %>%
  transmute(
    qtl,
    n_rep = as.integer(.data[[rep_col]])
  ) %>%
  count(qtl, n_rep, name = "n_features") %>%
  complete(qtl, n_rep = 0:max_rep, fill = list(n_features = 0)) %>%
  mutate(
    qtl   = factor(qtl, levels = c("eQTL", "sQTL")),
    n_rep = factor(n_rep, levels = 0:max_rep, ordered = TRUE)
  )

totals <- plot_df %>%
  group_by(qtl) %>%
  summarise(total = sum(n_features), .groups = "drop")
 
pal <- colorspace::sequential_hcl(n = length(levels(plot_df$n_rep)), palette = "Viridis")

#  base theme
base_theme <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.title.position = "plot",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    axis.text    = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.margin  = margin(8, 10, 8, 10)
  )

#  panel builder
make_panel <- function(which_qtl) {
  df_q  <- plot_df %>% filter(qtl == which_qtl)
  tot_q <- totals  %>% filter(qtl == which_qtl)
  
  ggplot(df_q, aes(x = qtl, y = n_features, fill = n_rep)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.4) +
    geom_text(
      data = tot_q,
      aes(x = qtl, y = total, label = total),
      inherit.aes = FALSE,
      vjust = -0.6,
      family = "Arial",
      size = 3.6
    ) +
    scale_fill_manual(values = pal, drop = FALSE, name = "Replications\n(# studies)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = NULL,
      y = "Number of TWAS-significant features",
      title = as.character(which_qtl)   # panel title
    ) +
    base_theme +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(
        hjust = 0.5, 
        face  = "bold",
        size  = 12
      )
    )
}

p_eqtl <- make_panel("eQTL")
p_sqtl <- make_panel("sQTL")

p <- (p_eqtl | p_sqtl) +
  plot_layout(guides = "collect") &  
  theme(legend.position = "right")

p <- p + plot_annotation(
  title = "Replication across independent studies",
  subtitle = paste0("Stacked by ", rep_col)
) & theme(
  plot.title = element_text(family = "Arial", face = "bold", size = 13),
  plot.subtitle = element_text(family = "Arial", size = 11)
)

print(p)

#  save
out_base <- paste0("replication_side_by_side_", rep_col)
ggsave(paste0(out_base, ".tiff"), p, width = 170, height = 95, units = "mm", dpi = 600, compression = "lzw")

# coloc
cis_kb <- 1000
gwas_n <- 31557
eqtl_n <- 150
gwas <- filename_gwas_primary
twas <- get_primary_twas("eqtl", "smultixcan", "brain_tissues")

# Load and clean GWAS
cols_keep <- c("chromosome","position","effect_size","standard_error","frequency")
dt_gwas <- fread(file.path(dir_gwas, gwas), select=cols_keep)
setnames(dt_gwas, "chromosome", "chr")
setnames(dt_gwas, old = c("position", "effect_size", "standard_error", "frequency"),
         new = c("pos", "beta", "se", "maf"))
dt_gwas[, chr := as.character(chr)]
dt_gwas[, chr := gsub("chr", "", chr)]
dt_gwas[, pos := as.numeric(pos)]
setkey(dt_gwas, chr, pos)

twas_sig_genes <- twas[twas$fdr_significant == TRUE,]
gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)
unique_chrs <- unique(twas_sig_genes$chromosome_name)

total_tasks <- nrow(twas_sig_genes) * length(gtex_brain_tissues)
results <- list()
gene_counter <- 0 # Tracks total gene-tissue combinations processed

for (chr_val in unique_chrs) {
  
  dt_gwas_chr <- dt_gwas[chr == chr_val]
  genes_on_chr <- twas_sig_genes[twas_sig_genes$chromosome_name == chr_val,]
  
  for (gtex_tissue in gtex_brain_tissues) {
    
    # LOAD & CLEAN eQTL DATA 
    gtex_tissue_file <- file.path(dir_gtex_eqtl,
                                  paste0("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_eQTL_all_associations_Brain_",gtex_tissue,".v8.EUR.allpairs.chr",chr_val,".parquet"))
    
    if (!file.exists(gtex_tissue_file)) { next }
    
    gtex_eqtl <- read_parquet(gtex_tissue_file)
    data.table::setDT(gtex_eqtl)
    
    gtex_eqtl[, c("chr", "pos", "ref", "alt", "build") := data.table::tstrsplit(variant_id, "_")]
    gtex_eqtl[, pos := as.numeric(pos)]
    gtex_eqtl[, chr := gsub("chr", "", chr)]
    gtex_eqtl <- gtex_eqtl[, .(phenotype_id, chr, pos, maf, pval_nominal, slope, slope_se)]
    gtex_eqtl[, chr := as.character(chr)]
    setkey(gtex_eqtl, chr, pos)
    
    for(i in 1:nrow(genes_on_chr)) {
      
      gene_counter <- gene_counter + 1
      gene_symbol <- genes_on_chr$gene[i]
      gene_name<- genes_on_chr$gene_name[i]
      
      gene_eqtl <- gtex_eqtl[phenotype_id == gene_symbol]
      if(nrow(gene_eqtl) == 0) { next }
      
      gene_start <- min(gene_eqtl$pos)
      gene_end <- max(gene_eqtl$pos)
      
      cis_start <- gene_start - cis_kb * 1000
      cis_end <- gene_end + cis_kb * 1000
      gwas_sub <- dt_gwas_chr[pos >= cis_start & pos <= cis_end]
      
      merged <- merge(gwas_sub, gene_eqtl, by = c("chr","pos"), suffixes = c(".x",".y"), all = FALSE)
      
      merged <- merged[!is.na(beta) & !is.na(se) & !is.na(slope) & !is.na(slope_se) & se != 0 & slope_se != 0]
      merged <- merged[!is.na(maf.x) & !is.na(maf.y) & maf.x > 0 & maf.y > 0 & maf.x < 1 & maf.y < 1] 
      merged[, snp := paste0(chr, ":", pos)]
      merged <- merged[, .SD[which.min(pval_nominal)], by = snp]
      merged <- merged[!is.infinite(beta) & !is.infinite(se) & !is.infinite(slope) & !is.infinite(slope_se)]
      
      if(nrow(merged) < 10) { next }
      
      dataset1 <- list(
        snp = merged$snp, beta = merged$beta, varbeta = merged$se^2, N = gwas_n,
        MAF = merged$maf.x, 
        type = "quant"
      )
      dataset2 <- list(
        snp = merged$snp, beta = merged$slope, varbeta = merged$slope_se^2, N = eqtl_n,
        MAF = merged$maf.y, 
        type = "quant"
      )
      coloc_res <- coloc.abf(dataset1, dataset2)
      
      res <- list(
        name = gene_name, gene = gene_symbol, chr = chr_val,
        start = gene_start, end = gene_end, tissue = gtex_tissue,
        summary = coloc_res$summary
      )
      
      gene_tissue <- paste0(gene_name,"_",gtex_tissue)
      results[[gene_tissue]] <- res
      
      pp4_prob <- format(coloc_res$summary["PP.H4.abf"], digits = 3, nsmall = 3)
      cat(paste0("[", gene_counter, "/", total_tasks, "] Chr ", chr_val, " | Gene: ", gene_name, " | Tissue: ", gtex_tissue," | PP4: ",pp4_prob,"\n"))
    }
    rm(gtex_eqtl)
    gc()
  }
  
  rm(dt_gwas_chr, genes_on_chr)
  gc()
}

coloc_summary_dt <- data.table::rbindlist(
  lapply(results, function(res_item) {
    data.table(
      gene_symbol = res_item$gene, gene_name = res_item$name, tissue = res_item$tissue,
      t(res_item$summary),
      chr = res_item$chr, start = res_item$start, end = res_item$end
    )
  }),
  fill = TRUE
)
coloc_summary_dt[, is_colocalized := PP.H4.abf > 0.8]
output_filename <- paste0("coloc",gwas,".csv")
data.table::fwrite(coloc_summary_dt, file = file.path(dir_results,"coloc",output_filename))

# coloc sqtl

cis_kb <- 1000
gwas_n <- 29400
sqtl_n <- 150
gwas <- "wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz"
twas <- df_twas_results[["wen_davatzikos_combined"]][["sqtl"]][["smultixcan"]][["brain_tissues"]]
dir_gtex_sqtl <- [[PLACEHOLDER]]
dir_gwas <- [[PLACEHOLDER]]
gtex_brain_tissues <- c(
    "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
    "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
    "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
    "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

# load and setup gwas
dt_gwas <- data.table::fread(
    file.path(dir_gwas, gwas),
    select = c("panel_variant_id","chromosome","position","pvalue","frequency","sample_size","imputation_status")
)
data.table::setnames(dt_gwas, old = c("chromosome","position","frequency"), new = c("chr_raw","pos","maf_gwas"))
dt_gwas[, chr := gsub("^chr", "", as.character(chr_raw))]
dt_gwas[, pos := as.numeric(pos)]
dt_gwas <- dt_gwas[!is.na(panel_variant_id) & !is.na(pos) & !is.na(pvalue) & !is.na(maf_gwas)]
data.table::setkey(dt_gwas, panel_variant_id)

# load twas
twas_sig <- twas[twas$fdr_significant == TRUE, , drop = FALSE]
data.table::setDT(twas_sig)

results <- list()

# for each chromosome
for (chr_val in unique(twas_sig$intron_chromosome_name))
{
    chr_clean <- gsub("^chr", "", as.character(chr_val))
    message("chromosome: ", chr_clean)
    
    message("subset gwas to this chromosome")
    dt_gwas_chr <- dt_gwas[chr == chr_clean]
    data.table::setkey(dt_gwas_chr, panel_variant_id)
    
    introns_on_chr <- twas_sig[intron_chromosome_name == chr_val]
    if (nrow(introns_on_chr) == 0) next
    
    # for each tissue
    for (gtex_tissue in gtex_brain_tissues) {
        message("tissue: ", gtex_tissue)
        
        gtex_file <- file.path(
            dir_gtex_sqtl,
            paste0(
                "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_sQTL_all_associations_Brain_",
                gtex_tissue,
                ".v8.EUR.sqtl_allpairs.chr",
                chr_clean,
                ".parquet"
            )
        )
        if (!file.exists(gtex_file)) next
        
        gtex_sqtl <- arrow::read_parquet(
            gtex_file,
            col_select = c("phenotype_id","variant_id","maf","pval_nominal","slope","slope_se")
        )
        data.table::setDT(gtex_sqtl)
        gtex_sqtl <- gtex_sqtl[!is.na(phenotype_id) & !is.na(variant_id) & !is.na(maf) & !is.na(pval_nominal)]
        
        # for each intron
        for (i in seq_len(nrow(introns_on_chr))) {
            
            intron_id <- as.character(introns_on_chr[["gene"]][i])
            intr_start <- introns_on_chr$intron_start[i]
            intr_end   <- introns_on_chr$intron_end[i]
            
            prefix <- paste0("chr", chr_clean, ":", intr_start, ":", intr_end, ":")
            
            gene_sqtl <- gtex_sqtl[startsWith(phenotype_id, prefix)]
            if (nrow(gene_sqtl) == 0L) next
            
            cis_start <- intr_start - cis_kb*1000
            cis_end   <- intr_end   + cis_kb*1000
            
            gene_sqtl[, pos := as.numeric(data.table::tstrsplit(variant_id, "_", fixed = TRUE, keep = 2L)[[1]])]
            gene_sqtl <- gene_sqtl[!is.na(pos)]
            gene_sqtl <- gene_sqtl[pos >= cis_start & pos <= cis_end]
            if (nrow(gene_sqtl) == 0L) next
            
            gene_sqtl[, panel_variant_id := variant_id]
            gene_sqtl[, maf_sqtl := maf]
            data.table::setkey(gene_sqtl, panel_variant_id)
            
            gwas_sub  <- dt_gwas_chr[pos >= cis_start & pos <= cis_end]
            if (nrow(gwas_sub) == 0L) next
            
            merged <- merge(gwas_sub, gene_sqtl, by = "panel_variant_id", all = FALSE)
            
            merged <- merged[
                !is.na(pvalue) & !is.na(maf_gwas) & maf_gwas > 0 & maf_gwas < 1 &
                    !is.na(slope) & !is.na(slope_se) & slope_se != 0 &
                    !is.na(maf_sqtl) & maf_sqtl > 0 & maf_sqtl < 1
                
                
                
            ]
            
            merged <- merged[, .SD[which.min(pval_nominal)], by = panel_variant_id]
            
            merged <- merged[!is.infinite(slope) & !is.infinite(slope_se)]
            
            merged[, pos2 := as.numeric(data.table::tstrsplit(panel_variant_id, "_", fixed = TRUE, keep = 2L)[[1]])]
            merged <- merged[!is.na(pos2)]
            if (nrow(merged) < 10) next
            
            lead_pos_gwas <- merged[which.min(pvalue)]$pos2[1]
            lead_pos_sqtl <- merged[which.min(pval_nominal)]$pos2[1]
            lead_pos <- if (is.finite(lead_pos_gwas)) lead_pos_gwas else lead_pos_sqtl
            if (!is.finite(lead_pos)) next
            
            win_kb <- 500
            win_start <- lead_pos - win_kb*1000
            win_end   <- lead_pos + win_kb*1000
            
            merged2 <- merged[pos2 >= win_start & pos2 <= win_end]
            if (nrow(merged2) < 50) {
                merged2 <- merged
            }
            if (nrow(merged) < 10) next
            
            dataset1 <- list(
                snp     = merged2$panel_variant_id,
                pvalues = merged2$pvalue,
                N       = 29400,
                MAF     = merged2$maf_gwas,
                type    = "quant"
            )
            
            dataset2 <- list(
                snp     = merged2$panel_variant_id,
                pvalues = merged2$pval_nominal,
                N       = sqtl_n,           
                MAF     = merged2$maf,      
                type    = "quant"
            )
            
            coloc_res <- coloc::coloc.abf(dataset1, dataset2)
            
            if (coloc_res$summary["PP.H4.abf"] > 0.8) {
                message("COLOC: ", gtex_tissue, " | chr", chr_clean, " | ", intron_id)
            }
            
            key_name <- paste0(intron_id, "_", gtex_tissue)
            results[[key_name]] <- list(
                intron_id = intron_id,
                chr = chr_clean,
                start = intr_start,
                end = intr_end,
                tissue = gtex_tissue,
                summary = coloc_res$summary
            )
        }
        
        rm(gtex_sqtl); gc()
    }
}

coloc_summary_dt <- data.table::rbindlist(
    lapply(results, function(res_item) {
        data.table::data.table(
            intron_id = res_item$intron_id,
            tissue = res_item$tissue,
            chr = res_item$chr,
            start = res_item$start,
            end = res_item$end,
            t(res_item$summary)
        )
    }),
    fill = TRUE
)

coloc_summary_dt[, is_colocalized := PP.H4.abf > 0.8]

dir.create(file.path(dir_results, "coloc"), recursive = TRUE, showWarnings = FALSE)
output_filename <- paste0("coloc_sqtl_", gwas, ".csv")
data.table::fwrite(coloc_summary_dt, file = file.path(dir_results, "coloc", output_filename))

# calculate max pp4 from coloc results

dir_coloc <- [[PLACEHOLDER]]
eqtl_file <- file.path(dir_coloc, "colocwen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz.csv")
sqtl_file <- file.path(dir_coloc, "coloc_sqtl_wen_davatzikos_combined_gwas_harmonised_imputed_merged_250kb.txt.gz.csv")

eqtl <- fread(file = eqtl_file)
sqtl <- fread(file = sqtl_file)

pick_id_col <- function(dt, kind = c("eqtl","sqtl")) {
  kind <- match.arg(kind)
  cand <- if (kind == "sqtl") {
    c("intron_id")
  } else {
    c("gene_name")
  }
}

summarise_coloc_A1A2 <- function(dt, kind = c("eqtl","sqtl")) {
  kind <- match.arg(kind)
  
  if (!("tissue" %in% names(dt))) stop(kind, ": missing tissue column")
  
  id_col  <- pick_id_col(dt, kind)
  pp4_col <- "PP.H4.abf"
  
  dt <- copy(dt)
  dt[, (pp4_col) := as.numeric(get(pp4_col))]
  dt[, id := as.character(get(id_col))]
  
  out <- dt[
    ,
    .(
      PP4_max = suppressWarnings(max(get(pp4_col), na.rm = TRUE)),
      tissue_max = tissue[which.max(get(pp4_col))][1],
      n_tissues_PP4_ge_0.8 = sum(get(pp4_col) >= 0.8, na.rm = TRUE),
      n_tissues_PP4_ge_0.5 = sum(get(pp4_col) >= 0.5, na.rm = TRUE)
    ),
    by = id
  ]
  
  out[is.infinite(PP4_max), `:=`(PP4_max = NA_real_, tissue_max = NA_character_)]
  
  setnames(out, "id", if (kind == "sqtl") "intron_id" else "gene")
  
  out[]
}

#  run A1/A2
eqtl_A1A2 <- summarise_coloc_A1A2(eqtl, "eqtl")
sqtl_A1A2 <- summarise_coloc_A1A2(sqtl, "sqtl")

#  write outputs
out_eqtl <- file.path(dir_coloc, "coloc_eqtl_A1A2_summary.csv")
out_sqtl <- file.path(dir_coloc, "coloc_sqtl_A1A2_summary.csv")

fwrite(eqtl_A1A2, file = out_eqtl)
fwrite(sqtl_A1A2, file = out_sqtl)

cat("Wrote:\n", out_eqtl, "\n", out_sqtl, "\n", sep = "")

# ternary 1
#install.packages("plotly")
library(plotly)

# Read coloc results
coloc <- read.csv([[PLACEHOLDER]], stringsAsFactors = FALSE)

# Compute ternary coordinates
coloc$PP.null_or_single <- coloc$PP.H0.abf + coloc$PP.H1.abf + coloc$PP.H2.abf
coloc$PP.independent    <- coloc$PP.H3.abf
coloc$PP.coloc          <- coloc$PP.H4.abf

coloc$hover <- paste0(
  "Gene: ", coloc$name,
  "<br>Tissue: ", coloc$tissue,
  "<br>PP.H4: ", signif(coloc$PP.H4.abf, 3)
)

# interactive ternary plot
plot_ly(
  data = coloc,
  type = 'scatterternary',
  a = ~PP.null_or_single,
  b = ~PP.independent,
  c = ~PP.coloc,
  mode = 'markers',
  marker = list(
    size = 7,
    color = ~PP.H4.abf,
    colorscale = 'RdYlBu',
    showscale = TRUE,
    opacity = 0.8
  ),
  text = ~hover,
  hoverinfo = 'text'
) %>%
  layout(
    ternary = list(
      sum = 1,
      aaxis = list(title = "H0 + H1 + H2 (no/single)"),
      baxis = list(title = "H3 (independent)"),
      caxis = list(title = "H4 (colocalised)")
    ),
    margin = list(l = 100, r = 100, b = 100, t = 100, pad = 20)
  )

# Interactive ternary plot coloured by tissue
plot_ly(
  data = coloc,
  type = 'scatterternary',
  a = ~PP.null_or_single,
  b = ~PP.independent,
  c = ~PP.coloc,
  color = ~is_colocalized,
  mode = 'markers',
  marker = list(size = 7, opacity = 0.5, symbol = 'circle'),
  text = ~hover,
  hoverinfo = 'text'
) %>%
  layout(
    ternary = list(
      sum = 1,
      aaxis = list(title = "H0 + H1 + H2 (no/single)"),
      baxis = list(title = "H3 (independent)"),
      caxis = list(title = "H4 (colocalised)")
    ),
    legend = list(title = list(text = "<b>Tissue</b>")),
    margin = list(l = 100, r = 100, b = 100, t = 100, pad = 20)
  )


# ternary 2
coloc <- read.csv(
  [[PLACEHOLDER]],
  stringsAsFactors = FALSE
)

# Compute ternary coordinates 
coloc <- coloc %>%
  mutate(
    PP.null_or_single = PP.H0.abf + PP.H1.abf + PP.H2.abf,
    PP.independent    = PP.H3.abf,
    PP.coloc          = PP.H4.abf,
    hover = paste0(
      "<b>Feature:</b> ", intron_id,
      "<br><b>Tissue:</b> ", tissue,
      "<br><b>PP.H4:</b> ", signif(PP.H4.abf, 3),
      "<br><b>H0+H1+H2:</b> ", signif(PP.null_or_single, 3),
      "<br><b>H3:</b> ", signif(PP.independent, 3)
    )
  )

NATURE_FONT   <- "Arial"
TITLE_SIZE    <- 16
AXIS_SIZE     <- 18
TICK_SIZE     <- 16

GRIDCOL       <- "rgba(0,0,0,0.12)"
AXISCOL       <- "rgba(0,0,0,0.55)"

PT_SIZE       <- 7
PT_OPACITY    <- 0.80
PT_LINE_COL   <- "rgba(0,0,0,0.35)"
PT_LINE_W     <- 0

H4_COLORSCALE <- "Viridis"

nature_ternary_layout <- function(title_text = NULL) {
  list(
    title = list(
      text = title_text,
      x = 0.02, xanchor = "left",
      font = list(family = NATURE_FONT, size = TITLE_SIZE, color = "black")
    ),
    paper_bgcolor = "white",
    plot_bgcolor  = "white",
    font = list(family = NATURE_FONT, size = TICK_SIZE, color = "black"),
    ternary = list(
      sum = 1,
      bgcolor = "white",
      aaxis = list(
        title = list(text = "H0 + H1 + H2 (no/single)", font = list(family = NATURE_FONT, size = AXIS_SIZE)),
        showgrid = TRUE, gridcolor = GRIDCOL, gridwidth = 1,
        showline = TRUE, linecolor = AXISCOL, linewidth = 1,
        ticks = "outside", ticklen = 4, tickcolor = AXISCOL,
        tickfont = list(family = NATURE_FONT, size = TICK_SIZE),
        min = 0
      ),
      baxis = list(
        title = list(text = "H3 (independent)", font = list(family = NATURE_FONT, size = AXIS_SIZE)),
        showgrid = TRUE, gridcolor = GRIDCOL, gridwidth = 1,
        showline = TRUE, linecolor = AXISCOL, linewidth = 1,
        ticks = "outside", ticklen = 4, tickcolor = AXISCOL,
        tickfont = list(family = NATURE_FONT, size = TICK_SIZE),
        min = 0
      ),
      caxis = list(
        title = list(text = "H4 (colocalised)", font = list(family = NATURE_FONT, size = AXIS_SIZE)),
        showgrid = TRUE, gridcolor = GRIDCOL, gridwidth = 1,
        showline = TRUE, linecolor = AXISCOL, linewidth = 1,
        ticks = "outside", ticklen = 4, tickcolor = AXISCOL,
        tickfont = list(family = NATURE_FONT, size = TICK_SIZE),
        min = 0
      )
    ),
    margin = list(l = 70, r = 30, b = 55, t = ifelse(is.null(title_text), 30, 60), pad = 4),
    legend = list(
      bgcolor = "rgba(255,255,255,0.85)",
      bordercolor = "rgba(0,0,0,0.15)",
      borderwidth = 1,
      font = list(family = NATURE_FONT, size = 11),
      orientation = "v"
    )
  )
}

# Plot 1: coloured by PP.H4 (continuous)
p_h4 <- plot_ly(
  data = coloc,
  type = "scatterternary",
  a = ~PP.null_or_single,
  b = ~PP.independent,
  c = ~PP.coloc,
  mode = "markers",
  marker = list(
    size = PT_SIZE,
    opacity = PT_OPACITY,
    color = ~PP.H4.abf,
    colorscale = H4_COLORSCALE,
    reversescale = FALSE,
    showscale = TRUE,
    line = list(color = PT_LINE_COL, width = PT_LINE_W),
    colorbar = list(
      title = list(text = "PP.H4", font = list(family = NATURE_FONT, size = 12)),
      tickfont = list(family = NATURE_FONT, size = 10),
      len = 0.75,
      outlinewidth = 0.5
    )
  )
) %>%
  layout(nature_ternary_layout("COLOC posterior mass (ternary): coloured by PP.H4"),
    ternary = list(
      sum = 1,
      aaxis = list(title = "H0 + H1 + H2 (no/single)"),
      baxis = list(title = "H3 (independent)"),
      caxis = list(title = "H4 (colocalised)")
    ),
    legend = list(title = list(text = "<b>Tissue</b>")),
    margin = list(l = 100, r = 100, b = 100, t = 100, pad = 20)
  )

p_h4


plotly::install_kaleido()


## tern 3

install.packages("ggtern")
library(dplyr)
library(ggplot2)
library(ggtern)
library(viridis)
library(ragg)

colocfile=[[PLACEHOLDER]]
coloc <- read.csv(
  colocfile,
  stringsAsFactors = FALSE
) %>%
  mutate(
    PP.null_or_single = PP.H0.abf + PP.H1.abf + PP.H2.abf,
    PP.independent    = PP.H3.abf,
    PP.coloc          = PP.H4.abf
  )

p_static <- ggtern(
  data = coloc,
  aes(x = PP.null_or_single, y = PP.independent, z = PP.coloc, colour = PP.H4.abf)
) +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_viridis_c(name = "PP.H4", option = "viridis", end = 0.98) +
  labs(
    title = "COLOC posterior mass (ternary): coloured by PP.H4",
    T = "H0 + H1 + H2 (no/single)",
    L = "H3 (independent)",
    R = "H4 (colocalised)"
  ) +
  theme_bw(base_family = "Arial", base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "plain", hjust = 0),
    tern.axis.title = element_text(size = 14),
    tern.axis.text  = element_text(size = 12, colour = "black"),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

# Save as TIFF
out_file <- [[PLACEHOLDER]]

ragg::agg_tiff(
  filename = out_file,
  width = 180, height = 160, units = "mm", res = 600, compression = "lzw"
)
print(p_static)
dev.off()

# coloc figures - pp4 heatmap

# Heatmap
#BiocManager::install("ComplexHeatmap")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

dir_coloc <- [[PLACEHOLDER]]

eqtl_coloc_file <- file.path(dir_coloc, "colocwen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz.csv")
sqtl_coloc_file <- file.path(dir_coloc, "coloc_sqtl_wen_davatzikos_combined_gwas_harmonised_imputed_merged_250kb.txt.gz.csv")

eqtl_A1A2_file  <- file.path(dir_coloc, "coloc_eqtl_A1A2_summary.csv")
sqtl_A1A2_file  <- file.path(dir_coloc, "coloc_sqtl_A1A2_summary.csv")

# GTEx brain tissues order 
gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

pick_pp4_col <- function(dt) {
  cand <- c("PP.H4.abf","PP.H4","PP4","pp4","PP4.abf")
  hit <- cand[cand %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  rx <- grep("^PP\\.H4(\\.abf)?$", names(dt), value = TRUE)
  if (length(rx)) return(rx[1])
  stop("Could not find PP4 column (expected PP.H4.abf / PP.H4).")
}

pick_id_col <- function(dt, kind = c("eqtl","sqtl")) {
  kind <- match.arg(kind)
  cand <- if (kind == "eqtl") c("gene_name") else c("intron_id")
  hit <- cand[cand %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  stop(kind, ": could not find an ID column. Columns: ", paste(names(dt), collapse=", "))
}

# - load --
eqtl <- fread(file = eqtl_coloc_file)
sqtl <- fread(file = sqtl_coloc_file)

eqtl_A1A2 <- fread(file = eqtl_A1A2_file)
sqtl_A1A2 <- fread(file = sqtl_A1A2_file)

#  build wide PP4 matrices -
make_pp4_matrix <- function(dt, kind = c("eqtl","sqtl"), tissue_levels) {
  kind <- match.arg(kind)
  id_col  <- pick_id_col(dt, kind)
  pp4_col <- pick_pp4_col(dt)
  
  x <- dt[, .(feature = as.character(get(id_col)),
              tissue  = as.character(tissue),
              pp4     = as.numeric(get(pp4_col)))]
  
  x <- x[tissue %in% tissue_levels]
  x[, tissue := factor(tissue, levels = tissue_levels)]
  
  x <- x[, .(pp4 = max(pp4, na.rm = TRUE)), by = .(feature, tissue)]
  
  wide <- dcast(x, feature ~ tissue, value.var = "pp4")
  mat <- as.matrix(wide[, -1, with = FALSE])
  rownames(mat) <- wide$feature
  mat
}

mat_eqtl <- make_pp4_matrix(eqtl, "eqtl", gtex_brain_tissues)
mat_sqtl <- make_pp4_matrix(sqtl, "sqtl", gtex_brain_tissues)

#order rows by PP4_max
order_by_pp4max <- function(mat, A1A2, kind = c("eqtl","sqtl")) {
  kind <- match.arg(kind)
  idcol <- if (kind == "eqtl") {
    if ("gene" %in% names(A1A2)) "gene" else names(A1A2)[1]
  } else {
    if ("intron_id" %in% names(A1A2)) "intron_id" else names(A1A2)[1]
  }
  
  A <- copy(A1A2)
  A[, feature := as.character(get(idcol))]
  A <- A[feature %in% rownames(mat)]
  setorder(A, -PP4_max)
  
  mat[match(A$feature, rownames(mat)), , drop = FALSE]
}

mat_eqtl_ord <- order_by_pp4max(mat_eqtl, eqtl_A1A2, "eqtl")
mat_sqtl_ord <- order_by_pp4max(mat_sqtl, sqtl_A1A2, "sqtl")

# - Panel A: Heatmaps 
ht_scale <- colorRamp2(c(0, 0.5, 0.8, 1.0), c("white", "#D9D9D9", "#6BAED6", "#08519C"))

save_single_heatmap <- function(mat, title, out_file) {
  
  col_fun <- colorRamp2(
    c(0, 0.5, 0.8, 1.0),
    c("white", "#D9D9D9", "#6BAED6", "#08519C")
  )
  
  ht <- Heatmap(
    mat,
    name = "PP4",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_title = title,
    column_title_gp = gpar(fontface = "bold"),
    heatmap_legend_param = list(
      title = "PP4",
      at = c(0, 0.5, 0.8, 1.0)
    )
  )
  
  tiff(
    filename = out_file,
    width = 180, height = 220, units = "mm",
    res = 300, compression = "lzw"
  )
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}

save_single_heatmap(
  mat_eqtl_ord,
  title = "eQTL colocalisation (PP4)",
  out_file = file.path(dir_coloc, "FigureX_PP4_heatmap_eQTL.tiff")
)

save_single_heatmap(
  mat_sqtl_ord,
  title = "sQTL colocalisation (PP4)",
  out_file = file.path(dir_coloc, "FigureX_PP4_heatmap_sQTL.tiff")
)


# SMR HEIDI

smr_bin  <- [[PLACEHOLDER]] 
bfile_ld_prefix <- [[PLACEHOLDER]]  
gwas_in  <- [[PLACEHOLDER]] 
besd_in_dir  <- [[PLACEHOLDER]]
twas_in <- [[PLACEHOLDER]] 
out_pref <- [[PLACEHOLDER]] 
threads  <- 8   
gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

#  Read S-PrediXcan TWAS results 
twas <- df_twas_results[["wen_davatzikos_combined"]][["eqtl"]][["smultixcan"]][["brain_tissues"]]
twas_sig <- twas[twas$fdr_significant == TRUE, , drop = FALSE]
data.table::setDT(twas_sig)
message(nrow(twas_sig), " TWAS genes selected (pvalue < 1e-5).")

# Prepare/validate GWAS summary for SMR
gwas <- fread(gwas_in, header = TRUE, nrows = 5)
print(names(gwas))
#  map  column names -> the SMR required ones 
cols_to_read <- c("variant_id","effect_allele","non_effect_allele",
                  "frequency","effect_size","standard_error","pvalue","n_cases")
full_gwas <- fread(gwas_in, select = cols_to_read)
# Filter and clean first
full_gwas2 <- full_gwas[
  grepl("^rs", variant_id, ignore.case = TRUE) &
    !is.na(effect_size) & !is.na(standard_error) &
    !is.na(n_cases) & !is.na(pvalue) &
    effect_allele %in% c("A","C","G","T") &
    non_effect_allele %in% c("A","C","G","T")]
# Fill missing frequency
full_gwas2[, frequency := fifelse(is.na(frequency), 0.5, frequency)]
# Sort by SNP and p-value, then keep first occurrence of each SNP
setorder(full_gwas2, variant_id, pvalue)
full_gwas2 <- unique(full_gwas2, by = "variant_id")
# Rename columns to SMR format
setnames(full_gwas2,
         c("variant_id","effect_allele","non_effect_allele","frequency",
           "effect_size","standard_error","pvalue","n_cases"),
         c("SNP","A1","A2","freq","b","se","p","n"))
gwas_formatted <- paste0(out_pref, "primary.gwas.formatted.txt") 
fwrite(full_gwas2, gwas_formatted, sep = "\t", quote = FALSE, na = "NA")

# run smr 13x
for (tissue in gtex_brain_tissues)
{
    besd_in <- file.path(besd_in_dir,paste0("Brain_",tissue,".lite"))

    smr_cmd <- paste0(
      shQuote(smr_bin),
      " --bfile ", shQuote(bfile_ld_prefix),
      " --gwas-summary ", shQuote(gwas_formatted),
      " --beqtl-summary ", shQuote(besd_in),
      " --peqtl-smr ", "1e-4",
      " --heidi-min-m ", "3",
      " --extract-probe ", shQuote([[PLACEHOLDER]]),
      " --out ", shQuote(paste0(out_pref,tissue)),
      " --thread-num ", threads
    )
    cmd_win <- chartr("/", "\\", smr_cmd)
    rc <- system(smr_cmd)
}

# SMR HEIDI SQTL

smr_bin  <- [[PLACEHOLDER]] 
bfile_ld_prefix <- [[PLACEHOLDER]]  
gwas_in  <- [[PLACEHOLDER]] 
besd_in_dir  <- [[PLACEHOLDER]]
twas_in <- [[PLACEHOLDER]] 
out_pref <- [[PLACEHOLDER]] 
threads  <- 8   
gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

#  Read S-PrediXcan TWAS results
twas <- df_twas_results[["wen_davatzikos_combined"]][["sqtl"]][["smultixcan"]][["brain_tissues"]]
twas_sig <- twas[twas$fdr_significant == TRUE, , drop = FALSE]
data.table::setDT(twas_sig)
message(nrow(twas_sig), " TWAS genes selected (pvalue < 1e-5).")

# Prepare/validate GWAS summary for SMR 
gwas <- fread(gwas_in, header = TRUE, nrows = 5)
print(names(gwas))
#  map  column names the SMR required ones 
cols_to_read <- c("variant_id","effect_allele","non_effect_allele",
                  "frequency","effect_size","standard_error","pvalue","n_cases")
full_gwas <- fread(gwas_in, select = cols_to_read)
# Filter and clean first
full_gwas2 <- full_gwas[
    grepl("^rs", variant_id, ignore.case = TRUE) &
        !is.na(effect_size) & !is.na(standard_error) &
        !is.na(n_cases) & !is.na(pvalue) &
        effect_allele %in% c("A","C","G","T") &
        non_effect_allele %in% c("A","C","G","T")]
# Fill missing frequency
full_gwas2[, frequency := fifelse(is.na(frequency), 0.5, frequency)]
# Sort by SNP and p-value, then keep first occurrence of each SNP
setorder(full_gwas2, variant_id, pvalue)
full_gwas2 <- unique(full_gwas2, by = "variant_id")
# Rename columns to SMR format
setnames(full_gwas2,
         c("variant_id","effect_allele","non_effect_allele","frequency",
           "effect_size","standard_error","pvalue","n_cases"),
         c("SNP","A1","A2","freq","b","se","p","n"))
gwas_formatted <- paste0(out_pref, "primary.gwas.formatted.txt") 
fwrite(full_gwas2, gwas_formatted, sep = "\t", quote = FALSE, na = "NA")

for (tissue in gtex_brain_tissues)
{
  besd_in <- file.path(besd_in_dir,paste0("sQTL_Brain_",tissue,".lite"))
 
  smr_cmd <- paste0(
    shQuote(smr_bin),
    " --bfile ", shQuote(bfile_ld_prefix),
    " --gwas-summary ", shQuote(gwas_formatted),
    " --beqtl-summary ", shQuote(besd_in),
    " --peqtl-smr ", "1e-4",
    " --heidi-min-m ", "3",
    " --extract-probe ", shQuote([[PLACEHOLDER]]),
    " --out ", shQuote(paste0(out_pref,tissue)),
    " --thread-num ", threads
  )
  cmd_win <- chartr("/", "\\", smr_cmd)
  
  rc <- system(smr_cmd)

}


# summarise SMR HEIDI

    library(data.table)
    library(dplyr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(purrr)

# Folder containing all SMR outputs 
dir_smr <- [[PLACEHOLDER]] 
dir_smr_sqtl <- [[PLACEHOLDER]] 

# Output folder
out_dir <- file.path(dir_smr, "summary_integrated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Thresholds
p_smr_thresh   <- 5e-8
p_heidi_thresh <- 0.01
do_fisher_combine <- TRUE
msg <- function(...) cat(sprintf(...), "\n")

get_tissue_from_filename <- function(f) {
    b <- basename(f)
    b <- sub("\\.smr$", "", b, ignore.case = TRUE)
    b <- sub("^sqtl", "", b, ignore.case = TRUE)
    b
}

read_smr_one <- function(path, qtl, tissue) {
    dt <- suppressWarnings(data.table::fread(path, sep = "\t", data.table = TRUE, showProgress = FALSE))
    
    required <- c("probeID","Gene","topSNP","p_SMR","p_HEIDI","b_SMR","se_SMR")
    miss <- setdiff(required, names(dt))
    
    dt[, `:=`(
        qtl = qtl,
        tissue = tissue,
        p_SMR   = as.numeric(p_SMR),
        p_HEIDI = as.numeric(p_HEIDI),
        b_SMR   = as.numeric(b_SMR),
        se_SMR  = as.numeric(se_SMR),
        probeID = as.character(probeID),
        Gene    = as.character(Gene),
        topSNP  = as.character(topSNP)
    )]
    
    dt[, probeID_noversion := sub("\\.\\d+$", "", probeID)]
    
    dt[, intron_key := {
        x <- probeID
        m <- str_match(x, "^(chr)?([0-9XYM]+):([0-9]+):([0-9]+):")
        ifelse(is.na(m[,1]), NA_character_, paste0("chr", m[,3], ":", m[,4], ":", m[,5]))
    }]
    
    dt
}

fisher_combine_p <- function(pvals) {
    p <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
    if (length(p) == 0) return(NA_real_)
    stat <- -2 * sum(log(p))
    df <- 2 * length(p)
    pchisq(stat, df = df, lower.tail = FALSE)
}

twas_eqtl <- df_twas_results[["wen_davatzikos_combined"]][["eqtl"]][["smultixcan"]][["brain_tissues"]]
twas_sqtl <- df_twas_results[["wen_davatzikos_combined"]][["sqtl"]][["smultixcan"]][["brain_tissues"]]

prep_twas_eqtl <- function(df) {
    df <- as.data.frame(df)
    stopifnot("fdr_significant" %in% names(df))
    df %>%
        filter(.data$fdr_significant == TRUE) %>%
        mutate(
            feature_id = as.character(.data$gene_noversion %||% sub("\\.\\d+$","", .data$gene))
        ) %>%
        dplyr::select(any_of(c("feature_id","gene","gene_name","gene_noversion","pvalue","pvalue_fdr","z_mean","n","n_indep"))) %>%
        distinct()
}

prep_twas_sqtl <- function(df) {
    df <- as.data.frame(df)
    stopifnot("fdr_significant" %in% names(df))
    
    if (!all(c("intron_chromosome_name","intron_start_position","intron_end_position") %in% names(df))) {
        stop("twas_sqtl is missing intron coordinate columns needed for joining to SMR.")
    }
    
    df %>%
        filter(.data$fdr_significant == TRUE) %>%
        mutate(
            intron_chr = as.character(.data$intron_chromosome_name),
            intron_chr = ifelse(str_detect(intron_chr, "^chr"), intron_chr, paste0("chr", intron_chr)),
            intron_start = as.integer(.data$intron_start_position),
            intron_end   = as.integer(.data$intron_end_position),
            feature_id   = as.character(.data$gene_noversion %||% sub("\\.\\d+$","", .data$gene)),  # gene for reporting
            intron_key   = paste0(intron_chr, ":", intron_start, ":", intron_end)
        ) %>%
        dplyr::select(any_of(c("intron_key","feature_id","gene","gene_name","gene_noversion",
                        "pvalue","pvalue_fdr","z_mean","n","n_indep",
                        "intron_chromosome_name","intron_start_position","intron_end_position"))) %>%
        distinct()
}

`%||%` <- function(a, b) if (!is.null(a) && !all(is.na(a))) a else b

twas_eqtl_sig <- prep_twas_eqtl(twas_eqtl)
twas_sqtl_sig <- prep_twas_sqtl(twas_sqtl)

msg("TWAS FDR-significant hits: eQTL=%d genes; sQTL=%d introns", nrow(twas_eqtl_sig), nrow(twas_sqtl_sig))

# LOAD SMR FILES

all_files <- list.files(dir_smr, pattern = "\\.smr$", full.names = TRUE)
if (length(all_files) == 0) stop(sprintf("No .smr files found in %s", dir_smr))
eqtl_files <- all_files[!str_detect(basename(all_files), "^sqtl")]


all_files_sqtl <- list.files(dir_smr_sqtl, pattern = "\\.smr$", full.names = TRUE)
if (length(all_files_sqtl) == 0) stop(sprintf("No .smr files found in %s", dir_smr))
sqtl_files <- all_files_sqtl[!str_detect(basename(all_files_sqtl), "^sqtl")]

if (length(eqtl_files) == 0 || length(sqtl_files) == 0) {
    stop("Did not detect both eQTL and sQTL SMR files. Expect eqtl: *.smr and sqtl: sqtl*.smr")
}

msg("Found SMR files: eQTL=%d, sQTL=%d", length(eqtl_files), length(sqtl_files))

smr_eqtl_all <- rbindlist(lapply(eqtl_files, function(f) {
    read_smr_one(f, qtl="eqtl", tissue=get_tissue_from_filename(f))
}), use.names = TRUE, fill = TRUE)

smr_sqtl_all <- rbindlist(lapply(sqtl_files, function(f) {
    read_smr_one(f, qtl="sqtl", tissue=get_tissue_from_filename(f))
}), use.names = TRUE, fill = TRUE)

msg("Loaded SMR rows: eQTL=%d; sQTL=%d", nrow(smr_eqtl_all), nrow(smr_sqtl_all))

# MERGE + SUMMARISE

summarise_eqtl <- function(smr_dt, twas_sig) {
    smr_hits <- smr_dt %>%
        as_tibble() %>%
        inner_join(twas_sig, by = c("probeID_noversion" = "feature_id"))
    
    msg("After merge (eQTL): %d rows (tissue x gene)", nrow(smr_hits))
    
    per_tissue <- smr_hits %>%
        group_by(tissue) %>%
        summarise(
            n_genes_with_result = n_distinct(probeID_noversion),
            n_rows = n(),
            n_smr_sig = sum(p_SMR < p_smr_thresh, na.rm = TRUE),
            n_heidi_pass = sum(p_HEIDI >= p_heidi_thresh, na.rm = TRUE),
            n_smr_sig_and_heidi_pass = sum((p_SMR < p_smr_thresh) & (p_HEIDI >= p_heidi_thresh), na.rm = TRUE),
            .groups = "drop"
        ) %>% arrange(desc(n_smr_sig_and_heidi_pass), desc(n_smr_sig))
    
    feature_summary <- smr_hits %>%
        group_by(probeID_noversion) %>%
        summarise(
            n_tissues_with_result = n_distinct(tissue),
            min_p_SMR = min(p_SMR, na.rm = TRUE),
            best_tissue = tissue[which.min(p_SMR)][1],
            best_p_HEIDI = p_HEIDI[which.min(p_SMR)][1],
            best_b_SMR = b_SMR[which.min(p_SMR)][1],
            best_se_SMR = se_SMR[which.min(p_SMR)][1],
            n_smr_sig = sum(p_SMR < p_smr_thresh, na.rm = TRUE),
            n_heidi_pass = sum(p_HEIDI >= p_heidi_thresh, na.rm = TRUE),
            n_smr_sig_and_heidi_pass = sum((p_SMR < p_smr_thresh) & (p_HEIDI >= p_heidi_thresh), na.rm = TRUE),
            fisher_p_SMR = if (do_fisher_combine) fisher_combine_p(p_SMR) else NA_real_,
            .groups = "drop"
        ) %>%
        left_join(twas_sig, by = c("probeID_noversion" = "feature_id")) %>%
        arrange(min_p_SMR) %>%
        rename(gene_id = probeID_noversion)
    
    list(long = smr_hits, per_tissue = per_tissue, feature = feature_summary)
}

summarise_sqtl <- function(smr_dt, twas_sig) {
    smr_hits <- smr_dt %>%
        as_tibble() %>%
        filter(!is.na(intron_key)) %>%
        inner_join(twas_sig, by = "intron_key")
    
    msg("After merge (sQTL): %d rows (tissue x intron)", nrow(smr_hits))
    
    per_tissue <- smr_hits %>%
        group_by(tissue) %>%
        summarise(
            n_introns_with_result = n_distinct(intron_key),
            n_rows = n(),
            n_smr_sig = sum(p_SMR < p_smr_thresh, na.rm = TRUE),
            n_heidi_pass = sum(p_HEIDI >= p_heidi_thresh, na.rm = TRUE),
            n_smr_sig_and_heidi_pass = sum((p_SMR < p_smr_thresh) & (p_HEIDI >= p_heidi_thresh), na.rm = TRUE),
            .groups = "drop"
        ) %>% arrange(desc(n_smr_sig_and_heidi_pass), desc(n_smr_sig))
    
    feature_summary <- smr_hits %>%
        group_by(intron_key) %>%
        summarise(
            n_tissues_with_result = n_distinct(tissue),
            min_p_SMR = min(p_SMR, na.rm = TRUE),
            best_tissue = tissue[which.min(p_SMR)][1],
            best_p_HEIDI = p_HEIDI[which.min(p_SMR)][1],
            best_b_SMR = b_SMR[which.min(p_SMR)][1],
            best_se_SMR = se_SMR[which.min(p_SMR)][1],
            n_smr_sig = sum(p_SMR < p_smr_thresh, na.rm = TRUE),
            n_heidi_pass = sum(p_HEIDI >= p_heidi_thresh, na.rm = TRUE),
            n_smr_sig_and_heidi_pass = sum((p_SMR < p_smr_thresh) & (p_HEIDI >= p_heidi_thresh), na.rm = TRUE),
            fisher_p_SMR = if (do_fisher_combine) fisher_combine_p(p_SMR) else NA_real_,
            .groups = "drop"
        ) %>%
        left_join(twas_sig, by = "intron_key") %>%
        arrange(min_p_SMR)
    
    list(long = smr_hits, per_tissue = per_tissue, feature = feature_summary)
}

res_eqtl <- summarise_eqtl(smr_eqtl_all, twas_eqtl_sig)
res_sqtl <- summarise_sqtl(smr_sqtl_all, twas_sqtl_sig)

# WRITE OUTPUTS

write_csv_safe <- function(df, path) {
    readr::write_csv(df, path)
    msg("Wrote: %s (%d rows)", path, nrow(df))
}

write_csv_safe(res_eqtl$long,       file.path(out_dir, "eqtl_smr_merged_long.csv"))
write_csv_safe(res_eqtl$per_tissue, file.path(out_dir, "eqtl_smr_per_tissue_summary.csv"))
write_csv_safe(res_eqtl$feature,    file.path(out_dir, "eqtl_smr_integrated_feature_summary.csv"))

write_csv_safe(res_sqtl$long,       file.path(out_dir, "sqtl_smr_merged_long.csv"))
write_csv_safe(res_sqtl$per_tissue, file.path(out_dir, "sqtl_smr_per_tissue_summary.csv"))
write_csv_safe(res_sqtl$feature,    file.path(out_dir, "sqtl_smr_integrated_feature_summary.csv"))

msg("DONE.")
msg("Thresholds: p_SMR < %g ; p_HEIDI >= %g ; fisher=%s",
    p_smr_thresh, p_heidi_thresh, ifelse(do_fisher_combine, "ON", "OFF"))



# SMR HEIDI volcano plots

    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(stringr)

# Paths
input_dir  <- [[PLACEHOLDER]] 
input_file <- file.path(input_dir, "sqtl_smr_merged_long.csv")
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
df <- fread(input_file)

# Basic checks
required_cols <- c("b_SMR", "p_SMR", "p_HEIDI", "Gene", "gene_name", "qtl", "tissue")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_plot <- df %>%
    mutate(
        b_SMR   = as.numeric(b_SMR),
        p_SMR   = as.numeric(p_SMR),
        p_HEIDI = as.numeric(p_HEIDI),
        label   = dplyr::coalesce(gene_name, Gene),
        heidi_pass = !is.na(p_HEIDI) & p_HEIDI >= 0.01,
        neglog10_p_SMR = -log10(p_SMR)
    ) %>%
    filter(
        !is.na(b_SMR),
        !is.na(p_SMR),
        p_SMR > 0,
        is.finite(neglog10_p_SMR)
    )

cat("Rows in input: ", nrow(df), "\n")
cat("Rows plotted:  ", nrow(df_plot), "\n")
cat("HEIDI pass:    ", sum(df_plot$heidi_pass, na.rm = TRUE), "\n")

# Significance threshold
p_thresh <- 0.05 / nrow(df_plot)
y_thresh <- -log10(p_thresh)

df_plot <- df_plot %>%
    mutate(
        smr_sig = p_SMR < p_thresh,
        smr_sig_heidi = smr_sig & heidi_pass
    )
n_label <- 12

df_labels <- df_plot %>%
    filter(heidi_pass) %>%
    arrange(p_SMR) %>%
    slice_head(n = n_label)

theme_nature <- function(base_size = 7, base_family = "sans") {
    theme_classic(base_size = base_size, base_family = base_family) +
        theme(
            axis.line = element_line(linewidth = 0.4, colour = "black"),
            axis.ticks = element_line(linewidth = 0.4, colour = "black"),
            axis.ticks.length = unit(1.5, "mm"),
            axis.title = element_text(size = base_size),
            axis.text = element_text(size = base_size, colour = "black"),
            plot.title = element_text(size = base_size, face = "plain", hjust = 0),
            plot.margin = margin(4, 4, 4, 4)
        )
}

# Plot
p <- ggplot(df_plot, aes(x = b_SMR, y = neglog10_p_SMR)) +
    
    geom_point(
        data = df_plot %>% filter(!smr_sig_heidi),
        colour = "grey70",
        size = 2,
        alpha = 0.75,
        stroke = 0
    ) +
    
    geom_point(
        data = df_plot %>% filter(smr_sig_heidi),
        colour = "#D62728",
        size = 2,
        alpha = 0.9,
        stroke = 0
    ) +
    
    geom_hline(
        yintercept = y_thresh,
        linetype = "dashed",
        linewidth = 0.35,
        colour = "black"
    ) +
    
    geom_vline(
        xintercept = 0,
        linetype = "solid",
        linewidth = 0.3,
        colour = "black"
    ) +
    
    labs(
        x = expression(italic("b")["SMR"]),
        y = expression(-log[10](italic("p")["SMR"]))
    ) +
    
    theme_nature(base_size = 7)

# Save outputs
plot_width_in  <- 3.35
plot_height_in <- 2.8

print(p)
ggsave(
    filename = file.path(output_dir, "smr_volcano_heidi_pass.tiff"),
    plot = p,
    width = plot_width_in,
    height = plot_height_in,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

# eqtl sqtl integration/overlap

    library(dplyr)
    library(GenomicRanges)
    library(IRanges)

eq <- df_twas_results[["wen_davatzikos_combined"]][["eqtl"]][["smultixcan"]][["brain_tissues"]]
sq <- df_twas_results[["wen_davatzikos_combined"]][["sqtl"]][["smultixcan"]][["brain_tissues"]]

use_significant_only <- TRUE 
cis_window_bp <- 1e6 

out_dir <- [[PLACEHOLDER]] 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (use_significant_only) {
    eq <- eq %>% filter(fdr_significant == TRUE)
    sq <- sq %>% filter(fdr_significant == TRUE)
}

message("eQTL rows used: ", nrow(eq))
message("sQTL rows used: ", nrow(sq))

# Build GRanges eqtl
gene_gr <- GRanges(
    seqnames = eq$chromosome_name,
    ranges   = IRanges(start = as.integer(eq$start_position), end = as.integer(eq$end_position)),
    strand   = "*"
)
mcols(gene_gr)$gene         <- eq$gene
mcols(gene_gr)$gene_name    <- eq$gene_name
mcols(gene_gr)$pvalue       <- eq$pvalue
mcols(gene_gr)$pvalue_fdr   <- eq$pvalue_fdr
mcols(gene_gr)$z_mean       <- eq$z_mean

# Build GRanges sqtl
intron_gr <- GRanges(
    seqnames = sq$intron_chromosome_name,
    ranges   = IRanges(start = as.integer(sq$intron_start_position), end = as.integer(sq$intron_end_position)),
    strand   = "*"
)
mcols(intron_gr)$intron       <- sq$gene_name 
mcols(intron_gr)$gene         <- sq$gene      
mcols(intron_gr)$pvalue       <- sq$pvalue
mcols(intron_gr)$pvalue_fdr   <- sq$pvalue_fdr
mcols(intron_gr)$z_mean       <- sq$z_mean

hits_body <- findOverlaps(gene_gr, intron_gr, ignore.strand = TRUE)

df_body <- tibble(
    overlap_type = "intron_overlaps_gene_body",
    eqtl_gene        = mcols(gene_gr)$gene[queryHits(hits_body)],
    eqtl_gene_name   = mcols(gene_gr)$gene_name[queryHits(hits_body)],
    eqtl_chr         = as.character(seqnames(gene_gr))[queryHits(hits_body)],
    eqtl_start       = start(gene_gr)[queryHits(hits_body)],
    eqtl_end         = end(gene_gr)[queryHits(hits_body)],
    eqtl_pvalue      = mcols(gene_gr)$pvalue[queryHits(hits_body)],
    eqtl_pvalue_fdr  = mcols(gene_gr)$pvalue_fdr[queryHits(hits_body)],
    eqtl_z_mean      = mcols(gene_gr)$z_mean[queryHits(hits_body)],
    
    sqtl_intron_id   = mcols(intron_gr)$intron[subjectHits(hits_body)],
    sqtl_chr         = as.character(seqnames(intron_gr))[subjectHits(hits_body)],
    sqtl_start       = start(intron_gr)[subjectHits(hits_body)],
    sqtl_end         = end(intron_gr)[subjectHits(hits_body)],
    sqtl_pvalue      = mcols(intron_gr)$pvalue[subjectHits(hits_body)],
    sqtl_pvalue_fdr  = mcols(intron_gr)$pvalue_fdr[subjectHits(hits_body)],
    sqtl_z_mean      = mcols(intron_gr)$z_mean[subjectHits(hits_body)]
) %>% distinct()

df_cis <- tibble()

if (!is.null(cis_window_bp) && is.finite(cis_window_bp) && cis_window_bp > 0) {
    tss <- start(gene_gr)  # proxy
    
    gene_tss_win <- GRanges(
        seqnames = seqnames(gene_gr),
        ranges   = IRanges(start = pmax(1L, tss - cis_window_bp), end = tss + cis_window_bp),
        strand   = "*"
    )
    mcols(gene_tss_win) <- mcols(gene_gr)
    
    hits_cis <- findOverlaps(gene_tss_win, intron_gr, ignore.strand = TRUE)
    
    df_cis <- tibble(
        overlap_type = paste0("intron_within_", cis_window_bp, "bp_of_gene_TSSproxy"),
        eqtl_gene        = mcols(gene_tss_win)$gene[queryHits(hits_cis)],
        eqtl_gene_name   = mcols(gene_tss_win)$gene_name[queryHits(hits_cis)],
        eqtl_chr         = as.character(seqnames(gene_tss_win))[queryHits(hits_cis)],
        eqtl_start       = start(gene_tss_win)[queryHits(hits_cis)],
        eqtl_end         = end(gene_tss_win)[queryHits(hits_cis)],
        eqtl_pvalue      = mcols(gene_tss_win)$pvalue[queryHits(hits_cis)],
        eqtl_pvalue_fdr  = mcols(gene_tss_win)$pvalue_fdr[queryHits(hits_cis)],
        eqtl_z_mean      = mcols(gene_tss_win)$z_mean[queryHits(hits_cis)],
        
        sqtl_intron_id   = mcols(intron_gr)$intron[subjectHits(hits_cis)],
        sqtl_chr         = as.character(seqnames(intron_gr))[subjectHits(hits_cis)],
        sqtl_start       = start(intron_gr)[subjectHits(hits_cis)],
        sqtl_end         = end(intron_gr)[subjectHits(hits_cis)],
        sqtl_pvalue      = mcols(intron_gr)$pvalue[subjectHits(hits_cis)],
        sqtl_pvalue_fdr  = mcols(intron_gr)$pvalue_fdr[subjectHits(hits_cis)],
        sqtl_z_mean      = mcols(intron_gr)$z_mean[subjectHits(hits_cis)]
    ) %>% distinct()
}

gene_class <- tibble(
    eqtl_gene = unique(eq$gene),
    eqtl_gene_name = eq$gene_name[match(unique(eq$gene), eq$gene)]
) %>%
    left_join(df_body %>% distinct(eqtl_gene) %>% mutate(has_body_overlap = TRUE), by = "eqtl_gene") %>%
    left_join(df_cis  %>% distinct(eqtl_gene) %>% mutate(has_cis_overlap  = TRUE), by = "eqtl_gene") %>%
    mutate(
        has_body_overlap = ifelse(is.na(has_body_overlap), FALSE, has_body_overlap),
        has_cis_overlap  = ifelse(is.na(has_cis_overlap),  FALSE, has_cis_overlap),
        overlap_class = case_when(
            has_body_overlap ~ "overlap_gene_body",
            !has_body_overlap & has_cis_overlap ~ "overlap_cis_window_only",
            TRUE ~ "eqtl_only"
        )
    ) %>%
    arrange(overlap_class, eqtl_gene_name)

overlap_summary <- tibble(
    n_eqtl_genes_used = length(unique(eq$gene)),
    n_sqtl_introns_used = length(unique(sq$gene_name)),
    n_pairs_gene_body = nrow(df_body),
    n_eqtl_genes_with_body_overlap = nrow(df_body %>% distinct(eqtl_gene)),
    n_sqtl_introns_with_body_overlap = nrow(df_body %>% distinct(sqtl_intron_id)),
    n_pairs_cis_window = ifelse(nrow(df_cis) == 0, NA_integer_, nrow(df_cis)),
    n_eqtl_genes_with_cis_overlap = ifelse(nrow(df_cis) == 0, NA_integer_, nrow(df_cis %>% distinct(eqtl_gene))),
    n_sqtl_introns_with_cis_overlap = ifelse(nrow(df_cis) == 0, NA_integer_, nrow(df_cis %>% distinct(sqtl_intron_id)))
)

print(overlap_summary)
print(gene_class %>% count(overlap_class))

write.csv(df_body, file.path(out_dir, "overlap_pairs_gene_body.csv"), row.names = FALSE)
if (nrow(df_cis) > 0) write.csv(df_cis, file.path(out_dir, "overlap_pairs_cis_window.csv"), row.names = FALSE)
write.csv(overlap_summary, file.path(out_dir, "overlap_summary.csv"), row.names = FALSE)
write.csv(gene_class, file.path(out_dir, "gene_overlap_classification.csv"), row.names = FALSE)

# merge evidence

  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)

# Paths
dir_results   <- [[PLACEHOLDER]] 
dir_coloc     <- file.path(dir_results, "coloc")
dir_smr       <- file.path(dir_results, "smr")
dir_twas      <- file.path(dir_results, "twas")
dir_repl      <- file.path(dir_results, "replication")
dir_out       <- file.path(dir_results, "evidence_tables")  # will create

if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

# Thresholds
twas_fdr_thresh   <- 0.05
coloc_pp4_strict  <- 0.80
coloc_pp4_relaxed <- 0.50

# SMR thresholds
smr_p_thresh      <- 5e-8
heidi_p_thresh    <- 0.01

# Helpers
read_twas_tab <- function(path) {
  readr::read_delim(path, delim = "\t", show_col_types = FALSE) %>%
    janitor::clean_names()
}

strip_ens_version <- function(x) {
  sub("\\..*$", "", x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x, na.rm = TRUE)
}

safe_which_min <- function(x) {
  if (all(!is.finite(x))) NA_integer_ else which.min(x)
}

f_twas_eqtl   <- file.path(dir_twas, "eqtl_smultixcan_wen_davatzikos_combined_brain_tissues.csv")
f_coloc_eqtl  <- file.path(dir_coloc, "coloc_eqtl_A1A2_summary.csv")
f_smr_eqtl    <- file.path(dir_smr,"summary_integrated", "eqtl_smr_merged_long.csv")
f_repl_eqtl   <- file.path(dir_repl, "replication_feature_level_eqtl_summary.csv")

twas_eqtl <- read_twas_tab(f_twas_eqtl) %>%
  mutate(
    gene_id        = strip_ens_version(gene),
    gene_name      = as.character(gene_name),
    twas_p         = as.numeric(pvalue),
    twas_fdr       = p.adjust(twas_p, method = "BH"),
    twas_sig_fdr   = twas_fdr < twas_fdr_thresh
  )

# Discovery set (FDR)
twas_eqtl_sig <- twas_eqtl %>%
  filter(twas_sig_fdr) %>%
  dplyr::select(
    gene_id, gene_name,
    twas_p, twas_fdr,
    z_mean, z_min, z_max, z_sd,
    p_i_best, t_i_best, p_i_worst, t_i_worst,
    n, n_indep, tmi, status
  ) %>%
  arrange(twas_p)

# Replication 
repl_eqtl <- readr::read_csv(f_repl_eqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(
    gene_id = as.character(feature),
    repl_z_primary = as.numeric(z_primary),
    repl_p_primary = as.numeric(p_primary),
    repl_q_primary = as.numeric(q_primary),
    n_rep_available = as.integer(n_rep_available),
    n_rep_nominal   = as.integer(n_rep_nominal),
    n_rep_consistent= as.integer(n_rep_consistent),
    any_rep_nominal = as.logical(any_rep_nominal),
    any_rep_consistent = as.logical(any_rep_consistent)
  )

# COLOC summary
coloc_eqtl <- readr::read_csv(f_coloc_eqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(
    gene_name = as.character(gene),
    coloc_pp4_max = as.numeric(pp4_max),
    coloc_tissue_max = as.character(tissue_max),
    coloc_n_tissues_pp4_ge_0_8 = as.integer(n_tissues_pp4_ge_0_8),
    coloc_n_tissues_pp4_ge_0_5 = as.integer(n_tissues_pp4_ge_0_5),
    coloc_support_strict  = coloc_pp4_max >= coloc_pp4_strict,
    coloc_support_relaxed = coloc_pp4_max >= coloc_pp4_relaxed
  )

# SMR+HEIDI 
smr_eqtl_long <- readr::read_csv(f_smr_eqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    gene_id   = as.character(gene_noversion),
    gene_name = as.character(gene_name),
    tissue    = as.character(tissue),
    p_smr     = as.numeric(p_smr),
    p_heidi   = as.numeric(p_heidi)
  )

smr_eqtl_sum <- smr_eqtl_long %>%
  group_by(gene_id, gene_name) %>%
  summarise(
    smr_min_p = safe_min(p_smr),
    smr_n_tests = sum(is.finite(p_smr), na.rm = TRUE),
    smr_n_sig = sum(is.finite(p_smr) & p_smr < smr_p_thresh, na.rm = TRUE),
    
    heidi_n_tested = sum(is.finite(p_heidi), na.rm = TRUE),
    smr_heidi_n_pass = sum(
      is.finite(p_smr) & p_smr < smr_p_thresh & is.finite(p_heidi) & (p_heidi > heidi_p_thresh),
      na.rm = TRUE
    ),
    smr_heidi_any_pass = any(
      is.finite(p_smr) & p_smr < smr_p_thresh &
        ( (is.finite(p_heidi) & p_heidi > heidi_p_thresh) ),
      na.rm = TRUE
    ),
    
    smr_best_tissue = {
      idx <- safe_which_min(p_smr)
      if (is.na(idx)) NA_character_ else tissue[idx]
    },
    smr_best_p = smr_min_p,
    smr_best_heidi_p = {
      idx <- safe_which_min(p_smr)
      if (is.na(idx)) NA_real_ else p_heidi[idx]
    },
    .groups = "drop"
  )

# Merge into one evidence table 
eqtl_evidence <- twas_eqtl_sig %>%
  left_join(repl_eqtl, by = "gene_id") %>%
  left_join(coloc_eqtl, by = "gene_name") %>%
  left_join(smr_eqtl_sum, by = c("gene_id","gene_name")) %>%
  mutate(
    high_confidence = (
      (isTRUE(any_rep_consistent)) &
        (isTRUE(coloc_support_strict) | isTRUE(smr_heidi_any_pass))
    )
  ) %>%
  arrange(twas_p)

out_eqtl <- file.path(dir_out, "high_confidence_eqtl_primary.csv")
readr::write_csv(eqtl_evidence, out_eqtl)
message("Wrote: ", out_eqtl)

# sQTL
f_twas_sqtl   <- file.path(dir_twas, "sqtl_smultixcan_wen_davatzikos_combined_brain_tissues.csv")
f_coloc_sqtl  <- file.path(dir_coloc, "coloc_sqtl_A1A2_summary.csv")
f_smr_sqtl    <- file.path(dir_smr, "summary_integrated", "sqtl_smr_merged_long.csv")
f_repl_sqtl   <- file.path(dir_repl, "replication_feature_level_sqtl_summary.csv")

twas_sqtl <- read_twas_tab(f_twas_sqtl) %>%
  mutate(
    intron_id      = as.character(gene),  
    twas_p         = as.numeric(pvalue),
    twas_fdr       = p.adjust(twas_p, method = "BH"),
    twas_sig_fdr   = twas_fdr < twas_fdr_thresh
  )

twas_sqtl_sig <- twas_sqtl %>%
  filter(twas_sig_fdr) %>%
  transmute(
    intron_id,
    twas_p, twas_fdr,
    z_mean, z_min, z_max, z_sd,
    p_i_best, t_i_best, p_i_worst, t_i_worst,
    n, n_indep, tmi, status
  ) %>%
  arrange(twas_p)

repl_sqtl <- readr::read_csv(f_repl_sqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(
    intron_id = as.character(feature),
    repl_z_primary = as.numeric(z_primary),
    repl_p_primary = as.numeric(p_primary),
    repl_q_primary = as.numeric(q_primary),
    n_rep_available = as.integer(n_rep_available),
    n_rep_nominal   = as.integer(n_rep_nominal),
    n_rep_consistent= as.integer(n_rep_consistent),
    any_rep_nominal = as.logical(any_rep_nominal),
    any_rep_consistent = as.logical(any_rep_consistent)
  )

coloc_sqtl <- readr::read_csv(f_coloc_sqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(
    intron_id = as.character(intron_id),
    coloc_pp4_max = as.numeric(pp4_max),
    coloc_tissue_max = as.character(tissue_max),
    coloc_n_tissues_pp4_ge_0_8 = as.integer(n_tissues_pp4_ge_0_8),
    coloc_n_tissues_pp4_ge_0_5 = as.integer(n_tissues_pp4_ge_0_5),
    coloc_support_strict  = coloc_pp4_max >= coloc_pp4_strict,
    coloc_support_relaxed = coloc_pp4_max >= coloc_pp4_relaxed
  )

smr_sqtl_long <- readr::read_csv(f_smr_sqtl, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    intron_id = as.character(feature_id),
    gene_symbol_from_smr = as.character(gene), 
    ensg_from_probe = as.character(str_replace(probe_id, ".*:(ENSG[0-9]+)\\..*$", "\\1")),
    tissue = as.character(tissue),
    p_smr   = as.numeric(p_smr),
    p_heidi = as.numeric(p_heidi)
  )

smr_sqtl_sum <- smr_sqtl_long %>%
  group_by(intron_id) %>%
  summarise(
    smr_min_p = safe_min(p_smr),
    smr_n_tests = sum(is.finite(p_smr), na.rm = TRUE),
    smr_n_sig = sum(is.finite(p_smr) & p_smr < smr_p_thresh, na.rm = TRUE),
    heidi_n_tested = sum(is.finite(p_heidi), na.rm = TRUE),
    smr_heidi_n_pass = sum(
      is.finite(p_smr) & p_smr < smr_p_thresh & is.finite(p_heidi) & (p_heidi > heidi_p_thresh),
      na.rm = TRUE
    ),
    smr_heidi_any_pass = any(
      is.finite(p_smr) & p_smr < smr_p_thresh & (is.finite(p_heidi) & p_heidi > heidi_p_thresh),
      na.rm = TRUE
    ),
    smr_best_tissue = {
      idx <- safe_which_min(p_smr)
      if (is.na(idx)) NA_character_ else tissue[idx]
    },
    smr_best_p = smr_min_p,
    smr_best_heidi_p = {
      idx <- safe_which_min(p_smr)
      if (is.na(idx)) NA_real_ else p_heidi[idx]
    },
    gene_symbol = dplyr::first(na.omit(gene_symbol_from_smr)),
    ensg_id     = dplyr::first(na.omit(ensg_from_probe)),
    .groups = "drop"
  )

sqtl_evidence <- twas_sqtl_sig %>%
  left_join(repl_sqtl, by = "intron_id") %>%
  left_join(coloc_sqtl, by = "intron_id") %>%
  left_join(smr_sqtl_sum, by = "intron_id") %>%
  mutate(
    high_confidence = (
      (isTRUE(any_rep_consistent)) &
        (isTRUE(coloc_support_strict) | isTRUE(smr_heidi_any_pass))
    )
  ) %>%
  arrange(twas_p)

out_sqtl <- file.path(dir_out, "high_confidence_sqtl_primary.csv")
readr::write_csv(sqtl_evidence, out_sqtl)
message("Wrote: ", out_sqtl)

readr::write_csv(eqtl_evidence %>% filter(high_confidence), file.path(dir_out, "high_confidence_eqtl_primary_ONLY.csv"))
readr::write_csv(sqtl_evidence %>% filter(high_confidence), file.path(dir_out, "high_confidence_sqtl_primary_ONLY.csv"))

# conditional

chr <- 17
from_bp <- 44630158
to_bp   <- 46949182
cond_snp <- "rs242561"
gwas_in  <- [[PLACEHOLDER]] 
gcta_exe <- [[PLACEHOLDER]] 
bfile_eur <- [[PLACEHOLDER]] 
dir_out <- [[PLACEHOLDER]] 
out_prefix <- file.path(dir_out, "primary_chr17_rs242561_cond")
gwas_region_cojo <- file.path(dir_out, "primary_chr17_44630158_46949182_for_cojo.txt")

# Read  GWAS 
gwas <- fread(gwas_in)
# map cols
setnames(gwas,
         old = c("variant_id","effect_allele","non_effect_allele","frequency",  "effect_size","standard_error","pvalue","n_cases"),
         new = c("SNP",       "A1",           "A2",               "freq",       "b",          "se",             "p","N"),
         skip_absent = TRUE)

# Remove missing
gwas <- gwas[is.finite(b) & is.finite(se) & is.finite(p) & !is.na(freq)]

# filter to region
gwas_reg <- gwas[chromosome == paste0("chr",as.character(chr)) & position >= as.character(from_bp) & position <= as.character(to_bp)]
fwrite(gwas_reg[, .(SNP,A1,A2,freq,b,se,p,N)], gwas_region_cojo, sep = "\t")

cond_file <- [[PLACEHOLDER]] 
writeLines("rs242561", cond_file)

# run
args <- c(
  "--bfile", shQuote(bfile_eur),
  "--cojo-file", shQuote(gwas_region_cojo),
  "--cojo-cond", shQuote(cond_file),
  "--out", shQuote(out_prefix)
)
res <- system2(gcta_exe, args = args, stdout = "", stderr = "")
cat(paste(res, collapse = "\n"))

# coloc on conditioned gwas
cojo_cma <- [[PLACEHOLDER]] 
gc <- fread(cojo_cma)  
gwas <- "wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz"
dir_gwas <- [[PLACEHOLDER]] 
gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

# load and setup gwas
dt_gwas <- data.table::fread(
  file.path(dir_gwas, gwas)
)
data.table::setnames(dt_gwas, old = c("chromosome","position","frequency"), new = c("chr_raw","pos","maf_gwas"))
dt_gwas[, chr := gsub("^chr", "", as.character(chr_raw))]
dt_gwas[, pos := as.numeric(pos)]
dt_gwas <- dt_gwas[!is.na(panel_variant_id) & !is.na(pos) & !is.na(pvalue) & !is.na(maf_gwas)]
data.table::setkey(dt_gwas, panel_variant_id)

rs_col <- "variant_id"

map <- unique(dt_gwas[, .(panel_variant_id, rsid = get(rs_col))])
setkey(map, rsid)

dt_gwas_cond <- gc[, .(rsid = SNP, pvalue_cond = pC)][map, on="rsid"]
dt_gwas_cond <- dt_gwas_cond[!is.na(panel_variant_id), .(panel_variant_id, pvalue_cond)]
setkey(dt_gwas_cond, panel_variant_id)

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(coloc)
})

cis_kb <- 1000
gwas_n <- 29400
sqtl_n <- 150

gwas <- "wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz"
twas <- df_twas_results[["wen_davatzikos_combined"]][["sqtl"]][["smultixcan"]][["brain_tissues"]]

dir_gtex_sqtl <- [[PLACEHOLDER]] 
dir_results   <- [[PLACEHOLDER]] 

gtex_brain_tissues <- c(
  "Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
  "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
  "Hippocampus", "Hypothalamus", "Nucleus_accumbens_basal_ganglia",
  "Putamen_basal_ganglia", "Spinal_cord_cervical_c-1", "Substantia_nigra"
)

use_conditioned_gwas <- TRUE

if (use_conditioned_gwas) {
  stopifnot(exists("dt_gwas_cond"))
  stopifnot(is.data.table(dt_gwas_cond))
  stopifnot(all(c("panel_variant_id", "pvalue_cond") %in% names(dt_gwas_cond)))
  if (!identical(key(dt_gwas_cond), "panel_variant_id")) setkey(dt_gwas_cond, panel_variant_id)
  
  dt_gwas <- dt_gwas[dt_gwas_cond, on = "panel_variant_id"]
}

stopifnot(all(c("chr", "pos", "pvalue", "maf_gwas") %in% names(dt_gwas)))

twas_sig <- twas[twas$fdr_significant == TRUE, , drop = FALSE]
setDT(twas_sig)

results <- list()
n_total <- nrow(twas_sig)
message("TWAS significant introns: ", n_total)

# Loop per chromosome
for (chr_val in unique(twas_sig$intron_chromosome_name)) {
  
  chr_clean <- gsub("^chr", "", as.character(chr_val))
  dt_gwas_chr <- dt_gwas[chr == chr_clean]
  setkey(dt_gwas_chr, panel_variant_id)
  
  introns_on_chr <- twas_sig[intron_chromosome_name == chr_val]
  if (nrow(introns_on_chr) == 0) next
    
  # Loop per tissue
  for (gtex_tissue in gtex_brain_tissues) {
    
    gtex_file <- file.path(
      dir_gtex_sqtl,
      paste0(
        "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_sQTL_all_associations_Brain_",
        gtex_tissue,
        ".v8.EUR.sqtl_allpairs.chr",
        chr_clean,
        ".parquet"
      )
    )
    
    if (!file.exists(gtex_file)) {
      message("  (skip) missing parquet: ", basename(gtex_file))
      next
    }
    
    gtex_sqtl <- arrow::read_parquet(
      gtex_file,
      col_select = c("phenotype_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
    )
    setDT(gtex_sqtl)
    
    gtex_sqtl <- gtex_sqtl[
      !is.na(phenotype_id) & !is.na(variant_id) &
        !is.na(maf) & !is.na(pval_nominal) &
        !is.na(slope) & !is.na(slope_se) & slope_se != 0
    ]
    if (nrow(gtex_sqtl) == 0L) {
      message("  (skip) empty gtex_sqtl after filtering")
      rm(gtex_sqtl); gc()
      next
    }
    
    for (i in seq_len(nrow(introns_on_chr))) {
      
      intron_id  <- as.character(introns_on_chr[["gene"]][i]) 
      intr_start <- as.numeric(introns_on_chr$intron_start[i])
      intr_end   <- as.numeric(introns_on_chr$intron_end[i])
      if (!is.finite(intr_start) || !is.finite(intr_end)) next
      
      # phenotype_id prefix: chr17:START:END:
      prefix <- paste0("chr", chr_clean, ":", intr_start, ":", intr_end, ":")
      
      gene_sqtl <- gtex_sqtl[startsWith(phenotype_id, prefix)]
      if (nrow(gene_sqtl) == 0L) next
      
      cis_start <- intr_start - cis_kb * 1000
      cis_end   <- intr_end   + cis_kb * 1000
      
      gene_sqtl[, pos := as.numeric(tstrsplit(variant_id, "_", fixed = TRUE, keep = 2L)[[1]])]
      gene_sqtl <- gene_sqtl[!is.na(pos) & pos >= cis_start & pos <= cis_end]
      if (nrow(gene_sqtl) == 0L) next
      
      gene_sqtl[, panel_variant_id := variant_id]
      gene_sqtl[, maf_sqtl := maf]
      setkey(gene_sqtl, panel_variant_id)
      
      gwas_sub <- dt_gwas_chr[pos >= cis_start & pos <= cis_end]
      if (nrow(gwas_sub) == 0L) next
      
      merged <- merge(gwas_sub, gene_sqtl, by = "panel_variant_id", all = FALSE)
      if (nrow(merged) < 10L) next
      
      if (use_conditioned_gwas) {
        merged <- merged[!is.na(pvalue_cond)]
        if (nrow(merged) < 10L) next
        merged[, pvalue_use := pvalue_cond]
      } else {
        merged[, pvalue_use := pvalue]
      }
      
      # QC filters
      merged <- merged[
        !is.na(pvalue_use) &
          !is.na(maf_gwas) & maf_gwas > 0 & maf_gwas < 1 &
          !is.na(slope) & !is.na(slope_se) & slope_se != 0 &
          !is.na(maf_sqtl) & maf_sqtl > 0 & maf_sqtl < 1
      ]
      if (nrow(merged) < 10L) next
      
      # one sQTL row per SNP
      merged <- merged[, .SD[which.min(pval_nominal)], by = panel_variant_id]
      if (nrow(merged) < 10L) next
      
      # ensure pos2 exists
      merged[, pos2 := as.numeric(tstrsplit(panel_variant_id, "_", fixed = TRUE, keep = 2L)[[1]])]
      merged <- merged[!is.na(pos2)]
      if (nrow(merged) < 10L) next
      
      lead_pos_gwas <- merged[which.min(pvalue_use)]$pos2[1]
      lead_pos_sqtl <- merged[which.min(pval_nominal)]$pos2[1]
      lead_pos <- if (is.finite(lead_pos_gwas)) lead_pos_gwas else lead_pos_sqtl
      if (!is.finite(lead_pos)) next
      
      win_kb <- 500
      win_start <- lead_pos - win_kb * 1000
      win_end   <- lead_pos + win_kb * 1000
      
      merged2 <- merged[pos2 >= win_start & pos2 <= win_end]
      if (nrow(merged2) < 50L) {
        merged2 <- merged
      }
      if (nrow(merged2) < 10L) next
      
      dataset1 <- list(
        snp     = merged2$panel_variant_id,
        pvalues = merged2$pvalue_use,
        N       = gwas_n,
        MAF     = merged2$maf_gwas,
        type    = "quant"
      )
      
      dataset2 <- list(
        snp     = merged2$panel_variant_id,
        pvalues = merged2$pval_nominal,
        N       = sqtl_n,
        MAF     = merged2$maf_sqtl,
        type    = "quant"
      )
      
      coloc_res <- coloc::coloc.abf(dataset1, dataset2)
      
      if (!is.null(coloc_res$summary) && coloc_res$summary["PP.H4.abf"] > 0.8) {
        message("  COLOC>0.8: ", gtex_tissue, " | chr", chr_clean, " | ", intron_id,
                " | PP4=", signif(coloc_res$summary["PP.H4.abf"], 3))
      }
      
      key_name <- paste0(intron_id, "_", gtex_tissue, "_", chr_clean, "_", intr_start, "_", intr_end)
      results[[key_name]] <- list(
        intron_id = intron_id,
        chr = chr_clean,
        start = intr_start,
        end = intr_end,
        tissue = gtex_tissue,
        summary = coloc_res$summary
      )
    }
    
    rm(gtex_sqtl); gc()
  }
}

coloc_summary_dt <- rbindlist(
  lapply(results, function(res_item) {
    data.table(
      intron_id = res_item$intron_id,
      tissue   = res_item$tissue,
      chr      = res_item$chr,
      start    = res_item$start,
      end      = res_item$end,
      t(res_item$summary)
    )
  }),
  fill = TRUE
)

if (nrow(coloc_summary_dt) == 0L) {
} else {
  coloc_summary_dt[, is_colocalized := PP.H4.abf > 0.8]
  
  out_dir <- file.path(dir_results, "gcta_coloc")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- if (use_conditioned_gwas) "COND_rs242561" else "RAW"
  output_filename <- paste0("coloc_sqtl_", tag, "_", gwas, ".csv")
  out_path <- file.path(out_dir, output_filename)
  
  fwrite(coloc_summary_dt, file = out_path)
}

