library(TwoSampleMR)
library(data.table)
library(ggplot2)
library(dplyr)

RUN_MODE <- "e2o"  # "e2o" or "o2e"

if (RUN_MODE == "e2o") {
  SOURCE_DIRS <- c("/mnt/f/result/post_ld/e2o/e2o_1-1200", 
                   "/mnt/f/result/post_ld/e2o/e2o_1201-2400", 
                   "/mnt/f/result/post_ld/e2o/e2o_2401-3935")
  OUTPUT_PATH <- "/mnt/f/result/mr_results_e2o"
  FILTER_COL <- "pval.exposure"
} else {
  SOURCE_DIRS <- c("/mnt/f/result/post_ld/o2e/o2e_1-1200", 
                   "/mnt/f/result/post_ld/o2e/o2e_1201-2400", 
                   "/mnt/f/result/post_ld/o2e/o2e_2401-3935")
  OUTPUT_PATH <- "/mnt/f/result/mr_results_o2e"
  FILTER_COL <- "pval.outcome"
}

SNP_OUTPUT_PATH <- file.path(OUTPUT_PATH, "significant_snps")
dir.create(SNP_OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

PVAL_THRESHOLD <- 5e-8
NOMINAL_THRESHOLD <- 0.05
MIN_SNPS <- 3
BATCH_SIZE <- 100

cat("Configuration loaded.\n")
cat(sprintf("Mode: %s\n", RUN_MODE))
cat(sprintf("Output path: %s\n", OUTPUT_PATH))


get_processed_files <- function() {
  status_file <- file.path(OUTPUT_PATH, paste0("processing_status_", RUN_MODE, ".csv"))
  if (!file.exists(status_file)) return(character(0))
  status <- fread(status_file)
  return(as.character(status[status == "success", source_file]))
}

update_status <- function(stat_file, source_file, status, reason = NA) {
  fwrite(data.table(
    source_file = source_file,
    status = status,
    reason = reason,
    timestamp = Sys.time()
  ), stat_file, append = TRUE)
}

cat("Resume functions loaded.\n")


create_mr_datasets <- function(dt, id_a, id_b, mode) {
  if (mode == "e2o") {
    exp_dat <- data.frame(
      SNP = dt$SNP,
      beta.exposure = dt$beta.exposure,
      se.exposure = dt$se.exposure,
      pval.exposure = dt$pval.exposure,
      effect_allele.exposure = dt$A1.exposure,
      other_allele.exposure = dt$A2.exposure,
      eaf.exposure = if("eaf.exposure" %in% names(dt)) dt$eaf.exposure else NA,
      samplesize.exposure = if("samplesize.exposure" %in% names(dt)) dt$samplesize.exposure else NA,
      exposure = id_a,
      id.exposure = id_a
    )
    
    out_dat <- data.frame(
      SNP = dt$SNP,
      beta.outcome = dt$beta.outcome,
      se.outcome = dt$se.outcome,
      pval.outcome = dt$pval.outcome,
      effect_allele.outcome = dt$A1.outcome,
      other_allele.outcome = dt$A2.outcome,
      eaf.outcome = if("eaf.outcome" %in% names(dt)) dt$eaf.outcome else NA,
      samplesize.outcome = if("samplesize.outcome" %in% names(dt)) dt$samplesize.outcome else NA,
      outcome = id_b,
      id.outcome = id_b
    )
  } else {
    exp_dat <- data.frame(
      SNP = dt$SNP,
      beta.exposure = dt$beta.outcome,
      se.exposure = dt$se.outcome,
      pval.exposure = dt$pval.outcome,
      effect_allele.exposure = dt$A1.outcome,
      other_allele.exposure = dt$A2.outcome,
      eaf.exposure = if("eaf.outcome" %in% names(dt)) dt$eaf.outcome else NA,
      samplesize.exposure = if("samplesize.outcome" %in% names(dt)) dt$samplesize.outcome else NA,
      exposure = id_a,
      id.exposure = id_a
    )
    
    out_dat <- data.frame(
      SNP = dt$SNP,
      beta.outcome = dt$beta.exposure,
      se.outcome = dt$se.exposure,
      pval.outcome = dt$pval.exposure,
      effect_allele.outcome = dt$A1.exposure,
      other_allele.outcome = dt$A2.exposure,
      eaf.outcome = if("eaf.exposure" %in% names(dt)) dt$eaf.exposure else NA,
      samplesize.outcome = if("samplesize.exposure" %in% names(dt)) dt$samplesize.exposure else NA,
      outcome = id_b,
      id.outcome = id_b
    )
  }
  
  exp_dat$effect_allele.exposure <- as.factor(as.character(exp_dat$effect_allele.exposure))
  exp_dat$other_allele.exposure <- as.factor(as.character(exp_dat$other_allele.exposure))
  out_dat$effect_allele.outcome <- as.factor(as.character(out_dat$effect_allele.outcome))
  out_dat$other_allele.outcome <- as.factor(as.character(out_dat$other_allele.outcome))
  
  return(list(exposure = exp_dat, outcome = out_dat))
}

cat("Data mapping function loaded.\n")


process_mr_file <- function(file_path, mode, filter_col) {
  tryCatch({
    dt <- fread(file_path)
    f_name <- basename(file_path)
    
    parts <- strsplit(f_name, "_to_")[[1]]
    if (length(parts) != 2) {
      return(list(status = "error", reason = "Invalid filename format"))
    }
    
    id_a <- gsub("\\.txt$", "", parts[1])
    id_b <- gsub("\\.txt\\.csv$|\\.csv$", "", parts[2])
    
    required_cols <- c("SNP", "beta.exposure", "se.exposure", "pval.exposure", 
                       "beta.outcome", "se.outcome", "pval.outcome",
                       "A1.exposure", "A2.exposure", "A1.outcome", "A2.outcome")
    missing_cols <- required_cols[!required_cols %in% names(dt)]
    if (length(missing_cols) > 0) {
      return(list(status = "error", reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))))
    }
    
    dt_filtered <- dt[get(filter_col) < PVAL_THRESHOLD, ]
    
    if (nrow(dt_filtered) < MIN_SNPS) {
      return(list(status = "insufficient_snps", 
                  reason = sprintf("Only %d SNPs with p < %.0e", nrow(dt_filtered), PVAL_THRESHOLD),
                  nsnp = nrow(dt_filtered)))
    }
    
    datasets <- create_mr_datasets(dt_filtered, id_a, id_b, mode)
    harmonized <- harmonise_data(datasets$exposure, datasets$outcome, action = 2)
    
    if (nrow(harmonized) < MIN_SNPS) {
      return(list(status = "insufficient_after_harmonization",
                  reason = paste("Only", nrow(harmonized), "SNPs after harmonization"),
                  nsnp = nrow(harmonized)))
    }
    
    mr_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    mr_res <- mr(harmonized, method_list = mr_methods)
    
    ivw_result <- mr_res[mr_res$method == "Inverse variance weighted", ]
    ivw_pval <- if(nrow(ivw_result) > 0) ivw_result$pval else NA
    
    results_row <- data.table(
      exposure = id_a,
      outcome = id_b,
      nsnp_original = nrow(dt),
      nsnp_filtered = nrow(dt_filtered),
      nsnp_harmonized = nrow(harmonized),
      ivw_pval = ivw_pval
    )
    
    for (i in 1:nrow(mr_res)) {
      method_clean <- gsub("mr_", "", mr_res$method[i])
      method_clean <- gsub(" ", "_", method_clean)
      results_row[[paste0("b_", method_clean)]] <- mr_res$b[i]
      results_row[[paste0("se_", method_clean)]] <- mr_res$se[i]
      results_row[[paste0("pval_", method_clean)]] <- mr_res$pval[i]
    }
    
    return(list(status = "success", results = results_row, snps = dt_filtered,
                harmonized = harmonized, source_file = f_name, nsnp = nrow(harmonized)))
    
  }, error = function(e) {
    return(list(status = "error", reason = e$message, source_file = basename(file_path)))
  })
}

cat("Main processing function loaded.\n")


cat("\n=== Starting Batch MR Analysis ===\n")

processed_files <- get_processed_files()
res_file <- file.path(OUTPUT_PATH, paste0("mr_results_intermediate_", RUN_MODE, ".csv"))
stat_file <- file.path(OUTPUT_PATH, paste0("processing_status_", RUN_MODE, ".csv"))

total_processed <- 0
total_success <- 0
total_errors <- 0
start_time <- Sys.time()

for (dir_path in SOURCE_DIRS) {
  if (!dir.exists(dir_path)) {
    cat(sprintf("Warning: Directory does not exist: %s\n", dir_path))
    next
  }
  
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("\n=== Processing Directory: %s ===\n", dir_path))
  cat(sprintf("Found %d files\n", length(files)))
  
  files_to_process <- files[!basename(files) %in% processed_files]
  cat(sprintf("Files to process: %d\n", length(files_to_process)))
  
  if (length(files_to_process) == 0) next
  
  for (i in seq_along(files_to_process)) {
    f <- files_to_process[i]
    
    if (i %% BATCH_SIZE == 0 || i == 1 || i == length(files_to_process)) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      cat(sprintf("  Progress: %d/%d files (%.1f%%) | Time: %.1f min | Success: %d | Errors: %d\n",
                  i, length(files_to_process), 100*i/length(files_to_process),
                  elapsed, total_success, total_errors))
    }
    
    res <- process_mr_file(f, RUN_MODE, FILTER_COL)
    total_processed <- total_processed + 1
    
    if (res$status == "success") {
      total_success <- total_success + 1
      fwrite(res$results, res_file, append = file.exists(res_file))
      
      if (!is.na(res$results$ivw_pval) && res$results$ivw_pval < NOMINAL_THRESHOLD) {
        fwrite(res$snps, file.path(SNP_OUTPUT_PATH, paste0(res$source_file, "_snps.csv")))
        fwrite(res$harmonized, file.path(SNP_OUTPUT_PATH, paste0(res$source_file, "_harmonized.csv")))
      }
    } else {
      total_errors <- total_errors + 1
      if (total_errors <= 5) {
        cat(sprintf("  Error on %s: %s\n", basename(f), res$reason))
      } else if (total_errors == 6) {
        cat("  ... suppressing further error messages ...\n")
      }
    }
    
    update_status(stat_file, basename(f), res$status, if(!is.null(res$reason)) res$reason else NA)
  }
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Total files processed: %d\n", total_processed))
cat(sprintf("Successful: %d\n", total_success))
cat(sprintf("Errors/Insufficient: %d\n", total_errors))
cat(sprintf("Success rate: %.1f%%\n", 100 * total_success / total_processed))
cat(sprintf("Total time: %.2f minutes\n", total_time))


cat("\n=== Cleaning Results and Computing FDR ===\n")

res_file <- file.path(OUTPUT_PATH, paste0("mr_results_intermediate_", RUN_MODE, ".csv"))

if (file.exists(res_file)) {
  raw_results <- fread(res_file, fill = TRUE, blank.lines.skip = TRUE)
  cat(sprintf("Raw results read: %d rows, %d columns\n", nrow(raw_results), ncol(raw_results)))
  
  if ("ivw_pval" %in% names(raw_results)) {
    complete_rows <- !is.na(raw_results$ivw_pval)
    valid_rows <- complete_rows & !is.na(raw_results$exposure) & !is.na(raw_results$outcome)
    clean_results <- raw_results[valid_rows, ]
    cat(sprintf("Cleaned results: %d rows\n", nrow(clean_results)))
    
    fwrite(clean_results, file.path(OUTPUT_PATH, "mr_results_cleaned.csv"))
    
    clean_results$pval_fdr <- p.adjust(clean_results$ivw_pval, method = "fdr")
    clean_results$is_significant_fdr <- clean_results$pval_fdr < NOMINAL_THRESHOLD
    
    fdr_sig <- clean_results[is_significant_fdr == TRUE]
    cat(sprintf("FDR-significant pairs: %d\n", nrow(fdr_sig)))
    
    if (nrow(fdr_sig) > 0) {
      fwrite(fdr_sig, file.path(OUTPUT_PATH, "fdr_significant_pairs.csv"))
      cat("✓ FDR-significant pairs saved\n")
    }
    
    fwrite(clean_results, file.path(OUTPUT_PATH, "mr_results_complete_with_fdr.csv"))
    cat("✓ Complete results with FDR saved\n")
  } else {
    cat("ERROR: ivw_pval column not found\n")
    clean_results <- NULL
  }
} else {
  cat("No results file found. Please run Part 5 first.\n")
  clean_results <- NULL
}


library(MRPRESSO)

cat("\n=== Running MR-PRESSO and Steiger Tests ===\n")

snp_files <- list.files(SNP_OUTPUT_PATH, pattern = "_snps.csv$", full.names = TRUE)

if (length(snp_files) > 0) {
  
  cat(sprintf("Found %d SNP files to process\n", length(snp_files)))
  
  presso_results <- list()
  steiger_results <- list()
  
  for (i in seq_along(snp_files)) {
    if (i %% 100 == 0) cat(sprintf("  Processing %d/%d...\n", i, length(snp_files)))
    
    snp_data <- fread(snp_files[i])
    
    fname <- basename(snp_files[i])
    fname <- gsub("_snps.csv$", "", fname)
    parts <- strsplit(fname, "_to_")[[1]]
    
    if (length(parts) == 2) {
      exposure_id <- as.character(parts[1])  # Force character
      outcome_id <- as.character(gsub("\\.txt$", "", parts[2]))  
      
      if (nrow(snp_data) >= 4) {
        presso_res <- tryCatch({
          mr_presso(BetaOutcome = "beta.outcome", 
                    BetaExposure = "beta.exposure", 
                    SdOutcome = "se.outcome", 
                    SdExposure = "se.exposure", 
                    OUTLIERtest = TRUE, 
                    DISTORTIONtest = TRUE, 
                    data = as.data.frame(snp_data),
                    NbDistribution = 1000, 
                    SignifThreshold = 0.05)
        }, error = function(e) NULL)
        
        if (!is.null(presso_res)) {
          global_p <- presso_res$`MR-PRESSO results`$`Global Test`$Pvalue
          distortion_p <- presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
          outliers <- !is.null(presso_res$`MR-PRESSO results`$`Outlier Test`$OutliersIndices)
        } else {
          global_p <- distortion_p <- NA_real_
          outliers <- FALSE
        }
      } else {
        global_p <- distortion_p <- NA_real_
        outliers <- FALSE
      }
      
      presso_results[[i]] <- data.table(
        exposure_id = exposure_id,
        outcome_id = outcome_id,
        presso_global_p = global_p,
        presso_distortion_p = distortion_p,
        outliers_detected = outliers
      )
      }
  }
  
  presso_df <- rbindlist(presso_results, fill = TRUE)
  steiger_df <- rbindlist(steiger_results, fill = TRUE)
  
  presso_df[, exposure_id := as.character(exposure_id)]
  presso_df[, outcome_id := as.character(outcome_id)]
  steiger_df[, exposure_id := as.character(exposure_id)]
  steiger_df[, outcome_id := as.character(outcome_id)]
  
  fwrite(presso_df, file.path(OUTPUT_PATH, "mr_presso_results.csv"))
  fwrite(steiger_df, file.path(OUTPUT_PATH, "steiger_results.csv"))
  cat("✓ MR-PRESSO and Steiger results saved\n")
  cat(sprintf("  Processed %d pairs\n", nrow(presso_df)))
  
} else {
  cat("No SNP files found. Run Part 5 first.\n")
  presso_df <- data.table()
  steiger_df <- data.table()
}

cat("\n=== Combining Robustness Criteria ===\n")

if (exists("clean_results") && !is.null(clean_results) && nrow(clean_results) > 0) {
  
  if ("exposure" %in% names(clean_results)) {
    clean_results[, exposure := as.character(exposure)]
    setnames(clean_results, "exposure", "exposure_id")
  }
  if ("outcome" %in% names(clean_results)) {
    clean_results[, outcome := as.character(outcome)]
    setnames(clean_results, "outcome", "outcome_id")
  }
  if ("ivw_pval" %in% names(clean_results)) {
    setnames(clean_results, "ivw_pval", "IVW_P")
  }
  if ("nsnp_harmonized" %in% names(clean_results)) {
    setnames(clean_results, "nsnp_harmonized", "SNPs")
  }
  
  if (!"IVW_Beta" %in% names(clean_results) && "b_Inverse_variance_weighted" %in% names(clean_results)) {
    setnames(clean_results, "b_Inverse_variance_weighted", "IVW_Beta")
  }
  
  clean_results[, exposure_id := as.character(exposure_id)]
  clean_results[, outcome_id := as.character(outcome_id)]
  
  if (exists("presso_df") && nrow(presso_df) > 0 && exists("steiger_df") && nrow(steiger_df) > 0) {
    
    presso_df[, exposure_id := as.character(exposure_id)]
    presso_df[, outcome_id := as.character(outcome_id)]
    steiger_df[, exposure_id := as.character(exposure_id)]
    steiger_df[, outcome_id := as.character(outcome_id)]
    
    combined <- merge(clean_results, presso_df, by = c("exposure_id", "outcome_id"), all.x = TRUE)
    combined <- merge(combined, steiger_df, by = c("exposure_id", "outcome_id"), all.x = TRUE)
    
    combined$pass_fdr <- combined$pval_fdr < NOMINAL_THRESHOLD
    combined$pass_presso <- is.na(combined$presso_global_p) | combined$presso_global_p > 0.05
    combined$pass_steiger <- combined$valid_direction == TRUE
    
    combined$is_robust <- combined$pass_fdr & combined$pass_presso & combined$pass_steiger
    
    cat("\n=== Results Summary by Robustness Criteria ===\n")
    cat(sprintf("Total pairs analyzed: %d\n", nrow(combined)))
    cat(sprintf("FDR-significant (q < 0.05): %d\n", sum(combined$pass_fdr, na.rm=TRUE)))
    cat(sprintf("FDR + PRESSO: %d\n", sum(combined$pass_fdr & combined$pass_presso, na.rm=TRUE)))
    cat(sprintf("FDR + PRESSO + Steiger (ROBUST): %d\n", sum(combined$is_robust, na.rm=TRUE)))
    
    fwrite(combined, file.path(OUTPUT_PATH, "master_results_with_robustness.csv"))
    cat("✓ Master results saved with robustness flags\n")
    
    robust_pairs <- combined[is_robust == TRUE]
    
    if (nrow(robust_pairs) > 0) {
      fwrite(robust_pairs, file.path(OUTPUT_PATH, "robust_pairs_gold_standard.csv"))
      cat(sprintf("✓ %d gold standard robust pairs saved\n", nrow(robust_pairs)))
      
      cat("\n========================================\n")
      cat("   GOLD STANDARD ROBUST PAIRS\n")
      cat("   (FDR + PRESSO + Steiger)\n")
      cat("========================================\n")
      print(robust_pairs[, .(
        exposure_id, 
        outcome_id, 
        IVW_P = sprintf("%.2e", IVW_P),
        FDR_Q = sprintf("%.2e", pval_fdr),
        PRESSO_P = sprintf("%.2e", presso_global_p),
        valid_direction,
        SNPs
      )])
    } else {
      cat("\n No robust pairs found meeting all criteria.\n")
    }
    
    # Use robust pairs for visualization, or fall back to FDR-significant
    if (nrow(robust_pairs) > 0) {
      sig_results <- robust_pairs
    } else {
      sig_results <- combined[pass_fdr == TRUE]
      cat("\nUsing FDR-significant pairs for visualization (no robust pairs found)\n")
    }
    
  } else {
    cat("\n MR-PRESSO or Steiger results not found.\n")
    cat("Using FDR-significant pairs only.\n")
    sig_results <- clean_results[is_significant_fdr == TRUE]
    
    if (nrow(sig_results) == 0) {
      sig_results <- clean_results[order(IVW_P)][1:min(100, nrow(clean_results))]
    }
  }
  
} else {
  cat("No results to summarize. Run Part 6 first.\n")
  sig_results <- NULL
}


if (exists("sig_results") && !is.null(sig_results) && nrow(sig_results) > 0) {
  
  cat(sprintf("\n Generating visualizations for %d associations\n", nrow(sig_results)))
  
  # Calculate confidence intervals
  sig_results[, se_approx := abs(IVW_Beta) / qnorm(1 - IVW_P/2)]
  sig_results[, lower := IVW_Beta - 1.96 * se_approx]
  sig_results[, upper := IVW_Beta + 1.96 * se_approx]
  
  top_20 <- sig_results[order(IVW_P)][1:min(20, nrow(sig_results))]
  top_table <- data.table(
    Rank = 1:nrow(top_20),
    Outcome = top_20$outcome_id,
    Beta = round(top_20$IVW_Beta, 4),
    CI_Lower = round(top_20$lower, 4),
    CI_Upper = round(top_20$upper, 4),
    P_value = sprintf("%.2e", top_20$IVW_P),
    FDR = sprintf("%.2e", top_20$pval_fdr),
    SNPs = top_20$SNPs
  )
  
  cat("\n=== TOP 20 ASSOCIATIONS ===\n")
  print(top_table)
  fwrite(top_table, file.path(OUTPUT_PATH, "top_20_associations.csv"))
  cat("✓ Top 20 table saved\n")
  
  effect_summary <- sig_results[, .(
    Mean_Beta = mean(IVW_Beta),
    Median_Beta = median(IVW_Beta),
    SD_Beta = sd(IVW_Beta),
    Min_Beta = min(IVW_Beta),
    Max_Beta = max(IVW_Beta),
    N_Positive = sum(IVW_Beta > 0),
    N_Negative = sum(IVW_Beta < 0)
  )]
  
  cat("\n=== Effect Size Distribution ===\n")
  print(effect_summary)
  
  excel_output <- sig_results[, .(
    Exposure_ID = exposure_id,
    Outcome_ID = outcome_id,
    Beta_IVW = round(IVW_Beta, 4),
    CI_Lower = round(lower, 4),
    CI_Upper = round(upper, 4),
    P_Value = sprintf("%.2e", IVW_P),
    FDR_Q = sprintf("%.2e", pval_fdr),
    N_SNPs = SNPs
  )]
  
  cat("\n✓ Excel output created with columns:", paste(names(excel_output), collapse=", "), "\n")
  
  fwrite(excel_output, file.path(OUTPUT_PATH, "results_for_excel.csv"))
  cat("✓ Excel-friendly output saved\n")
  
  COLOR_SCHEME <- c("positive" = "#008080", "negative" = "#E64B35FF")
  
  p1 <- ggplot(sig_results, aes(x = IVW_Beta, fill = IVW_Beta > 0)) +
    geom_histogram(bins = 50, alpha = 0.7, color = "black") +
    scale_fill_manual(values = c("TRUE" = COLOR_SCHEME["positive"], 
                                 "FALSE" = COLOR_SCHEME["negative"]),
                      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                      name = "Effect Direction") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(title = "Distribution of Causal Effect Sizes (e2o Direction)",
         subtitle = paste0("N = ", nrow(sig_results), " FDR-significant associations"),
         x = "IVW Beta Estimate", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(file.path(OUTPUT_PATH, "effect_size_distribution.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(OUTPUT_PATH, "effect_size_distribution.png"), p1, width = 8, height = 6, dpi = 300)
  cat("✓ Effect size distribution saved\n")
  
  if (nrow(top_20) >= 5) {
    top_20[, label := paste0(outcome_id, " (", SNPs, " SNPs)")]
    
    p2 <- ggplot(top_20, aes(x = reorder(label, IVW_P), y = IVW_Beta)) +
      geom_pointrange(aes(ymin = lower, ymax = upper), 
                      color = "#008080", size = 0.6, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
      coord_flip() +
      labs(title = "Top 20 Most Significant Associations (e2o Direction)",
           x = "Outcome (Brain Trait)", 
           y = "Causal Effect Size (Beta)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text.y = element_text(size = 8))
    
    ggsave(file.path(OUTPUT_PATH, "forest_plot_top_20.pdf"), p2, width = 8, height = 10)
    cat("✓ Forest plot saved\n")
  }
  
  p3 <- ggplot(sig_results, aes(x = IVW_P)) +
    geom_histogram(bins = 50, fill = "#008080", color = "black", alpha = 0.7) +
    scale_x_log10() +
    labs(title = "P-value Distribution (e2o FDR-significant pairs)",
         x = "IVW P-value (log10 scale)", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  ggsave(file.path(OUTPUT_PATH, "pvalue_distribution_e2o.pdf"), p3, width = 7, height = 5)
  cat("✓ P-value distribution saved\n")
  
} else {
  cat("No results to visualize.\n")
}

cat("\n============================================================\n")
cat("           ANALYSIS COMPLETE - e2o DIRECTION\n")
cat("============================================================\n")
cat(sprintf("\n FINAL SUMMARY (e2o: Traits → Brain):\n"))
cat(sprintf("  • Total pairs analyzed: %d\n", nrow(combined)))
cat(sprintf("  • FDR-significant pairs: %d\n", sum(combined$pass_fdr, na.rm=TRUE)))
cat(sprintf("  • FDR + PRESSO: %d\n", sum(combined$pass_fdr & combined$pass_presso, na.rm=TRUE)))
cat(sprintf("  • FDR + PRESSO + Steiger (ROBUST): %d\n", sum(combined$is_robust, na.rm=TRUE)))
cat("\n INTERPRETATION:\n")
cat("  • o2e direction (Brain → Traits): ~9 robust pairs\n")
cat("  • e2o direction (Traits → Brain): 0 robust pairs\n")
cat("  • This suggests causal direction is primarily from brain to traits\n")
cat("\n All analyses completed successfully!\n")
cat("============================================================\n")


if (exists("sig_results") && nrow(sig_results) > 0) {
  
  COLOR_SCHEME <- c("positive" = "#008080", "negative" = "#E64B35FF")
  
  p1 <- ggplot(sig_results, aes(x = IVW_Beta, fill = IVW_Beta > 0)) +
    geom_histogram(bins = 50, alpha = 0.7, color = "black") +
    scale_fill_manual(values = c("TRUE" = COLOR_SCHEME["positive"], 
                                 "FALSE" = COLOR_SCHEME["negative"]),
                      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                      name = "Effect Direction") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(title = "Distribution of Causal Effect Sizes",
         subtitle = paste0("N = ", nrow(sig_results), " significant associations"),
         x = "IVW Beta Estimate", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "bottom")
  
  ggsave(file.path(OUTPUT_PATH, "effect_size_distribution.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(OUTPUT_PATH, "effect_size_distribution.png"), p1, width = 8, height = 6, dpi = 300)
  cat("✓ Effect size distribution saved\n")
  
  p2 <- ggplot(sig_results, aes(x = IVW_P)) +
    geom_histogram(bins = 50, fill = "#008080", color = "black", alpha = 0.7) +
    scale_x_log10() +
    labs(title = "P-value Distribution",
         x = "IVW P-value (log10 scale)", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  ggsave(file.path(OUTPUT_PATH, "pvalue_distribution.pdf"), p2, width = 7, height = 5)
  cat("✓ P-value distribution saved\n")
  
  if (nrow(top_20) >= 5) {
    top_20[, label := paste0(outcome_id, " (", SNPs, " SNPs)")]
    
    p3 <- ggplot(top_20, aes(x = reorder(label, IVW_P), y = IVW_Beta)) +
      geom_pointrange(aes(ymin = lower, ymax = upper), 
                      color = "#008080", size = 0.6, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
      coord_flip() +
      labs(title = "Top 20 Most Significant Associations",
           x = "Outcome (Brain Trait)", 
           y = "Causal Effect Size (Beta)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text.y = element_text(size = 8))
    
    ggsave(file.path(OUTPUT_PATH, "forest_plot_top_20.pdf"), p3, width = 8, height = 10)
    cat("✓ Forest plot saved\n")
  }
  
  excel_output <- sig_results[, .(
    Outcome_ID = as.character(outcome_id),
    Beta_IVW = round(IVW_Beta, 4),
    CI_Lower = round(lower, 4),
    CI_Upper = round(upper, 4),
    P_Value = sprintf("%.2e", IVW_P),
    FDR_Q = sprintf("%.2e", pval_fdr),
    N_SNPs = SNPs
  )]
  
  fwrite(excel_output, file.path(OUTPUT_PATH, "results_for_excel.csv"))
  cat("✓ Excel-friendly output saved\n")
  
} else {
  cat("No results to plot. Run Part 7 first.\n")
}


cat("\n============================================================\n")
cat("                    PIPELINE COMPLETE\n")
cat("============================================================\n")
cat(sprintf("\n All outputs saved to: %s\n", OUTPUT_PATH))
cat("\nGenerated files:\n")
cat("  • mr_results_cleaned.csv - Cleaned results\n")
cat("  • mr_results_complete_with_fdr.csv - Results with FDR\n")
cat("  • fdr_significant_pairs.csv - FDR-significant pairs\n")
cat("  • top_20_significant.csv - Top 20 table\n")
cat("  • results_for_excel.csv - Excel-friendly format\n")
cat("  • effect_size_distribution.pdf/png - Histogram\n")
cat("  • pvalue_distribution.pdf - P-value distribution\n")
cat("  • forest_plot_top_20.pdf - Forest plot\n")
cat("============================================================\n")






library(data.table)
library(ggplot2)

E2O_FILE <- "/mnt/f/result/mr_results_e2o/master_results_with_robustness.csv"
O2E_FILE <- "/mnt/f/result/mr_results_o2e/master_results_with_robustness.csv"

e2o_results <- fread(E2O_FILE)
o2e_results <- fread(O2E_FILE)

cat("=== LOADED DATA ===\n")
cat(sprintf("e2o (Traits → Brain): %d pairs\n", nrow(e2o_results)))
cat(sprintf("o2e (Brain → Traits): %d pairs\n\n", nrow(o2e_results)))


comparison_summary <- data.table(
  Metric = c(
    "Total pairs analyzed",
    "",
    "FDR-significant (q < 0.05)",
    "FDR + PRESSO (p > 0.05)",
    "FDR + Steiger (valid direction)",
    "FDR + PRESSO + Steiger (ROBUST)",
    "",
    "Steiger valid direction (all pairs)",
    "PRESSO global test p < 0.05",
    "PRESSO outliers detected"
  ),
  e2o_Traits_to_Brain = c(
    sprintf("%d", nrow(e2o_results)),
    "",
    sprintf("%d (%.1f%%)", sum(e2o_results$pass_fdr, na.rm = TRUE), 
            100 * sum(e2o_results$pass_fdr, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$pass_fdr & e2o_results$pass_presso, na.rm = TRUE), 
            100 * sum(e2o_results$pass_fdr & e2o_results$pass_presso, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$pass_fdr & e2o_results$pass_steiger, na.rm = TRUE), 
            100 * sum(e2o_results$pass_fdr & e2o_results$pass_steiger, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$is_robust, na.rm = TRUE), 
            100 * sum(e2o_results$is_robust, na.rm = TRUE) / nrow(e2o_results)),
    "",
    sprintf("%d (%.1f%%)", sum(e2o_results$valid_direction == TRUE, na.rm = TRUE), 
            100 * sum(e2o_results$valid_direction == TRUE, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$presso_global_p < 0.05, na.rm = TRUE), 
            100 * sum(e2o_results$presso_global_p < 0.05, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$outliers_detected == TRUE, na.rm = TRUE), 
            100 * sum(e2o_results$outliers_detected == TRUE, na.rm = TRUE) / nrow(e2o_results))
  ),
  o2e_Brain_to_Traits = c(
    sprintf("%d", nrow(o2e_results)),
    "",
    sprintf("%d (%.1f%%)", sum(o2e_results$pass_fdr, na.rm = TRUE), 
            100 * sum(o2e_results$pass_fdr, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$pass_fdr & o2e_results$pass_presso, na.rm = TRUE), 
            100 * sum(o2e_results$pass_fdr & o2e_results$pass_presso, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$pass_fdr & o2e_results$pass_steiger, na.rm = TRUE), 
            100 * sum(o2e_results$pass_fdr & o2e_results$pass_steiger, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$is_robust, na.rm = TRUE), 
            100 * sum(o2e_results$is_robust, na.rm = TRUE) / nrow(o2e_results)),
    "",
    sprintf("%d (%.1f%%)", sum(o2e_results$valid_direction == TRUE, na.rm = TRUE), 
            100 * sum(o2e_results$valid_direction == TRUE, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$presso_global_p < 0.05, na.rm = TRUE), 
            100 * sum(o2e_results$presso_global_p < 0.05, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$outliers_detected == TRUE, na.rm = TRUE), 
            100 * sum(o2e_results$outliers_detected == TRUE, na.rm = TRUE) / nrow(o2e_results))
  )
)

cat("\n=== COMPARISON SUMMARY: e2o vs o2e ===\n")
print(comparison_summary, row.names = FALSE)

fwrite(comparison_summary, "/mnt/f/result/comparison_summary_e2o_vs_o2e.csv")
cat("\n✓ Comparison summary saved\n")


create_detailed_summary <- function(results, direction_name) {
  
  cat(sprintf("\n=== %s ===\n", direction_name))
  
  cat(sprintf("\nTotal pairs: %d\n", nrow(results)))
  cat(sprintf("FDR-significant: %d (%.1f%%)\n", 
              sum(results$pass_fdr, na.rm = TRUE), 
              100 * sum(results$pass_fdr, na.rm = TRUE) / nrow(results)))
  cat(sprintf("Robust (FDR + PRESSO + Steiger): %d (%.1f%%)\n", 
              sum(results$is_robust, na.rm = TRUE), 
              100 * sum(results$is_robust, na.rm = TRUE) / nrow(results)))
  
  cat(sprintf("PRESSO global test p < 0.05: %d (%.1f%%)\n", 
              sum(results$presso_global_p < 0.05, na.rm = TRUE), 
              100 * sum(results$presso_global_p < 0.05, na.rm = TRUE) / nrow(results)))
  cat(sprintf("PRESSO outliers detected: %d (%.1f%%)\n", 
              sum(results$outliers_detected == TRUE, na.rm = TRUE), 
              100 * sum(results$outliers_detected == TRUE, na.rm = TRUE) / nrow(results)))
  
  cat(sprintf("\nEffect sizes (IVW Beta):\n"))
  cat(sprintf("  Mean: %.4f\n", mean(results$IVW_Beta, na.rm = TRUE)))
  cat(sprintf("  Median: %.4f\n", median(results$IVW_Beta, na.rm = TRUE)))
  cat(sprintf("  Positive effects: %d (%.1f%%)\n", 
              sum(results$IVW_Beta > 0, na.rm = TRUE), 
              100 * sum(results$IVW_Beta > 0, na.rm = TRUE) / nrow(results)))
  
  summary_table <- data.table(
    Metric = c(
      "Total pairs",
      "FDR-significant",
      "FDR + PRESSO",
      "FDR + Steiger",
      "ROBUST (FDR + PRESSO + Steiger)",
      "Steiger valid (all pairs)",
      "PRESSO significant (p < 0.05)",
      "Outliers detected",
      "Mean IVW Beta",
      "Median IVW Beta",
      "Positive effects (%)"
    ),
    Count = c(
      sprintf("%d", nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$pass_fdr, na.rm = TRUE), 100 * sum(results$pass_fdr, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$pass_fdr & results$pass_presso, na.rm = TRUE), 100 * sum(results$pass_fdr & results$pass_presso, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$pass_fdr & results$pass_steiger, na.rm = TRUE), 100 * sum(results$pass_fdr & results$pass_steiger, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$is_robust, na.rm = TRUE), 100 * sum(results$is_robust, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$valid_direction == TRUE, na.rm = TRUE), 100 * sum(results$valid_direction == TRUE, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$presso_global_p < 0.05, na.rm = TRUE), 100 * sum(results$presso_global_p < 0.05, na.rm = TRUE) / nrow(results)),
      sprintf("%d (%.1f%%)", sum(results$outliers_detected == TRUE, na.rm = TRUE), 100 * sum(results$outliers_detected == TRUE, na.rm = TRUE) / nrow(results)),
      sprintf("%.4f", mean(results$IVW_Beta, na.rm = TRUE)),
      sprintf("%.4f", median(results$IVW_Beta, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(results$IVW_Beta > 0, na.rm = TRUE), 100 * sum(results$IVW_Beta > 0, na.rm = TRUE) / nrow(results))
    )
  )
  
  return(summary_table)
}

# Create detailed summaries
e2o_detailed <- create_detailed_summary(e2o_results, "e2o (Traits → Brain)")
o2e_detailed <- create_detailed_summary(o2e_results, "o2e (Brain → Traits)")

# Save detailed summaries
fwrite(e2o_detailed, "/mnt/f/result/mr_results_e2o/detailed_summary.csv")
fwrite(o2e_detailed, "/mnt/f/result/mr_results_o2e/detailed_summary.csv")

# ============================================================
# 3. ROBUST PAIRS FOR EACH DIRECTION
# ============================================================

# e2o robust pairs
e2o_robust <- e2o_results[is_robust == TRUE]
if (nrow(e2o_robust) > 0) {
  cat("\n=== e2o ROBUST PAIRS ===\n")
  print(e2o_robust[, .(exposure_id, outcome_id, IVW_P = sprintf("%.2e", IVW_P), 
                       pval_fdr = sprintf("%.2e", pval_fdr), SNPs)])
  fwrite(e2o_robust, "/mnt/f/result/mr_results_e2o/robust_pairs_list.csv")
} else {
  cat("\n=== e2o: No robust pairs found ===\n")
}

# o2e robust pairs
o2e_robust <- o2e_results[is_robust == TRUE]
if (nrow(o2e_robust) > 0) {
  cat("\n=== o2e ROBUST PAIRS ===\n")
  print(o2e_robust[, .(exposure_id, outcome_id, IVW_P = sprintf("%.2e", IVW_P), 
                       pval_fdr = sprintf("%.2e", pval_fdr), SNPs)])
  fwrite(o2e_robust, "/mnt/f/result/mr_results_o2e/robust_pairs_list.csv")
} else {
  cat("\n=== o2e: No robust pairs found ===\n")
}

# ============================================================
# 4. P-VALUE DISTRIBUTION COMPARISON
# ============================================================

# Create p-value distribution table
pval_distribution <- data.table(
  P_value_Range = c("< 1e-10", "1e-10 to 1e-8", "1e-8 to 1e-6", 
                    "1e-6 to 1e-4", "1e-4 to 0.001", "0.001 to 0.01", 
                    "0.01 to 0.05", "> 0.05"),
  e2o_Count = c(
    sum(e2o_results$IVW_P < 1e-10, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 1e-10 & e2o_results$IVW_P < 1e-8, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 1e-8 & e2o_results$IVW_P < 1e-6, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 1e-6 & e2o_results$IVW_P < 1e-4, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 1e-4 & e2o_results$IVW_P < 0.001, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 0.001 & e2o_results$IVW_P < 0.01, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 0.01 & e2o_results$IVW_P < 0.05, na.rm = TRUE),
    sum(e2o_results$IVW_P >= 0.05, na.rm = TRUE)
  ),
  o2e_Count = c(
    sum(o2e_results$IVW_P < 1e-10, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 1e-10 & o2e_results$IVW_P < 1e-8, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 1e-8 & o2e_results$IVW_P < 1e-6, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 1e-6 & o2e_results$IVW_P < 1e-4, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 1e-4 & o2e_results$IVW_P < 0.001, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 0.001 & o2e_results$IVW_P < 0.01, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 0.01 & o2e_results$IVW_P < 0.05, na.rm = TRUE),
    sum(o2e_results$IVW_P >= 0.05, na.rm = TRUE)
  )
)

cat("\n=== P-VALUE DISTRIBUTION COMPARISON ===\n")
print(pval_distribution)
fwrite(pval_distribution, "/mnt/f/result/pvalue_distribution_comparison.csv")

# ============================================================
# 5. EFFECT SIZE COMPARISON
# ============================================================

effect_comparison <- data.table(
  Metric = c("Mean Beta", "SD Beta", "Median Beta", "Min Beta", "Max Beta", 
             "Positive (%)", "Negative (%)"),
  e2o = c(
    sprintf("%.4f", mean(e2o_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", sd(e2o_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", median(e2o_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", min(e2o_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", max(e2o_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%d (%.1f%%)", sum(e2o_results$IVW_Beta > 0, na.rm = TRUE), 
            100 * sum(e2o_results$IVW_Beta > 0, na.rm = TRUE) / nrow(e2o_results)),
    sprintf("%d (%.1f%%)", sum(e2o_results$IVW_Beta < 0, na.rm = TRUE), 
            100 * sum(e2o_results$IVW_Beta < 0, na.rm = TRUE) / nrow(e2o_results))
  ),
  o2e = c(
    sprintf("%.4f", mean(o2e_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", sd(o2e_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", median(o2e_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", min(o2e_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%.4f", max(o2e_results$IVW_Beta, na.rm = TRUE)),
    sprintf("%d (%.1f%%)", sum(o2e_results$IVW_Beta > 0, na.rm = TRUE), 
            100 * sum(o2e_results$IVW_Beta > 0, na.rm = TRUE) / nrow(o2e_results)),
    sprintf("%d (%.1f%%)", sum(o2e_results$IVW_Beta < 0, na.rm = TRUE), 
            100 * sum(o2e_results$IVW_Beta < 0, na.rm = TRUE) / nrow(o2e_results))
  )
)

cat("\n=== EFFECT SIZE COMPARISON ===\n")
print(effect_comparison)
fwrite(effect_comparison, "/mnt/f/result/effect_size_comparison.csv")

# ============================================================
# 6. FINAL SUMMARY REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("           FINAL COMPARATIVE SUMMARY\n")
cat("============================================================\n")
cat("\n📊 DIRECTIONAL ANALYSIS SUMMARY:\n")
cat(sprintf("\n  e2o (Traits → Brain):\n"))
cat(sprintf("    • Total pairs: %d\n", nrow(e2o_results)))
cat(sprintf("    • FDR-significant: %d (%.1f%%)\n", 
            sum(e2o_results$pass_fdr, na.rm = TRUE), 
            100 * sum(e2o_results$pass_fdr, na.rm = TRUE) / nrow(e2o_results)))
cat(sprintf("    • ROBUST pairs: %d (%.1f%%)\n", 
            sum(e2o_results$is_robust, na.rm = TRUE), 
            100 * sum(e2o_results$is_robust, na.rm = TRUE) / nrow(e2o_results)))
cat(sprintf("    • Steiger valid: %d (%.1f%%)\n", 
            sum(e2o_results$valid_direction == TRUE, na.rm = TRUE), 
            100 * sum(e2o_results$valid_direction == TRUE, na.rm = TRUE) / nrow(e2o_results)))

cat(sprintf("\n  o2e (Brain → Traits):\n"))
cat(sprintf("    • Total pairs: %d\n", nrow(o2e_results)))
cat(sprintf("    • FDR-significant: %d (%.1f%%)\n", 
            sum(o2e_results$pass_fdr, na.rm = TRUE), 
            100 * sum(o2e_results$pass_fdr, na.rm = TRUE) / nrow(o2e_results)))
cat(sprintf("    • ROBUST pairs: %d (%.1f%%)\n", 
            sum(o2e_results$is_robust, na.rm = TRUE), 
            100 * sum(o2e_results$is_robust, na.rm = TRUE) / nrow(o2e_results)))
cat(sprintf("    • Steiger valid: %d (%.1f%%)\n", 
            sum(o2e_results$valid_direction == TRUE, na.rm = TRUE), 
            100 * sum(o2e_results$valid_direction == TRUE, na.rm = TRUE) / nrow(o2e_results)))

cat("\n🔍 KEY FINDINGS:\n")
if (sum(o2e_results$is_robust, na.rm = TRUE) > sum(e2o_results$is_robust, na.rm = TRUE)) {
  cat("  • More robust associations found in o2e direction (Brain → Traits)\n")
  cat("  • This suggests causal direction is primarily from brain to traits\n")
} else if (sum(e2o_results$is_robust, na.rm = TRUE) > sum(o2e_results$is_robust, na.rm = TRUE)) {
  cat("  • More robust associations found in e2o direction (Traits → Brain)\n")
} else {
  cat("  • Similar number of robust associations in both directions\n")
}

cat(sprintf("\n  • Ratio of robust pairs (o2e/e2o): %.2f\n", 
            ifelse(sum(e2o_results$is_robust, na.rm = TRUE) > 0,
                   sum(o2e_results$is_robust, na.rm = TRUE) / sum(e2o_results$is_robust, na.rm = TRUE),
                   NA)))

cat("\n✅ All summary files saved to /mnt/f/result/\n")
cat("   • comparison_summary_e2o_vs_o2e.csv\n")
cat("   • pvalue_distribution_comparison.csv\n")
cat("   • effect_size_comparison.csv\n")
cat("   • mr_results_e2o/detailed_summary.csv\n")
cat("   • mr_results_o2e/detailed_summary.csv\n")
cat("============================================================\n")