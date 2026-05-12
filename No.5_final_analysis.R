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
  cat("No results file found. Please run part 5 first.\n")
  clean_results <- NULL
}


library(MRPRESSO)

cat("\n=== Running MR-PRESSO ===\n")

snp_files <- list.files(SNP_OUTPUT_PATH, pattern = "_snps.csv$", full.names = TRUE)

if (length(snp_files) > 0) {
  
  cat(sprintf("Found %d SNP files to process\n", length(snp_files)))
  
  presso_results <- list()
  
  for (i in seq_along(snp_files)) {
    if (i %% 100 == 0) cat(sprintf("  Processing %d/%d...\n", i, length(snp_files)))
    
    snp_data <- fread(snp_files[i])
    
    fname <- basename(snp_files[i])
    fname <- gsub("_snps.csv$", "", fname)
    parts <- strsplit(fname, "_to_")[[1]]
    
    if (length(parts) == 2) {
      exposure_id <- as.character(parts[1])
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
  
  presso_df[, exposure_id := as.character(exposure_id)]
  presso_df[, outcome_id := as.character(outcome_id)]
  
  fwrite(presso_df, file.path(OUTPUT_PATH, "mr_presso_results.csv"))
  cat("✓ MR-PRESSO results saved\n")
  cat(sprintf("  Processed %d pairs\n", nrow(presso_df)))
  
} else {
  cat("No SNP files found. Run part 5 first.\n")
  presso_df <- data.table()
}

cat("\n=== Combining Results ===\n")

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
  
  if (exists("presso_df") && nrow(presso_df) > 0) {
    
    presso_df[, exposure_id := as.character(exposure_id)]
    presso_df[, outcome_id := as.character(outcome_id)]
    
    combined <- merge(clean_results, presso_df, by = c("exposure_id", "outcome_id"), all.x = TRUE)
    
    combined$pass_fdr <- combined$pval_fdr < NOMINAL_THRESHOLD
    combined$pass_presso <- is.na(combined$presso_global_p) | combined$presso_global_p > 0.05
    combined$is_robust <- combined$pass_fdr & combined$pass_presso
    
    cat("\n=== Results Summary ===\n")
    cat(sprintf("Total pairs analyzed: %d\n", nrow(combined)))
    cat(sprintf("FDR-significant (q < 0.05): %d\n", sum(combined$pass_fdr, na.rm=TRUE)))
    cat(sprintf("FDR + PRESSO (ROBUST): %d\n", sum(combined$is_robust, na.rm=TRUE)))
    
    fwrite(combined, file.path(OUTPUT_PATH, "master_results_with_robustness.csv"))
    cat("✓ Master results saved\n")
    
    robust_pairs <- combined[is_robust == TRUE]
    
    if (nrow(robust_pairs) > 0) {
      fwrite(robust_pairs, file.path(OUTPUT_PATH, "robust_pairs_gold_standard.csv"))
      cat(sprintf("✓ %d robust pairs saved\n", nrow(robust_pairs)))
      
      cat("\n========================================\n")
      cat("   ROBUST PAIRS (FDR + PRESSO)\n")
      cat("========================================\n")
      print(robust_pairs[, .(
        exposure_id, 
        outcome_id, 
        IVW_P = sprintf("%.2e", IVW_P),
        FDR_Q = sprintf("%.2e", pval_fdr),
        PRESSO_P = sprintf("%.2e", presso_global_p),
        SNPs
      )])
    } else {
      cat("\n No robust pairs found meeting all criteria.\n")
    }
    
    if (nrow(robust_pairs) > 0) {
      sig_results <- robust_pairs
    } else {
      sig_results <- combined[pass_fdr == TRUE]
      cat("\nUsing FDR-significant pairs for visualization (no robust pairs found)\n")
    }
    
  } else {
    cat("\n MR-PRESSO results not found.\n")
    cat("Using FDR-significant pairs only.\n")
    sig_results <- clean_results[is_significant_fdr == TRUE]
    
    if (nrow(sig_results) == 0) {
      sig_results <- clean_results[order(IVW_P)][1:min(100, nrow(clean_results))]
    }
  }
  
} else {
  cat("No results to summarize. Run part 6 first.\n")
  sig_results <- NULL
}


if (exists("sig_results") && !is.null(sig_results) && nrow(sig_results) > 0) {
  
  cat(sprintf("\n Generating visualizations for %d associations\n", nrow(sig_results)))
  
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
  
  cat("\n✓ Excel output created\n")
  
  fwrite(excel_output, file.path(OUTPUT_PATH, "results_for_excel.csv"))
  cat("✓ Excel-friendly output saved\n")
  
} else {
  cat("No results to visualize.\n")
}

cat("\n============================================================\n")
cat("           ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("\n All outputs saved to: %s\n", OUTPUT_PATH))
cat("\nGenerated files:\n")
cat("  • mr_results_cleaned.csv\n")
cat("  • mr_results_complete_with_fdr.csv\n")
cat("  • fdr_significant_pairs.csv\n")
cat("  • mr_presso_results.csv\n")
cat("  • master_results_with_robustness.csv\n")
cat("  • robust_pairs_gold_standard.csv\n")
cat("  • top_20_associations.csv\n")
cat("  • results_for_excel.csv\n")
cat("============================================================\n")