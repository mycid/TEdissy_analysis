# scripts/process_plate_run.R
required_packages <- c(
  "readxl", "dplyr", "tidyr", "tidyverse", "purrr", "zoo",
  "ggplot2", "lubridate", "stringr", "tools", "plotly", 
  "here", "broom", "easystats", "performance", "ggrepel", 
  "effectsize", "rlang", "boot"
)
#
#
#
for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
#
#
#
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(zoo)
library(readxl)
library(performance)
library(broom)
library(easystats)
library(ggplot2)
library(lubridate)
library(stringr)
library(tools)
library(plotly)
library(tidyverse)
library(here)
library(ggrepel)
library(ggpubr)
#
#
load_excel_files <- function(input_folder) {
  
  input_path <- here::here(input_folder)
  
  excel_files <- list.files(path = input_path, pattern = "\\.xlsx?$", full.names = TRUE)
  
  if (length(excel_files) == 0) {
    warning("No Excel files found in ", input_path)
    return(invisible(NULL))
  }
  
  for (file in excel_files) {
    name <- tools::file_path_sans_ext(basename(file))
    
    sheet_names <- readxl::excel_sheets(file)
    sheet_to_read <- if (length(sheet_names) >= 2) sheet_names[2] else sheet_names[1]
    
    data <- suppressMessages(suppressWarnings(
      readxl::read_excel(file, sheet = sheet_to_read)
    ))
    
    assign(name, data, envir = .GlobalEnv)
  }
  
  message("Finished loading ", length(excel_files), " Excel files from ", input_path)
}
#
#
#
#
#
#
#
tidy_and_correct <- function(input_df, blanks, df_name = NULL) {
  # Check input
  if (!is.data.frame(input_df)) stop("Input must be a dataframe.")
  
  # Convert all "NA" strings to real NA
  input_df[input_df == "NA"] <- NA
  
  # Find first "A" in first two columns
  col_A_match <- which(input_df[[1]] == "A")
  if (length(col_A_match) == 0 && ncol(input_df) >= 2) {
    col_A_match <- which(input_df[[2]] == "A")
    if (length(col_A_match) > 0) {
      row_start <- col_A_match[1]
      row_col <- 2
    } else {
      stop("Could not find the first data row (looking for 'A' in first two columns).")
    }
  } else {
    row_start <- col_A_match[1]
    row_col <- 1
  }
  
  df <- input_df[row_start:nrow(input_df), ]
  
  if (row_col != 1) {
    df <- df[, c(row_col, setdiff(seq_along(df), row_col))]
  }
  
  df[[1]] <- as.character(df[[1]])
  df[[1]] <- zoo::na.locf(df[[1]])
  
  wl_column <- df[[ncol(df)]]
  wl_vector <- unique(na.omit(wl_column)) %>% stringr::str_trim() %>% paste0("X", .)
  num_wl <- length(wl_vector)
  
  df <- df[, 1:13]
  colnames(df)[1:13] <- c("Row", as.character(1:12))
  df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], function(x) as.numeric(as.character(x)))
  
  tidy_df <- df %>%
    tidyr::pivot_longer(cols = -Row, names_to = "Column", values_to = "Value",
                        values_transform = list(Value = as.numeric)) %>%
    dplyr::mutate(
      Column = as.integer(Column),
      Cell_ID = sprintf("%s%02d", Row, Column)
    ) %>%
    dplyr::select(Cell_ID, Value)
  
  if (num_wl > 1) {
    if (nrow(tidy_df) %% num_wl != 0) {
      stop("Data rows are not divisible by number of wavelengths. Check format.")
    }
    tidy_df <- tidy_df %>%
      dplyr::group_by(Cell_ID) %>%
      dplyr::mutate(Wavelength = rep(wl_vector, times = n() / num_wl)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = Wavelength, values_from = Value)
  } else {
    colnames(tidy_df)[2] <- wl_vector
  }
  
  blanks_df <- tidy_df %>% dplyr::filter(Cell_ID %in% blanks)
  blank_means <- blanks_df %>% dplyr::select(-Cell_ID) %>%
    dplyr::summarise(dplyr::across(everything(), ~mean(.x, na.rm = TRUE)))
  
  tidy_corrected <- tidy_df %>%
    dplyr::mutate(dplyr::across(-Cell_ID, ~ .x - unlist(blank_means)))
  
  spectral_cols <- names(tidy_corrected)[-1]
  negative_flags <- apply(tidy_corrected[, spectral_cols], 1, function(row) any(row < 0, na.rm = TRUE))
  
  removed_rows <- tidy_corrected[negative_flags, ]
  if (nrow(removed_rows) > 0) {
    message("Removed rows due to negative values:")
    print(removed_rows$Cell_ID)
  }
  
  tidy_filtered <- tidy_corrected[!negative_flags, ]
  attr(tidy_filtered, "removed_rows_spectral") <- removed_rows
  
  # Save to global environment if df_name provided
  if (!is.null(df_name)) {
    assign(paste0(df_name, "_tidy"), tidy_filtered, envir = .GlobalEnv)
  }
  
  return(tidy_filtered)
}

#
#
# 
tidy_all <- function(df_list, blanks_list) {
  if (length(df_list) != length(blanks_list)) {
    stop("Length of dataframe list and blanks list must match.")
  }
  
  for (i in seq_along(df_list)) {
    df_name <- names(df_list)[i]
    blanks <- blanks_list[[i]]
    df <- df_list[[i]]
    
    cat("Processing:", df_name, "with blanks:", toString(blanks), "\n")
    
    # Call tidy_and_correct, capture output
    tidy_df <- tidy_and_correct(df, blanks = blanks)
    
    # Assign result to global env with _tidy appended
    assign(paste0(df_name, "_tidy"), tidy_df, envir = .GlobalEnv)
    cat("Saved:", paste0(df_name, "_tidy"), "to global environment.\n")
  }
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#

  #
  #
  #
calculate_pe_and_filter <- function(tidy_df, 
                                    output_basename,                  # <- NEW: name for files and objects
                                    extract_volume_mL = 1.5, 
                                    sample_weight_col = "sample weights",
                                    sample_id_col = "Cell_ID",
                                    simulate_lower_pct = 0) {
  
  # --- Set up directories ---
  filtered_dir <- file.path("output PE", "export data", "pe filtered")
  removed_dir  <- file.path("output PE", "export data", "pe removed")
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(removed_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- File paths based on provided output_basename ---
  filtered_file <- file.path(filtered_dir, paste0(output_basename, ".csv"))
  removed_file  <- file.path(removed_dir, paste0(output_basename, "_removed_rows.csv"))
  
  # --- Step 0: Optional simulation ---
  if(simulate_lower_pct > 0) {
    tidy_df <- tidy_df %>%
      mutate(
        X455 = X455 * (1 - simulate_lower_pct / 100),
        X592 = X592 * (1 - simulate_lower_pct / 100)
      )
    message(paste("Simulated", simulate_lower_pct, "% lower values for wavelengths 455 and 592."))
  }
  
  # --- Step 1: Calculate PE in mg/mL ---
  tidy_df <- tidy_df %>%
    mutate(
      PE_mg_per_mL = ((`X564` - `X592`) - ((`X455` - `X592`) * 0.20)) * 0.12
    )
  
  # --- Step 2: Remove rows with PE < 0 ---
  to_remove_pe <- which(tidy_df$PE_mg_per_mL < 0)
  removed_rows <- data.frame()
  
  if(length(to_remove_pe) > 0) {
    removed_rows <- tidy_df[to_remove_pe, ]
    removed_rows$Removal_Reason <- "removed because PE < 0"
    
    write.csv(removed_rows, removed_file, row.names = FALSE)
    assign(paste0(output_basename, "_removed"), removed_rows, envir = .GlobalEnv)
    
    message("Removed observations due to negative PE values:")
    print(removed_rows %>% select(all_of(sample_id_col), PE_mg_per_mL, Removal_Reason))
  } else {
    message("No observations removed due to negative PE values.")
  }
  
  # --- Step 3–4: Calculate total PE and PE per g sample ---
  tidy_final <- tidy_df[-to_remove_pe, ] %>%
    mutate(
      total_PE_mg = PE_mg_per_mL * extract_volume_mL,
      PE_mg_per_g_sample = total_PE_mg / (.data[[sample_weight_col]] * 0.001)
    )
  
  # --- Save filtered output and assign globally ---
  write.csv(tidy_final, filtered_file, row.names = FALSE)
  assign(output_basename, tidy_final, envir = .GlobalEnv)
  
  message("Filtered data saved to: ", filtered_file)
  
  # Attach removed rows as attribute and return
  attr(tidy_final, "removed_rows_pe") <- removed_rows
  return(tidy_final)
}



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
joindf_by_id <- function(df1, df2, 
                         output_name,          # full path to joined CSV
                         unmatched_out,        # full path to unmatched CSV
                         key_df1, 
                         key_df2, 
                         assign_to_global = TRUE) {
  
  # Clean column names
  colnames(df1) <- trimws(colnames(df1))
  colnames(df2) <- trimws(colnames(df2))
  
  # Validate keys
  if (!(key_df1 %in% colnames(df1))) stop(paste0("Error: '", key_df1, "' not found in df1."))
  if (!(key_df2 %in% colnames(df2))) stop(paste0("Error: '", key_df2, "' not found in df2."))
  
  # Rename join columns
  df1_renamed <- df1 %>% dplyr::rename(join_id = all_of(key_df1))
  df2_renamed <- df2 %>% dplyr::rename(join_id = all_of(key_df2))
  
  # Join and identify unmatched
  if (nrow(df1_renamed) <= nrow(df2_renamed)) {
    result <- dplyr::left_join(df1_renamed, df2_renamed, by = "join_id")
    base_ids <- df1_renamed$join_id
    compare_ids <- df2_renamed$join_id
    unmatched_df <- dplyr::filter(df2_renamed, !(join_id %in% base_ids))
  } else {
    result <- dplyr::left_join(df2_renamed, df1_renamed, by = "join_id")
    base_ids <- df2_renamed$join_id
    compare_ids <- df1_renamed$join_id
    unmatched_df <- dplyr::filter(df1_renamed, !(join_id %in% base_ids))
  }
  
  # Save results
  write.csv(result, output_name, row.names = FALSE)
  write.csv(unmatched_df, unmatched_out, row.names = FALSE)
  
  # Assign to global environment
  object_name <- tools::file_path_sans_ext(basename(output_name))
  if (assign_to_global) {
    assign(object_name, result, envir = .GlobalEnv)
  }
  
  # Print unmatched join_ids
  repeat_list <- unique(unmatched_df$join_id)
  cat("\nValues needing repetition:\n")
  cat(paste0(repeat_list, " needs to be repeated\n"), sep = "")
  
  # Summary report
  report <- list(
    result_df = result,
    merged_cells = sum(compare_ids %in% base_ids),
    unmatched_cells = nrow(unmatched_df),
    unmatched_saved_to = unmatched_out,
    assigned_to_global = assign_to_global
  )
  
  message("Join complete. Output saved to: ", output_name)
  return(report)
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# Helper function to select the best 3 replicates minimizing CV
select_best_three <- function(df) {
  n <- nrow(df)
  rows <- seq_len(n)
  
  # If 3 or fewer replicates, include all, exclude none (NA)
  if (n <= 3) {
    return(list(
      included_rows = rows,
      excluded_rows = NA_character_
    ))
  }
  
  # For more than 3, find best 3 replicates minimizing CV
  combos <- combn(rows, 3)
  cvs <- apply(combos, 2, function(idx) {
    vals <- df$value[idx]
    vals <- vals[!is.na(vals)]
    # Avoid division by zero or insufficient data
    if (length(vals) < 3 || mean(vals, na.rm = TRUE) == 0) return(Inf)
    sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE)
  })
  
  best_idx <- which.min(cvs)
  best_rows <- combos[, best_idx]
  excluded_rows <- setdiff(rows, best_rows)
  
  return(list(
    included_rows = best_rows,
    excluded_rows = excluded_rows
  ))
}
#
#
#
#
#
#
#
#
#
#
#
#
analyze_replicates <- function(data,
                               id_col = "ID",
                               join_col = "join_id",
                               weight_col = "sample weight",
                               date_col = "date",
                               output_prefix = "replicate_analysis",
                               choose_best_3 = FALSE) {
  
  # Ensure proper types
  data[[id_col]] <- as.factor(data[[id_col]])
  data[[date_col]] <- as.Date(data[[date_col]])
  
  # Identify numeric columns to analyze (exclude metadata)
  numeric_cols <- data %>%
    select(where(is.numeric)) %>%
    select(-any_of(c(weight_col))) %>%  # exclude sample weight
    colnames()
  
  # Pivot to long format for analysis
  long_data <- data %>%
    pivot_longer(cols = all_of(numeric_cols), names_to = "measure", values_to = "value")
  
  # Summary stats with best 3 selection if enabled
  summary_long <- long_data %>%
    group_by(.data[[id_col]], measure) %>%
    group_modify(~ {
      df_sub <- .x
      n <- nrow(df_sub)
      
      if (choose_best_3 && n > 3) {
        selected <- select_best_three(df_sub)
        included_rows <- selected$included_rows
        excluded_rows <- selected$excluded_rows
        
        # Filter to best 3 for stats
        df_best <- df_sub[included_rows, , drop = FALSE]
        
        # Summarize stats on best 3
        mean_val <- mean(df_best$value, na.rm = TRUE)
        sd_val <- sd(df_best$value, na.rm = TRUE)
        se_val <- sd_val / sqrt(nrow(df_best))
        cv_val <- sd_val / mean_val
        max_dev_pct_val <- max(abs(df_best$value - mean_val) / mean_val * 100, na.rm = TRUE)
        
        # Included/excluded rows as comma-separated strings (or NA)
        included_str <- paste(included_rows, collapse = ",")
        excluded_str <- if (all(is.na(excluded_rows))) NA_character_ else paste(excluded_rows, collapse = ",")
        
        tibble(
          n = nrow(df_best),
          mean = mean_val,
          sd = sd_val,
          se = se_val,
          cv = cv_val,
          max_dev_pct = max_dev_pct_val,
          included_rows = included_str,
          excluded_rows = excluded_str
        )
      } else {
        # Normal stats on all replicates (no best 3)
        mean_val <- mean(df_sub$value, na.rm = TRUE)
        sd_val <- sd(df_sub$value, na.rm = TRUE)
        se_val <- sd_val / sqrt(n)
        cv_val <- sd_val / mean_val
        max_dev_pct_val <- max(abs(df_sub$value - mean_val) / mean_val * 100, na.rm = TRUE)
        
        tibble(
          n = n,
          mean = mean_val,
          sd = sd_val,
          se = se_val,
          cv = cv_val,
          max_dev_pct = max_dev_pct_val
        )
      }
    }) %>%
    ungroup()
  
  # Pivot stats to wide format for readability
  if (choose_best_3) {
    summary_wide <- summary_long %>%
      pivot_wider(
        names_from = measure,
        values_from = c(mean, sd, se, cv, max_dev_pct, n, included_rows, excluded_rows),
        names_glue = "{measure}_{.value}"
      )
  } else {
    summary_wide <- summary_long %>%
      pivot_wider(
        names_from = measure,
        values_from = c(mean, sd, se, cv, max_dev_pct, n),
        names_glue = "{measure}_{.value}"
      )
  }
 
  
  # Add join_id list, average sample weight, date, and replicate count
  replicate_info <- data %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      join_ids = paste(unique(.data[[join_col]]), collapse = ", "),
      avg_sample_weight = mean(.data[[weight_col]], na.rm = TRUE),
      replicate_date = ifelse(length(unique(.data[[date_col]])) == 1,
                              as.character(unique(.data[[date_col]])),
                              paste0("Mixed (", paste(unique(.data[[date_col]]), collapse = "; "), ")")),
      replicate_count = n(),
      .groups = "drop"
    )
  
  # Final output
  final_summary <- summary_wide %>%
    left_join(replicate_info, by = id_col)
  
  # Write summary
  summary_csv <- paste0(output_prefix, "_summary.csv")
  write.csv(final_summary, summary_csv, row.names = FALSE)
  message("Summary saved to: ", summary_csv)
  
  return(final_summary)
}
#
#
#
#
#
#
#
#
#
#
graph_replicates_custom_error <- function(data, 
                                          id_col = "ID", 
                                          value_col, 
                                          se_col, 
                                          output_prefix = "replicate_analysis") {
  plot_dir <- paste0(output_prefix, "_plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Rename columns to standard names for plotting
  df_plot <- data %>%
    select(all_of(c(id_col, value_col, se_col))) %>%
    drop_na() %>%
    rename(ID = !!sym(id_col),
           value = !!sym(value_col),
           se = !!sym(se_col))
  
  # Create the plot
  p <- ggplot(df_plot, aes(x = ID, y = value)) +
    geom_col(fill = "steelblue", width = 0.6) +
    geom_errorbar(aes(ymin = value - se, ymax = value + se), 
                  width = 0.2, color = "black") +
    labs(title = paste("Replicate Summary:", value_col),
         x = id_col, y = value_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # Save static plot
  plot_filename <- file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_]", "_", value_col), "_replicates.png"))
  ggsave(plot_filename, p)
  
  message("Plot saved to: ", plot_filename)
  
  # Return interactive version
  return(plotly::ggplotly(p))
}
#
#
#
#
#
#
#
#
#
#
#
#
#
graph_histograms_with_error <- function(
    data,
    variables,
    id_col = "ID",
    output_prefix = "histogram_analysis"
) {
  plot_dir <- paste0(output_prefix, "_plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  plots <- purrr::map(variables, function(var) {
    if (!all(c(id_col, var) %in% colnames(data))) {
      warning("Missing required columns for variable: ", var)
      return(NULL)
    }
    
    # Summary statistics: mean and SE
    summary_df <- data %>%
      group_by(.data[[id_col]]) %>%
      summarise(
        mean_val = mean(.data[[var]], na.rm = TRUE),
        se_val = sd(.data[[var]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[var]]))),
        .groups = "drop"
      )
    
    # Create plot
    p <- ggplot(summary_df, aes_string(x = id_col, y = "mean_val")) +
      geom_col(fill = "steelblue", width = 0.6) +
      geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                    width = 0.1, color = "black") +
      geom_jitter(data = data, aes_string(x = id_col, y = var),
                  width = 0.15, height = 0, alpha = 0.5, color = "black", size = 1) +
      labs(title = paste("Mean ± SE:", var),
           x = id_col, y = var) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    
    # Save static plot
    ggsave(file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_]", "_", var), "_histogram.png")), p)
    return(ggplotly(p))
  })
  
  names(plots) <- variables
  message("Plots saved to folder: ", plot_dir)
  return(plots)
}

# Function to extract annotation
get_stats <- function(model, label) {
  mod_summary <- summary(model)
  
  # Safely extract R² and p-value
  r2 <- round(mod_summary$r.squared, 3)
  
  # Catch models with missing or weird coefficients
  if (nrow(mod_summary$coefficients) < 2) {
    p <- NA
  } else {
    p <- round(mod_summary$coefficients[2, 4], 3)  # p-value for x term
  }
  
  paste0(label, ": R²=", r2, ", p=", p)
}

compare_groups <- function(data, response_var, group_var, facet_var = NULL) {
  library(dplyr)
  library(ggpubr)
  library(rlang)
  
  # Ensure variables are symbols for tidy evaluation
  response_sym <- sym(response_var)
  group_sym <- sym(group_var)
  
  if (!is.null(facet_var)) {
    facet_sym <- sym(facet_var)
  }
  
  # Convert group to factor and response to numeric
  data <- data %>%
    mutate(
      !!group_sym := as.factor(!!group_sym),
      !!response_sym := as.numeric(!!response_sym)
    )
  
  # Number of groups
  num_groups <- nlevels(data[[group_var]])
  
  # Run parametric and non-parametric tests
  if (num_groups == 2) {
    param_test <- t.test(as.formula(paste(response_var, "~", group_var)), data = data)
    nonparam_test <- wilcox.test(as.formula(paste(response_var, "~", group_var)), data = data)
    
    message("Parametric test result (t-test):")
    print(param_test)
    
    message("Non-parametric test result (Wilcoxon test):")
    print(nonparam_test)
    
    compare_method <- "t.test"
    compare_method_np <- "wilcox.test"
    
  } else if (num_groups > 2) {
    param_test <- aov(as.formula(paste(response_var, "~", group_var)), data = data)
    nonparam_test <- kruskal.test(as.formula(paste(response_var, "~", group_var)), data = data)
    
    message("Parametric test result (ANOVA):")
    print(summary(param_test))
    
    message("Non-parametric test result (Kruskal-Wallis test):")
    print(nonparam_test)
    
    compare_method <- "anova"
    compare_method_np <- "kruskal.test"
    
  } else {
    stop("Group variable must have at least 2 levels.")
  }
  
  # Create boxplot with ggpubr
  p <- ggboxplot(data, x = group_var, y = response_var,
                 color = group_var, palette = "jco",
                 add = "jitter", add.params = list(size = 1.2, alpha = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    stat_compare_means(method = compare_method, label = "p.signif") +
    stat_compare_means(method = compare_method_np, label = "p.signif", label.y = max(data[[response_var]], na.rm = TRUE)*1.1)
  
  # Add faceting if requested
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }
  
  return(p)
}
#
#
#
#
#
#
#
# Function to calculate protein concentration from Lowry model
calculate_protein_from_lowry <- function(tidy_df, model, absorbance_column = "X600") {
  if (!absorbance_column %in% names(tidy_df)) {
    stop(paste("Column", absorbance_column, "not found in tidy_df"))
  }
  tidy_df %>%
    mutate(Protein_mg_per_mL = (.[[absorbance_column]] - coef(model)[1]) / coef(model)[2])
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
compare_groups_bootstrap <- function(data,
                                     response_var,
                                     group_var,
                                     facet_var     = NULL,
                                     n_boot         = 5000,
                                     conf_level     = 0.95,
                                     assumptions_test = FALSE) {
  
  # Convert string inputs to symbols
  response_sym <- sym(response_var)
  group_sym    <- sym(group_var)
  
  if (!is.null(facet_var)) {
    facet_sym <- sym(facet_var)
  }
  
  # Ensure group is factor and response is numeric; drop rows with NA in key columns
  df <- data %>%
    mutate(
      !!group_sym    := as.factor(!!group_sym),
      !!response_sym := as.numeric(!!response_sym)
    ) %>%
    filter(!is.na(!!response_sym), !is.na(!!group_sym))
  
  # If facet_var is provided, split dataset by each facet level; else treat all as one block
  if (!is.null(facet_var)) {
    df_split     <- df %>% group_split(!!facet_sym)
    facet_levels <- df %>% pull(!!facet_sym) %>% as.character() %>% unique()
  } else {
    df_split     <- list(df)
    facet_levels <- "ALL"
  }
  
  # Utility: bootstrap‐statistic function for a single group’s mean
  boot_mean_fn <- function(dat, indices) {
    d_sub <- dat[indices]
    return(mean(d_sub))
  }
  
  # Utility: compute pairwise effect‐size bootstraps (raw mean diff & Cohen’s d)
  compute_effsize_boot <- function(x, y, n_boot, conf_level) {
    pooled  <- c(x, y)
    grp_ind <- c(rep(1, length(x)), rep(2, length(y)))
    dat_for_boot <- data.frame(val = pooled, grp = grp_ind)
    
    # (1) Raw mean‐difference bootstrap
    stat_mean_diff <- function(d, i) {
      d2 <- d[i, ]
      m1 <- mean(d2$val[d2$grp == 1])
      m2 <- mean(d2$val[d2$grp == 2])
      return(m1 - m2)
    }
    boot_mean_diff <- boot(data = dat_for_boot, statistic = stat_mean_diff, R = n_boot)
    mean_diff_ci  <- boot.ci(boot_mean_diff, conf = conf_level, type = "perc")
    
    mean_diff_est   <- mean(x) - mean(y)
    mean_diff_lower <- mean_diff_ci$percent[4]
    mean_diff_upper <- mean_diff_ci$percent[5]
    
    # (2) Cohen’s d bootstrap (pooled‐SD version)
    cohen_d_stat <- function(d, i) {
      d2 <- d[i, ]
      x2 <- d2$val[d2$grp == 1]
      y2 <- d2$val[d2$grp == 2]
      coh <- cohens_d(x2, y2, pooled_sd = TRUE)
      return(coh$Cohens_d)
    }
    boot_cohend   <- boot(data = dat_for_boot, statistic = cohen_d_stat, R = n_boot)
    cohen_d_ci    <- boot.ci(boot_cohend, conf = conf_level, type = "perc")
    
    cohens_d_est   <- cohens_d(x, y, pooled_sd = TRUE)$Cohens_d
    cohens_d_lower <- cohen_d_ci$percent[4]
    cohens_d_upper <- cohen_d_ci$percent[5]
    
    return(tibble(
      mean_diff       = mean_diff_est,
      mean_diff_lower = mean_diff_lower,
      mean_diff_upper = mean_diff_upper,
      cohens_d        = cohens_d_est,
      cohens_d_lower  = cohens_d_lower,
      cohens_d_upper  = cohens_d_upper
    ))
  }
  
  # Prepare storage for summaries and plots
  all_summaries <- list()
  all_plots     <- list()
  
  # Iterate over each facet (or single “ALL” block)
  for (i in seq_along(df_split)) {
    sub_df     <- df_split[[i]]
    facet_name <- facet_levels[i]
    
    # Identify group levels within this facet
    grp_levels <- levels(sub_df[[group_var]])
    
    # --------------------------------------
    # 1) Bootstrap each group’s mean + keep the boot objects
    # --------------------------------------
    group_boot_list <- map(grp_levels, function(g) {
      vals     <- sub_df %>% filter((!!group_sym) == g) %>% pull(!!response_sym)
      boot_res <- boot(data = vals, statistic = boot_mean_fn, R = n_boot)
      
      # Extract CI
      ci <- boot.ci(boot_res, conf = conf_level, type = "perc")
      tibble(
        group     = g,
        mean_est  = mean(vals),
        ci_lower  = ci$percent[4],
        ci_upper  = ci$percent[5],
        n         = length(vals),
        boot_obj  = list(boot_res)   # store the boot object
      )
    }) %>%
      bind_rows()
    
    # Extract a clean summary of means ± CI (dropping the boot objects themselves)
    group_boot_res <- group_boot_list %>%
      select(group, mean_est, ci_lower, ci_upper, n)
    
    # --------------------------------------
    # 2) Pairwise effect‐size bootstraps
    # --------------------------------------
    pairwise <- combn(grp_levels, 2, simplify = FALSE) %>%
      set_names(map_chr(., ~ paste0(.x[1], "_vs_", .x[2]))) %>%
      map_dfr(function(pair) {
        x <- sub_df %>% filter((!!group_sym) == pair[1]) %>% pull(!!response_sym)
        y <- sub_df %>% filter((!!group_sym) == pair[2]) %>% pull(!!response_sym)
        
        compute_effsize_boot(x, y, n_boot, conf_level) %>%
          mutate(
            group1 = pair[1],
            group2 = pair[2]
          )
      })
    
    # --------------------------------------
    # 3) (Optional) Assumptions tests & bootstrap‐distribution plot
    # --------------------------------------
    assumption_plot <- NULL
    bartlett_test   <- NULL
    
    if (assumptions_test) {
      # 3a) Combine bootstrap replicates of group means into one data frame
      boot_dists_df <- group_boot_list %>%
        select(group, boot_obj) %>%
        unnest_longer(boot_obj) %>% 
        # boot_obj is a list; each entry is a 'boot' object whose $t is a matrix (n_boot × 1)
        mutate(boot_mean = map_dbl(boot_obj, ~ .x$t[, 1])) %>%
        select(group, boot_mean)
      
      # 3b) Plot the bootstrap distributions of group means (faceted by group)
      assumption_plot <- ggplot(boot_dists_df, aes(x = boot_mean, fill = group)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~ group, scales = "free") +
        theme_minimal() +
        labs(
          title = paste0("Bootstrap Distribution of Group Means (", response_var, ") [", facet_name, "]"),
          x     = paste("Bootstrapped means of", response_var),
          y     = "Density"
        ) +
        theme(legend.position = "none")
      
      # 3c) Run Bartlett’s test for homogeneity of variances (exploratory)
      bartlett_test <- tryCatch(
        {
          bartlett.test(
            formula = as.formula(paste(response_var, "~", group_var)),
            data    = sub_df
          )
        },
        error = function(e) {
          warning("Bartlett's test failed: ", e$message)
          return(NULL)
        }
      )
    }
    
    # --------------------------------------
    # 4) Assemble this facet’s summary
    # --------------------------------------
    summary_this <- list(
      facet                 = facet_name,
      group_bootstrap_means = group_boot_res,
      pairwise_effectsizes  = pairwise,
      bartlett_test         = bartlett_test
    )
    
    # --------------------------------------
    # 5) Create the main plots
    # --------------------------------------
    # 5a) Plot 1: Bootstrapped group means ± CI
    p1 <- ggplot(group_boot_res, aes(x = group, y = mean_est)) +
      geom_point(size = 3, color = "steelblue") +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.25, color = "steelblue") +
      theme_minimal() +
      labs(
        title = paste0("Bootstrapped Means ± ", conf_level*100, "% CI (", facet_name, ")"),
        x     = group_var,
        y     = paste0("Mean of ", response_var)
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 5b) Plot 2: Pairwise effect sizes (mean_diff & Cohen’s d) ± CI
    pairwise_plot_data <- pairwise %>%
      pivot_longer(cols = c(mean_diff, cohens_d),
                   names_to   = "metric",
                   values_to  = "estimate") %>%
      mutate(
        lower      = if_else(metric == "mean_diff", mean_diff_lower, cohens_d_lower),
        upper      = if_else(metric == "mean_diff", mean_diff_upper, cohens_d_upper),
        comparison = paste0(group1, " vs ", group2)
      )
    
    p2 <- ggplot(pairwise_plot_data, aes(x = comparison, y = estimate, color = metric)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    position = position_dodge(width = 0.5),
                    width    = 0.2) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste0("Pairwise Effect Sizes ± ", conf_level*100, "% CI (", facet_name, ")"),
        x     = "Comparison",
        y     = "Estimate"
      ) +
      scale_color_manual(
        values = c("mean_diff" = "darkgreen", "cohens_d" = "firebrick"),
        labels = c("Mean Difference", "Cohen's d"),
        name   = NULL
      )
    
    # 5c) If assumptions_test = TRUE, include the bootstrap‐distribution plot as well
    if (assumptions_test) {
      all_plots[[facet_name]] <- list(
        boot_means_plot       = p1,
        pairwise_effects_plot = p2,
        boot_dist_plot        = assumption_plot
      )
    } else {
      all_plots[[facet_name]] <- list(
        boot_means_plot       = p1,
        pairwise_effects_plot = p2,
        boot_dist_plot        = NULL
      )
    }
    
    # 6) Store this facet’s summary
    all_summaries[[facet_name]] <- summary_this
  }
  
  # Return a list containing:
  #  - $summary:    a named list of summaries (one per facet), including group means, pairwise effects, and Bartlett’s test object (if run)
  #  - $plots:      a named list of lists (one per facet) containing:
  #                  • boot_means_plot
  #                  • pairwise_effects_plot
  #                  • boot_dist_plot (NULL unless assumptions_test = TRUE)
  return(list(
    summary = all_summaries,
    plots   = all_plots
  ))
}


#
#
#
#
#
#
#
#
#
#
#
# ── (Optional) Function to also add standard‐error-of‐fit using the clean Run 2 model ────────

#' add_PE_with_se
#'
#' Given a dataframe containing a column of fluorescence values (Xred),
#' adds two new columns:
#'   1) Predicted PE (mg/g) based on the cleaned Run 2 model
#'   2) Standard error of the fitted mean (SE) for each prediction
#'
#' @param df        A data frame.
#' @param fluor_col String: name of the column in df containing the fluorescence (Xred) values.
#' @param pred_col  String: name of the new predicted PE column.
#' @param se_col    String: name of the new standard‐error column.
#' @return          A new data frame with additional columns `pred_col` and `se_col`.
add_PE_with_se <- function(df, fluor_col = "Xred",
                           pred_col = "PE_pred_run2_mg_per_g",
                           se_col   = "PE_se_run2") {
  # Build a newdata frame for prediction
  newdata <- df %>%
    transmute(!!sym(fluor_col) := .data[[fluor_col]])
  
  # Use predict.lm with se.fit = TRUE
  preds <- predict(model_clean_run2, newdata = newdata, se.fit = TRUE)
  
  # Attach predicted values and SE to the original df
  df[[pred_col]] <- preds$fit
  df[[se_col]]   <- preds$se.fit
  return(df)
}

# Example: apply to the master dataset
df_combined <- add_PE_with_se(df_combined,
                              fluor_col = "Xred",
                              pred_col  = "PE_pred_run2_mg_per_g",
                              se_col    = "PE_se_run2")

# View the first few rows to confirm
df_combined %>%
  select(join_id, run, Xred, PE_mg_per_g_sample,
         PE_pred_run2_mg_per_g, PE_se_run2) %>%
  head(10) %>%
  print()








