# scripts/process_plate_run.R
required_packages <- c(
  "readxl", "dplyr", "tidyr", "tidyverse", "purrr", "zoo",
  "ggplot2", "lubridate", "stringr", "tools", "plotly", 
  "here", "broom", "easystats", "performance", "ggrepel", 
  "effectsize", "rlang", "boot", "emmeans", "multcompView", "fs", "knitr")
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
##' Tidy and Correct Microspectrometer Plate Data
##'
##' This function takes a raw data.frame organized like a microwell plate,
##' pivots it to long format by Cell_ID, applies blank subtraction per wavelength,
##' and returns a tidy data.frame of corrected spectral values.
##'
##' @param input_df A data.frame containing raw plate readings. One row per well,
##'   with the first column (or second) marking Row labels (A-H), multiple columns of readings,
##'   and a final column containing wavelength identifiers.
##' @param blanks A character vector of Cell_IDs marking which wells are blanks.
##' @param df_name Optional name. If provided, will assign the result to
##'   the global environment as <df_name>_tidy.
##' @return A data.frame in wide format: columns = Cell_ID + one column per wavelength,
##'   with blank-corrected values. Wells with any negative values are removed and
##'   their IDs stored in an attribute 'removed_rows_spectral'.
##' @examples
##' corrected <- tidy_and_correct(raw_plate_df, blanks = c("A01", "H12"))

# Primary function
tidy_and_correct <- function(input_df, df_name = NULL) {
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
  
  # Remove rows with negative values (if any)
  spectral_cols <- names(tidy_df)[-1]
  negative_flags <- apply(tidy_df[, spectral_cols], 1, function(row) any(row < 0, na.rm = TRUE))
  
  removed_rows <- tidy_df[negative_flags, ]
  if (nrow(removed_rows) > 0) {
    message("Removed rows due to negative values:")
    print(removed_rows$Cell_ID)
  }
  
  tidy_filtered <- tidy_df[!negative_flags, ]
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
#
#
#
tidy_all <- function(df_list, blanks_list) {
  if (length(df_list) != length(blanks_list)) {
    stop("Length of df_list and blanks_list must match.")
  }
  if (is.null(names(df_list))) {
    stop("df_list must be a named list.")
  }
  
  for (i in seq_along(df_list)) {
    name   <- names(df_list)[i]
    raw_df <- df_list[[i]]
    blanks <- blanks_list[[i]]
    
    message("Tidying ", name)
    tidy_df <- tidy_and_correct(raw_df)
    assign(paste0(name, "_start"), tidy_df, envir = .GlobalEnv)
    message("Correcting ", name, " using blank(s): ", paste(blanks, collapse = ", "))
    
    # Check Cell_ID column exists
    if (!"Cell_ID" %in% colnames(tidy_df)) {
      stop(paste("Tidy dataframe", name, "must contain a 'Cell_ID' column."))
    }
    
    # Calculate column-wise mean of blanks
    blank_means <- tidy_df %>%
      dplyr::filter(Cell_ID %in% blanks) %>%
      dplyr::select(where(is.numeric)) %>%
      dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
      unlist()
    
    # Subtract from all numeric columns
    corrected_df <- tidy_df %>%
      dplyr::mutate(dplyr::across(
        where(is.numeric),
        ~ .x - blank_means[cur_column()]
      ))
    
    assign(paste0(name, "_tidy"), corrected_df, envir = .GlobalEnv)
    message("Saved ", name, "_final to global environment.")
  }
  
  invisible(NULL)
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
correct_all <- function(df_list, blanks_list) {
  if (length(df_list) != length(blanks_list)) {
    stop("Length of df_list and blanks_list must match.")
  }
  if (is.null(names(df_list))) {
    stop("df_list must be a named list.")
  }
  
  for (i in seq_along(df_list)) {
    name    <- names(df_list)[i]
    df      <- df_list[[i]]
    blanks  <- blanks_list[[i]]
    
    message("Correcting ", name, " using blank(s): ", paste(blanks, collapse = ", "))
    
    # Check that Cell_ID column exists
    if (!"Cell_ID" %in% colnames(df)) {
      stop(paste("Dataframe", name, "must contain a 'Cell_ID' column."))
    }
    
    # Extract blank rows and calculate mean for each numeric column
    blank_means <- df %>%
      dplyr::filter(Cell_ID %in% blanks) %>%
      dplyr::select(where(is.numeric)) %>%
      dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
    
    # Subtract blank means from numeric columns
    df_corrected <- df %>%
      dplyr::mutate(across(where(is.numeric), ~ .x - unlist(blank_means)))
    
    # Save to global environment
    assign(paste0(name, "_corrected"), df_corrected, envir = .GlobalEnv)
    message("Saved ", name, "_corrected to global environment.")
  }
  
  invisible(NULL)
}


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
# Enhanced analyze_replicates: choose best three, and if exactly three replicates, choose best two
analyze_replicates <- function(data,
                               id_col = "ID",
                               join_col = "join_id",
                               weight_col = "sample weight",
                               date_col = "date",
                               output_prefix = "replicate_analysis",
                               choose_best_3 = FALSE,
                               dir = "output PE/export data") {
  
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
  
  # Summary stats with best selection if enabled
  summary_long <- long_data %>%
    group_by(.data[[id_col]], measure) %>%
    group_modify(~ {
      df_sub <- .x
      n <- nrow(df_sub)
      
      if (choose_best_3 && n > 3) {
        # More than three: select best three
        selected <- select_best_three(df_sub)
        keep_rows <- selected$included_rows
      } else if (choose_best_3 && n == 3) {
        # Exactly three: select best two
        selected <- select_best_two(df_sub)
        keep_rows <- selected$included_rows
      } else {
        # Default: keep all
        keep_rows <- seq_len(n)
      }
      
      # Ensure keep_rows is a simple integer vector (not matrix)
      keep_rows <- as.vector(keep_rows)
      
      # Subset for stats
      df_best <- df_sub[keep_rows, , drop = FALSE]
      m <- nrow(df_best)
      mean_val <- mean(df_best$value, na.rm = TRUE)
      sd_val <- sd(df_best$value, na.rm = TRUE)
      se_val <- sd_val / sqrt(m)
      cv_val <- sd_val / mean_val
      max_dev_pct_val <- max(abs(df_best$value - mean_val) / mean_val * 100, na.rm = TRUE)
      
      # Build result tibble
      out <- tibble(
        n = m,
        mean = mean_val,
        sd = sd_val,
        se = se_val,
        cv = cv_val,
        max_dev_pct = max_dev_pct_val
      )
      
      # If using best selection, record included/excluded rows
      if (choose_best_3 && n >= 3) {
        included_str <- paste(keep_rows, collapse = ",")
        excluded <- setdiff(seq_len(n), keep_rows)
        excluded_str <- if (length(excluded) == 0) NA_character_ else paste(excluded, collapse = ",")
        out <- out %>%
          mutate(included_rows = included_str,
                 excluded_rows = excluded_str)
      }
      
      out
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
  summary_csv <- paste0(dir, output_prefix, "_summary.csv")
  write.csv(final_summary, summary_csv, row.names = FALSE)
  message("Summary saved to: ", summary_csv)
  
  return(final_summary)
}

# Helper: select_best_two (chooses the pair with smallest difference)
select_best_two <- function(df) {
  vals <- df$value
  pairs <- combn(seq_along(vals), 2)
  diffs <- apply(pairs, 2, function(idx) abs(vals[idx[1]] - vals[idx[2]]))
  best_pair <- pairs[, which.min(diffs)]
  list(included_rows = best_pair)
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
compare_groups <- function(data,
                           response_var,
                           group_var,
                           facet_var  = NULL,
                           subfolder_name = NULL) {
  # Load required packages
  if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
  library(ggpubr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(rlang)
  
  # Use global params from Rmd for base dirs
  base_plot_dir   <- plot_dir
  base_report_dir <- report_dir
  base_data_dir   <- data_dir
  
  # If subfolder name is provided, append it to each base dir
  if (!is.null(subfolder_name)) {
    base_plot_dir   <- file.path(base_plot_dir, subfolder_name)
    base_report_dir <- file.path(base_report_dir, subfolder_name)
    base_data_dir   <- file.path(base_data_dir, subfolder_name)
    dir.create(base_plot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(base_report_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(base_data_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Clean name helper
  clean_name <- function(x) {
    x %>% str_replace_all("[^a-zA-Z]+", "_") %>% str_remove_all("_$") %>% str_remove_all("^_")
  }
  
  # Compact letter display
  generate_group_letters <- function(posthoc_df) {
    comps <- str_split_fixed(posthoc_df$comparison, "-", 2)
    pvals <- setNames(posthoc_df$`p adj`, paste(comps[,1], comps[,2], sep = "-"))
    cld <- multcompLetters(pvals)
    tibble(group = names(cld$Letters), letter = cld$Letters)
  }
  
  response_sym <- sym(response_var)
  group_sym    <- sym(group_var)
  if (!is.null(facet_var)) facet_sym <- sym(facet_var)
  
  df <- data %>%
    mutate(
      !!group_sym    := factor(!!group_sym, levels = unique(data[[group_var]])),
      !!response_sym := as.numeric(!!response_sym)
    ) %>%
    filter(!is.na(!!response_sym), !is.na(!!group_sym))
  
  if (!is.null(facet_var)) {
    df_split     <- df %>% group_split(!!facet_sym)
    facet_levels <- df %>% pull(!!facet_sym) %>% as.character() %>% unique()
  } else {
    df_split     <- list(df)
    facet_levels <- "ALL"
  }
  
  clean_response <- clean_name(response_var)
  clean_group    <- clean_name(group_var)
  
  all_results <- list()
  
  for (i in seq_along(df_split)) {
    sub_df     <- df_split[[i]]
    facet_name <- facet_levels[i]
    
    summary_tbl <- sub_df %>%
      group_by(!!group_sym) %>%
      summarise(
        n = n(),
        mean = mean(!!response_sym, na.rm = TRUE),
        sd   = sd(!!response_sym, na.rm = TRUE),
        median = median(!!response_sym, na.rm = TRUE),
        IQR    = IQR(!!response_sym, na.rm = TRUE),
        .groups = "drop"
      )
    save_object(summary_tbl, paste0("summary_", facet_name), directory = "data", subdir = subfolder_name, format = "csv")
    
    ng <- nlevels(sub_df[[group_var]])
    test_results <- NULL
    posthoc_results <- NULL
    group_letters <- NULL
    
    if (ng == 2) {
      fmla <- as.formula(paste(response_var, "~", group_var))
      t_res <- t.test(fmla, data = sub_df)
      w_res <- wilcox.test(fmla, data = sub_df)
      d     <- cohens_d(fmla, data = sub_df)
      
      test_results <- tibble(
        test        = c("t.test", "wilcox.test"),
        statistic   = c(t_res$statistic, w_res$statistic),
        df          = c(t_res$parameter, NA),
        p_value     = c(t_res$p.value, w_res$p.value),
        effect_size = c(d$Cohens_d, NA),
        CI_lower    = c(t_res$conf.int[1], NA),
        CI_upper    = c(t_res$conf.int[2], NA)
      )
    } else {
      fmla <- as.formula(paste(response_var, "~", group_var))
      aov_res <- aov(fmla, data = sub_df)
      kw_res  <- kruskal.test(fmla, data = sub_df)
      
      eta2  <- eta_squared(aov_res)
      r_eta <- rank_eta_squared(fmla, data = sub_df)
      
      test_results <- tibble(
        test        = c("ANOVA", "Kruskal-Wallis"),
        statistic   = c(summary(aov_res)[[1]]$`F value`[1], kw_res$statistic),
        df          = c(paste(summary(aov_res)[[1]]$Df[1], summary(aov_res)[[1]]$Df[2], sep = ", "), kw_res$parameter),
        p_value     = c(summary(aov_res)[[1]]$`Pr(>F)`[1], kw_res$p.value),
        effect_size = c(as.numeric(eta2)[1], r_eta$Rank_Eta2[1]),
        CI_lower    = NA,
        CI_upper    = NA
      )
      
      if (summary(aov_res)[[1]]$`Pr(>F)`[1] < 0.05) {
        tukey <- TukeyHSD(aov_res)
        posthoc_results <- as.data.frame(tukey[[1]])
        posthoc_results$comparison <- rownames(posthoc_results)
        rownames(posthoc_results) <- NULL
        group_letters <- generate_group_letters(posthoc_results)
      } else if (kw_res$p.value < 0.05) {
        ph <- pairwise.wilcox.test(sub_df[[response_var]], sub_df[[group_var]], p.adjust.method = "BH")
        posthoc_results <- as.data.frame(as.table(ph$p.value)) %>%
          filter(!is.na(Freq)) %>%
          rename(comparison_1 = Var1, comparison_2 = Var2, p_value = Freq)
      }
    }
    
    save_object(test_results, paste0("tests_", facet_name), directory = "data", subdir = subfolder_name, format = "csv")
    if (!is.null(posthoc_results)) save_object(posthoc_results, paste0("posthoc_", facet_name), directory = "data", subdir = subfolder_name, format = "csv")
    if (!is.null(group_letters))    save_object(group_letters, paste0("group_letters_", facet_name), directory = "data", subdir = subfolder_name, format = "csv")
    
    if (!is.null(group_letters)) {
      label_positions <- sub_df %>%
        mutate(group = as.character(!!group_sym)) %>%
        group_by(group) %>%
        summarise(y_pos = max(!!response_sym, na.rm = TRUE) * 1.05, .groups = "drop")
      
      group_letters <- group_letters %>% mutate(group = as.character(group))
      label_df <- left_join(group_letters, label_positions, by = "group")
    }
    
    sub_df[[group_var]] <- factor(sub_df[[group_var]], levels = levels(df[[group_var]]))
    
    p <- ggboxplot(
      sub_df,
      x        = group_var,
      y        = response_var,
      color    = group_var,
      palette  = "jco",
      add      = "jitter",
      add.params = list(size = 2, alpha = 0.7, width = 0.3),
      outlier.shape = NA
    ) +
      coord_cartesian(ylim = c(0, max(df[[response_var]], na.rm = TRUE) * 1.2)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)))
    }
    
    if (exists("label_df") && nrow(label_df) > 0 && all(!is.na(label_df$y_pos))) {
      p <- p + geom_text(
        data = label_df,
        aes(x = group, y = y_pos, label = letter),
        vjust = -0.5,
        fontface = "bold",
        size = 2.5,
        inherit.aes = FALSE
      )
    }
    
    save_object(p,
                filename  = paste0("plot_", clean_response, "_by_", clean_group, "_", facet_name),
                directory = "plots",
                subdir    = subfolder_name,
                width     = 8,
                height    = 6,
                dpi       = 300)
    
    print(p)
    
    all_results[[facet_name]] <- list(
      summary   = summary_tbl,
      tests     = test_results,
      posthoc   = posthoc_results,
      cld       = group_letters
    )
  }
  
  invisible(all_results)
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
#
compare_groups_bootstrap <- function(data,
                                     response_var,
                                     group_var,
                                     facet_var       = NULL,
                                     n_boot           = 5000,
                                     conf_level       = 0.95,
                                     assumptions_test = FALSE,
                                     report_dir       = "reports",
                                     plot_dir         = "plots") {
  library(dplyr); library(purrr); library(boot); library(ggplot2)
  library(tibble); library(tidyr); library(effsize); library(stringr)
  
  clean_name <- function(x) {
    x %>%
      str_replace_all("[^a-zA-Z]+", "_") %>%
      str_remove_all("_$") %>%
      str_remove_all("^_")
  }
  
  response_sym <- sym(response_var)
  group_sym    <- sym(group_var)
  if (!is.null(facet_var)) facet_sym <- sym(facet_var)
  
  df <- data %>%
    mutate(
      !!group_sym    := factor(!!group_sym, levels = unique(data[[group_var]])),
      !!response_sym := as.numeric(!!response_sym)
    ) %>%
    filter(!is.na(!!response_sym), !is.na(!!group_sym))
  
  if (!is.null(facet_var)) {
    df_split     <- df %>% group_split(!!facet_sym)
    facet_levels <- df %>% pull(!!facet_sym) %>% as.character() %>% unique()
  } else {
    df_split     <- list(df)
    facet_levels <- "ALL"
  }
  
  clean_response <- clean_name(response_var)
  clean_group    <- clean_name(group_var)
  sub_report_dir <- file.path(report_dir, paste0(clean_response, "_by_", clean_group))
  sub_plot_dir   <- file.path(plot_dir, paste0(clean_response, "_by_", clean_group))
  dir.create(sub_report_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(sub_plot_dir,   showWarnings = FALSE, recursive = TRUE)
  
  boot_mean_fn <- function(dat, i) mean(dat[i])
  
  compute_effsize_boot <- function(x, y, n_boot, conf_level) {
    pooled <- data.frame(val = c(x, y), grp = rep(c(1,2), c(length(x), length(y))))
    
    md_stat <- function(d, ix) {
      d2 <- d[ix,]
      mean(d2$val[d2$grp==1]) - mean(d2$val[d2$grp==2])
    }
    cd_stat <- function(d, ix) {
      d2 <- d[ix,]
      effsize::cohen.d(d2$val[d2$grp==1], d2$val[d2$grp==2], pooled = TRUE)$estimate
    }
    
    b_md <- boot(pooled, md_stat, R = n_boot)
    ci_md <- tryCatch({
      boot.ci(b_md, conf = conf_level, type = "perc")$percent[4:5]
    }, error = function(e) c(NA, NA))
    
    b_cd <- boot(pooled, cd_stat, R = n_boot)
    ci_cd <- tryCatch({
      boot.ci(b_cd, conf = conf_level, type = "perc")$percent[4:5]
    }, error = function(e) c(NA, NA))
    
    cohens_d_point <- effsize::cohen.d(x, y, pooled = TRUE)$estimate
    
    tibble(
      mean_diff       = mean(x) - mean(y),
      mean_diff_lower = ci_md[1],
      mean_diff_upper = ci_md[2],
      cohens_d        = cohens_d_point,
      cohens_d_lower  = ci_cd[1],
      cohens_d_upper  = ci_cd[2]
    )
  }
  
  assign_overlap_letters <- function(df) {
    letters <- LETTERS
    current_letter <- 1
    df$letter <- NA
    for (i in seq_len(nrow(df))) {
      if (is.na(df$letter[i])) {
        df$letter[i] <- letters[current_letter]
        for (j in seq((i + 1), nrow(df))) {
          if (!anyNA(c(df$ci_lower[j], df$ci_upper[j], df$ci_lower[i], df$ci_upper[i])) &&
              df$ci_lower[j] <= df$ci_upper[i] && df$ci_upper[j] >= df$ci_lower[i]) {
            df$letter[j] <- letters[current_letter]
          }
        }
        current_letter <- current_letter + 1
      }
    }
    return(df)
  }
  
  
  all_summaries <- list()
  
  for (i in seq_along(df_split)) {
    sub_df     <- df_split[[i]]
    facet_name <- facet_levels[i]
    grp_levels <- levels(sub_df[[group_var]])
    
    if (length(grp_levels) < 2) {
      warning("Facet '", facet_name, "' skipped: not enough groups for comparison.")
      next
    }
    
    group_boot_list <- map(grp_levels, function(g) {
      vals <- sub_df %>% filter((!!group_sym)==g) %>% pull(!!response_sym)
      b    <- boot(vals, boot_mean_fn, R = n_boot)
      ci   <- boot.ci(b, conf = conf_level, type = "perc")$percent[4:5]
      tibble(
        group    = g,
        mean_est = mean(vals),
        ci_lower = ci[1],
        ci_upper = ci[2],
        n        = length(vals)
      )
    }) %>% bind_rows()
    
    pairwise <- combn(grp_levels, 2, simplify = FALSE) %>%
      map_dfr(function(p) {
        compute_effsize_boot(
          x = sub_df %>% filter((!!group_sym)==p[1]) %>% pull(!!response_sym),
          y = sub_df %>% filter((!!group_sym)==p[2]) %>% pull(!!response_sym),
          n_boot, conf_level
        ) %>%
          mutate(group1 = p[1], group2 = p[2])
      })
    
    bartlett_test <- NULL
    if (assumptions_test) {
      bartlett_test <- tryCatch(
        bartlett.test(as.formula(paste(response_var, "~", group_var)), data = sub_df),
        error = function(e) NULL
      )
    }
    
    group_boot_list <- group_boot_list %>%
      arrange(desc(mean_est)) %>%
      mutate(group = as.character(group))
    
    group_boot_list <- assign_overlap_letters(group_boot_list)
    
    write.csv(group_boot_list, file = file.path(sub_report_dir, paste0("means_", facet_name, ".csv")), row.names = FALSE)
    write.csv(pairwise,        file = file.path(sub_report_dir, paste0("effects_", facet_name, ".csv")), row.names = FALSE)
    write.csv(group_boot_list %>% select(group, letter),
              file = file.path(sub_report_dir, paste0("group_letters_", facet_name, ".csv")), row.names = FALSE)
    
    all_summaries[[facet_name]] <- list(
      group_bootstrap_means = group_boot_list,
      pairwise_effectsizes  = pairwise,
      bartlett_test         = bartlett_test,
      group_letters         = group_boot_list %>% select(group, letter)
    )
    
    p1 <- ggplot(group_boot_list, aes(x = group, y = mean_est)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
      geom_text(aes(label = letter, y = ci_upper * 1.05), vjust = 0, fontface = "bold", size = 3) +
      theme_minimal() +
      labs(
        title = paste0("Bootstrapped Means ± ", conf_level * 100, "% CI [", facet_name, "]"),
        x     = group_var, y = response_var
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p1)
    ggsave(file.path(sub_plot_dir, paste0("means_", facet_name, ".png")), p1, width = 8, height = 6)
    
    df2 <- pairwise %>%
      pivot_longer(c(mean_diff, cohens_d), names_to = "metric", values_to = "est") %>%
      mutate(
        lower = ifelse(metric=="mean_diff", mean_diff_lower, cohens_d_lower),
        upper = ifelse(metric=="mean_diff", mean_diff_upper, cohens_d_upper),
        comp  = paste(group1, "vs", group2)
      )
    p2 <- ggplot(df2, aes(comp, est, color = metric)) +
      geom_point(position = position_dodge(0.5), size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    position = position_dodge(0.5), width = 0.2) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste0("Pairwise Effects ± ", conf_level * 100, "% CI [", facet_name, "]"),
        x = "", y = "Estimate"
      )
    print(p2)
    ggsave(file.path(sub_plot_dir, paste0("effects_", facet_name, ".png")), p2, width = 9, height = 6)
    
    if (assumptions_test) {
      boot_dists_df <- map2_df(group_boot_list$group, group_boot_list$n, function(grp, n) {
        tibble(
          group     = grp,
          boot_mean = boot(data = sub_df %>% filter((!!group_sym)==grp) %>% pull(!!response_sym),
                           statistic = boot_mean_fn, R = n_boot)$t[,1]
        )
      })
      p3 <- ggplot(boot_dists_df, aes(boot_mean)) +
        geom_density(fill = "steelblue", alpha = 0.4) +
        facet_wrap(~group, scales = "free") +
        theme_minimal() +
        labs(
          title = paste0("Bootstrap Density [", facet_name, "]"),
          x     = paste("Bootstrapped", response_var)
        )
      print(p3)
      ggsave(file.path(sub_plot_dir, paste0("density_", facet_name, ".png")), p3, width = 9, height = 6)
    }
  }
  
  message("Reports saved to: ", normalizePath(sub_report_dir))
  message("Plots saved to: ", normalizePath(sub_plot_dir))
  return(invisible(all_summaries))
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
calc_PE_from_Xred <- function(df, fluor_col = "Xred", new_col = "PE_predicted_mg_per_g") {
  intercept <- coef_intercept
  slope     <- coef_slope
  
  # Ensure the fluorescence column exists
  if (!fluor_col %in% names(df)) {
    stop(glue::glue("Column '{fluor_col}' not found in the input data frame."))
  }
  
  # Compute predicted PE and add as a new column
  df[[new_col]] <- intercept + slope * df[[fluor_col]]
  return(df)
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
generate_group_letters <- function(posthoc_df) {
  library(multcompView)
  library(stringr)
  
  # Split the "comparison" column into two groups
  comps <- str_split_fixed(posthoc_df$comparison, "-", 2)
  
  # Create a named vector of p-values
  pvals <- setNames(posthoc_df$`p adj`, paste(comps[,1], comps[,2], sep = "-"))
  
  # Run multcompView
  cld <- multcompLetters(pvals)
  
  # Return a data frame: group and its assigned letter
  tibble::tibble(
    group = names(cld$Letters),
    letter = cld$Letters
  )
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
save_object <- function(object,
                        filename,
                        directory = c("plots", "data", "reports"),
                        format    = NULL,
                        subdir    = NULL,
                        ...) {
  
  directory <- match.arg(directory)
  base_dir  <- switch(directory,
                      plots   = plot_dir,
                      data    = data_dir,
                      reports = report_dir)
  
  # Build full target directory (with optional subdir)
  target_dir <- if (!is.null(subdir)) file.path(base_dir, subdir) else base_dir
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  # Determine file extension
  if (is.null(format)) {
    format <- if (inherits(object, "ggplot")) "png"
    else if (is.data.frame(object)) "csv"
    else "rds"
  }
  
  out_path <- file.path(target_dir, paste0(filename, ".", format))
  
  # Save the object
  if (inherits(object, "ggplot") && format %in% c("png", "pdf", "svg", "jpeg")) {
    ggsave(filename = out_path, plot = object, ...)
  } else if (is.data.frame(object) && format == "csv") {
    write.csv(object, out_path, row.names = FALSE, ...)
  } else if (format == "rds") {
    saveRDS(object, out_path, ...)
  } else {
    stop("Unsupported object type or format.")
  }
  
  invisible(out_path)
}





