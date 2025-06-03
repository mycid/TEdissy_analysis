# scripts/process_plate_run.R
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)
library(ggplot2)
library(lubridate)
library(stringr)
library(tools)
library(plotly)
library(tidyverse)

required_packages <- c(
  "readxl", "dplyr", "tidyr", "tidyverse", "purrr", "zoo",
  "ggplot2", "lubridate", "stringr", "tools", "plotly"
)

for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(zoo)
library(readxl)
library(performance)
library(parameters)
library(broom)
library(easystats)

tidy_and_correct <- function(input_path, blanks, sheet = 2) {
  
  # Load data
  df_raw <- readxl::read_excel(input_path, sheet = sheet, col_names = FALSE)
  
  # Remove decorative header row
  df <- df_raw[-1, ]
  
  # Fill down row labels (first column: A, B, C...)
  df[[1]] <- zoo::na.locf(df[[1]])
  
  # Extract all wavelength values from the last column
  wl_column <- df[[ncol(df)]]
  unique_wavelengths <- unique(na.omit(wl_column))
  wl_vector <- as.character(unique_wavelengths)
  wl_vector <- stringr::str_trim(wl_vector)
  wl_vector <- paste0("X", wl_vector)  # Add X prefix
  num_wl <- length(wl_vector)
  
  # Drop the wavelength column, keep only the plate data
  df <- df[, 1:13]  # First column is Row, next 12 are plate values
  
  # Assign numeric column names to plate columns
  colnames(df)[1:13] <- c("Row", as.character(1:12))
  
  # Pivot to long format
  tidy_df <- df %>%
    tidyr::pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
    dplyr::mutate(
      Column = as.integer(Column),
      Cell_ID = sprintf("%s%02d", Row, Column),
      Value = as.numeric(Value)
    ) %>%
    dplyr::select(Cell_ID, Value)
  
  # Handle multiple wavelengths
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
  
  # Blank correction
  blanks_df <- tidy_df %>% dplyr::filter(Cell_ID %in% blanks)
  blank_means <- blanks_df %>% dplyr::select(-Cell_ID) %>%
    dplyr::summarise(dplyr::across(everything(), ~mean(.x, na.rm = TRUE)))
  
  tidy_corrected <- tidy_df %>%
    dplyr::mutate(dplyr::across(-Cell_ID, ~ .x - unlist(blank_means)))
  
  # Remove rows with negative values
  spectral_cols <- names(tidy_corrected)[-1]
  negative_flags <- apply(tidy_corrected[, spectral_cols], 1, function(row) any(row < 0, na.rm = TRUE))
  
  removed_rows <- tidy_corrected[negative_flags, ]
  if (nrow(removed_rows) > 0) {
    message("Removed rows due to negative values:")
    print(removed_rows$Cell_ID)
  }
  
  tidy_filtered <- tidy_corrected[!negative_flags, ]
  
  # Attach removed rows as attribute
  attr(tidy_filtered, "removed_rows_spectral") <- removed_rows
  
  return(tidy_filtered)
}

        



calculate_pe_and_filter <- function(tidy_df, 
                                    extract_volume_mL = 1.5, 
                                    sample_weight_col = "sample weights",
                                    sample_id_col = "Cell_ID") {
  library(dplyr)
  
  # Calculate PE (Beer & Eshel 1985)
  tidy_df <- tidy_df %>%
    mutate(
      PE_ug_per_g = ((`X564` - `X592`) - ((`X455` - `X592`) * 0.20)) * 0.12
    )
  
  # Identify rows with negative PE
  to_remove_pe <- which(tidy_df$PE_ug_per_g < 0)
  
  removed_rows <- data.frame()
  if(length(to_remove_pe) > 0) {
    removed_rows <- tidy_df[to_remove_pe, ]
    removed_rows$Removal_Reason <- "removed because PE < 0"
    message("Removed observations due to negative PE values:")
    print(removed_rows %>% select(all_of(sample_id_col), PE_ug_per_g, Removal_Reason))
  } else {
    message("No observations removed due to negative PE values.")
  }
  
  tidy_final <- tidy_df[-to_remove_pe, ]
  
  # Convert PE from µg/mL to mg/g dry weight using sample weight column
  tidy_final <- tidy_final %>%
    mutate(
      PE_mg_per_g_sample = (PE_ug_per_g * (extract_volume_mL / .data[[sample_weight_col]])) / 1000
    )
  
  # Attach removed rows as attribute
  attr(tidy_final, "removed_rows_pe") <- removed_rows
  
  return(tidy_final)
}

  
joindf_by_id <- function(df1, df2, output_name, key_df1, key_df2) {
    
    # Trim column names
    colnames(df1) <- trimws(colnames(df1))
    colnames(df2) <- trimws(colnames(df2))
    
    # Check that the specified keys exist in their respective dataframes
    if (!(key_df1 %in% colnames(df1))) {
      stop(paste0("Error: '", key_df1, "' column not found in df1."))
    }
    if (!(key_df2 %in% colnames(df2))) {
      stop(paste0("Error: '", key_df2, "' column not found in df2."))
    }
    
    # Rename join keys to a common name "join_id" for joining
    df1_renamed <- df1 %>% rename(join_id = all_of(key_df1))
    df2_renamed <- df2 %>% rename(join_id = all_of(key_df2))
    
    # Join logic: left join from smaller to larger dataframe
    if (nrow(df1_renamed) <= nrow(df2_renamed)) {
      result <- df1_renamed %>% left_join(df2_renamed, by = "join_id")
      base_ids <- df1_renamed$join_id
      compare_ids <- df2_renamed$join_id
      unmatched_df <- df2_renamed %>% filter(!(join_id %in% base_ids))
    } else {
      result <- df2_renamed %>% left_join(df1_renamed, by = "join_id")
      base_ids <- df2_renamed$join_id
      compare_ids <- df1_renamed$join_id
      unmatched_df <- df1_renamed %>% filter(!(join_id %in% base_ids))
    }
    
    # Save joined result
    write.csv(result, output_name, row.names = FALSE)
    
    # Assign to global environment
    object_name <- file_path_sans_ext(basename(output_name))
    assign(object_name, result, envir = .GlobalEnv)
    
    # Save unmatched rows to CSV
    unmatched_csv <- paste0("unmatched_rows_", object_name, ".csv")
    write.csv(unmatched_df, unmatched_csv, row.names = FALSE)
    
    # Print "needs to be repeated" messages
    repeat_list <- unique(unmatched_df[[1]])  # first column of unmatched
    cat("\nValues needing repetition:\n")
    cat(paste0(repeat_list, " needs to be repeated\n"), sep = "")
    
    # Report
    matched <- sum(compare_ids %in% base_ids)
    unmatched_count <- nrow(unmatched_df)
    
    report <- list(
      merged_cells = matched,
      unmatched_cells = unmatched_count,
      unmatched_saved_to = unmatched_csv
    )
    
    message("Join complete. Output saved to: ", output_name)
    return(report)
}


analyze_replicates <- function(data,
                               id_col = "ID",
                               join_col = "join_id",
                               weight_col = "sample weight",
                               date_col = "date",
                               output_prefix = "replicate_analysis") {
  
  # Ensure proper types
  data[[id_col]] <- as.factor(data[[id_col]])
  data[[date_col]] <- as.Date(data[[date_col]])
  
  # Identify numeric columns to analyze (exclude metadata)
  numeric_cols <- data %>%
    select(where(is.numeric)) %>%
    select(-any_of(c(weight_col))) %>%  # exclude sample weight from numeric cols
    colnames()
  
  # Pivot to long format for analysis
  long_data <- data %>%
    pivot_longer(cols = all_of(numeric_cols), names_to = "measure", values_to = "value")
  
  # Summary stats
  summary_long <- long_data %>%
    group_by(.data[[id_col]], measure) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n),
      cv = sd / mean,
      max_dev_pct = max(abs(value - mean) / mean * 100, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Pivot stats to wide format for readability
  summary_wide <- summary_long %>%
    pivot_wider(
      names_from = measure,
      values_from = c(mean, sd, se, cv, max_dev_pct),
      names_glue = "{measure}_{.value}"
    )
  
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
  final_summary <<- summary_wide %>%
    left_join(replicate_info, by = id_col)
  
  # Write summary
  summary_csv <- paste0(output_prefix, "_summary.csv")
  write.csv(final_summary, summary_csv, row.names = FALSE)
  message("Summary saved to: ", summary_csv)
  
  # Plotting
  plot_dir <- paste0(output_prefix, "_plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  plots <- long_data %>%
    group_by(measure) %>%
    group_split() %>%
    map(~ {
      p <- ggplot(.x, aes_string(x = id_col, y = "value")) +
        stat_summary(fun = mean, geom = "bar", fill = "steelblue", width = 0.6) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
        labs(title = paste("Replicate Fit:", unique(.x$measure)),
             x = "Sample ID", y = "Value") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # ensures better y-axis spacing
      ggsave(file.path(plot_dir, paste0(unique(.x$measure), "_replicates.png")), p)
      return(ggplotly(p))
    })
  
  
  message("Plots saved to folder: ", plot_dir)
  
  return(list(summary = final_summary, plots = plots))
}

get_stats <- function(model, label) {
  stats <- glance(model)
  r2 <- round(stats$r.squared, 3)
  p <- round(tidy(model)$p.value[2], 3)  # second row = x term
  paste0(label, ": R²=", r2, ", p=", p)
}




