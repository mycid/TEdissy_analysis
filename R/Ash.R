---
  title: "Ash Content Analysis"
author: "Trevor Eakes"
date: "`r Sys.Date()`"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

## Load Required Packages

```{r packages}
required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "plotly", "broom",
  "ggrepel", "ggpubr", "tibble", "purrr", "zoo", "lubridate", "tools"
)

for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
```

## Load the Data

```{r load-data}
ash_data <- read_excel("ash_final.xlsx", sheet = 1)
head(ash_data)
```

## Clean and Prepare Data

```{r clean-data}
ash_clean <- ash_data %>%
  select(same_id, rep_id, ash_per) %>%
  mutate(ash_per = as.numeric(ash_per)) %>%
  drop_na()
```

## Summary Statistics by Replicate

```{r summary-stats}
ash_summary <- ash_clean %>%
  group_by(same_id) %>%
  summarise(
    mean_ash_percent = mean(ash_per, na.rm = TRUE),
    sd_ash_percent = sd(ash_per, na.rm = TRUE),
    n_reps = n()
  ) %>%
  arrange(desc(n_reps))

ash_summary
```

## Bar Plot with Error Bars

```{r plot-bar}
ash_plot <- ggplot(ash_summary, aes(x = reorder(same_id, -mean_ash_percent), y = mean_ash_percent)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = mean_ash_percent - sd_ash_percent, ymax = mean_ash_percent + sd_ash_percent), width = 0.2) +
  labs(title = "Ash Content by Replicate Group", x = "Replicate Group", y = "Mean Ash %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ash_plot
```

## Interactive Plot

```{r plotly-bar, results='asis'}
plotly::ggplotly(ash_plot)
```

## Save Summary to CSV

```{r export-summary, eval=FALSE}
write.csv(ash_summary, "ash_summary_by_replicate.csv", row.names = FALSE)