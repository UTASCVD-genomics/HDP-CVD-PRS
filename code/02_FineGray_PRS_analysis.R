# =============================================================================
# Fine-Gray competing-risk analysis of HDP-related polygenic risk scores
#
# Manuscript:
# Polygenic Risk Scores for Hypertensive Disorders of Pregnancy and
# Later-Life Cardiovascular Disease in Women
#
# Code author: Gayathry Krishnamurthy
#
# Associated manuscript authors:
# Gayathry Krishnamurthy, James A. Temple, Mustafa A. Ahmed, Jennie Hui,
# Peter J. Meikle, Corey Giles, Hoang T. Phan, Shaun P. Brennecke,
# Eric K. Moses, and Phillip E. Melton
#
# Repository: https://github.com/UTASCVD-genomics/PE-CVD-PRS
#
# This script implements the primary competing-risk analysis reported in the
# manuscript and generates the corresponding result tables and forest plots.
#
# Polygenic risk scores:
#   PGS003586 - preeclampsia
#   PGS003587 - gestational hypertension
#   PGS004593 - preeclampsia-related score
#   PGSHDP    - broad hypertensive disorders of pregnancy score
#
# For each first cardiovascular disease outcome, that outcome is coded as the
# event of interest, other first cardiovascular disease outcomes are coded as
# competing events, and women without an incident cardiovascular event are
# censored. Polygenic risk scores are standardised to mean 0 and standard
# deviation 1 in the analytical cohort.
#
# Models include the PRS, HDP history, a PRS-by-HDP-history term, age, age
# squared, systolic blood pressure, waist circumference, total cholesterol to
# HDL-cholesterol ratio, diabetes, lipid-lowering medication use, ever-smoking
# status, and genetic principal components 1-10. Group-specific PRS effects
# are estimated for women with and without a history of HDP and are reported
# as subdistribution hazard ratios with 95% confidence intervals and p-values.
# The PRS-by-HDP-history interaction coefficient is not reported in the
# manuscript.
#
# Peripheral vascular disease and cardiac arrest were excluded from the PRS
# association analyses because of sparse event numbers among women with a
# history of HDP (n = 1 and n = 0, respectively).
#
# Individual-level Busselton Health Study genotype, phenotype, pregnancy-history
# and linked health data are subject to data-access, governance and privacy
# restrictions and are not distributed with this repository.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(cmprsk)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

set.seed(123)

# ============================================================
# 1) SETTINGS
# ============================================================

# Restricted BHS analysis dataset (not distributed publicly)
infile <- "data/BHS_analysis_data.csv"

# Analysis outputs
out_base <- "results"

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

prs_list <- c("PGS003586", "PGS003587", "PGS004593", "PGSHDP")

# FULL model only
covar_sets <- c("FULL")

min_total_interest_events <- 5
min_group_interest_events <- 1

# Outcomes excluded from PRS association analyses because of sparse event
# numbers among women with a history of HDP.
exclude_from_analysis <- c("Peripheral_vascular", "Cardiac_arrest")

# No additional outcome exclusions are required at the plotting stage because
# sparse outcomes are removed before model fitting.
exclude_from_forest <- character(0)

COL_NOHDP <- "#1F78B4"
COL_HDP   <- "#D62728"

BASE_FONT_FAMILY <- "Arial"

PUB <- list(
  wrap_width = 22,
  row_off = 0.13,
  base_font = 10,
  lab_outcome = 3.15,
  lab_row = 3.00,
  head_size = 3.25,
  title_size = 12.2,
  point_size = 2.6,
  err_lwd = 0.80,
  vline_lwd = 0.50,
  out_width = 9.4,
  out_height = 5.35,
  out_dpi = 1000
)

PANEL <- list(
  width = 9.4,
  two_panel_height = 11.0,
  dpi = 1000
)


results_by_key <- list()

# ============================================================
# 2) READ DATA
# ============================================================

df <- read_csv(infile, show_col_types = FALSE)

cat("Original N:", nrow(df), "\n")
cat("Column names:\n")
print(names(df))

# ============================================================
# 3) REQUIRED COLUMN CHECK
# ============================================================

required_cols <- c(
  "ID",
  prs_list,
  "age",
  "age2",
  "SBP",
  "waist",
  "chol_hdl_ratio",
  "diabetes",
  "smokever",
  "meds",
  "CVD_status",
  "HDP",
  paste0("PC", 1:10),
  "AgeFirstCVD_or_Censor",
  "CVD_bucket",
  "follow_up_years"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# 4) HELPER FUNCTIONS
# ============================================================

fmt_num_scalar <- function(x, d = 2) {
  if (is.na(x) || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = d)
}

fmt_p_scalar <- function(p) {
  if (is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

fmt_ci_scalar <- function(hr, lcl, ucl, d = 2) {
  paste0(
    fmt_num_scalar(hr, d),
    " (",
    fmt_num_scalar(lcl, d),
    "\u2013",
    fmt_num_scalar(ucl, d),
    ")"
  )
}

pretty_subtype <- function(x) {
  dplyr::recode(
    x,
    "AF_flutter"                 = "Atrial fibrillation and flutter",
    "Cardiac_arrest"             = "Cardiac arrest",
    "Cerebrovascular"            = "Cerebrovascular diseases",
    "Peripheral_vascular"        = "Peripheral vascular diseases",
    "Venous_lymphatic"           = "Venous and lymphatic diseases",
    "Heart_failure"              = "Heart failure",
    "Hypertensive_CVD"           = "Hypertensive diseases",
    "IHD"                        = "Ischaemic heart disease",
    "Ill_defined_or_unspecified" = "Ill-defined or unspecified circulatory diseases",
    "Other_CVD"                  = "Other cardiovascular diseases",
    "Unknown_CVD"                = "Unknown cardiovascular diseases",
    .default = gsub("_", " ", x)
  )
}

pretty_prs <- function(x) {
  dplyr::recode(
    x,
    "PGS003586" = "PGS003586 (PE PRS)",
    "PGS004593" = "PGS004593 (PE PRS)",
    "PGS003587" = "PGS003587 (GH PRS)",
    "PGSHDP"    = "Broad HDP PRS",
    .default = x
  )
}

recode_hdp <- function(hdp_vec) {
  
  hdp_num <- suppressWarnings(as.numeric(hdp_vec))
  hdp_chr <- tolower(trimws(as.character(hdp_vec)))
  
  vals <- sort(unique(hdp_num[!is.na(hdp_num)]))
  
  if (length(vals) > 0 && all(vals %in% c(0, 1))) {
    return(
      ifelse(
        hdp_num == 1, 1L,
        ifelse(hdp_num == 0, 0L, NA_integer_)
      )
    )
  }
  
  if (length(vals) > 0 && all(vals %in% c(1, 2))) {
    return(
      ifelse(
        hdp_num == 2, 1L,
        ifelse(hdp_num == 1, 0L, NA_integer_)
      )
    )
  }
  
  dplyr::case_when(
    hdp_chr %in% c("hdp", "yes", "y", "true", "case", "cases") ~ 1L,
    hdp_chr %in% c("no hdp", "no_hdp", "no", "n", "false", "control", "controls") ~ 0L,
    TRUE ~ NA_integer_
  )
}

save_pub_figure <- function(fig, outfile_base, width, height, dpi = 1000) {
  
  pdf_file  <- paste0(outfile_base, ".pdf")
  eps_file  <- paste0(outfile_base, ".eps")
  png_file  <- paste0(outfile_base, ".png")
  tiff_file <- paste0(outfile_base, ".tiff")
  
  if (capabilities("cairo")) {
    ggsave(
      filename = pdf_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      family = BASE_FONT_FAMILY,
      limitsize = FALSE,
      bg = "white"
    )
    
    ggsave(
      filename = eps_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_ps,
      family = BASE_FONT_FAMILY,
      onefile = FALSE,
      fallback_resolution = dpi,
      limitsize = FALSE,
      bg = "white"
    )
  } else {
    ggsave(
      filename = pdf_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      device = "pdf",
      family = BASE_FONT_FAMILY,
      limitsize = FALSE,
      bg = "white"
    )
    
    ggsave(
      filename = eps_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      device = "eps",
      family = BASE_FONT_FAMILY,
      limitsize = FALSE,
      bg = "white"
    )
  }
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      filename = png_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_png,
      limitsize = FALSE,
      bg = "white"
    )
    
    ggsave(
      filename = tiff_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw",
      limitsize = FALSE,
      bg = "white"
    )
  } else {
    ggsave(
      filename = png_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE,
      bg = "white"
    )
    
    ggsave(
      filename = tiff_file,
      plot = fig,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      limitsize = FALSE,
      bg = "white"
    )
  }
  
  cat("\nSaved figure files:\n")
  cat(pdf_file, "\n")
  cat(eps_file, "\n")
  cat(tiff_file, "\n")
  cat(png_file, "\n")
}

extract_prs_by_group <- function(beta, V, prs_term) {
  nm <- names(beta)

  prs_candidates <- unique(c(
    prs_term,
    paste0("`", prs_term, "`")
  ))

  int_candidates <- unique(c(
    paste0(prs_term, ":HDP_hist"),
    paste0("HDP_hist:", prs_term),
    paste0("`", prs_term, "`:HDP_hist"),
    paste0("HDP_hist:`", prs_term, "`")
  ))

  i_prs <- which(nm %in% prs_candidates)
  i_int <- which(nm %in% int_candidates)

  if (length(i_prs) == 0 || length(i_int) == 0) {
    return(data.frame(
      beta_NoHDP = NA_real_,
      HR_NoHDP = NA_real_,
      LCL_NoHDP = NA_real_,
      UCL_NoHDP = NA_real_,
      p_NoHDP = NA_real_,
      beta_HDP = NA_real_,
      HR_HDP = NA_real_,
      LCL_HDP = NA_real_,
      UCL_HDP = NA_real_,
      p_HDP = NA_real_
    ))
  }

  i_prs <- i_prs[1]
  i_int <- i_int[1]

  b_prs <- beta[i_prs]
  var_prs <- V[i_prs, i_prs]
  se_prs <- ifelse(is.finite(var_prs) && var_prs > 0, sqrt(var_prs), NA_real_)

  beta_nohdp <- b_prs
  HR_nohdp <- exp(beta_nohdp)
  LCL_nohdp <- exp(beta_nohdp - 1.96 * se_prs)
  UCL_nohdp <- exp(beta_nohdp + 1.96 * se_prs)
  p_nohdp <- 2 * pnorm(abs(beta_nohdp / se_prs), lower.tail = FALSE)

  b_int <- beta[i_int]
  var_int <- V[i_int, i_int]
  cov_prs_int <- V[i_prs, i_int]

  beta_hdp <- b_prs + b_int
  var_hdp <- var_prs + var_int + 2 * cov_prs_int
  se_hdp <- ifelse(is.finite(var_hdp) && var_hdp > 0, sqrt(var_hdp), NA_real_)

  HR_hdp <- exp(beta_hdp)
  LCL_hdp <- exp(beta_hdp - 1.96 * se_hdp)
  UCL_hdp <- exp(beta_hdp + 1.96 * se_hdp)
  p_hdp <- 2 * pnorm(abs(beta_hdp / se_hdp), lower.tail = FALSE)

  data.frame(
    beta_NoHDP = beta_nohdp,
    HR_NoHDP = HR_nohdp,
    LCL_NoHDP = LCL_nohdp,
    UCL_NoHDP = UCL_nohdp,
    p_NoHDP = p_nohdp,
    beta_HDP = beta_hdp,
    HR_HDP = HR_hdp,
    LCL_HDP = LCL_hdp,
    UCL_HDP = UCL_hdp,
    p_HDP = p_hdp
  )
}

# ============================================================
# 5) PREPARE DATA
# ============================================================

cat("\nObserved HDP values:\n")
print(sort(unique(df$HDP)))

cat("\nObserved CVD_status values:\n")
print(sort(unique(df$CVD_status)))

cat("\nCVD_bucket distribution before filtering:\n")
print(table(df$CVD_bucket, useNA = "ifany"))

df2 <- df %>%
  mutate(
    ID = as.character(ID),
    fu_time = as.numeric(follow_up_years),
    age_fg = as.numeric(age),
    age2_fg = as.numeric(age2),
    SBP = as.numeric(SBP),
    waist = as.numeric(waist),
    chol_hdl_ratio = as.numeric(chol_hdl_ratio),
    diabetes = as.numeric(diabetes),
    smokever = as.numeric(smokever),
    meds = as.numeric(meds),
    across(all_of(paste0("PC", 1:10)), as.numeric),
    event_allCVD = as.integer(as.numeric(CVD_status) == 1),
    HDP_hist = recode_hdp(HDP),
    CVD_bucket = as.character(CVD_bucket),
    CVD_bucket = case_when(
      event_allCVD == 0 ~ "No_CVD",
      event_allCVD == 1 & (is.na(CVD_bucket) | CVD_bucket == "") ~ "Unknown_CVD",
      TRUE ~ CVD_bucket
    )
  ) %>%
  filter(
    !is.na(fu_time),
    fu_time > 0,
    !is.na(event_allCVD),
    !is.na(HDP_hist),
    !is.na(age_fg),
    !is.na(age2_fg)
  )

cat("\nN after Fine-Gray filtering:", nrow(df2), "\n")
cat("All CVD events:", sum(df2$event_allCVD == 1, na.rm = TRUE), "\n")
cat("Controls/non-events:", sum(df2$event_allCVD == 0, na.rm = TRUE), "\n")
cat("No-HDP N:", sum(df2$HDP_hist == 0, na.rm = TRUE), "\n")
cat("HDP N:", sum(df2$HDP_hist == 1, na.rm = TRUE), "\n")
cat("No-HDP CVD events:", sum(df2$event_allCVD == 1 & df2$HDP_hist == 0, na.rm = TRUE), "\n")
cat("HDP CVD events:", sum(df2$event_allCVD == 1 & df2$HDP_hist == 1, na.rm = TRUE), "\n")

cat("\nCheck HDP recoding:\n")
print(table(df2$HDP, df2$HDP_hist, useNA = "ifany"))

cat("\nCVD_bucket distribution after filtering:\n")
print(table(df2$CVD_bucket, useNA = "ifany"))

# ============================================================
# 6) STANDARDISE PRS
# ============================================================

df2 <- df2 %>%
  mutate(
    across(
      all_of(prs_list),
      ~ as.numeric(scale(as.numeric(.x))),
      .names = "{.col}_z"
    )
  )

prs_z_list <- paste0(prs_list, "_z")

cat("\nPRS standardisation check:\n")
print(
  df2 %>%
    summarise(
      across(
        all_of(prs_z_list),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE)
        )
      )
    )
)

# IMPORTANT: participant-level analysis data are not written to disk
# in the public repository version.

# ============================================================
# 7) FOREST PLOT AND TABLE FUNCTION
# ============================================================

make_fg_pub_forest <- function(results_tbl_with_counts,
                               prs_name,
                               covar_set,
                               outfile_base = NULL,
                               exclude_subtypes = character(0),
                               save_now = TRUE,
                               show_legend = TRUE,
                               show_x_title = TRUE,
                               panel_label = NULL,
                               title_size = PUB$title_size) {
  
  base <- results_tbl_with_counts %>%
    filter(!CVD_subtype %in% exclude_subtypes) %>%
    mutate(
      subtype_pretty = pretty_subtype(CVD_subtype)
    ) %>%
    arrange(desc(events_total), subtype_pretty)
  
  if (nrow(base) == 0) return(NULL)
  
  long <- base %>%
    select(
      CVD_subtype,
      subtype_pretty,
      events_NoHDP,
      events_HDP,
      HR_NoHDP,
      LCL_NoHDP,
      UCL_NoHDP,
      p_NoHDP,
      HR_HDP,
      LCL_HDP,
      UCL_HDP,
      p_HDP
    ) %>%
    pivot_longer(
      cols = c(
        HR_NoHDP, LCL_NoHDP, UCL_NoHDP, p_NoHDP,
        HR_HDP,   LCL_HDP,   UCL_HDP,   p_HDP
      ),
      names_to = c(".value", "group"),
      names_pattern = "(HR|LCL|UCL|p)_(NoHDP|HDP)"
    ) %>%
    rename(p_value = p) %>%
    mutate(
      group = ifelse(group == "NoHDP", "No HDP", "HDP"),
      group = factor(group, levels = c("HDP", "No HDP")),
      event_n = ifelse(group == "No HDP", events_NoHDP, events_HDP),
      group_label = paste0(as.character(group), " (n=", event_n, ")"),
      hr_ci = mapply(fmt_ci_scalar, HR, LCL, UCL),
      p_label = vapply(p_value, fmt_p_scalar, character(1)),
      subtype_wrapped = str_wrap(subtype_pretty, width = PUB$wrap_width)
    ) %>%
    filter(is.finite(HR), is.finite(LCL), is.finite(UCL))
  
  if (nrow(long) == 0) return(NULL)
  
  subtype_levels <- rev(unique(long$subtype_wrapped))
  
  y_map <- data.frame(
    subtype_wrapped = subtype_levels,
    y_base = seq_along(subtype_levels)
  )
  
  long <- long %>%
    left_join(y_map, by = "subtype_wrapped") %>%
    mutate(
      y_row = y_base + ifelse(group == "No HDP", PUB$row_off, -PUB$row_off)
    )
  
  subtype_df <- y_map %>%
    mutate(label = subtype_wrapped)
  
  y_min <- min(long$y_row) - 0.50
  y_max <- max(long$y_row) + 0.95
  y_header <- y_max - 0.22
  
  x_outcome <- 0.00
  x_group   <- 1.55
  x_hr      <- 2.85
  x_p       <- 4.05
  
  x_min_tbl <- -0.02
  x_max_tbl <- 4.45
  
  p_table <- ggplot() +
    geom_rect(
      data = y_map %>% mutate(band = y_base %% 2),
      aes(
        xmin = x_min_tbl,
        xmax = x_max_tbl,
        ymin = y_base - 0.52,
        ymax = y_base + 0.52,
        fill = factor(band)
      ),
      alpha = 0.10,
      colour = NA
    ) +
    scale_fill_manual(
      values = c("0" = "grey85", "1" = "white"),
      guide = "none"
    ) +
    geom_text(
      data = subtype_df,
      aes(x = x_outcome, y = y_base, label = label),
      hjust = 0,
      size = PUB$lab_outcome,
      family = BASE_FONT_FAMILY
    ) +
    geom_text(
      data = long,
      aes(x = x_group, y = y_row, label = group_label),
      hjust = 0,
      size = PUB$lab_row,
      family = BASE_FONT_FAMILY
    ) +
    geom_text(
      data = long,
      aes(x = x_hr, y = y_row, label = hr_ci),
      hjust = 0,
      size = PUB$lab_row,
      family = BASE_FONT_FAMILY
    ) +
    geom_text(
      data = long,
      aes(x = x_p, y = y_row, label = p_label),
      hjust = 0,
      size = PUB$lab_row,
      family = BASE_FONT_FAMILY
    ) +
    annotate(
      "text",
      x = x_outcome,
      y = y_header,
      label = "First CVD event",
      hjust = 0,
      fontface = "bold",
      size = PUB$head_size,
      family = BASE_FONT_FAMILY
    ) +
    annotate(
      "text",
      x = x_group,
      y = y_header,
      label = "Group (event n)",
      hjust = 0,
      fontface = "bold",
      size = PUB$head_size,
      family = BASE_FONT_FAMILY
    ) +
    annotate(
      "text",
      x = x_hr,
      y = y_header,
      label = "sHR (95% CI)",
      hjust = 0,
      fontface = "bold",
      size = PUB$head_size,
      family = BASE_FONT_FAMILY
    ) +
    annotate(
      "text",
      x = x_p,
      y = y_header,
      label = "P-value",
      hjust = 0,
      fontface = "bold",
      size = PUB$head_size,
      family = BASE_FONT_FAMILY
    ) +
    scale_x_continuous(
      limits = c(x_min_tbl, x_max_tbl),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_map$y_base,
      labels = NULL,
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = PUB$base_font, base_family = BASE_FONT_FAMILY) +
    theme(
      plot.margin = margin(t = 8, r = 2, b = 8, l = 6)
    )
  # Fixed x-axis range for comparability across panels
  x_min <- 0.25
  x_max <- 8.0
  
  axis_breaks <- c(0.3, 1.0, 3.0, 7.0)

  adjust_lab <- ifelse(
    covar_set == "FULL",
    "fully adjusted",
    "minimally adjusted"
  )
  
  x_title <- paste0(
    "sHR per 1-SD increase in PRS\n",
    "(", adjust_lab, "; log scale)"
  )
  
  p_forest <- ggplot(long, aes(x = HR, y = y_row)) +
    geom_vline(
      xintercept = 1,
      linetype = 2,
      linewidth = PUB$vline_lwd,
      colour = "black"
    ) +
    geom_segment(
      aes(
        x = LCL,
        xend = UCL,
        y = y_row,
        yend = y_row,
        colour = group
      ),
      linewidth = PUB$err_lwd,
      lineend = "round"
    ) +
    geom_point(
      aes(shape = group, colour = group),
      size = PUB$point_size
    ) +
    scale_colour_manual(
      values = c("No HDP" = COL_NOHDP, "HDP" = COL_HDP),
      breaks = c("HDP", "No HDP")
    ) +
    scale_shape_manual(
      values = c("No HDP" = 17, "HDP" = 16),
      breaks = c("HDP", "No HDP")
    ) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = axis_breaks,
      labels = format(axis_breaks, nsmall = 1),
      expand = expansion(mult = c(0.04, 0.10))
    ) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_map$y_base,
      labels = NULL,
      expand = expansion(mult = 0)
    ) +
    labs(
      x = ifelse(show_x_title, x_title, ""),
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = PUB$base_font, base_family = BASE_FONT_FAMILY) +
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      legend.justification = "center",
      legend.text = element_text(size = 9),
      legend.margin = margin(b = 2),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = if (show_x_title) {
        element_text(
          face = "bold",
          colour = "black",
          size = 9.5,
          margin = margin(t = 8)
        )
      } else {
        element_blank()
      },
      axis.text.x = element_text(size = 8.5),
      plot.margin = margin(t = 8, r = 10, b = 12, l = 2)
    )
  
  fig_title <- pretty_prs(prs_name)
  
  if (!is.null(panel_label)) {
    fig_title <- paste0(panel_label, ". ", fig_title)
  }
  
  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = fig_title,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = title_size / ggplot2::.pt,
      family = BASE_FONT_FAMILY,
      colour = "black"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(t = 2, r = 8, b = 0, l = 8)
    )
  
  body_plot <- (p_table | p_forest) +
    plot_layout(widths = c(0.95, 1.00))
  
  fig <- title_plot / body_plot +
    plot_layout(heights = c(0.10, 1.00)) &
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(t = 4, r = 10, b = 6, l = 8)
    )
  
  if (isTRUE(save_now) && !is.null(outfile_base)) {
    save_pub_figure(
      fig = fig,
      outfile_base = outfile_base,
      width = PUB$out_width,
      height = PUB$out_height,
      dpi = PUB$out_dpi
    )
  }
  
  return(fig)
}

# ============================================================
# 8) RUN FINE-GRAY MODELS FOR ONE PRS
# ============================================================

run_fg <- function(prs_name, covar_set = c("FULL")) {
  
  covar_set <- match.arg(covar_set)
  prs_z_name <- paste0(prs_name, "_z")
  
  out_dir <- file.path(out_base, paste0("FG_", prs_name, "__", covar_set))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  counts_tbl <- df2 %>%
    filter(event_allCVD == 1) %>%
    count(CVD_bucket, HDP_hist, name = "n_events") %>%
    pivot_wider(
      names_from = HDP_hist,
      values_from = n_events,
      values_fill = 0
    )
  
  if (!"0" %in% names(counts_tbl)) counts_tbl$`0` <- 0
  if (!"1" %in% names(counts_tbl)) counts_tbl$`1` <- 0
  
  counts_tbl <- counts_tbl %>%
    rename(
      CVD_subtype = CVD_bucket,
      events_NoHDP = `0`,
      events_HDP = `1`
    ) %>%
    mutate(
      events_total = events_NoHDP + events_HDP
    ) %>%
    filter(
      !CVD_subtype %in% c(
        "Controls",
        "Control",
        "controls",
        "control",
        "No_CVD",
        "No CVD",
        "Non_CVD",
        "Unknown",
        ""
      )
    )
  
  subtypes_valid <- counts_tbl %>%
    filter(
      !CVD_subtype %in% exclude_from_analysis,
      events_NoHDP >= min_group_interest_events,
      events_HDP >= min_group_interest_events,
      events_total >= min_total_interest_events
    ) %>%
    pull(CVD_subtype)
  
  if (length(subtypes_valid) == 0) {
    message("No valid subtypes for ", prs_name, " | ", covar_set)
    return(NULL)
  }
  
  covar_vec <- c(
    "age_fg",
    "age2_fg",
    "waist",
    "chol_hdl_ratio",
    "SBP",
    "diabetes",
    "smokever",
    "meds",
    paste0("PC", 1:10)
  )
  
  covar_vec <- intersect(covar_vec, names(df2))
  
  results_list <- list()
  
  for (st in subtypes_valid) {
    
    dat <- df2 %>%
      mutate(
        cr_status = case_when(
          event_allCVD == 1 & CVD_bucket == st ~ 1L,
          event_allCVD == 1 & CVD_bucket != st ~ 2L,
          TRUE ~ 0L
        )
      )
    
    fmla <- as.formula(
      paste0(
        "~ `",
        prs_z_name,
        "` * HDP_hist + ",
        paste(covar_vec, collapse = " + ")
      )
    )
    
    X <- model.matrix(fmla, data = dat)[, -1, drop = FALSE]
    
    keep <- complete.cases(dat$fu_time, dat$cr_status, X)
    
    dat_keep <- dat[keep, ]
    X_keep <- X[keep, , drop = FALSE]
    
    n_used <- nrow(dat_keep)
    n_interest <- sum(dat_keep$cr_status == 1, na.rm = TRUE)
    n_competing <- sum(dat_keep$cr_status == 2, na.rm = TRUE)
    n_censored <- sum(dat_keep$cr_status == 0, na.rm = TRUE)
    
    if (n_interest < min_total_interest_events) next
    
    fg <- tryCatch(
      {
        crr(
          ftime = dat_keep$fu_time,
          fstatus = dat_keep$cr_status,
          cov1 = X_keep,
          failcode = 1,
          cencode = 0
        )
      },
      error = function(e) {
        message(
          "Fine-Gray failed for ",
          prs_name,
          " | ",
          covar_set,
          " | ",
          st,
          ": ",
          e$message
        )
        return(NULL)
      }
    )
    
    if (is.null(fg)) next
    
    eff <- extract_prs_by_group(fg$coef, fg$var, prs_z_name)
    
    eff$PRS <- prs_name
    eff$PRS_used <- prs_z_name
    eff$Covariates <- covar_set
    eff$CVD_subtype <- st
    eff$N_used <- n_used
    eff$Interest_events <- n_interest
    eff$Competing_events <- n_competing
    eff$Censored <- n_censored
    
    results_list[[length(results_list) + 1]] <- eff
  }
  
  results_tbl <- bind_rows(results_list)
  
  if (nrow(results_tbl) == 0) {
    message("No valid Fine-Gray results for ", prs_name, " | ", covar_set)
    return(NULL)
  }
  
  results_tbl_with_counts <- results_tbl %>%
    left_join(counts_tbl, by = "CVD_subtype") %>%
    arrange(desc(events_total), CVD_subtype)
  
  public_results <- results_tbl_with_counts %>%
    select(
      PRS, PRS_used, Covariates, CVD_subtype,
      N_used, Interest_events, Competing_events, Censored,
      events_total, events_NoHDP, events_HDP,
      beta_NoHDP, HR_NoHDP, LCL_NoHDP, UCL_NoHDP, p_NoHDP,
      beta_HDP, HR_HDP, LCL_HDP, UCL_HDP, p_HDP
    )

  write_csv(
    public_results,
    file.path(
      out_dir,
      paste0("FineGray_results_", prs_name, "_", covar_set, ".csv")
    )
  )
  
  results_for_paper <- results_tbl_with_counts %>%
    mutate(
      CVD_event = pretty_subtype(CVD_subtype),
      PRS_label = pretty_prs(PRS),
      NoHDP_sHR_95CI = mapply(fmt_ci_scalar, HR_NoHDP, LCL_NoHDP, UCL_NoHDP),
      HDP_sHR_95CI = mapply(fmt_ci_scalar, HR_HDP, LCL_HDP, UCL_HDP),
      p_NoHDP_fmt = vapply(p_NoHDP, fmt_p_scalar, character(1)),
      p_HDP_fmt = vapply(p_HDP, fmt_p_scalar, character(1))
    ) %>%
    select(
      PRS_label,
      Covariates,
      CVD_event,
      events_NoHDP,
      events_HDP,
      NoHDP_sHR_95CI,
      p_NoHDP_fmt,
      HDP_sHR_95CI,
      p_HDP_fmt
    )
  
  write_csv(
    results_for_paper,
    file.path(
      out_dir,
      paste0("Paper_table_", prs_name, "_", covar_set, ".csv")
    )
  )
  
  key <- paste(prs_name, covar_set, sep = "_")
  results_by_key[[key]] <<- results_tbl_with_counts
  
  make_fg_pub_forest(
    results_tbl_with_counts = results_tbl_with_counts,
    prs_name = prs_name,
    covar_set = covar_set,
    outfile_base = file.path(
      out_dir,
      paste0("FineGray_FOREST_TABLE_", prs_name, "_", covar_set)
    ),
    exclude_subtypes = exclude_from_forest,
    save_now = TRUE,
    show_legend = TRUE,
    show_x_title = TRUE,
    panel_label = NULL
  )
  
  return(results_tbl_with_counts)
}

# ============================================================
# 9) RUN ALL FULLY ADJUSTED MODELS
# ============================================================

master_rows <- list()
k <- 0

for (prs in prs_list) {
  
  cat("\nRunning:", prs, "| FULL\n")
  
  res <- run_fg(prs, "FULL")
  
  if (!is.null(res) && nrow(res) > 0) {
    k <- k + 1
    master_rows[[k]] <- res
  }
}

# ============================================================
# 10) SAVE MASTER TABLES
# ============================================================

if (length(master_rows) > 0) {
  
  master_fg <- bind_rows(master_rows)
  
  master_fg_clean <- master_fg %>%
    select(
      PRS,
      PRS_used,
      Covariates,
      CVD_subtype,
      N_used,
      Interest_events,
      Competing_events,
      Censored,
      events_total,
      events_NoHDP,
      events_HDP,
      beta_NoHDP,
      HR_NoHDP,
      LCL_NoHDP,
      UCL_NoHDP,
      p_NoHDP,
      beta_HDP,
      HR_HDP,
      LCL_HDP,
      UCL_HDP,
      p_HDP
    )
  
  write_csv(
    master_fg_clean,
    file.path(out_base, "MASTER_FineGray_PRS_effects_CLEAN_FULL.csv")
  )
  
  master_paper_table <- master_fg_clean %>%
    mutate(
      PRS_label = pretty_prs(PRS),
      CVD_event = pretty_subtype(CVD_subtype),
      NoHDP_sHR_95CI = mapply(fmt_ci_scalar, HR_NoHDP, LCL_NoHDP, UCL_NoHDP),
      HDP_sHR_95CI = mapply(fmt_ci_scalar, HR_HDP, LCL_HDP, UCL_HDP),
      p_NoHDP_fmt = vapply(p_NoHDP, fmt_p_scalar, character(1)),
      p_HDP_fmt = vapply(p_HDP, fmt_p_scalar, character(1))
    ) %>%
    select(
      PRS_label,
      Covariates,
      CVD_event,
      events_NoHDP,
      events_HDP,
      NoHDP_sHR_95CI,
      p_NoHDP_fmt,
      HDP_sHR_95CI,
      p_HDP_fmt
    )
  
  write_csv(
    master_paper_table,
    file.path(out_base, "MASTER_PAPER_TABLE_formatted_FULL.csv")
  )
  
  cat("\nFine-Gray master tables and individual figures saved in:\n")
  cat(out_base, "\n")
  
} else {
  cat("\nNo valid Fine-Gray models ran. Check event counts and missingness.\n")
}

# ============================================================
# 11) CREATE COMBINED PANEL FIGURES
# ============================================================

make_combined_prs_panel <- function(prs_order,
                                    outfile_name,
                                    covar_set = "FULL",
                                    show_only_first_legend = TRUE) {
  
  panel_keys <- paste(prs_order, covar_set, sep = "_")
  available_keys <- panel_keys[panel_keys %in% names(results_by_key)]
  
  if (length(available_keys) == 0) {
    message("No panels available for ", outfile_name)
    return(NULL)
  }
  
  panel_plots <- list()
  
  for (i in seq_along(available_keys)) {
    
    key <- available_keys[i]
    prs_name <- sub(paste0("_", covar_set, "$"), "", key)
    
    show_leg <- TRUE
    if (isTRUE(show_only_first_legend)) {
      show_leg <- i == 1
    }
    
    show_x <- i == length(available_keys)
    
    panel_plots[[i]] <- make_fg_pub_forest(
      results_tbl_with_counts = results_by_key[[key]],
      prs_name = prs_name,
      covar_set = covar_set,
      outfile_base = NULL,
      exclude_subtypes = exclude_from_forest,
      save_now = FALSE,
      show_legend = show_leg,
      show_x_title = show_x,
      panel_label = LETTERS[i],
      title_size = 10.8
    )
  }
  
  panel_plots <- panel_plots[!sapply(panel_plots, is.null)]
  
  if (length(panel_plots) == 0) {
    message("No valid plots to combine for ", outfile_name)
    return(NULL)
  }
  
  combined_fig <- wrap_plots(
    plotlist = panel_plots,
    ncol = 1
  ) &
    theme(
      plot.margin = margin(t = 6, r = 8, b = 6, l = 8),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  save_pub_figure(
    fig = combined_fig,
    outfile_base = file.path(out_base, outfile_name),
    width = PANEL$width,
    height = PANEL$two_panel_height,
    dpi = PANEL$dpi
  )
  
  return(combined_fig)
}

# Figure 2: two preeclampsia PRS
make_combined_prs_panel(
  prs_order = c("PGS003586", "PGS004593"),
  outfile_name = "Figure_2_Two_Preeclampsia_PRS_Forest_Table_FULL"
)

# Figure 3: gestational hypertension + broad HDP PRS
make_combined_prs_panel(
  prs_order = c("PGS003587", "PGSHDP"),
  outfile_name = "Figure_3_GH_Broad_HDP_PRS_Forest_Table_FULL"
)

# Save the R session information used for the analysis.
capture.output(
  sessionInfo(),
  file = file.path(out_base, "sessionInfo.txt")
)

cat("\nCombined figures saved.\n")
cat("\nMain combined files:\n")
cat(file.path(out_base, "Figure_2_Two_Preeclampsia_PRS_Forest_Table_FULL.pdf"), "\n")
cat(file.path(out_base, "Figure_3_GH_Broad_HDP_PRS_Forest_Table_FULL.pdf"), "\n")
cat("\nAnalysis complete.\n")
