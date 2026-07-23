# Survival analysis helpers for bulk RNA-seq templates.

# Validate survival data frame.
validate_surv_df <- function(surv_df, time_col = "time", status_col = "status") {
  if (!all(c(time_col, status_col) %in% colnames(surv_df))) {
    stop("surv_df must contain columns: ", paste(c(time_col, status_col), collapse = ", "))
  }
  if (any(surv_df[[time_col]] <= 0, na.rm = TRUE)) {
    stop("Survival time must be > 0")
  }
  if (!all(surv_df[[status_col]] %in% c(0, 1), na.rm = TRUE)) {
    stop("Survival status must be 0 (censored) or 1 (event)")
  }
  invisible(TRUE)
}

# Run univariate Cox regression for one or more variables.
run_univariate_cox <- function(surv_df, vars, time_col = "time", status_col = "status") {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  validate_surv_df(surv_df, time_col, status_col)
  rows <- list()
  for (var in vars) {
    if (!var %in% colnames(surv_df)) {
      warning("Variable not found: ", var)
      next
    }
    sub <- surv_df[!is.na(surv_df[[var]]), c(time_col, status_col, var), drop = FALSE]
    if (nrow(sub) < 5) {
      warning("Too few non-missing observations for: ", var)
      next
    }
    formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", status_col, ") ~ ", var))
    fit <- tryCatch(survival::coxph(formula, data = sub), error = function(e) NULL)
    if (is.null(fit)) {
      warning("Univariate Cox failed for: ", var)
      next
    }
    sm <- summary(fit)
    coef_row <- sm$coefficients
    conf <- sm$conf.int
    rows[[var]] <- data.frame(
      variable = var,
      n = sm$n,
      hazard_ratio = conf[, "exp(coef)"],
      HR_lower_95 = conf[, "lower .95"],
      HR_upper_95 = conf[, "upper .95"],
      beta = coef_row[, "coef"],
      se = coef_row[, "se(coef)"],
      z = coef_row[, "z"],
      pvalue = coef_row[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# Run multivariate Cox regression with selected variables.
run_multivariate_cox <- function(surv_df, vars, time_col = "time", status_col = "status") {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  validate_surv_df(surv_df, time_col, status_col)
  missing_vars <- setdiff(vars, colnames(surv_df))
  if (length(missing_vars) > 0) {
    stop("Variables not found: ", paste(missing_vars, collapse = ", "))
  }
  sub <- surv_df[, c(time_col, status_col, vars), drop = FALSE]
  sub <- sub[stats::complete.cases(sub), ]
  if (nrow(sub) < 5) stop("Too few complete cases for multivariate Cox.")

  rhs <- paste(vars, collapse = " + ")
  formula <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", status_col, ") ~ ", rhs))
  fit <- survival::coxph(formula, data = sub)
  sm <- summary(fit)
  coef_row <- sm$coefficients
  conf <- sm$conf.int

  data.frame(
    variable = rownames(coef_row),
    n = sm$n,
    hazard_ratio = conf[, "exp(coef)"],
    HR_lower_95 = conf[, "lower .95"],
    HR_upper_95 = conf[, "upper .95"],
    beta = coef_row[, "coef"],
    se = coef_row[, "se(coef)"],
    z = coef_row[, "z"],
    pvalue = coef_row[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

# Plot forest plot for Cox regression results.
plot_cox_forest_pdf <- function(cox_df, filename,
                                 title = "Cox Regression Forest Plot",
                                 width = 8, height = 6) {
  if (is.null(cox_df) || nrow(cox_df) == 0) return(invisible(NULL))
  cox_df <- cox_df |>
    dplyr::mutate(
      label = paste0(
        "HR=", sprintf("%.2f", .data$hazard_ratio),
        " (", sprintf("%.2f", .data$HR_lower_95), "-", sprintf("%.2f", .data$HR_upper_95), ")"
      ),
      sig = ifelse(.data$pvalue < 0.05, "p < 0.05", "ns")
    ) |>
    dplyr::arrange(.data$hazard_ratio)
  cox_df$variable <- factor(cox_df$variable, levels = cox_df$variable)

  p <- ggplot2::ggplot(cox_df, ggplot2::aes(x = .data$hazard_ratio, y = .data$variable)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#555555", linewidth = 0.6) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$HR_lower_95, xmax = .data$HR_upper_95), height = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(color = .data$sig), size = 3.5) +
    ggplot2::scale_color_manual(values = c("p < 0.05" = "#C8473E", "ns" = "#888888")) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), hjust = -0.15, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.35))) +
    ggplot2::labs(x = "Hazard Ratio (95% CI)", y = NULL, title = title, color = "Significance") +
    theme_publication(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Plot Kaplan-Meier by a categorical grouping variable.
plot_km_by_group_pdf <- function(surv_df, group_col, filename,
                                 title = NULL, ylab = "Survival probability",
                                 xlab = "Time in months", time_unit = "month",
                                 xlim = NULL, ylim = c(0, 1)) {
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("survminer", quietly = TRUE)) {
    stop("Packages 'survival' and 'survminer' are required for KM plots.")
  }
  validate_surv_df(surv_df)
  if (!group_col %in% colnames(surv_df)) {
    stop("Group column not found: ", group_col)
  }

  plot_data <- surv_df[!is.na(surv_df[[group_col]]), ]
  if (nrow(plot_data) == 0) return(invisible(NULL))
  plot_data[[group_col]] <- factor(plot_data[[group_col]])

  if (time_unit == "month") {
    plot_data$time <- plot_data$time / 30.4375
  } else if (time_unit == "year") {
    plot_data$time <- plot_data$time / 365.25
  }

  # Inline the formula into the survfit call so ggsurvplot() can re-evaluate it
  # later; passing a symbol (e.g. `surv_formula`) leaves fit$call$formula as that
  # symbol, which breaks with "object of type 'symbol' is not subsettable".
  surv_formula <- stats::as.formula(paste0("survival::Surv(time, status) ~ ", group_col))
  fit <- eval(substitute(survival::survfit(F, data = plot_data), list(F = surv_formula)))
  p <- survminer::ggsurvplot(
    fit,
    data = plot_data,
    title = title %||% paste(group_col, "survival"),
    xlab = xlab,
    ylab = ylab,
    legend.title = group_col,
    pval = TRUE,
    pval.size = 5,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    surv.scale = "percent",
    ggtheme = ggplot2::theme_classic(base_size = 12),
    xlim = xlim,
    ylim = ylim
  )

  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = 7, height = 8)
  print(p, newpage = FALSE)
  grDevices::dev.off()
  invisible(p)
}

# Split a continuous variable into quantile-based groups for KM plots.
stratify_by_quantile <- function(x, n_groups = 2, labels = NULL) {
  if (is.null(labels)) labels <- if (n_groups == 2) c("Low", "High") else paste0("Q", seq_len(n_groups))
  if (length(unique(x)) < n_groups) {
    stop("Not enough unique values for ", n_groups, " quantile groups")
  }
  cut(x, breaks = stats::quantile(x, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE),
      labels = labels, include.lowest = TRUE)
}
