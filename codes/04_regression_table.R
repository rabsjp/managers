# ============================================================
# 04_regression_table.R
#
# Exports fit1, fit1_1, fit2, fit2_1 (from 03_regressions.R) as a single
# LaTeX table. All group-level fixed effects (factor(group_uid)) are
# omitted from the printed table; a note below the table says they are
# included.
# ============================================================

repo_root <- {
  cand <- c(getwd(), "/home/claude/managers_v6", "/home/claude/managers_v5", "/home/claude/managers")
  cand[which(sapply(cand, function(x) dir.exists(file.path(x, "data"))))[1]]
}
source(file.path(repo_root, "codes", "03_regressions.R"))

out_dir <- file.path(repo_root, "output", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

models <- list(fit1 = fit1, fit1_1 = fit1_1, fit2 = fit2, fit2_1 = fit2_1)

fe_prefix <- "factor(group_uid)"

# display label for each coefficient (order = row order in the table)
var_labels <- c(
  "player.genderWoman" = "Woman",
  "rel_belief_inv"     = "Investment belief (log diff)",
  "rel_belief_price"   = "Price belief (log diff)"
)

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}
num  <- function(x) sprintf("%.3f", x)

# ---- pull coefficient/SE/p for each model, keyed by variable name --------
model_coef <- lapply(models, function(f) {
  co <- summary(f)$coefficients
  keep <- !startsWith(rownames(co), fe_prefix)
  co[keep, , drop = FALSE]
})

model_n    <- sapply(models, function(f) length(summary(f)$residuals))
model_r2   <- sapply(models, function(f) summary(f)$r.squared)
model_adjr2<- sapply(models, function(f) summary(f)$adj.r.squared)
model_nfe  <- sapply(models, function(f) {
  co <- summary(f)$coefficients
  sum(startsWith(rownames(co), fe_prefix))
})

# ---- build the LaTeX body ------------------------------------------------
n_models <- length(models)

coef_lines <- character(0)
for (vn in names(var_labels)) {
  est_cells <- character(n_models)
  se_cells  <- character(n_models)
  for (j in seq_len(n_models)) {
    co <- model_coef[[j]]
    if (vn %in% rownames(co)) {
      est <- co[vn, "Estimate"]; se <- co[vn, "Std. Error"]; p <- co[vn, "Pr(>|t|)"]
      est_cells[j] <- paste0(num(est), stars(p))
      se_cells[j]  <- paste0("(", num(se), ")")
    } else {
      est_cells[j] <- ""
      se_cells[j]  <- ""
    }
  }
  coef_lines <- c(coef_lines,
                   paste0(var_labels[vn], " & ", paste(est_cells, collapse = " & "), " \\\\"),
                   paste0(" & ", paste(se_cells, collapse = " & "), " \\\\"))
}

tex_lines <- c(
  "  \\begin{tabular}{l*{4}{c}}",
  "    \\toprule",
  "     & \\multicolumn{2}{c}{Exp $-$ No} & \\multicolumn{2}{c}{ M $-$ W} \\\\",
  "    \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "     & (1) & (2) & (3) & (4) \\\\",
  "    \\midrule",
  coef_lines,
  "    \\midrule",
  paste0("    Group fixed effects & ", paste(rep("Yes", n_models), collapse = " & "), " \\\\"),
  paste0("    $N$ & ", paste(model_n, collapse = " & "), " \\\\"),
  paste0("    $R^2$ & ", paste(num(model_r2), collapse = " & "), " \\\\"),
  paste0("    Adj. $R^2$ & ", paste(num(model_adjr2), collapse = " & "), " \\\\"),
  "    \\bottomrule",
  "  \\end{tabular}"
)

writeLines(tex_lines, file.path(out_dir, "regression_table.tex"))
message("Regression LaTeX table written to ", file.path(out_dir, "regression_table.tex"))
