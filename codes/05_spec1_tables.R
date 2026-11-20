# ============================================================
# 05_spec1_tables.R
#
# Rebuilds SPEC 1 (relative bid EXP - NO, rounds 1-10) as a 5-column
# progressive table, per the latest instructions:
#
#   spec1: rel_bid ~ rel_belief_inv
#   spec2: rel_bid ~ rel_belief_inv + Woman + Woman x rel_belief_inv
#   spec3: rel_bid ~ rel_belief_inv + rel_belief_price                      (spec1 + price beliefs)
#   spec4: rel_bid ~ rel_belief_inv + rel_belief_price + rel_belief_others_inv
#   spec5: rel_bid ~ rel_belief_inv + rel_belief_price + rel_belief_others_inv + rel_belief_others_price
#
# All five include group-level fixed effects (factor(group_uid)), which
# are fit but omitted from the printed table (note added instead).
#
# "rel_X" = X(exp) - X(no), i.e. a LEVEL difference (no logs) for the main
# table -- this fixes the earlier version of 03_regressions.R, which took
# log(exp) - log(no) for bid/belief_inv/belief_price.
#
# Three tables are produced, all sharing the same 5-spec structure:
#   1. output/tables/spec1_levels_all_rounds.tex   - levels, bid averaged over rounds 1-10
#   2. output/tables/spec1_levels_last5.tex        - levels, bid averaged over rounds 6-10 only
#   3. output/tables/spec1_logdiff_all_rounds.tex  - log differences (log(exp)-log(no)) for
#                                                     bid and all belief variables, rounds 1-10
# ============================================================

repo_root <- {
  cand <- c(getwd(), "/home/claude/managers_work", "/home/claude/managers_v6",
            "/home/claude/managers_v5", "/home/claude/managers")
  cand[which(sapply(cand, function(x) dir.exists(file.path(x, "data"))))[1]]
}
source(file.path(repo_root, "codes", "01_load_data.R"))

out_dir <- file.path(repo_root, "output", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

covars <- unique(dat[, c("participant.code", "session.code", "player.gender", "player.treatment")])
belief_bases_en <- c("player.belief_inv", "player.belief_price",
                      "player.belief_others_inv", "player.belief_others_price")

# ------------------------------------------------------------------
# Build relative (exp - no) belief differences -- single-shot, round 1.
# use_log = TRUE  -> log(exp) - log(no)  (0s recoded to 1 before logging)
# use_log = FALSE -> exp - no            (levels)
# ------------------------------------------------------------------
build_beliefs <- function(use_log) {
  bel <- belief_exp_no_long[, c("participant.code", "manager", belief_bases_en)]
  if (use_log) {
    for (v in belief_bases_en) bel[[v]][bel[[v]] == 0] <- 1
  }
  w <- reshape(bel, idvar = "participant.code", timevar = "manager", direction = "wide")
  out <- data.frame(participant.code = w$participant.code)
  for (v in belief_bases_en) {
    expcol <- paste0(v, ".exp"); nocol <- paste0(v, ".no")
    nm <- paste0("rel_", sub("^player\\.", "", v))
    if (use_log) {
      out[[nm]] <- log(w[[expcol]]) - log(w[[nocol]])
    } else {
      out[[nm]] <- w[[expcol]] - w[[nocol]]
    }
  }
  out
}

# ------------------------------------------------------------------
# Build relative (exp - no) bid, averaged over the given round set first
# (one number per participant per manager), then differenced.
# ------------------------------------------------------------------
build_bid <- function(use_log, rounds) {
  bd <- bidtrade_exp_no_long[bidtrade_exp_no_long$subsession.round_number %in% rounds, ]
  bp <- aggregate(player.bid ~ participant.code + manager + group_uid, data = bd, FUN = mean)
  bw <- reshape(bp, idvar = c("participant.code", "group_uid"), timevar = "manager", direction = "wide")
  names(bw) <- c("participant.code", "group_uid", "player.bid.no", "player.bid.exp")
  if (use_log) {
    bw$rel_bid <- log(bw$player.bid.exp) - log(bw$player.bid.no)
  } else {
    bw$rel_bid <- bw$player.bid.exp - bw$player.bid.no
  }
  bw[, c("participant.code", "group_uid", "rel_bid")]
}

build_spec1_dat <- function(use_log, rounds) {
  bid <- build_bid(use_log, rounds)
  bel <- build_beliefs(use_log)
  d <- merge(bid, bel, by = "participant.code")
  d <- merge(d, covars, by = "participant.code")
  d$player.gender <- factor(d$player.gender, levels = c("Man", "Woman"))
  d$group_uid <- factor(d$group_uid)
  needed <- c("rel_bid", "rel_belief_inv", "rel_belief_price",
              "rel_belief_others_inv", "rel_belief_others_price", "player.gender")
  d <- d[complete.cases(d[, needed]), ]
  d$Woman <- as.integer(d$player.gender == "Woman")
  d$Woman_x_inv <- d$Woman * d$rel_belief_inv
  d
}

fit_specs <- function(d) {
  list(
    spec1 = lm(rel_bid ~ rel_belief_inv + factor(group_uid), data = d),
    spec2 = lm(rel_bid ~ rel_belief_inv + Woman + Woman_x_inv + factor(group_uid), data = d),
    spec3 = lm(rel_bid ~ rel_belief_inv + rel_belief_price + factor(group_uid), data = d),
    spec4 = lm(rel_bid ~ rel_belief_others_inv + factor(group_uid), data = d),
    spec5 = lm(rel_bid ~ + rel_belief_others_inv + rel_belief_others_price +
                 factor(group_uid), data = d)
  )
}

# ------------------------------------------------------------------
# Generic LaTeX exporter: 5-column table, group-FE dummies omitted.
# ------------------------------------------------------------------
stars <- function(p) if (is.na(p)) "" else if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else ""
num   <- function(x) sprintf("%.3f", x)

export_table <- function(fits, var_labels, dep_label, caption, file_name, fe_prefix = "factor(group_uid)") {
  n_models <- length(fits)
  model_coef <- lapply(fits, function(f) {
    co <- summary(f)$coefficients
    co[!startsWith(rownames(co), fe_prefix), , drop = FALSE]
  })
  model_n     <- sapply(fits, function(f) length(summary(f)$residuals))
  model_r2    <- sapply(fits, function(f) summary(f)$r.squared)
  model_adjr2 <- sapply(fits, function(f) summary(f)$adj.r.squared)

  coef_lines <- character(0)
  for (vn in names(var_labels)) {
    est_cells <- character(n_models); se_cells <- character(n_models)
    for (j in seq_len(n_models)) {
      co <- model_coef[[j]]
      if (vn %in% rownames(co)) {
        est <- co[vn, "Estimate"]; se <- co[vn, "Std. Error"]; p <- co[vn, "Pr(>|t|)"]
        est_cells[j] <- paste0(num(est), stars(p))
        se_cells[j]  <- paste0("(", num(se), ")")
      } else {
        est_cells[j] <- ""; se_cells[j] <- ""
      }
    }
    coef_lines <- c(coef_lines,
                     paste0(var_labels[vn], " & ", paste(est_cells, collapse = " & "), " \\\\"),
                     paste0(" & ", paste(se_cells, collapse = " & "), " \\\\"))
  }

  tex_lines <- c(
    "\\begin{table}[!htbp]",
    "  \\centering",
    sprintf("  \\caption{%s}", caption),
    sprintf("  \\begin{tabular}{l*{%d}{c}}", n_models),
    "    \\toprule",
    sprintf("     & \\multicolumn{%d}{c}{%s} \\\\", n_models, dep_label),
    sprintf("    \\cmidrule(lr){2-%d}", n_models + 1),
    paste0("     & ", paste(sprintf("(%d)", seq_len(n_models)), collapse = " & "), " \\\\"),
    "    \\midrule",
    coef_lines,
    "    \\midrule",
    paste0("    Group fixed effects & ", paste(rep("Yes", n_models), collapse = " & "), " \\\\"),
    paste0("    $N$ & ", paste(model_n, collapse = " & "), " \\\\"),
    paste0("    $R^2$ & ", paste(num(model_r2), collapse = " & "), " \\\\"),
    paste0("    Adj. $R^2$ & ", paste(num(model_adjr2), collapse = " & "), " \\\\"),
    "    \\bottomrule",
    "  \\end{tabular}",
    "  \\begin{tablenotes}",
    "  \\footnotesize",
    "  \\item Note: all specifications include group-level fixed effects (group\\_uid); coefficients omitted from the table. Standard errors in parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.",
    "  \\end{tablenotes}",
    "\\end{table}"
  )
  writeLines(tex_lines, file.path(out_dir, file_name))
  message("Wrote ", file.path(out_dir, file_name))
}

var_labels_levels <- c(
  rel_belief_inv          = "Relative belief: investment (exp $-$ no)",
  Woman                   = "Woman (participant)",
  Woman_x_inv             = "Woman $\\times$ relative belief: investment",
  rel_belief_price        = "Relative belief: price (exp $-$ no)",
  rel_belief_others_inv   = "Relative belief: others' investment, 2nd order (exp $-$ no)",
  rel_belief_others_price = "Relative belief: others' price, 2nd order (exp $-$ no)"
)
var_labels_log <- c(
  rel_belief_inv          = "Relative belief: investment (log exp $-$ log no)",
  Woman                   = "Woman (participant)",
  Woman_x_inv             = "Woman $\\times$ relative belief: investment (log)",
  rel_belief_price        = "Relative belief: price (log exp $-$ log no)",
  rel_belief_others_inv   = "Relative belief: others' investment, 2nd order (log exp $-$ log no)",
  rel_belief_others_price = "Relative belief: others' price, 2nd order (log exp $-$ log no)"
)

# ==================================================================
# Table 1: LEVELS, all rounds 1-10
# ==================================================================
dat1 <- build_spec1_dat(use_log = FALSE, rounds = 1:10)
message("Table 1 (levels, rounds 1-10): N = ", nrow(dat1))
fits1 <- fit_specs(dat1)
export_table(fits1, var_labels_levels,
             dep_label = "Relative bid (Exp $-$ No), rounds 1-10",
             caption   = "Relative bid, EXP $-$ NO: level differences, all rounds (1-10)",
             file_name = "spec1_levels_all_rounds.tex")

# ==================================================================
# Table 2: LEVELS, last 5 rounds only (6-10)
# ==================================================================
dat2 <- build_spec1_dat(use_log = FALSE, rounds = 6:10)
message("Table 2 (levels, rounds 6-10): N = ", nrow(dat2))
fits2 <- fit_specs(dat2)
export_table(fits2, var_labels_levels,
             dep_label = "Relative bid (Exp $-$ No), last 5 rounds (6-10)",
             caption   = "Relative bid, EXP $-$ NO: level differences, last 5 rounds only (6-10)",
             file_name = "spec1_levels_last5.tex")

# ==================================================================
# Table 3: LOG DIFFERENCES, all rounds 1-10
# ==================================================================
dat3 <- build_spec1_dat(use_log = TRUE, rounds = 1:10)
message("Table 3 (log differences, rounds 1-10): N = ", nrow(dat3))
fits3 <- fit_specs(dat3)
export_table(fits3, var_labels_log,
             dep_label = "Relative bid (log Exp $-$ log No), rounds 1-10",
             caption   = "Relative bid, EXP $-$ NO: log differences, all rounds (1-10)",
             file_name = "spec1_logdiff_all_rounds.tex")

message("\nDone.")
