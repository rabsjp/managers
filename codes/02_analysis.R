# ============================================================
# 02_analysis.R
#
# Compares bids/beliefs/allocations ACROSS MANAGERS.
#
#   - "exp" vs "no"                         (task 1 allocation; round-1
#                                             abstract beliefs/bids; rounds
#                                             1-10 repeated bid/trade market)
#   - M_exp / M_no / W_exp / W_no           (round-11 beliefs about the
#                                             actual manager shown; rounds
#                                             11-20 repeated bid/trade market)
#
# Bids are cleared at the GROUP level every round (group.price_exp,
# group.price_no, group.price_l, group.price_r are shared by every member
# of a group-round), so bid observations within the same group are NOT
# independent. Every bid comparison below is therefore fit as a linear
# mixed-effects model with random intercepts for group (and participant
# nested within group, since each participant also contributes multiple
# rounds), using the pre-installed `nlme` package.
#
# Allocation, beliefs, and trade outcomes are individual, one-shot (or,
# for trade, averaged across a participant's rounds) elicitations with no
# such group-level dependence, so those are compared with plain paired
# t-tests / Wilcoxon signed-rank tests (exp vs no, matched within
# participant) or one-way ANOVA + Kruskal-Wallis with pairwise post-hoc
# (the 4-level manager-identity comparisons).
# ============================================================

library(nlme)
# Base R only -- no dplyr dependency (kept out so the script runs in a
# plain R install with no internet access to CRAN).
repo_root <- {
  cand <- c(getwd(), "/home/claude/managers")
  cand[which(sapply(cand, function(x) dir.exists(file.path(x, "data"))))[1]]
}
source(file.path(repo_root, "codes", "01_load_data.R"))

out_dir  <- file.path(repo_root, "output")
plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

report <- character(0)
say <- function(...) { txt <- paste0(...); report <<- c(report, txt); cat(txt, "\n") }
hr <- function(title) { say(""); say(strrep("=", 72)); say(title); say(strrep("=", 72)) }
fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")

# ------------------------------------------------------------------
# Descriptives by manager level
# ------------------------------------------------------------------
describe_by_manager <- function(df, var, label = var, group_col = "manager") {
  tab <- do.call(rbind, lapply(split(df[[var]], df[[group_col]]), function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(data.frame(n = 0, mean = NA, sd = NA, median = NA, min = NA, max = NA))
    data.frame(n = length(v), mean = mean(v), sd = sd(v), median = median(v), min = min(v), max = max(v))
  }))
  tab <- cbind(setNames(list(rownames(tab)), group_col), tab); rownames(tab) <- NULL
  tab <- tab[tab$n > 0, , drop = FALSE]
  say(sprintf("\n-- %s: descriptives by %s --", label, group_col))
  say(paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"))
  tab
}

# ------------------------------------------------------------------
# Paired within-participant test (exp vs no), no group dependence
# ------------------------------------------------------------------
paired_manager_test <- function(df, var, label = var, id_col = "participant.code") {
  w <- reshape(df[, c(id_col, "manager", var)], idvar = id_col, timevar = "manager", direction = "wide")
  xcol <- paste0(var, ".exp"); ycol <- paste0(var, ".no")
  if (!all(c(xcol, ycol) %in% names(w))) { say(sprintf("\n-- %s: exp/no not both present, skipped --", label)); return(invisible(NULL)) }
  ok <- complete.cases(w[[xcol]], w[[ycol]])
  x <- w[[xcol]][ok]; y <- w[[ycol]][ok]
  say(sprintf("\n-- %s: EXP vs NO (paired, n = %d participants) --", label, length(x)))
  say(sprintf("   mean(exp) = %s, mean(no) = %s, mean diff (exp-no) = %s", fmt(mean(x)), fmt(mean(y)), fmt(mean(x - y))))
  tt <- t.test(x, y, paired = TRUE)
  wt <- tryCatch(wilcox.test(x, y, paired = TRUE), error = function(e) NULL)
  say(sprintf("   paired t-test:       t(%d) = %s, p = %s, 95%% CI = [%s, %s]",
              round(tt$parameter), fmt(tt$statistic), fmt(tt$p.value, 4), fmt(tt$conf.int[1]), fmt(tt$conf.int[2])))
  if (!is.null(wt)) say(sprintf("   Wilcoxon signed-rank: V = %s, p = %s", fmt(wt$statistic, 0), fmt(wt$p.value, 4)))
  invisible(list(t = tt, w = wt))
}

# ------------------------------------------------------------------
# One-sample test against a benchmark (e.g. 50 = midpoint of 0-100 scale)
# ------------------------------------------------------------------
one_sample_vs <- function(x, mu, label) {
  x <- x[!is.na(x)]
  say(sprintf("\n-- %s (n = %d) --", label, length(x)))
  say(sprintf("   mean = %s, sd = %s", fmt(mean(x)), fmt(sd(x))))
  tt <- t.test(x, mu = mu)
  wt <- tryCatch(wilcox.test(x, mu = mu), error = function(e) NULL)
  say(sprintf("   one-sample t-test:   t(%d) = %s, p = %s, 95%% CI = [%s, %s]",
              round(tt$parameter), fmt(tt$statistic), fmt(tt$p.value, 4), fmt(tt$conf.int[1]), fmt(tt$conf.int[2])))
  if (!is.null(wt)) say(sprintf("   Wilcoxon signed-rank: V = %s, p = %s", fmt(wt$statistic, 0), fmt(wt$p.value, 4)))
  invisible(tt)
}

# ------------------------------------------------------------------
# Plain (no clustering) comparison across manager levels: one-way
# ANOVA + Kruskal-Wallis, with pairwise post-hoc when the omnibus test
# is suggestive. Used for beliefs (always) and for trade (per
# instructions -- no group dependence to model).
# ------------------------------------------------------------------
compare_across_manager_plain <- function(df, var, label = var) {
  sub <- df[!is.na(df[[var]]) & !is.na(df$manager), c(var, "manager")]
  names(sub) <- c("y", "manager")
  sub$manager <- droplevels(factor(sub$manager))
  n_groups <- nlevels(sub$manager)
  say(sprintf("\n-- %s across managers (%d manager level(s) present) --", label, n_groups))
  tab <- describe_by_manager(sub, "y", label = label)
  if (n_groups < 2) { say("   only one manager level has data; no comparison possible."); return(invisible(list(tab = tab))) }
  
  aov_fit <- aov(y ~ manager, data = sub)
  aov_s <- summary(aov_fit)[[1]]
  say(sprintf("   One-way ANOVA: F(%d,%d) = %s, p = %s",
              aov_s$Df[1], aov_s$Df[2], fmt(aov_s$`F value`[1]), fmt(aov_s$`Pr(>F)`[1], 4)))
  kw <- kruskal.test(y ~ manager, data = sub)
  say(sprintf("   Kruskal-Wallis: chi-sq(%d) = %s, p = %s", kw$parameter, fmt(kw$statistic), fmt(kw$p.value, 4)))
  
  if (n_groups > 2 && (aov_s$`Pr(>F)`[1] < 0.10 || kw$p.value < 0.10)) {
    say("   -> pairwise comparisons (Welch t-tests, BH-adjusted p-values):")
    pw <- pairwise.t.test(sub$y, sub$manager, p.adjust.method = "BH", pool.sd = FALSE)
    say(paste(capture.output(print(pw)), collapse = "\n"))
  }
  invisible(list(tab = tab, aov = aov_s, kw = kw))
}

# ------------------------------------------------------------------
# Group-clustered comparison across manager levels for BID: linear
# mixed model with random intercepts for group_uid, and participant
# nested within group_uid (a participant contributes several rounds).
# ------------------------------------------------------------------
compare_across_manager_clustered <- function(df, var, label = var) {
  sub <- df[!is.na(df[[var]]) & !is.na(df$manager), c(var, "manager", "group_uid", "participant.code")]
  names(sub) <- c("y", "manager", "group_uid", "pid")
  sub$manager   <- droplevels(factor(sub$manager))
  sub$group_uid <- factor(sub$group_uid)
  sub$pid       <- factor(sub$pid)
  n_groups <- nlevels(sub$manager)
  say(sprintf("\n-- %s across managers, GROUP-CLUSTERED (%d manager level(s), n = %d rounds x participants, %d groups) --",
              label, n_groups, nrow(sub), nlevels(sub$group_uid)))
  tab <- describe_by_manager(sub, "y", label = paste0(label, " (raw, unclustered descriptives)"))
  if (n_groups < 2) { say("   only one manager level has data; no comparison possible."); return(invisible(list(tab = tab))) }
  
  fit <- tryCatch(
    lme(y ~ manager, random = ~1 | group_uid/pid, data = sub, method = "REML",
        control = lmeControl(opt = "optim", msMaxIter = 200)),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    say(sprintf("   mixed model failed to fit (%s); falling back to plain ANOVA (NOT cluster-adjusted).", conditionMessage(fit)))
    return(compare_across_manager_plain(df, var, label))
  }
  a <- anova(fit)
  say(sprintf("   Mixed model (random intercepts: group_uid / participant): F(%d,%d) = %s, p = %s",
              a$numDF[2], a$denDF[2], fmt(a$`F-value`[2]), fmt(a[["p-value"]][2], 4)))
  vc <- VarCorr(fit)
  say(paste("   Variance components:", paste(capture.output(print(vc)), collapse = " | ")))
  
  levs <- levels(sub$manager)
  if (n_groups > 2) {
    say("   -> pairwise contrasts from the mixed model (each row refits with a different reference level):")
    pairs <- combn(levs, 2, simplify = FALSE)
    rows <- list()
    for (p in pairs) {
      sub2 <- sub[sub$manager %in% p, ]
      sub2$manager <- factor(sub2$manager, levels = p)
      sub2$group_uid <- droplevels(sub2$group_uid)
      f2 <- tryCatch(lme(y ~ manager, random = ~1 | group_uid/pid, data = sub2, method = "REML"), error = function(e) NULL)
      if (is.null(f2)) next
      s <- summary(f2)$tTable
      rows[[length(rows) + 1]] <- data.frame(
        contrast = paste0(p[2], " - ", p[1]),
        estimate = s[2, "Value"], se = s[2, "Std.Error"], p_raw = s[2, "p-value"]
      )
    }
    if (length(rows)) {
      pt <- do.call(rbind, rows)
      pt$p_BH <- p.adjust(pt$p_raw, method = "BH")
      pt$estimate <- fmt(pt$estimate); pt$se <- fmt(pt$se); pt$p_raw <- fmt(pt$p_raw, 4); pt$p_BH <- fmt(pt$p_BH, 4)
      say(paste(capture.output(print(pt, row.names = FALSE)), collapse = "\n"))
    }
  } else {
    s <- summary(fit)$tTable
    say(sprintf("   contrast %s - %s: estimate = %s, SE = %s, p = %s",
                levs[2], levs[1], fmt(s[2, "Value"]), fmt(s[2, "Std.Error"]), fmt(s[2, "p-value"], 4)))
  }
  invisible(list(tab = tab, fit = fit, anova = a))
}

# ------------------------------------------------------------------
# GENDER comparisons: does the PARTICIPANT's own gender (player.gender,
# Man/Woman) predict a different allocation / belief / bid, checked
# separately WITHIN each manager level so a participant is never
# compared against themselves and each row is one participant's single
# observation for that manager (no repeated-measure pseudo-replication).
# Plain version (independent two-sample t-test + Wilcoxon rank-sum) for
# allocation / beliefs / trade; clustered version (group_uid/participant
# mixed model) for bid, for the same group-dependence reason as before.
# ------------------------------------------------------------------
gender_by_manager_plain <- function(df, var, label = var) {
  say(sprintf("\n-- %s by participant GENDER, within each manager level --", label))
  levs <- levels(droplevels(factor(df$manager)))
  rows <- list()
  for (l in levs) {
    sub <- df[df$manager == l & !is.na(df[[var]]) & !is.na(df$player.gender), ]
    if (nrow(sub) == 0) next
    x <- sub[[var]][sub$player.gender == "Man"]
    y <- sub[[var]][sub$player.gender == "Woman"]
    if (length(x) < 2 || length(y) < 2) next
    tt <- t.test(x, y)
    wt <- tryCatch(wilcox.test(x, y), error = function(e) NULL)
    rows[[length(rows) + 1]] <- data.frame(
      manager = l, n_man = length(x), n_woman = length(y),
      mean_man = mean(x), mean_woman = mean(y), diff_man_minus_woman = mean(x) - mean(y),
      t = tt$statistic, p_t = tt$p.value, p_wilcox = if (!is.null(wt)) wt$p.value else NA
    )
  }
  if (!length(rows)) { say("   no manager level had both genders present with n>=2."); return(invisible(NULL)) }
  tab <- do.call(rbind, rows)
  disp <- tab
  for (cc in c("mean_man", "mean_woman", "diff_man_minus_woman", "t")) disp[[cc]] <- fmt(disp[[cc]])
  for (cc in c("p_t", "p_wilcox")) disp[[cc]] <- fmt(disp[[cc]], 4)
  say(paste(capture.output(print(disp, row.names = FALSE)), collapse = "\n"))
  invisible(tab)
}

gender_by_manager_clustered <- function(df, var, label = var) {
  say(sprintf("\n-- %s by participant GENDER, GROUP-CLUSTERED, within each manager level --", label))
  levs <- levels(droplevels(factor(df$manager)))
  rows <- list()
  for (l in levs) {
    sub <- df[df$manager == l & !is.na(df[[var]]) & !is.na(df$player.gender), ]
    if (nrow(sub) == 0) next
    sub$gender    <- factor(sub$player.gender, levels = c("Man", "Woman"))
    sub$group_uid <- factor(sub$group_uid)
    sub$pid       <- factor(sub$participant.code)
    if (nlevels(droplevels(sub$gender)) < 2) next
    fit <- tryCatch(
      lme(as.formula(paste(var, "~ gender")), random = ~1 | group_uid/pid, data = sub, method = "REML",
          control = lmeControl(opt = "optim", msMaxIter = 200)),
      error = function(e) e
    )
    if (inherits(fit, "error")) { say(sprintf("   [%s] mixed model failed (%s)", l, conditionMessage(fit))); next }
    s <- summary(fit)$tTable
    rows[[length(rows) + 1]] <- data.frame(
      manager = l, n_obs = nrow(sub),
      n_man = sum(sub$gender == "Man"), n_woman = sum(sub$gender == "Woman"),
      estimate_woman_minus_man = -s[2, "Value"], se = s[2, "Std.Error"], p = s[2, "p-value"]
    )
  }
  if (!length(rows)) { say("   no manager level could be fit."); return(invisible(NULL)) }
  tab <- do.call(rbind, rows)
  disp <- tab
  for (cc in c("estimate_woman_minus_man", "se")) disp[[cc]] <- fmt(disp[[cc]])
  disp$p <- fmt(disp$p, 4)
  say(paste(capture.output(print(disp, row.names = FALSE)), collapse = "\n"))
  invisible(tab)
}

# ==================================================================
# PART 1 - Task-1 allocation (player.allocation), manager in {exp, no}
# ==================================================================
hr("PART 1 - Allocation to the manager (task 1): EXP vs NO")
say("player.allocation is on the 0-100 scale of player.mgr_p1; the two manager")
say("allocations sum to 100 for each participant. A value of 50 = an unbiased 50/50 split.")

exp_share <- alloc_long$player.allocation[alloc_long$manager == "exp"]
one_sample_vs(exp_share, 50, "player.allocation to the EXPERIENCED manager vs. 50 (=.5 share)")
paired_manager_test(alloc_long, "player.allocation", "player.allocation")

exp_MW_ee_share <- alloc_long$player.allocation[alloc_long$manager == "M_exp" & 
                                                  alloc_long$player.treatment == "M_exp_W_exp"]
exp_MW_nn_share <- alloc_long$player.allocation[alloc_long$manager == "M_no" & 
                                                  alloc_long$player.treatment == "M_no_W_no"]

exp_MW_ne_share <- alloc_long$player.allocation[alloc_long$manager == "M_no" & 
                                                  alloc_long$player.treatment == "M_no_W_exp"]

exp_MW_en_share <- alloc_long$player.allocation[alloc_long$manager == "M_exp" & 
                                                  alloc_long$player.treatment == "M_exp_W_no"]

plot_levels <- c("exp_no", "M_exp_W_exp", "M_exp_W_no", "W_exp_M_no", "M_no_W_no")
alloc_long$manager_plot <- with(alloc_long, ifelse(
  manager == "exp", "exp_no",
  ifelse(manager == "W_exp" & player.treatment == "M_no_W_exp", "W_exp_M_no",
         ifelse(manager %in% c("M_exp", "M_no"), as.character(player.treatment), NA_character_))
))
alloc_long$manager_plot <- factor(alloc_long$manager_plot, levels = plot_levels)


one_sample_vs(exp_MW_ee_share, 50, "player.allocation to the M_exp_W_exp manager vs. 50 (=.5 share)")
one_sample_vs(exp_MW_nn_share, 50, "player.allocation to the M_exp_W_exp manager vs. 50 (=.5 share)")
one_sample_vs(exp_MW_ne_share, 50, "player.allocation to the M_exp_W_exp manager vs. 50 (=.5 share)")
one_sample_vs(exp_MW_en_share, 50, "player.allocation to the M_exp_W_exp manager vs. 50 (=.5 share)")

# ==================================================================
# PART 2 - Round-1 beliefs/bid about the abstract EXP vs NO manager
# ==================================================================
hr("PART 2 - Round-1 beliefs: EXP vs NO manager (no group dependence)")
for (b in belief_bases) paired_manager_test(belief_exp_no_long, b, b)

hr("PART 3 - Rounds 1-10 repeated market: BID (group-clustered) and TRADE (plain), EXP vs NO")
compare_across_manager_clustered(bidtrade_exp_no_long, "player.bid", "player.bid (rounds 1-10)")

# trade: average across a participant's 10 rounds first (one obs per
# participant per manager), then a plain paired test -- no group model needed.
trade_agg_en <- aggregate(player.trade ~ participant.code + manager, data = bidtrade_exp_no_long, FUN = mean)
paired_manager_test(trade_agg_en, "player.trade", "player.trade (participant-level mean over rounds 1-10)")

# ==================================================================
# PART 4 - Round-11 beliefs about the ACTUAL manager shown
#           (M_exp / M_no / W_exp / W_no)
# ==================================================================
hr("PART 4 - Round-11 beliefs across the 4 manager identities (no group dependence)")
for (b in belief_bases) compare_across_manager_plain(belief_lr_long, b, b)

hr("PART 5 - Rounds 11-20 repeated market: BID (group-clustered) and TRADE (plain) across manager identities")
compare_across_manager_clustered(bidtrade_lr_long, "player.bid", "player.bid (rounds 11-20)")

trade_agg_lr <- aggregate(player.trade ~ participant.code + manager, data = bidtrade_lr_long, FUN = mean)
compare_across_manager_plain(trade_agg_lr, "player.trade", "player.trade (participant-level mean over rounds 11-20)")

# ==================================================================
# PART 6 - Gender (participant.gender: Man vs Woman) and ALLOCATION
# ==================================================================
hr("PART 6 - Do Man vs Woman participants allocate differently?")
say(sprintf("Sample: %d Man, %d Woman (of %d participants with a survey match).",
            sum(alloc_long$player.gender[!duplicated(alloc_long$participant.code)] == "Man", na.rm = TRUE),
            sum(alloc_long$player.gender[!duplicated(alloc_long$participant.code)] == "Woman", na.rm = TRUE),
            length(unique(alloc_long$participant.code))))
gender_by_manager_plain(alloc_long, "player.allocation", "player.allocation")

# ==================================================================
# PART 7 - Gender and BELIEFS (round 1 exp/no, round 11 identity)
# ==================================================================
hr("PART 7 - Do Man vs Woman participants hold different beliefs?")
say("-- round-1 beliefs about the abstract EXP/NO manager --")
for (b in belief_bases) gender_by_manager_plain(belief_exp_no_long, b, b)
say("\n-- round-11 beliefs about the actual manager shown (M_exp/M_no/W_exp/W_no) --")
for (b in belief_bases) gender_by_manager_plain(belief_lr_long, b, b)

# ==================================================================
# PART 8 - Gender and BID / TRADE (rounds 1-10 exp/no; rounds 11-20 identity)
# ==================================================================
hr("PART 8 - Do Man vs Woman participants bid / trade differently?")
say("-- BID, group-clustered --")
gender_by_manager_clustered(bidtrade_exp_no_long, "player.bid", "player.bid (rounds 1-10, exp/no)")
gender_by_manager_clustered(bidtrade_lr_long,     "player.bid", "player.bid (rounds 11-20, identity)")

say("\n-- TRADE, plain (participant-level mean over rounds; no group dependence) --")
trade_agg_en_g <- merge(trade_agg_en, unique(bidtrade_exp_no_long[, c("participant.code", "player.gender")]), by = "participant.code")
trade_agg_lr_g <- merge(trade_agg_lr, unique(bidtrade_lr_long[, c("participant.code", "player.gender")]), by = "participant.code")
gender_by_manager_plain(trade_agg_en_g, "player.trade", "player.trade (rounds 1-10, exp/no)")
gender_by_manager_plain(trade_agg_lr_g, "player.trade", "player.trade (rounds 11-20, identity)")

# ==================================================================
# PART 9 - Fill in the "Summary of choices" table (paper_oct_2026/results.tex)
#
# Layout: 5 comparisons (1st stage = exp vs no, pooled across all 4
# treatments; T1..T4 = the two managers actually shown in that treatment),
# and for each comparison, one value per manager plus one p-value testing
# whether the two managers' values are equal in that comparison.
#
#   Allocation, Investment (1st/2nd order), Price (1st/2nd order):
#     INDIVIDUAL level. Per participant, compute manager_A - manager_B,
#     test that this difference is 0 (equivalent to a paired t-test).
#
#   Bid: the market clears at the GROUP level each round, so the group
#     dependence has to be accounted for by aggregating to one number per
#     group (mean over participants and rounds within that group) BEFORE
#     testing manager_A - manager_B = 0 across groups. The displayed mean
#     is still the raw individual-level mean bid, for comparability with
#     the rest of the report; only the p-value uses the group aggregate.
#
#   Price: already a group-level variable (group.price_exp/no/l/r are
#     shared by every participant in a group-round); aggregated to one
#     number per group the same way as Bid.
# ==================================================================
hr("PART 9 - Summary-of-choices table (Allocation / Bid / Price / Beliefs)")

table_comparisons <- list(
  list(name = "1st stage", pair = c("exp", "no"),     treatment = NA),
  list(name = "T1",        pair = c("M_exp", "W_exp"), treatment = "M_exp_W_exp"),
  list(name = "T2",        pair = c("M_exp", "W_no"),  treatment = "M_exp_W_no"),
  list(name = "T3",        pair = c("W_exp", "M_no"),  treatment = "M_no_W_exp"),
  list(name = "T4",        pair = c("M_no", "W_no"),   treatment = "M_no_W_no")
)

# individual-level paired test: per-participant difference tested against 0
table_test_individual <- function(df, var, pair, treatment = NA, id_col = "participant.code") {
  if (!is.na(treatment)) df <- df[df$player.treatment == treatment, ]
  df <- df[df$manager %in% pair & !is.na(df[[var]]), ]
  w <- reshape(df[, c(id_col, "manager", var)], idvar = id_col, timevar = "manager", direction = "wide")
  colA <- paste0(var, ".", pair[1]); colB <- paste0(var, ".", pair[2])
  if (!all(c(colA, colB) %in% names(w))) return(list(meanA = NA, meanB = NA, n = 0, p = NA))
  ok <- complete.cases(w[[colA]], w[[colB]])
  a <- w[[colA]][ok]; b <- w[[colB]][ok]
  if (length(a) < 2) return(list(meanA = mean(a), meanB = mean(b), n = length(a), p = NA))
  tt <- t.test(a - b)
  list(meanA = mean(a), meanB = mean(b), n = length(a), p = tt$p.value)
}

# group-level paired test: aggregate to one value per group_uid per manager
# (mean over participants/rounds within the group) THEN test the group-level
# difference against 0 -- this is what "accounts for the group dependence".
table_test_group <- function(df, var, pair, treatment = NA) {
  raw <- df
  if (!is.na(treatment)) raw <- raw[raw$player.treatment == treatment, ]
  raw <- raw[raw$manager %in% pair & !is.na(raw[[var]]), ]
  meanA_raw <- mean(raw[[var]][raw$manager == pair[1]], na.rm = TRUE)
  meanB_raw <- mean(raw[[var]][raw$manager == pair[2]], na.rm = TRUE)
  
  agg <- aggregate(raw[[var]], by = list(group_uid = raw$group_uid, manager = raw$manager), FUN = mean)
  names(agg)[3] <- "val"
  w <- reshape(agg, idvar = "group_uid", timevar = "manager", direction = "wide")
  colA <- paste0("val.", pair[1]); colB <- paste0("val.", pair[2])
  if (!all(c(colA, colB) %in% names(w))) return(list(meanA = meanA_raw, meanB = meanB_raw, n_groups = 0, p = NA))
  ok <- complete.cases(w[[colA]], w[[colB]])
  a <- w[[colA]][ok]; b <- w[[colB]][ok]
  if (length(a) < 2) return(list(meanA = meanA_raw, meanB = meanB_raw, n_groups = length(a), p = NA))
  tt <- t.test(a - b)
  list(meanA = meanA_raw, meanB = meanB_raw, n_groups = length(a), p = tt$p.value)
}



table_rows <- list(
  list(name = "Allocation",              kind = "individual", var = "player.allocation",              ds1 = "alloc_long",           dslr = "alloc_long"),
  list(name = "Bid",                     kind = "group",      var = "player.bid",                      ds1 = "bidtrade_exp_no_long", dslr = "bidtrade_lr_long"),
  list(name = "Price",                   kind = "group",      var = "player.price",                    ds1 = "price_exp_no_long",    dslr = "price_lr_long"),
  list(name = "Investment (1st order)",  kind = "individual", var = "player.belief_inv",               ds1 = "belief_exp_no_long",   dslr = "belief_lr_long"),
  list(name = "Price (1st order)",       kind = "individual", var = "player.belief_price",             ds1 = "belief_exp_no_long",   dslr = "belief_lr_long"),
  list(name = "Investment (2nd order)",  kind = "individual", var = "player.belief_others_inv",        ds1 = "belief_exp_no_long",   dslr = "belief_lr_long"),
  list(name = "Price (2nd order)",       kind = "individual", var = "player.belief_others_price",      ds1 = "belief_exp_no_long",   dslr = "belief_lr_long")
)

table_results <- list()
for (r in table_rows) {
  cells <- list()
  for (cmp in table_comparisons) {
    ds <- get(if (cmp$name == "1st stage") r$ds1 else r$dslr)
    res <- if (r$kind == "individual") {
      table_test_individual(ds, r$var, cmp$pair, cmp$treatment)
    } else {
      table_test_group(ds, r$var, cmp$pair, cmp$treatment)
    }
    cells[[cmp$name]] <- res
  }
  table_results[[r$name]] <- cells
}

# ---- print a plain-text version of the table into the report ----------
say("\nSummary-of-choices table (mean for each manager, and the p-value that")
say("the two managers are equal in that comparison):\n")
for (r in table_rows) {
  say(sprintf("-- %s --", r$name))
  for (cmp in table_comparisons) {
    res <- table_results[[r$name]][[cmp$name]]
    say(sprintf("   %-10s %-6s = %-8s %-6s = %-8s   p = %s",
                cmp$name, cmp$pair[1], fmt(res$meanA), cmp$pair[2], fmt(res$meanB), fmt(res$p, 3)))
  }
}

# ---- write the filled LaTeX table (matches paper_oct_2026/results.tex) --
tex_num  <- function(x) if (is.na(x)) "" else sprintf("%.0f", x)
tex_pval <- function(x) if (is.na(x)) "" else sprintf("%.3f", x)

value_row <- function(row_name, display_name = row_name) {
  vals <- character(0)
  for (cmp in table_comparisons) {
    res <- table_results[[row_name]][[cmp$name]]
    vals <- c(vals, tex_num(res$meanA), tex_num(res$meanB))
  }
  paste0(display_name, " & ", paste(vals, collapse = " & "), " \\\\")
}

pvalue_row <- function(row_name) {
  vals <- vapply(table_comparisons, function(cmp) {
    res <- table_results[[row_name]][[cmp$name]]
    sprintf("\\multicolumn{2}{c|}{%s}", tex_pval(res$p))
  }, character(1))
  paste0("$p$-value & ", paste(vals, collapse = " & "), "\\\\")
}

tex_lines <- c(
  "\\begin{tabular}{lll|ll|ll|ll|ll|}",
  "  \\toprule",
  "  & \\multicolumn{2}{c}{1st stage} &  \\multicolumn{2}{c}{T1} &  \\multicolumn{2}{c}{T2} &  \\multicolumn{2}{c}{T3} &  \\multicolumn{2}{c}{T4}\\\\",
  "  & Exp & No & M\\_exp & W\\_exp & M\\_exp & W\\_no & W\\_exp & M\\_no & M\\_no & W\\_no   \\\\",
  "  \\midrule",
  "\\emph{Choices} & & & & &  & & & & &\\\\",
  value_row("Allocation"),
  pvalue_row("Allocation"),
  value_row("Bid"),
  pvalue_row("Bid"),
  value_row("Price"),
  pvalue_row("Price"),
  "  \\midrule",
  "\\emph{First order beliefs}  & & & & &  & & & & &\\\\",
  value_row("Investment (1st order)", "Investment"),
  pvalue_row("Investment (1st order)"),
  value_row("Price (1st order)", "Price"),
  pvalue_row("Price (1st order)"),
  "  \\midrule",
  "\\emph{Second order beliefs}  & & & & &  & & & & &\\\\",
  value_row("Investment (2nd order)", "Investment"),
  pvalue_row("Investment (2nd order)"),
  value_row("Price (2nd order)", "Price"),
  pvalue_row("Price (2nd order)"),
  "  \\midrule",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(paste(out_dir,"/tables",sep=""), "table_comparison.tex"))
say(sprintf("\nFilled LaTeX table written to %s", file.path(out_dir, "table_comparison.tex")))

# ==================================================================
# Plots
# ==================================================================
png(file.path(plot_dir, "allocation_exp_no.png"), width = 700, height = 550, res = 120)
par(mar = c(4, 4, 3, 1))
boxplot(player.allocation ~ manager_plot, data = alloc_long,
        main = "Allocation to the first manager", ylab = "allocation (0-100)", xlab = "manager first_second",cex.axis = 0.6)
abline(h = 50, lty = 2, col = "red")
dev.off()

png(file.path(plot_dir, "belief_exp_no.png"), width = 1100, height = 800, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (b in belief_bases) boxplot(as.formula(paste(b, "~ manager")), data = belief_exp_no_long, main = sub("^player\\.", "",b), xlab = "", ylab = "")
dev.off()

png(file.path(plot_dir, "belief_lr_by_identity.png"), width = 1100, height = 800, res = 120)
par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))
for (b in belief_bases) boxplot(as.formula(paste(b, "~ manager")), data = belief_lr_long, main = b, xlab = "", ylab = "", las = 2)
dev.off()


png(file.path(plot_dir, "bid_exp_no_by_round.png"), width = 900, height = 600, res = 120)
m <- aggregate(player.bid ~ subsession.round_number + manager, data = bidtrade_exp_no_long, FUN = mean)
plot(m$subsession.round_number[m$manager == "exp"], m$player.bid[m$manager == "exp"], type = "b", col = "steelblue",
     ylim = range(m$player.bid), xlab = "round", ylab = "mean bid", main = "Mean bid by round: EXP vs NO")
lines(m$subsession.round_number[m$manager == "no"], m$player.bid[m$manager == "no"], type = "b", col = "firebrick")
legend("topright", legend = c("exp", "no"), col = c("steelblue", "firebrick"), lty = 1, pch = 1)
dev.off()

png(file.path(plot_dir, "bid_lr_by_identity.png"), width = 900, height = 600, res = 120)
boxplot(player.bid ~ manager, data = bidtrade_lr_long, main = "player.bid (rounds 11-20) by manager identity",
        ylab = "player.bid", xlab = "")
dev.off()

png(file.path(plot_dir, "bid_gender_by_round.png"), width = 900, height = 600, res = 120)
m <- aggregate(player.bid ~ subsession.round_number + manager, data = bidtrade_lr_long, FUN = mean)
plot(m$subsession.round_number[m$manager == "M_exp"], m$player.bid[m$manager == "M_exp"], type = "b", col = "steelblue",
     ylim = range(m$player.bid), xlab = "round", ylab = "mean bid", main = "Mean bid by round")
lines(m$subsession.round_number[m$manager == "M_no"], m$player.bid[m$manager == "M_no"], type = "b", col = "firebrick")
lines(m$subsession.round_number[m$manager == "W_no"], m$player.bid[m$manager == "W_no"], type = "b", col = "firebrick",pch=2)
lines(m$subsession.round_number[m$manager == "W_exp"], m$player.bid[m$manager == "W_exp"], type = "b", col = "steelblue",pch=2)

legend("topright", legend = c("M_exp", "M_no","W_no","W_exp"), col = c("steelblue", "firebrick", "firebrick","steelblue"), lty = 1, pch = c (1,1,2,2))
dev.off()


png(file.path(plot_dir, "allocation_by_gender.png"), width = 900, height = 600, res = 120)
par(mar = c(4, 4, 3, 1))
alloc_expno <- alloc_long[alloc_long$manager %in% c("exp", "no"), ]
alloc_expno$manager <- droplevels(alloc_expno$manager)
alloc_expno$player.gender <- factor(alloc_expno$player.gender)
boxplot(player.allocation ~ player.gender + manager, data = alloc_expno,
        main = "Task-1 allocation by participant gender", ylab = "player.allocation (0-100)", xlab = "gender.manager",
        col = c("lightblue", "lightpink"), las = 2)
abline(h = 50, lty = 2, col = "red")
dev.off()

say(sprintf("\nPlots written to %s", plot_dir))

# ==================================================================
# PART 10 - firm_data.csv: invest.other by gender / experience.
#
# firm_data.csv is an independent data set (real managers' investment
# choices), not part of the portfolio/survey pipeline above, so it is
# read directly here. Each of the 4 bar-pairs below is compared with an
# unpaired two-sample t-test (Welch); a significance bracket connects
# each pair of bars with the p-value printed above the bracket.
# ==================================================================
hr("PART 10 - firm_data.csv (invest.other) by gender / experience")

firm_dat <- read.csv(file.path(repo_root, "data", "firm_data.csv"), stringsAsFactors = FALSE)
say(sprintf("firm_data.csv: %d rows (gender: %s; experience: %s)",
            nrow(firm_dat), paste(table(firm_dat$gender), collapse = "/"), paste(table(firm_dat$experience), collapse = "/")))

ttest_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  t.test(x, y)$p.value
}

bar_specs <- list(
  list(bar_label = "Woman", pair_label = "Gender",      value = mean(firm_dat$invest.other[firm_dat$gender == "W"], na.rm = TRUE)),
  list(bar_label = "Man",   pair_label = "Gender",      value = mean(firm_dat$invest.other[firm_dat$gender == "M"], na.rm = TRUE)),
  list(bar_label = "Exp",   pair_label = "Experience",  value = mean(firm_dat$invest.other[firm_dat$experience == "exp"], na.rm = TRUE)),
  list(bar_label = "No",    pair_label = "Experience",  value = mean(firm_dat$invest.other[firm_dat$experience == "no"], na.rm = TRUE)),
  list(bar_label = "Woman", pair_label = "Exp: gender", value = mean(firm_dat$invest.other[firm_dat$experience == "exp" & firm_dat$gender == "W"], na.rm = TRUE)),
  list(bar_label = "Man",   pair_label = "Exp: gender", value = mean(firm_dat$invest.other[firm_dat$experience == "exp" & firm_dat$gender == "M"], na.rm = TRUE)),
  list(bar_label = "Woman", pair_label = "No Exp: gender",  value = mean(firm_dat$invest.other[firm_dat$experience == "no" & firm_dat$gender == "W"], na.rm = TRUE)),
  list(bar_label = "Man",   pair_label = "No Exp: gender",  value = mean(firm_dat$invest.other[firm_dat$experience == "no" & firm_dat$gender == "M"], na.rm = TRUE))
)
bar_vals   <- sapply(bar_specs, `[[`, "value")
bar_labels <- sapply(bar_specs, `[[`, "bar_label")
pair_names <- sapply(bar_specs, `[[`, "pair_label")[c(1, 3, 5, 7)]

pair_pvals <- c(
  gender        = ttest_p(firm_dat$invest.other[firm_dat$gender == "W"], firm_dat$invest.other[firm_dat$gender == "M"]),
  experience    = ttest_p(firm_dat$invest.other[firm_dat$experience == "exp"], firm_dat$invest.other[firm_dat$experience == "no"]),
  exp_by_gender = ttest_p(firm_dat$invest.other[firm_dat$experience == "exp" & firm_dat$gender == "W"],
                          firm_dat$invest.other[firm_dat$experience == "exp" & firm_dat$gender == "M"]),
  no_by_gender  = ttest_p(firm_dat$invest.other[firm_dat$experience == "no" & firm_dat$gender == "W"],
                          firm_dat$invest.other[firm_dat$experience == "no" & firm_dat$gender == "M"])
)

say("\n-- bar values and pairwise t-test p-values --")
for (i in seq_along(pair_pvals)) {
  b1 <- bar_specs[[2 * i - 1]]; b2 <- bar_specs[[2 * i]]
  say(sprintf("   %-14s %-6s = %-8s %-6s = %-8s   p = %s",
              pair_names[i], b1$bar_label, fmt(b1$value), b2$bar_label, fmt(b2$value), fmt(pair_pvals[i], 4)))
}

png(file.path(plot_dir, "invest_other_gender_experience.png"), width = 1000, height = 700, res = 120)
par(mar = c(6, 4, 3, 1))
bar_cols <- rep(c("lightpink", "lightblue"), 4)
y_max   <- max(bar_vals, na.rm = TRUE)
bp <- barplot(bar_vals, names.arg = bar_labels, col = bar_cols,
              ylim = c(0, y_max * 1.35),
              main = "Investment choice (managers)",
              ylab = "mean")

tick_h <- y_max * 0.02   # length of the small vertical ticks dropping to each bar
for (i in seq_along(pair_pvals)) {
  idx <- c(2 * i - 1, 2 * i)
  x1 <- bp[idx[1]]; x2 <- bp[idx[2]]
  y_bar_top <- max(bar_vals[idx], na.rm = TRUE)
  y_bracket <- y_bar_top + y_max * 0.10   # height of the horizontal bracket bar
  
  # bracket: two vertical ticks (one per bar) joined by a horizontal line
  segments(x1, bar_vals[idx[1]] + tick_h * 0.3, x1, y_bracket)
  segments(x2, bar_vals[idx[2]] + tick_h * 0.3, x2, y_bracket)
  segments(x1, y_bracket, x2, y_bracket)
  
  # p-value text sits just above the bracket, pointed to by the two ticks
  text(mean(c(x1, x2)), y_bracket + y_max * 0.03, labels = sprintf("p = %.3f", pair_pvals[i]), cex = 0.75)
  mtext(pair_names[i], side = 1, line = 4.2, at = mean(c(x1, x2)), cex = 0.7)
}
legend("topleft", legend = c("Woman / Exp", "Man / No"), fill = c("lightpink", "lightblue"), cex = 0.7, bty = "n")
dev.off()

say(sprintf("\ninvest_other_gender_experience.png written to %s", plot_dir))

writeLines(report, file.path(out_dir, "analysis_report.txt"))
message("\nDone. Report written to ", file.path(out_dir, "analysis_report.txt"))
