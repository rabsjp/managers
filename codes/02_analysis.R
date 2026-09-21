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

writeLines(report, file.path(out_dir, "analysis_report.txt"))
message("\nDone. Report written to ", file.path(out_dir, "analysis_report.txt"))
