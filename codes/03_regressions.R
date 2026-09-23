# ============================================================
# 03_regressions.R
#
# Two OLS specifications with RELATIVE BID as the dependent variable,
# both with session-level fixed effects (factor(session.code)).
#
# Unit of observation: one row per PARTICIPANT. Bids are elicited every
# round (1-10 for the abstract exp/no market, 11-20 for the actual
# M/W market), so each participant's bid for a given manager is first
# averaged across their rounds, and the "relative bid" is the
# difference of those participant-level means. This keeps one
# observation per participant (matching the single-shot relative
# belief variable in spec 1) and avoids treating repeated rounds from
# the same participant as independent observations.
#
# Spec 1: relative bid EXP vs NO (rounds 1-10)
#   rel_bid_en = mean(bid | manager=="exp") - mean(bid | manager=="no")
#   rel_belief_inv = belief_inv(exp) - belief_inv(no)      [round-1, single-shot]
#   rel_bid_en ~ player.gender + rel_belief_inv + factor(session.code)
#
# Spec 2: relative bid MAN vs WOMAN (rounds 11-20)
#   rel_bid_mw = mean(bid | manager in {M_exp,M_no}) - mean(bid | manager in {W_exp,W_no})
#   T1 = 1{player.treatment == "M_exp_W_exp"}
#   T2 = 1{player.treatment == "M_exp_W_no"}
#   T3 = 1{player.treatment == "M_no_W_exp"}
#   (T4 = M_no_W_no is the omitted baseline, captured by the intercept)
#   rel_bid_mw ~ T1 + T2 + T3 + player.gender + factor(session.code)
# ============================================================

repo_root <- {
  cand <- c(getwd(), "/home/claude/managers_v5", "/home/claude/managers")
  cand[which(sapply(cand, function(x) dir.exists(file.path(x, "data"))))[1]]
}
source(file.path(repo_root, "codes", "01_load_data.R"))

out_dir <- file.path(repo_root, "output")
report <- character(0)
say <- function(...) { txt <- paste0(...); report <<- c(report, txt); cat(txt, "\n") }
hr <- function(title) { say(""); say(strrep("=", 72)); say(title); say(strrep("=", 72)) }

print_lm <- function(fit, fe_prefix = "factor(session.code)") {
  s <- summary(fit)
  co <- s$coefficients
  keep <- !startsWith(rownames(co), fe_prefix)
  co_show <- co[keep, , drop = FALSE]
  say(paste(capture.output(print(round(co_show, 4))), collapse = "\n"))
  say(sprintf("Session fixed effects: yes (%d session dummies, omitted from display above)",
              sum(startsWith(rownames(co), fe_prefix))))
  say(sprintf("N = %d, R-squared = %.4f, Adj. R-squared = %.4f",
              length(s$residuals), s$r.squared, s$adj.r.squared))
}

# ------------------------------------------------------------------
# SPEC 1: relative bid EXP vs NO
# ------------------------------------------------------------------
hr("SPEC 1 - Relative bid (EXP - NO), rounds 1-10")

bid_en_p <- aggregate(log(player.bid) ~ participant.code + manager + group_uid, data = bidtrade_exp_no_long, FUN = mean)
bid_en_w <- reshape(bid_en_p, idvar = c("participant.code","group_uid"), timevar = "manager", direction = "wide")
names(bid_en_w)<-c("participant.code","group_uid","player.bid.no","player.bid.exp")
bid_en_w$rel_bid_en <- bid_en_w$player.bid.exp - bid_en_w$player.bid.no

belief_inv_p <- belief_exp_no_long[, c("participant.code", "manager", "player.belief_inv","player.belief_price")]
belief_inv_p$player.belief_inv[belief_inv_p$player.belief_inv==0]<-1
belief_inv_w <- reshape(belief_inv_p, idvar = "participant.code", timevar = "manager", direction = "wide")
belief_inv_w$rel_belief_inv <- log(belief_inv_w$player.belief_inv.exp) - log(belief_inv_w$player.belief_inv.no)
belief_inv_w$rel_belief_price <- log(belief_inv_w$player.belief_price.exp) - log(belief_inv_w$player.belief_price.no)

covars <- unique(dat[, c("participant.code", "session.code", "player.gender", "player.treatment")])

spec1_dat <- merge(bid_en_w[, c("participant.code", "group_uid","rel_bid_en")],
                    belief_inv_w[, c("participant.code", "rel_belief_inv","rel_belief_price")],
                    by = "participant.code")
spec1_dat <- merge(spec1_dat, covars, by = "participant.code")
spec1_dat <- spec1_dat[complete.cases(spec1_dat[, c("rel_bid_en", "rel_belief_inv","rel_belief_price", "player.gender", "session.code")]), ]
spec1_dat$player.gender <- factor(spec1_dat$player.gender, levels = c("Man", "Woman"))
spec1_dat$session.code  <- factor(spec1_dat$session.code)
spec1_dat$group_uid  <- factor(spec1_dat$group_uid)

say(sprintf("N = %d participants with non-missing rel_bid_en, rel_belief_inv, gender", nrow(spec1_dat)))
fit1 <- lm(rel_bid_en ~ player.gender + rel_belief_inv + factor(group_uid), data = spec1_dat)
print_lm(fit1)
fit1_1 <- lm(rel_bid_en ~ player.gender + rel_belief_inv + rel_belief_price + factor(group_uid), data = spec1_dat)
print_lm(fit1_1)

# ------------------------------------------------------------------
# SPEC 2: relative bid MAN vs WOMAN, treatment dummies (baseline = T4/M_no_W_no)
# ------------------------------------------------------------------
hr("SPEC 2 - Relative bid (M - W), rounds 11-20, treatment dummies (baseline = M_no_W_no)")

bidtrade_lr_long$mw <- ifelse(bidtrade_lr_long$manager %in% c("M_exp", "M_no"), "M",
                        ifelse(bidtrade_lr_long$manager %in% c("W_exp", "W_no"), "W", NA))
bid_mw_p <- aggregate(log(player.bid) ~ participant.code + group_uid + mw, data = bidtrade_lr_long, FUN = mean)
bid_mw_w <- reshape(bid_mw_p, idvar = c("participant.code","group_uid"), timevar = "mw", direction = "wide")

names(bid_mw_w)<-c("participant.code","group_uid","player.bid.M","player.bid.W")

bid_mw_w$rel_bid_mw <- bid_mw_w$player.bid.M - bid_mw_w$player.bid.W


belief_lr_long$mw <- ifelse(belief_lr_long$manager %in% c("M_exp", "M_no"), "M",
                              ifelse(belief_lr_long$manager %in% c("W_exp", "W_no"), "W", NA))

belief_lr_inv_p <- belief_lr_long[, c("participant.code", "group_uid","mw", "player.belief_inv","player.belief_price")]

belief_lr_inv_p$player.belief_inv[belief_lr_inv_p$player.belief_inv==0]<-1
belief_lr_inv_w <- reshape(belief_lr_inv_p, idvar = "participant.code", timevar = "mw", direction = "wide")
belief_lr_inv_w$rel_belief_inv <- log(belief_lr_inv_w$player.belief_inv.M) - log(belief_lr_inv_w$player.belief_inv.W)
belief_lr_inv_w$rel_belief_price <- log(belief_lr_inv_w$player.belief_price.M) - log(belief_lr_inv_w$player.belief_price.W)


spec2_dat <- merge(bid_mw_w[, c("participant.code", "group_uid", "rel_bid_mw")],
                   belief_lr_inv_w[, c("participant.code", "rel_belief_inv","rel_belief_price")],
                   by = "participant.code")

spec2_dat <- merge(spec2_dat, covars, by = "participant.code")
spec2_dat$T1 <- as.integer(spec2_dat$player.treatment == "M_exp_W_exp")
spec2_dat$T2 <- as.integer(spec2_dat$player.treatment == "M_exp_W_no")
spec2_dat$T3 <- as.integer(spec2_dat$player.treatment == "M_no_W_exp")
# (T4 = M_no_W_no is the omitted baseline)

spec2_dat <- spec2_dat[complete.cases(spec2_dat[, c("rel_bid_mw","rel_belief_inv","rel_belief_price" ,"T1", "T2", "T3", "player.gender", "session.code")]), ]
spec2_dat$player.gender <- factor(spec2_dat$player.gender, levels = c("Man", "Woman"))
spec2_dat$session.code  <- factor(spec2_dat$session.code)

say(sprintf("N = %d participants with non-missing rel_bid_mw, gender", nrow(spec2_dat)))
say(sprintf("Treatment counts: T1(M_exp_W_exp)=%d, T2(M_exp_W_no)=%d, T3(M_no_W_exp)=%d, baseline T4(M_no_W_no)=%d",
            sum(spec2_dat$T1), sum(spec2_dat$T2), sum(spec2_dat$T3),
            sum(spec2_dat$T1 == 0 & spec2_dat$T2 == 0 & spec2_dat$T3 == 0)))

spec2_dat$invT3 <- spec2_dat$rel_belief_inv*spec2_dat$T3
spec2_dat$invT2 <- spec2_dat$rel_belief_inv*spec2_dat$T2
spec2_dat$invT1 <- spec2_dat$rel_belief_inv*spec2_dat$T1

spec2_dat$priceT3 <- spec2_dat$rel_belief_price*spec2_dat$T3
spec2_dat$priceT2 <- spec2_dat$rel_belief_price*spec2_dat$T2
spec2_dat$priceT1 <- spec2_dat$rel_belief_price*spec2_dat$T1



fit2 <- lm(rel_bid_mw ~ player.gender +  rel_belief_inv + factor(group_uid), data = spec2_dat)
fit2_1 <- lm(rel_bid_mw ~ player.gender +  rel_belief_inv + rel_belief_price + factor(group_uid), data = spec2_dat)
fit2_2 <- lm(rel_bid_mw ~ player.gender +  rel_belief_inv + invT1 + invT2 + invT3 + factor(group_uid), data = spec2_dat)
fit2_3 <- lm(rel_bid_mw ~ player.gender +  rel_belief_inv + invT1 + invT2 + invT3 + rel_belief_price + factor(group_uid), data = spec2_dat)
fit2_4 <- lm(rel_bid_mw ~ player.gender +  rel_belief_inv + invT1 + invT2 + invT3 + rel_belief_price + priceT1 + priceT2 + priceT3 +
               factor(group_uid), data = spec2_dat)

print_lm(fit2)
print_lm(fit2_1)
print_lm(fit2_2)
print_lm(fit2_3)
print_lm(fit2_4)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
writeLines(report, file.path(out_dir, "regression_report.txt"))
saveRDS(list(fit1 = fit1, fit2 = fit2), file.path(out_dir, "regression_fits.rds"))
message("\nDone. Report written to ", file.path(out_dir, "regression_report.txt"))
