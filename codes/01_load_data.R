# ============================================================
# 01_load_data.R
#
# Reads every portfolio*.csv file found in the subfolders of data/,
# stacks them into one data frame, drops the fixed set of housekeeping
# variables, and reshapes the manager-specific variables from WIDE
# (one column per manager) to LONG (one row per participant x manager)
# so that bids/beliefs can be compared *across* managers directly.
#
# Five long-format data sets are produced:
#
#   alloc_long            round 1 task-1 allocation, manager in {exp, no}
#                          -> player.allocation
#
#   belief_exp_no_long    round 1 beliefs about an abstract "experienced"
#                          vs "inexperienced" manager, manager in {exp, no}
#                          -> player.belief_inv, player.belief_price,
#                             player.belief_others_inv, player.belief_others_price
#
#   bidtrade_exp_no_long  rounds 1-10 repeated bidding market for the
#                          abstract exp/no manager, manager in {exp, no}
#                          -> player.bid, player.trade
#                          (kept at the round level; bids are cleared at
#                          the GROUP level each round -> group_uid lets
#                          later analysis cluster/model that dependence)
#
#   belief_lr_long        round 11 beliefs about the ACTUAL manager shown
#                          (left/right), manager in {M_exp,M_no,W_exp,W_no}
#                          as implied by player.treatment + group.p2_left_is_mgr_a
#                          -> player.belief_inv, player.belief_price,
#                             player.belief_others_inv, player.belief_others_price
#
#   bidtrade_lr_long      rounds 11-20 repeated bidding market for the
#                          ACTUAL manager shown, manager in
#                          {M_exp,M_no,W_exp,W_no}
#                          -> player.bid, player.trade (round level, group_uid kept)
#
# Source this file to get `dat` plus the five long data frames in the
# workspace, or run it directly to also write them to output/*.csv.
# ============================================================

# Base R only (no external package dependencies required for this file;
# 02_analysis.R uses the pre-installed `nlme` package for the group-
# clustered bid models).

# ---- locate the repo root robustly (works whether sourced or run) -------
find_repo_root <- function() {
  candidates <- c(getwd(), "/home/claude/managers")
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "data"))) return(cand)
  }
  stop("Could not locate the repo root (a folder containing 'data/').")
}
repo_root <- find_repo_root()
data_dir  <- file.path(repo_root, "data")

# ---- 1. find every portfolio**.csv file in every subfolder --------------
portfolio_files <- list.files(
  path = data_dir, pattern = "^portfolio.*\\.csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(portfolio_files) == 0) stop("No portfolio*.csv files found under ", data_dir)

message("Found ", length(portfolio_files), " portfolio file(s):")
message(paste(" -", portfolio_files, collapse = "\n"))

# ---- 2. read & stack, tagging each row with its source subfolder --------
read_one <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$source_folder <- basename(dirname(f))
  d$source_file   <- basename(f)
  d
}
raw_list <- lapply(portfolio_files, read_one)

all_cols <- unique(unlist(lapply(raw_list, names)))
raw_list <- lapply(raw_list, function(d) {
  missing <- setdiff(all_cols, names(d))
  for (m in missing) d[[m]] <- NA
  d[all_cols]
})

dat <- do.call(rbind, raw_list)
rownames(dat) <- NULL

# ---- 3. drop the variables we were told never to include ----------------
drop_vars <- c(
  "participant._is_bot", "participant._index_in_pages", "participant._max_page_index",
  "participant._current_app_name", "participant._current_page_name",
  "participant.time_started_utc", "participant.visited",
  "participant.mturk_worker_id", "participant.mturk_assignment_id",
  "session.mturk_HITId", "session.mturk_HITGroupId",
  "session.comment", "session.is_demo"
)
dat <- dat[, !(names(dat) %in% drop_vars)]
message("Dropped ", sum(drop_vars %in% all_cols), " excluded variables; ",
        ncol(dat), " columns remain (incl. source_folder/source_file).")

# a group identifier that is unique across sessions (group.id_in_subsession
# restarts at 1 in every subfolder, and source_file is just a date stamp
# that repeats across different subfolders run on the same day, so the
# unique key must use source_folder, not source_file)
dat$group_uid <- paste(dat$source_folder, dat$group.id_in_subsession, sep = "__g")

treat_levels <- c("M_no_W_no", "M_no_W_exp", "M_exp_W_no", "M_exp_W_exp")
dat$player.treatment <- factor(dat$player.treatment, levels = treat_levels)

# ==========================================================================
# 4. Manager identity implied by player.treatment + group.p2_left_is_mgr_a
#    (needed for the round-11-onwards "_l"/"_r" variables).
#    e.g. treatment "M_exp_W_no" -> manager A = "M_exp", manager B = "W_no"
#    group.p2_left_is_mgr_a == 1 -> A is on the left, B is on the right
#    group.p2_left_is_mgr_a == 0 -> B is on the left, A is on the right
# ==========================================================================
treat_parts <- strsplit(as.character(dat$player.treatment), "_")
mgr_a <- vapply(treat_parts, function(p) if (length(p) >= 4) paste(p[1], p[2], sep = "_") else NA_character_, character(1))
mgr_b <- vapply(treat_parts, function(p) if (length(p) >= 4) paste(p[3], p[4], sep = "_") else NA_character_, character(1))

dat$left_mgr_label  <- ifelse(dat$group.p2_left_is_mgr_a == 1, mgr_a, mgr_b)
dat$right_mgr_label <- ifelse(dat$group.p2_left_is_mgr_a == 1, mgr_b, mgr_a)

id_cols <- c("participant.code", "player.treatment", "source_file", "source_folder", "group_uid")

# ==========================================================================
# 5. alloc_long : task-1 allocation (round 1), manager in {exp, no}
#
#    player.mgr_p1 = allocation (0-100) to the manager shown on the RIGHT.
#    group.p1_left_is_exp == 0 -> LEFT is inexperienced -> RIGHT (mgr_p1) is
#                                  the experienced manager
#    group.p1_left_is_exp == 1 -> LEFT is experienced    -> RIGHT (mgr_p1) is
#                                  the inexperienced manager
#    The two allocations sum to 100, so the inexperienced-manager share is
#    always 100 - (experienced-manager share).
# ==========================================================================
r1 <- dat[dat$subsession.round_number == 1, ]

alloc_exp <- ifelse(r1$group.p1_left_is_exp == 0, r1$player.mgr_p1, 100 - r1$player.mgr_p1)
alloc_no  <- 100 - alloc_exp

# ==========================================================================
# 5.1 alloc_long : task-2 allocation (round 11), manager in left_mgr_label or right
r11 <- dat[dat$subsession.round_number == 11, ]

alloc_left<- cbind(r11[id_cols], manager = r11$left_mgr_label, player.allocation = 100-r11$player.mgr_p2)
alloc_right<- cbind(r11[id_cols], manager = r11$right_mgr_label, player.allocation = r11$player.mgr_p2)

#alloc_left$manager <- factor(alloc_left$manager, levels = c("M_exp", "M_no", "W_exp", "W_no"))


alloc_long <- rbind(
  cbind(r1[id_cols], manager = "exp", player.allocation = alloc_exp),
  cbind(r1[id_cols], manager = "no",  player.allocation = alloc_no),
  alloc_left,
  alloc_right
)
alloc_long <- alloc_long[!is.na(alloc_long$player.allocation), ]
alloc_long$manager <- factor(alloc_long$manager, levels = c("no", "exp","M_exp", "M_no", "W_exp", "W_no"))




# ==========================================================================
# 6. belief_exp_no_long : round-1 beliefs about the abstract exp/no manager
# ==========================================================================
belief_bases <- c("player.belief_inv", "player.belief_price",
                  "player.belief_others_inv", "player.belief_others_price")

build_exp_no_long <- function(df, bases) {
  exp_part <- df[id_cols]
  no_part  <- df[id_cols]
  for (b in bases) {
    newname <- b
    exp_part[[newname]] <- df[[paste0(b, "_exp")]]
    no_part[[newname]]  <- df[[paste0(b, "_no")]]
  }
  exp_part$manager <- "exp"
  no_part$manager  <- "no"
  out <- rbind(exp_part, no_part)
  out$manager <- factor(out$manager, levels = c("no", "exp"))
  out
}

belief_exp_no_long <- build_exp_no_long(r1, belief_bases)
keep_rows <- rowSums(!is.na(belief_exp_no_long[belief_bases])) > 0
belief_exp_no_long <- belief_exp_no_long[keep_rows, ]

# ==========================================================================
# 7. bidtrade_exp_no_long : rounds 1-10 repeated market, manager in {exp,no}
#    kept at the round level (not aggregated) because bids are cleared at
#    the group level each round -> group_uid + round identify the cluster.
# ==========================================================================
r1_10 <- dat[dat$subsession.round_number %in% 1:10, ]
rt_id_cols <- c(id_cols, "subsession.round_number")

bidtrade_exp_no_long <- rbind(
  cbind(r1_10[rt_id_cols], manager = "exp",
        player.bid = r1_10$player.bid_exp, player.trade = r1_10$player.trade_exp),
  cbind(r1_10[rt_id_cols], manager = "no",
        player.bid = r1_10$player.bid_no,  player.trade = r1_10$player.trade_no)
)
bidtrade_exp_no_long <- bidtrade_exp_no_long[!is.na(bidtrade_exp_no_long$player.bid), ]
bidtrade_exp_no_long$manager <- factor(bidtrade_exp_no_long$manager, levels = c("no", "exp"))

# ==========================================================================
# 8. belief_lr_long : round-11 beliefs about the ACTUAL manager shown,
#    identified as one of M_exp / M_no / W_exp / W_no
# ==========================================================================
r11 <- dat[dat$subsession.round_number == 11, ]

build_lr_long <- function(df, bases) {
  left_part  <- df[id_cols]
  right_part <- df[id_cols]
  for (b in bases) {
    left_part[[b]]  <- df[[paste0(b, "_l")]]
    right_part[[b]] <- df[[paste0(b, "_r")]]
  }
  left_part$manager  <- df$left_mgr_label
  right_part$manager <- df$right_mgr_label
  out <- rbind(left_part, right_part)
  out$manager <- factor(out$manager, levels = c("M_exp", "M_no", "W_exp", "W_no"))
  out
}

belief_lr_long <- build_lr_long(r11, belief_bases)
keep_rows <- rowSums(!is.na(belief_lr_long[belief_bases])) > 0
belief_lr_long <- belief_lr_long[keep_rows, ]

# ==========================================================================
# 9. bidtrade_lr_long : rounds 11-20 repeated market, manager identified as
#    M_exp / M_no / W_exp / W_no. Kept at the round level (group_uid + round
#    identify the market-clearing cluster for the "bid" variable).
# ==========================================================================
r11_20 <- dat[dat$subsession.round_number %in% 11:20, ]

bidtrade_lr_long <- rbind(
  cbind(r11_20[rt_id_cols], manager = r11_20$left_mgr_label,
        player.bid = r11_20$player.bid_l, player.trade = r11_20$player.trade_l),
  cbind(r11_20[rt_id_cols], manager = r11_20$right_mgr_label,
        player.bid = r11_20$player.bid_r, player.trade = r11_20$player.trade_r)
)
bidtrade_lr_long <- bidtrade_lr_long[!is.na(bidtrade_lr_long$player.bid), ]
bidtrade_lr_long$manager <- factor(bidtrade_lr_long$manager, levels = c("M_exp", "M_no", "W_exp", "W_no"))

# ---- 10. save everything ----------------------------------------------
out_dir <- file.path(repo_root, "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  write.csv(dat,                   file.path(out_dir, "master_data.csv"),           row.names = FALSE)
  write.csv(alloc_long,            file.path(out_dir, "alloc_long.csv"),            row.names = FALSE)
  write.csv(belief_exp_no_long,    file.path(out_dir, "belief_exp_no_long.csv"),    row.names = FALSE)
  write.csv(bidtrade_exp_no_long,  file.path(out_dir, "bidtrade_exp_no_long.csv"),  row.names = FALSE)
  write.csv(belief_lr_long,        file.path(out_dir, "belief_lr_long.csv"),        row.names = FALSE)
  write.csv(bidtrade_lr_long,      file.path(out_dir, "bidtrade_lr_long.csv"),      row.names = FALSE)
  message("Wrote master_data.csv and the 5 long-format data sets to ", out_dir)
}
