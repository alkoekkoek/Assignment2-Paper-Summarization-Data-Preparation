#!/usr/bin/env Rscript
# =============================================================================
# Reproduction harness for
#   Stekhoven & Buhlmann (2012), "MissForest - non-parametric missing value
#   imputation for mixed-type data", Bioinformatics 28(1):112-118.
#
# missForest itself is NOT re-implemented here: it comes from the authors' R
# package (CRAN `missForest`, GitHub stekhoven/missForest). The package has no
# scripts for the paper's experiments, so this file supplies them:
#   * the experimental protocol (MCAR deletion at 10/20/30 %, repeated runs,
#     NRMSE / PFC computed on the deleted entries only, paired Wilcoxon tests)
#   * the comparison methods: cross-validated KNNimpute (Algorithm 2 of the
#     paper, with its standardisation / {-1,1} dummy coding) and MICE
#   * the out-of-bag (OOB) error check (paper Fig. 4) and the ntree x mtry
#     runtime table (paper Table 2)
#
# Usage:   source() this file, then call run_main(opt) with opt <- DEFAULTS (see run_experiments.R)
# =============================================================================

suppressPackageStartupMessages({
  library(missForest)
  library(mice)
  library(parallel)
})

# ----------------------------------------------------------------------------
# 0. Options (copy DEFAULTS, change what you need, pass it to run_main() / run_table2())
# ----------------------------------------------------------------------------
DEFAULTS <- list(
  mode        = "main",          # main | table2
  datasets    = "paper",         # paper | demo | comma-separated names
  datadir     = "data",          # folder holding the UCI files
  outdir      = "results",
  rates       = "0.1,0.2,0.3",   # fraction of entries deleted (MCAR)
  nsim        = 50,              # simulations per dataset and rate (paper: 50)
  cores       = 1,
  seed        = 1,
  methods     = "auto",          # auto | comma-separated subset of missforest,knn,knn_variables,knn_samples,mice
  # missForest (paper: ntree = 100, mtry = floor(sqrt(p)); package default maxiter = 10)
  ntree       = 100,
  maxiter     = 10,
  mtry        = 0,               # 0 -> floor(sqrt(p))
  backend     = "randomForest",  # paper era; package default today is "ranger"
  # KNNimpute (paper Algorithm 2: k in 1..15 chosen by cross-validation)
  kmax        = 15,
  cv_sets     = 3,               # number of validation sets l (not given in the paper)
  knn_mode    = "variables",     # neighbours are variables (paper) | samples
  # MICE
  mice_maxit  = 5,
  mice_timeout = 600,            # seconds; a timeout is recorded as a failed run
  mice_remove_collinear = FALSE, # TRUE = current mice default: silently skips collinear variables and leaves NAs (e.g. Parkinson's)
  # Table 2 mode
  t2_ntree    = "10,50,100,250,500",
  t2_mtry     = "1,2,4,8,16",
  t2_rate     = 0.1
)

csv_chr <- function(s) trimws(strsplit(as.character(s), ",")[[1]])
csv_num <- function(s) as.numeric(csv_chr(s))

# ----------------------------------------------------------------------------
# 1. Data
# ----------------------------------------------------------------------------
PAPER_DATASETS <- c("parkinsons", "musk_cont", "spect", "promoter", "lymphography", "musk_mixed")
DEMO_DATASETS  <- c("demo_iris_cont", "demo_iris_mixed", "demo_mtcars_cat")

# Dimensions reported in the paper (n x p); used only to warn about wrong files
EXPECTED_DIMS <- list(parkinsons = c(195, 22), musk_cont = c(476, 166), spect = c(267, 22),
                      promoter = c(106, 57), lymphography = c(148, 19),
                      musk_mixed = c(476, 167))

# File (regular expression, searched recursively inside --datadir) needed per dataset
DATA_FILES <- list(parkinsons = "^parkinsons\\.data$", musk_cont = "^clean1\\.data(\\.Z)?$",
                   musk_mixed = "^clean1\\.data(\\.Z)?$", spect = "^SPECT\\.train$",
                   promoter = "^promoters\\.data$",
                   lymphography = "^lymphography\\.data$")

find_data_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
  if (length(f)) f[1] else NA_character_
}

need_file <- function(dir, pattern, what) {
  f <- find_data_file(dir, pattern)
  if (is.na(f)) stop("Cannot find a file matching '", pattern, "' for dataset '", what,
                     "' under '", dir, "'. See README.md for the download links.", call. = FALSE)
  f
}

check_dims <- function(name, x) {
  e <- EXPECTED_DIMS[[name]]
  if (!is.null(e) && !identical(as.integer(dim(x)), as.integer(e)))
    warning(sprintf("%s: loaded %d x %d, but the paper reports %d x %d - check the file/version.",
                    name, nrow(x), ncol(x), e[1], e[2]), call. = FALSE)
  x
}

as_factors <- function(d) as.data.frame(lapply(d, function(v) factor(as.character(v))))

musk1_raw <- function(datadir) {
  # UCI "Musk (Version 1)": clean1.data = molecule name, conformation name, 166 features, class
  f <- need_file(datadir, DATA_FILES$musk_cont, "musk")
  if (grepl("\\.Z$", f)) {                    # UCI ships this file in Unix `compress` format; gzip can unpack it
    out <- sub("\\.Z$", "", f)
    if (!file.exists(out) || file.size(out) == 0)
      suppressWarnings(system2("gzip", c("-dc", shQuote(f)), stdout = out, stderr = FALSE))
    if (!file.exists(out) || file.size(out) == 0)
      stop("Could not unpack '", f, "'. Unpack it manually (7-Zip on Windows, `gzip -d` elsewhere) so that clean1.data sits next to it.", call. = FALSE)
    f <- out
  }
  d <- read.csv(f, header = FALSE, stringsAsFactors = FALSE)
  p <- ncol(d) - 3
  feats <- as.data.frame(lapply(d[, 3:(2 + p), drop = FALSE], as.numeric))
  names(feats) <- paste0("f", seq_len(p))
  list(feats = feats, class = factor(d[[ncol(d)]]))
}

promoter_raw <- function(datadir) {
  # UCI "Promoter Gene Sequences": lines like  +,S10,<tabs>tactagcaat...  (class, name, 57 nucleotides)
  l <- readLines(need_file(datadir, DATA_FILES$promoter, "promoter"), warn = FALSE)
  l <- l[nzchar(trimws(l))]
  parts <- strsplit(l, ",")
  seqs <- vapply(parts, function(p) gsub("[[:space:]]", "", p[3]), "")
  m <- do.call(rbind, strsplit(tolower(seqs), ""))
  d <- as.data.frame(m, stringsAsFactors = FALSE)
  names(d) <- paste0("pos", seq_len(ncol(d)))
  list(x = as.data.frame(lapply(d, function(v) factor(v, levels = c("a", "c", "g", "t")))),
       class = factor(vapply(parts, function(p) trimws(p[1]), "")))
}

load_dataset <- function(name, datadir) {
  x <- switch(name,
    parkinsons = {
      d <- read.csv(need_file(datadir, DATA_FILES$parkinsons, name), stringsAsFactors = FALSE)
      d$name <- NULL; d$status <- NULL            # id and the health-status response are dropped
      as.data.frame(lapply(d, as.numeric))
    },
    musk_cont  = musk1_raw(datadir)$feats,
    musk_mixed = { m <- musk1_raw(datadir); cbind(m$feats, class = m$class) },
    spect = {
      f1 <- need_file(datadir, "^SPECT\\.train$", name)
      f2 <- need_file(datadir, "^SPECT\\.test$", name)
      d <- rbind(read.csv(f1, header = FALSE), read.csv(f2, header = FALSE))
      d <- d[, -1, drop = FALSE]                  # first column is the diagnosis
      names(d) <- paste0("F", seq_len(ncol(d)))
      as.data.frame(lapply(d, function(v) factor(v, levels = c(0, 1))))
    },
    promoter = promoter_raw(datadir)$x,                       # 57 nucleotides, as in the paper's Table 1
    lymphography = {
      d <- read.csv(need_file(datadir, DATA_FILES$lymphography, name), header = FALSE)
      names(d) <- c("class", paste0("attr", seq_len(ncol(d) - 1)))
      as_factors(d)                               # all 19 variables are nominal in the paper
    },
    demo_iris_cont  = iris[, 1:4],
    demo_iris_mixed = iris,
    demo_mtcars_cat = as_factors(mtcars[, c("cyl", "vs", "am", "gear", "carb")]),
    stop("Unknown dataset '", name, "'. Known: ", paste(c(PAPER_DATASETS, DEMO_DATASETS), collapse = ", "))
  )
  x <- check_dims(name, as.data.frame(x))
  isf <- vapply(x, is.factor, TRUE)
  type <- if (!any(isf)) "continuous" else if (all(isf)) "categorical" else "mixed"
  list(name = name, x = x, type = type)
}

resolve_datasets <- function(spec, datadir) {
  if (spec == "demo") return(DEMO_DATASETS)
  if (spec == "paper") {
    have <- PAPER_DATASETS[vapply(PAPER_DATASETS, function(n) !is.na(find_data_file(datadir, DATA_FILES[[n]])), TRUE)]
    if (!length(have)) stop("No UCI files found under '", datadir, "'. See README.md, or try --datasets demo.", call. = FALSE)
    message("Datasets found: ", paste(have, collapse = ", "))
    return(have)
  }
  csv_chr(spec)
}

# MCAR deletion with the authors' own helper; redraw if a variable is (almost) wiped out
make_missing <- function(x, rate, max_tries = 100) {
  for (i in seq_len(max_tries)) {
    xm <- prodNA(x, noNA = rate)
    if (any(is.na(xm)) && all(colSums(!is.na(xm)) >= 2)) return(xm)
  }
  stop("Could not draw a usable missingness pattern at rate ", rate)
}

# ----------------------------------------------------------------------------
# 2. Error measures (paper Section 2) - computed on the deleted entries only
# ----------------------------------------------------------------------------
imputation_error <- function(ximp, xmis, xtrue) {
  if (anyNA(ximp)) return(c(NRMSE = NA_real_, PFC = NA_real_))
  isf <- vapply(xtrue, is.factor, TRUE)
  mis <- is.na(xmis)
  nr <- NA_real_; pf <- NA_real_
  if (any(!isf) && any(mis[, !isf, drop = FALSE])) {
    m  <- mis[, !isf, drop = FALSE]
    xi <- as.matrix(ximp[, !isf, drop = FALSE]); xt <- as.matrix(xtrue[, !isf, drop = FALSE])
    nr <- sqrt(mean((xi[m] - xt[m])^2) / var(xt[m]))       # pooled over all deleted continuous entries
  }
  if (any(isf) && any(mis[, isf, drop = FALSE])) {
    wrong <- 0; tot <- 0
    for (j in which(isf)) {
      mj <- mis[, j]
      if (any(mj)) {
        wrong <- wrong + sum(as.character(ximp[[j]][mj]) != as.character(xtrue[[j]][mj]))
        tot   <- tot + sum(mj)
      }
    }
    pf <- wrong / tot
  }
  c(NRMSE = nr, PFC = pf)
}

# ----------------------------------------------------------------------------
# 3. missForest (authors' package)
# ----------------------------------------------------------------------------
run_missforest <- function(xmis, opt, ntree = opt$ntree, mtry = NULL) {
  if (is.null(mtry)) mtry <- if (opt$mtry > 0) opt$mtry else floor(sqrt(ncol(xmis)))
  log <- character()
  el <- system.time(
    log <- capture.output(
      mf <- missForest(xmis, maxiter = opt$maxiter, ntree = ntree, mtry = mtry,
                       backend = opt$backend, verbose = TRUE)
    )
  )[["elapsed"]]
  oob <- mf$OOBerror
  list(ximp = mf$ximp, time = el,
       oob_NRMSE = if ("NRMSE" %in% names(oob)) oob[["NRMSE"]] else NA_real_,
       oob_PFC   = if ("PFC"   %in% names(oob)) oob[["PFC"]]   else NA_real_,
       n_iter = sum(grepl("in progress", log)))   # includes the final iteration that triggers the stop
}

# ----------------------------------------------------------------------------
# 4. KNNimpute with cross-validated k (paper Algorithm 2 and Section 3)
#    Pre-processing as in the paper: categorical variables -> {-1,1} dummies,
#    all columns centred and scaled to unit sd, imputation, back-transformation.
# ----------------------------------------------------------------------------
mode_impute <- function(x) {                 # mean / most frequent level ("initial imputation")
  for (j in seq_along(x)) {
    v <- x[[j]]; na <- is.na(v)
    if (!any(na)) next
    if (is.factor(v)) { tb <- table(v); v[na] <- names(tb)[which.max(tb)] } else v[na] <- mean(v, na.rm = TRUE)
    x[[j]] <- v
  }
  x
}

encode_pm1 <- function(x) {
  n <- nrow(x); blocks <- vector("list", ncol(x)); info <- vector("list", ncol(x))
  for (j in seq_along(x)) {
    v <- x[[j]]
    if (is.factor(v)) {
      lv <- levels(v)
      B <- vapply(lv, function(l) ifelse(is.na(v), NA_real_, ifelse(v == l, 1, -1)), numeric(n))
      B <- matrix(B, nrow = n, dimnames = list(NULL, paste0(names(x)[j], "=", lv)))
      info[[j]] <- list(factor = TRUE, levels = lv)
    } else {
      B <- matrix(as.numeric(v), nrow = n, dimnames = list(NULL, names(x)[j]))
      info[[j]] <- list(factor = FALSE)
    }
    blocks[[j]] <- B
  }
  list(M = do.call(cbind, blocks), info = info, ncols = vapply(blocks, ncol, 1L))
}

decode_pm1 <- function(M, enc, template) {
  out <- template; start <- 1
  for (j in seq_along(enc$info)) {
    cols <- start:(start + enc$ncols[j] - 1); start <- start + enc$ncols[j]
    if (enc$info[[j]]$factor) {
      idx <- max.col(M[, cols, drop = FALSE], ties.method = "first")
      out[[j]] <- factor(enc$info[[j]]$levels[idx], levels = enc$info[[j]]$levels)
    } else out[[j]] <- M[, cols]
  }
  out
}

standardize <- function(M) {
  mu <- colMeans(M, na.rm = TRUE)
  s  <- apply(M, 2, sd, na.rm = TRUE); s[!is.finite(s) | s == 0] <- 1
  list(Z = sweep(sweep(M, 2, mu, "-"), 2, s, "/"), mu = mu, s = s)
}
unstandardize <- function(Z, mu, s) sweep(sweep(Z, 2, s, "*"), 2, mu, "+")

# Troyanskaya et al. (2001): each row of M is a "unit" (here: a variable). A missing entry (u, f) is
# the inverse-distance-weighted mean of the entries (v, f) of the k nearest units v that are observed
# at f. Distances are Euclidean over the features observed in both units (rescaled by their number).
# Returns one imputed matrix per value in `ks`; distances are computed once for all k.
knn_impute_units <- function(M, ks, eps = 1e-8) {
  obs  <- !is.na(M)
  outs <- replicate(length(ks), M, simplify = FALSE)
  if (all(obs)) return(outs)
  O <- obs * 1; M0 <- M; M0[!obs] <- 0
  shared <- O %*% t(O)
  sq <- (M0^2) %*% t(O)
  D2 <- pmax(sq + t(sq) - 2 * (M0 %*% t(M0)), 0)
  D  <- sqrt(D2 / pmax(shared, 1))
  D[shared == 0] <- Inf
  diag(D) <- Inf
  for (u in which(rowSums(!obs) > 0)) {
    miss <- which(!obs[u, ])
    ord  <- order(D[u, ]); ord <- ord[is.finite(D[u, ord])]
    if (!length(ord)) { for (a in seq_along(ks)) outs[[a]][u, miss] <- 0; next }
    w  <- 1 / (D[u, ord] + eps)
    Ob <- obs[ord, miss, drop = FALSE] * 1
    cs <- Ob;       cs[]  <- apply(Ob, 2, cumsum)               # rank among observed candidates
    cW <- Ob * w;   cWV <- cW * M0[ord, miss, drop = FALSE]
    cW[] <- apply(cW, 2, cumsum); cWV[] <- apply(cWV, 2, cumsum)
    nr <- length(ord); cols <- seq_along(miss)
    for (a in seq_along(ks)) {
      pos <- pmin(colSums(cs < ks[a]) + 1L, nr)                  # position of the k-th observed neighbour
      den <- cW[cbind(pos, cols)]; num <- cWV[cbind(pos, cols)]
      outs[[a]][u, miss] <- ifelse(den > 0, num / den, 0)        # 0 = column mean after centring
    }
  }
  outs
}

knn_on_standardized <- function(Z, ks, mode) {
  if (mode == "variables") lapply(knn_impute_units(t(Z), ks), t) else knn_impute_units(Z, ks)
}

knn_cv_impute <- function(xmis, kmax = 15, cv_sets = 3, mode = "variables", cv_frac = NULL) {
  ks <- seq_len(kmax)
  frac <- if (is.null(cv_frac)) mean(is.na(xmis)) else cv_frac     # delete as many entries as are missing
  xinit <- mode_impute(xmis)                                       # Algorithm 2, step 1
  obs_idx <- which(!is.na(xmis))
  err <- matrix(NA_real_, kmax, cv_sets)
  for (t in seq_len(cv_sets)) {                                    # steps 2-8
    ij  <- arrayInd(sample(obs_idx, max(1, round(frac * length(obs_idx)))), dim(xmis))
    xcv <- xinit
    for (j in unique(ij[, 2])) xcv[[j]][ij[ij[, 2] == j, 1]] <- NA
    st    <- standardize(encode_pm1(xcv)$M)
    Ztrue <- sweep(sweep(encode_pm1(xinit)$M, 2, st$mu, "-"), 2, st$s, "/")
    cells <- is.na(st$Z)
    imps  <- knn_on_standardized(st$Z, ks, mode)
    err[, t] <- vapply(imps, function(Zi) mean((Zi[cells] - Ztrue[cells])^2), 0)
  }
  kbest <- ks[which.min(rowMeans(err))]                            # step 9
  enc <- encode_pm1(xmis); st <- standardize(enc$M)                # step 10
  Zi  <- knn_on_standardized(st$Z, kbest, mode)[[1]]
  ximp <- decode_pm1(unstandardize(Zi, st$mu, st$s), enc, xmis)
  for (j in seq_along(xmis)) { o <- !is.na(xmis[[j]]); ximp[[j]][o] <- xmis[[j]][o] }
  list(ximp = ximp, kbest = kbest)
}

# ----------------------------------------------------------------------------
# 5. MICE (paper: normal-linear for continuous, logistic for binary, polytomous for >2 levels)
# ----------------------------------------------------------------------------
mice_impute <- function(xmis, maxit = 5, timeout = Inf, remove_collinear = FALSE) {
  setTimeLimit(elapsed = timeout, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  imp <- suppressWarnings(mice(xmis, m = 1, maxit = maxit, printFlag = FALSE, remove.collinear = remove_collinear,
                               defaultMethod = c("norm", "logreg", "polyreg", "polr")))
  mice::complete(imp, 1)
}

# ----------------------------------------------------------------------------
# 6. One simulation = one deletion pattern, all methods on the same pattern
# ----------------------------------------------------------------------------
run_method <- function(method, xmis, xtrue, opt) {
  switch(method,
    missforest = run_missforest(xmis, opt),
    knn = , knn_variables = , knn_samples = {     # "knn" uses opt$knn_mode; the other two fix the neighbour type
      mode <- switch(method, knn_variables = "variables", knn_samples = "samples", opt$knn_mode)
      el <- system.time(r <- knn_cv_impute(xmis, opt$kmax, opt$cv_sets, mode))[["elapsed"]]
      list(ximp = r$ximp, time = el, k_best = r$kbest)
    },
    mice = {
      el <- system.time(xi <- mice_impute(xmis, opt$mice_maxit, opt$mice_timeout, opt$mice_remove_collinear))[["elapsed"]]
      list(ximp = xi, time = el)
    },
    stop("unknown method ", method))
}

run_job <- function(ds, rate, sim, seed, methods, opt) {
  set.seed(seed)
  xmis <- make_missing(ds$x, rate)
  rows <- lapply(methods, function(m) {
    r <- tryCatch(run_method(m, xmis, ds$x, opt), error = function(e) e)
    if (inherits(r, "error")) {
      e <- c(NRMSE = NA_real_, PFC = NA_real_); status <- paste("error:", conditionMessage(r)); r <- list(time = NA_real_)
    } else {
      e <- imputation_error(r$ximp, xmis, ds$x)
      status <- if (anyNA(r$ximp)) "incomplete (NAs left)" else "ok"
    }
    g <- function(nm) if (is.null(r[[nm]])) NA_real_ else r[[nm]]
    data.frame(dataset = ds$name, rate = rate, sim = sim, method = m, NRMSE = e[["NRMSE"]], PFC = e[["PFC"]],
               time_s = g("time"), oob_NRMSE = g("oob_NRMSE"), oob_PFC = g("oob_PFC"),
               n_iter = g("n_iter"), k_best = g("k_best"), status = status, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

run_parallel <- function(X, FUN, cores) {
  if (cores > 1 && .Platform$OS.type == "windows") { message("mclapply is unavailable on Windows: using 1 core"); cores <- 1 }
  res <- if (cores > 1) mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE) else lapply(X, FUN)
  bad <- !vapply(res, is.data.frame, TRUE)
  if (any(bad)) warning(sum(bad), " job(s) crashed and were dropped")
  res[!bad]
}

default_methods <- function(type) if (type == "continuous") c("missforest", "knn") else c("missforest", "knn", "mice")

# ----------------------------------------------------------------------------
# 7. Summaries: tables, Wilcoxon tests, OOB check, figure
# ----------------------------------------------------------------------------
se <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

summarise_results <- function(ok) {
  parts <- split(ok, list(ok$dataset, ok$rate, ok$method), drop = TRUE)
  do.call(rbind, lapply(parts, function(d) data.frame(
    dataset = d$dataset[1], rate = d$rate[1], method = d$method[1], n = nrow(d),
    NRMSE = mean(d$NRMSE, na.rm = TRUE), NRMSE_se = se(d$NRMSE),
    PFC = mean(d$PFC, na.rm = TRUE), PFC_se = se(d$PFC),
    time_s = mean(d$time_s, na.rm = TRUE), n_iter = mean(d$n_iter, na.rm = TRUE),
    k_best = mean(d$k_best, na.rm = TRUE), row.names = NULL)))
}

compare_to_missforest <- function(ok) {
  out <- list()
  for (ds in unique(ok$dataset)) for (r in sort(unique(ok$rate[ok$dataset == ds]))) {
    d  <- ok[ok$dataset == ds & ok$rate == r, ]
    mf <- d[d$method == "missforest", ]
    for (m in setdiff(unique(d$method), "missforest")) {
      both <- merge(mf, d[d$method == m, ], by = "sim", suffixes = c(".mf", ".o"))
      for (metric in c("NRMSE", "PFC")) {
        a <- both[[paste0(metric, ".mf")]]; b <- both[[paste0(metric, ".o")]]
        if (length(a) < 2 || all(is.na(a)) || all(is.na(b))) next
        p <- tryCatch(suppressWarnings(wilcox.test(b, a, paired = TRUE, exact = FALSE)$p.value), error = function(e) NA_real_)
        out[[length(out) + 1]] <- data.frame(
          dataset = ds, rate = r, metric = metric, versus = m, n_pairs = length(a),
          mean_missforest = mean(a), mean_other = mean(b),
          reduction_pct = 100 * (mean(b) - mean(a)) / mean(b),
          p_value = p, stars = stars(p), favours = if (mean(a) < mean(b)) "missForest" else m)
      }
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

oob_vs_true <- function(ok) {                # paper Fig. 4: err_true - err_OOB
  mf <- ok[ok$method == "missforest", ]
  out <- list()
  for (ds in unique(mf$dataset)) for (r in sort(unique(mf$rate[mf$dataset == ds]))) for (metric in c("NRMSE", "PFC")) {
    d <- mf[mf$dataset == ds & mf$rate == r, ]
    tr <- d[[metric]]; es <- d[[paste0("oob_", metric)]]
    if (all(is.na(tr)) || all(is.na(es))) next
    out[[length(out) + 1]] <- data.frame(dataset = ds, rate = r, metric = metric, mean_true = mean(tr, na.rm = TRUE),
      mean_oob = mean(es, na.rm = TRUE), mean_signed_diff = mean(tr - es, na.rm = TRUE),
      mean_abs_diff = mean(abs(tr - es), na.rm = TRUE))
  }
  if (length(out)) do.call(rbind, out) else NULL
}

plot_results <- function(smry, file, tests = NULL) {  # grey = KNNimpute, white = MICE, black = missForest (as in the paper)
  # `tests` (from compare_to_missforest) adds significance stars above the competitor's bar:
  # paired Wilcoxon test vs. missForest. Black stars = missForest better, red = the competitor better.
  panels <- do.call(rbind, lapply(c("NRMSE", "PFC"), function(mt) {
    ds <- unique(smry$dataset[!is.na(smry[[mt]])])
    if (length(ds)) data.frame(dataset = ds, metric = mt, stringsAsFactors = FALSE)
  }))
  if (is.null(panels) || !nrow(panels)) return(invisible(NULL))
  cols <- c(knn = "grey70", knn_variables = "grey70", knn_samples = "grey40", mice = "white", missforest = "black")
  nr <- ceiling(nrow(panels) / 3)
  png(file, width = 1600, height = 520 * nr, res = 150)
  par(mfrow = c(nr, 3), mar = c(4, 4, 3, 1))
  for (i in seq_len(nrow(panels))) {
    d <- smry[smry$dataset == panels$dataset[i], ]
    meths <- intersect(c("knn", "knn_variables", "knn_samples", "mice", "missforest"), unique(d$method))
    rates <- sort(unique(d$rate))
    mat <- sapply(rates, function(r) sapply(meths, function(m) {
      v <- d[[panels$metric[i]]][d$rate == r & d$method == m]; if (length(v)) v else NA_real_ }))
    mat <- matrix(mat, nrow = length(meths), dimnames = list(meths, NULL))
    bp <- barplot(mat, beside = TRUE, col = cols[meths], names.arg = paste0(100 * rates, "%"),
            ylim = c(0, 1.3 * max(mat, na.rm = TRUE)),                  # headroom for the legend
            main = paste(panels$dataset[i], "-", panels$metric[i]), ylab = panels$metric[i])
    if (!is.null(tests) && nrow(tests)) for (mi in seq_along(meths)) for (ri in seq_along(rates)) {
      tt <- tests[tests$dataset == panels$dataset[i] & tests$metric == panels$metric[i] &
                  tests$versus == meths[mi] & tests$rate == rates[ri], ]
      if (nrow(tt) == 1 && nzchar(tt$stars) && is.finite(mat[mi, ri]))
        text(bp[mi, ri], mat[mi, ri], tt$stars, pos = 3, offset = 0.2, cex = 1.1,
             col = if (tt$favours == "missForest") "black" else "red", xpd = NA)
    }
    if (!is.null(tests)) mtext("stars: paired Wilcoxon vs missForest (* p<.05, ** p<.01, *** p<.001); red = competitor better",
                               side = 1, line = 2.8, cex = 0.45)
    legend("topleft", legend = meths, fill = cols[meths], bty = "n", cex = 0.8, horiz = TRUE)
  }
  dev.off()
}

write_provenance <- function(opt, outdir) {
  write.csv(data.frame(option = names(opt), value = vapply(opt, function(v) paste(v, collapse = ","), "")),
            file.path(outdir, "parameters_used.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "session_info.txt"))
}

# ----------------------------------------------------------------------------
# 8. Modes
# ----------------------------------------------------------------------------
run_main <- function(opt) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  write_provenance(opt, opt$outdir)
  rates <- csv_num(opt$rates)
  all_res <- list()
  for (name in resolve_datasets(opt$datasets, opt$datadir)) {
    ds <- load_dataset(name, opt$datadir)
    methods <- if (opt$methods == "auto") default_methods(ds$type) else csv_chr(opt$methods)
    message(sprintf("== %s (%s, n=%d, p=%d): %s", ds$name, ds$type, nrow(ds$x), ncol(ds$x), paste(methods, collapse = ", ")))
    for (ri in seq_along(rates)) {
      t0 <- Sys.time()
      jobs <- run_parallel(seq_len(opt$nsim), function(s)
        run_job(ds, rates[ri], s, opt$seed + 100000 * ri + s, methods, opt), opt$cores)
      all_res[[length(all_res) + 1]] <- do.call(rbind, jobs)
      message(sprintf("   rate %.0f%%: %d simulations in %.1f min", 100 * rates[ri], length(jobs),
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      write.csv(do.call(rbind, all_res), file.path(opt$outdir, "results_raw.csv"), row.names = FALSE)  # save as we go
    }
  }
  res <- do.call(rbind, all_res)
  failed <- res[res$status != "ok", ]
  if (nrow(failed)) {
    message("\nNot all runs succeeded. Failed runs per dataset / rate / method / message (column 'sim' = count):")
    print(aggregate(sim ~ dataset + rate + method + status, data = failed, FUN = length), row.names = FALSE)
  }
  ok <- res[res$status == "ok", ]
  smry <- summarise_results(ok); tests <- compare_to_missforest(ok); oob <- oob_vs_true(ok)
  write.csv(smry, file.path(opt$outdir, "results_summary.csv"), row.names = FALSE)
  if (!is.null(tests)) write.csv(tests, file.path(opt$outdir, "results_tests.csv"), row.names = FALSE)
  if (!is.null(oob))   write.csv(oob,   file.path(opt$outdir, "results_oob.csv"),   row.names = FALSE)
  plot_results(smry, file.path(opt$outdir, "results_figure.png"), tests)
  cat("\n--- mean errors (SE in results_summary.csv) ---\n"); print(smry[, c("dataset", "rate", "method", "n", "NRMSE", "PFC", "time_s")], row.names = FALSE, digits = 3)
  if (!is.null(tests)) { cat("\n--- missForest vs. other (paired Wilcoxon) ---\n"); print(tests[, c("dataset", "rate", "metric", "versus", "reduction_pct", "p_value", "stars")], row.names = FALSE, digits = 3) }
  if (!is.null(oob))   { cat("\n--- OOB estimate vs. true error ---\n"); print(oob, row.names = FALSE, digits = 3) }
  cat("\nFiles written to", normalizePath(opt$outdir), "\n")
}

run_table2 <- function(opt) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  write_provenance(opt, opt$outdir)
  for (name in resolve_datasets(opt$datasets, opt$datadir)) {
    ds <- load_dataset(name, opt$datadir)
    mtrys <- csv_num(opt$t2_mtry); mtrys <- mtrys[mtrys <= ncol(ds$x) - 1]
    grid <- expand.grid(mtry = mtrys, ntree = csv_num(opt$t2_ntree))
    one <- function(s) {
      set.seed(opt$seed + s)
      xmis <- make_missing(ds$x, opt$t2_rate)               # same deletion pattern for every setting
      do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
        r <- run_missforest(xmis, opt, ntree = grid$ntree[g], mtry = grid$mtry[g])
        e <- imputation_error(r$ximp, xmis, ds$x)
        data.frame(sim = s, mtry = grid$mtry[g], ntree = grid$ntree[g], NRMSE = e[["NRMSE"]], PFC = e[["PFC"]], time_s = r$time)
      }))
    }
    res <- do.call(rbind, run_parallel(seq_len(opt$nsim), one, opt$cores))
    agg <- do.call(rbind, lapply(split(res, list(res$mtry, res$ntree)), function(d) data.frame(
      mtry = d$mtry[1], ntree = d$ntree[1], NRMSE_pct = 100 * mean(d$NRMSE, na.rm = TRUE),
      PFC_pct = 100 * mean(d$PFC, na.rm = TRUE), time_s = mean(d$time_s), row.names = NULL)))
    write.csv(agg, file.path(opt$outdir, paste0("table2_", name, ".csv")), row.names = FALSE)
    cat("\nTable-2 style result for", name, "(error in %, runtime in s; compare with paper Table 2)\n")
    print(agg, row.names = FALSE, digits = 3)
  }
}
