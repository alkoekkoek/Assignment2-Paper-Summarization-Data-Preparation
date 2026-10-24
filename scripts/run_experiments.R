# =============================================================================
# One-click runner. Steps:
#   1. In RStudio: Session > Set Working Directory > To Source File Location (or setwd() to the scripts folder).
#   2. Choose SESSION below (1, 2 or 3) and press "Source".
#   3. To run in parallel, open a NEW RStudio session (Session > New Session), open this file there,
#      set a different SESSION number and press Source.
# Folder layout expected:  repro/scripts (this file + missforest_repro.R),  repro/data (UCI files),  repro/results (created).
# =============================================================================

# ---------------------------- SETTINGS ---------------------------------------

NSIM    <- 10       # simulations per dataset and missing rate 

RUNS <- list(
  small    = list(datasets = "parkinsons,spect,lymphography", methods = "missforest,knn_variables,mice")
  #,promoter = list(datasets = "promoter", methods = "missforest,knn_variables,mice")
  #,musk_cont_mice = list(datasets = "musk_cont", methods = "mice", nsim = 5)
  #, musk_mixed = list(datasets = "musk_mixed", methods = "missforest,knn_variables", nsim = 3)
)
# -----------------------------------------------------------------------------

here    <- getwd()
root    <- dirname(here)
DATADIR <- file.path(root, "data")

source(file.path(here, "missforest_repro.R"))

for (run_name in names(RUNS)) {
  r <- RUNS[[run_name]]
  opt <- DEFAULTS
  opt$datadir  <- DATADIR
  opt$outdir   <- file.path(root, "results", run_name)
  opt$datasets <- r$datasets
  opt$methods  <- r$methods
  opt$nsim     <- if (is.null(r$nsim)) NSIM else r$nsim
  cat(sprintf("\n########## %s: %s (nsim = %d) ##########\n", run_name, opt$datasets, opt$nsim))
  t0 <- Sys.time()
  tryCatch(run_main(opt), error = function(e) cat("!! Run '", run_name, "' failed: ", conditionMessage(e), "\n", sep = ""))
  cat(sprintf("Run '%s' took %.1f min\n", run_name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\nDone. Results are in: ", file.path(root, "results"), "\n")
