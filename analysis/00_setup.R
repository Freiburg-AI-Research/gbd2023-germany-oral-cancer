# 00_setup.R
# One-time setup: install Nordpred + supporting packages.

cat("\n=== Checking and installing required packages ===\n\n")

cran_pkgs <- c("data.table", "ggplot2", "dplyr", "tidyr", "readr", "scales",
               "patchwork", "remotes")
new_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(new_cran) > 0) {
  cat("Installing CRAN packages:", paste(new_cran, collapse = ", "), "\n")
  install.packages(new_cran, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

# nordpred from GitHub (the original author's fork)
if (!"nordpred" %in% rownames(installed.packages())) {
  cat("\nInstalling nordpred from haraldwf/nordpred GitHub...\n")
  ok <- FALSE
  tryCatch({
    remotes::install_github("haraldwf/nordpred", upgrade = "never", quiet = TRUE)
    ok <- TRUE
  }, error = function(e) {
    cat("haraldwf install failed:", conditionMessage(e), "\n")
  })
  if (!ok) {
    cat("Trying cancercentrum/nordpred fallback...\n")
    tryCatch({
      remotes::install_github("cancercentrum/nordpred", upgrade = "never", quiet = TRUE)
      ok <- TRUE
    }, error = function(e) {
      cat("cancercentrum install failed:", conditionMessage(e), "\n")
    })
  }
}

cat("\n=== Verification ===\n")
ok <- TRUE
for (p in c(cran_pkgs, "nordpred")) {
  available <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-12s %s\n", p, ifelse(available, "OK", "MISSING")))
  if (!available) ok <- FALSE
}

if (ok) {
  cat("\nAll packages installed successfully.\n")
  cat("\n=== nordpred package interface ===\n")
  cat("Exported functions:\n")
  print(ls("package:nordpred"))
  if (exists("nordpred.estimate", where = asNamespace("nordpred"))) {
    cat("\nnordpred.estimate signature:\n")
    print(args(nordpred::nordpred.estimate))
  }
} else {
  cat("\nSome packages failed to install. See messages above.\n")
}
