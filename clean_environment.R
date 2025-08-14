## --- Clean environment and detach all non-base packages ---

# 1. Remove all objects from the workspace
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

# 2. Detach all non-base packages
base_pkgs <- c("package:stats", "package:graphics", "package:grDevices",
               "package:utils", "package:datasets", "package:methods", "package:base")
loaded_pkgs <- search()[grepl("^package:", search())]
for (pkg in setdiff(loaded_pkgs, base_pkgs)) {
  detach(pkg, character.only = TRUE, unload = TRUE)
}

# 3. Clear plots (if in interactive mode)
if (dev.cur() > 1) dev.off()

# 4. Clear console (works in RStudio)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  rstudioapi::sendToConsole("\f", execute = TRUE)
}

# 5. Run garbage collection to free up memory
gc()
