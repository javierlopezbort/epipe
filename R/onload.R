.onAttach <- function(libname, pkgname) {
  # Path to _targets.R in the package
  targets_file <- system.file("_targets.R", package = pkgname)

  # Path to the user's working directory
  user_targets_file <- file.path(getwd(), "_targets.R")

  # Copy _targets.R if it doesn't exist in the user's directory
  if (!file.exists(user_targets_file)) {
    file.copy(targets_file, user_targets_file)
    packageStartupMessage(
      "_targets.R file has been copied to your current working directory."
    )
  } else {
    packageStartupMessage("_targets.R file already exists in your working directory.")
  }
}
