is_slow_enabled <- function() {
  run_slow <- Sys.getenv("RUN_SLOW_TESTS", unset = "false")
  run_giant <- Sys.getenv("RUN_GIANT_TESTS", unset = "false")
  tolower(run_slow) %in% c("true", "1", "yes") ||
    tolower(run_giant) %in% c("true", "1", "yes")
}

is_giant_enabled <- function() {
  run_giant <- Sys.getenv("RUN_GIANT_TESTS", unset = "false")
  tolower(run_giant) %in% c("true", "1", "yes")
}
