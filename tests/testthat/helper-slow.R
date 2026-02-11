is_slow_enabled <- function() {
  run_slow <- Sys.getenv("RUN_SLOW_TESTS", unset = "false")
  tolower(run_slow) %in% c("true", "1", "yes")
}
