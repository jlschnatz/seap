.onAttach <- function(libname, pkgname) {
  string <- "
           _                __        __    _
     _____(_)___ ___  ___  / /_____ _/ /_  (_)___ ______
    / ___/ / __ `__ \\/ _ \\/ __/ __ `/ __ \\/ / __ `/ ___/
   (__  ) / / / / / /  __/ /_/ /_/ / /_/ / / /_/ (__  )
  /____/_/_/ /_/ /_/\\___/\\__/\\__,_/_.___/_/\\__,_/____/
  "
  rlang::inform(cli::format_inline(string), class = "packageStartupMessage")
  #rlang::inform(cli::cli_text(cli::col_br_magenta("Simulation-Based Estimation and Adjustment of Publication Bias")), class = "packageStartupMessage")
}




