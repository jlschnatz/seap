.onAttach <- function(libname, pkgname) {
  if (requireNamespace("cli", quietly = TRUE) && requireNamespace("rlang", quietly = TRUE)) {
    packageStartupMessage(
      cli::cli({
        cli::cli_rule()
        #rlang::inform(cli::format_inline(string))
        cli::cli_text(cli::col_br_magenta("Simulation-Based Estimation and Adjustment of Publication Bias"))
        cli::cli_text('Loading {cli::col_br_magenta("simetabias")} version: {packageVersion("simetabias")}')
        cli::cli_rule()
        })
    )
  }
}



