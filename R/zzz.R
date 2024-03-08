# .onAttach <- function(libname, pkgname) {
#   if (requireNamespace("cli", quietly = TRUE) && requireNamespace("rlang", quietly = TRUE)) {
#     packageStartupMessage(
#       cli::cli({
#         cli::cli_rule()
#         cli::cli_text(cli::col_br_magenta("Simulation-Based Estimation Publication Bias Estimation and Effect Size Correction"))
#         cli::cli_text('Loading {cli::col_br_magenta("speec")} version: {utils::packageVersion("speec")}')
#         cli::cli_rule()
#         })
#     )
#   }
# }



