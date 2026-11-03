#' Launch the NTSminer Shiny Interface
#'
#' This function starts the interactive web-based dashboard for NTSminer.
#' It provides a user-friendly interface for project management, PFAS screening, 
#' MS2 purification, and molecular network analysis.
#'
#' @return No return value, opens a web browser.
#' @export
#' @importFrom shiny runApp
launch_ntsminer <- function() {
  app_dir <- system.file("shiny_app", package = "NTSminer")
  if (app_dir == "") {
    stop("Could not find shiny_app directory. Try re-installing `NTSminer`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal", launch.browser = TRUE)
}
