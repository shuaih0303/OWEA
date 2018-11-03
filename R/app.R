#' Shiny App for \code{design} function
#' 
#' A function to launch graphical interface to \code{design} function.
#' @export
#' @examples  
#' \dontrun{
#' design_app() # lauching the app.}

design_app <- function(){
  shiny::runApp(system.file('shiny', package = 'OWEA'), display.mode = 'normal')
}
