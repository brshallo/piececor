plot_splits_setup <- function(cors_mods, show_points){

  mod <- cors_mods[[1]]
  cors_dataframe <- cors_mods[[2]]

  data <- cors_dataframe %>%
    select(data) %>%
    tidyr::unnest(data)

  data <- bind_cols(data, predict(mod, data))

  col_names <- names(data)

  p <- data %>%
    ggplot2::ggplot(aes(x = .data[[col_names[[2]]]], y = .data[[col_names[[1]]]]))

  if(show_points){
    p <- p + ggplot2::geom_point(alpha = 0.3)
  }

  p <- p +
    ggplot2::geom_line(aes(y = .pred), colour = "red") +
    ggplot2::geom_vline(xintercept = cors_dataframe$gtoe[-1],
                        color = "blue",
                        alpha = 0.5)

  if(nrow(cors_dataframe) == 1){
    message(col_names[[2]], ": No splits, i.e. no local extrema in smoothed fit.")
  }

  p

}

#' Plot Splits and Smoothed Model Fit
#'
#' Plot data with cut points for segments and smoothed model overlaid on
#' observations.
#'
#' @param cors_mods List of lists of `mods` and `cors` of the form returned by
#'   running `piececor::piecewise_cors()`.
#' @param ... Column names of variables of interest (unquoted), or a tidy
#'   selection specification (see:
#'   https://dplyr.tidyverse.org/reference/select.html). Only specify for those
#'   results you *wish to plot*.
#' @param show_points Logical. Default is `TRUE`. Whether to include points in
#'   the graph or not. If you have lots of observations, you may want to set to
#'   `FALSE`.
#'
#' @return Plots of smoothing functions and splits between segments of
#'   observations (vertical lines) on observations for each variable specified
#'   in `...` (on x axis, `.target` variable will be on y).
#'
#' @export
#'
#' @examples
#' # See README
plot_splits <- function(cors_mods, ..., show_points = TRUE){

  data_for_plots <- cors_mods %>%
    purrr::transpose() %>%
    as_tibble() %>%
    select(...)

  if(ncol(data_for_plots) == 0) message("No columns selected.")

  map(data_for_plots, plot_splits_setup, show_points = show_points)
}
