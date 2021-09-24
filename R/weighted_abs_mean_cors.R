#' Weighted Absolute Mean of Correlations
#'
#' Returns the weighted mean of the absolute value of the correlations (with
#' Fisher's Z transformation applied if `fisher_z_transform` is left as `TRUE`).
#'
#' @param cors_mods List of lists of `mods` and `cors` of the form returned by running
#'   `piececor::piecewise_cors()`.
#' @param fisher_z_transform Logical. Default is `TRUE`. If `FALSE`, simply
#'   calculates the means of `cor` weighted by `n_obs` for each dataframe in
#'   `cors`.
#'
#'   If `TRUE` applies Fisher's z transformation on correlations of
#'   individual segments prior to calculating weighted mean (then applies
#'   inverse).
#'
#' @section Warning: See https://en.wikipedia.org/wiki/Fisher_transformation and
#'   https://www.researchgate.net/post/How-to-calculate-the-weighted-average-of-correlations-per-outcome
#'    for information on Fisher's z transformation. Appropriate methods for
#'   calculating weighted correlations should be further investigated.
#'
#' @return Dataframe with variables of interest in `name` column and weighted
#'   mean correlation (across segments) in `value` for each dataframe in `cors_mods$cors`.
#'
#' @export
#'
#' @examples
#' # See README
weighted_abs_mean_cors <- function(cors_mods,
                                   fisher_z_transform = TRUE){

  # Could improve performance by avoiding weighting schemes when don't have any
  # splits... but doesn't matter as isn't very costly anyways...
  cors <- cors_mods %>%
    purrr::pluck("cors")

  if(fisher_z_transform){
    cors <- cors %>%
      map_dbl(
        ~mutate(.x, cor = atanh(cor)) %>%
          with(weighted.mean(abs(cor), n_obs)) %>%
          tanh()
      )
  } else{
    cors <- cors %>%
      map_dbl(~with(.x, weighted.mean(abs(cor), n_obs)))
  }

  tibble::enframe(cors)

}
