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


# Copied from bhklab/survcomp:
# https://github.com/bhklab/survcomp/blob/master/R/combine.test.R Copied rather
# than add to suggests in to avoid various modeling dependencies that are
# unnecessary and just wanted this function.
combine_test <- function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {

  # Only addition not in survcomp
  if(method == "fisher"){ message("Fisher's method weights all observations equally and ignores `weights`.")}

  if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }

  method <- match.arg(method)
  na.ix <- is.na(p)
  if(any(na.ix) && !na.rm) { stop("missing values are present!") }
  if(all(na.ix)) { return(NA) } ## all p-values are missing
  p <- p[!na.ix]
  k <- length(p)
  if(k == 1) { return(p) }
  if(missing(weight)) { weight <- rep(1, k); }
  switch(method,
         "fisher"={
           cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
         },
         "z.transform"={
           z <- qnorm(p, lower.tail=FALSE)
           cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
         },
         "logit"={
           tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
           cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
         })
  return(cp)
}
