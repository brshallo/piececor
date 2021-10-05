fit_mod <- function(data,
                    model_spec,
                    fit_formula){
  model_spec %>%
    parsnip::fit(fit_formula, data = data)
}

signs_same <- function(x1, x2){
  sign(x1) == sign(x2)
}

find_cuts <- function(data,
                      mod,
                      x_name,
                      cor_function){

  data <- select(data, 1, .data[[x_name]])

  deriv <- gratia::fderiv(pluck(mod, "fit"), terms = x_name)

  extrema <-
    bind_cols(tibble(deriv = as.numeric(deriv$derivatives[[1]]$deriv)),
              deriv$eval) %>%
    filter(!signs_same(deriv, lead(deriv)))

  if(nrow(extrema) == 0){

    cor <- tibble(gtoe = -Inf, lt = Inf, data = list(data)) %>%
      mutate(n_obs = map_dbl(data, nrow),
             cor = map_dbl(data, formula(cor_function)))
    return(cor)
  }

  split_points <- tibble(lt = c(extrema[[x_name]], Inf)) %>%
    mutate(gtoe = lag(lt),
           gtoe = ifelse(is.na(gtoe), -Inf, gtoe)) %>%
    relocate(gtoe)

  # Split data by and calculate individual Spearman correlations
  data_splits <- map2(
    split_points$lt,
    split_points$gtoe,
    ~ filter(data, .data[[x_name]] < .x, .data[[x_name]] >= .y)
  )

  cors_splits <- split_points %>%
    bind_cols(tibble(data = data_splits)) %>%
    mutate(n_obs = map_int(data, nrow),
           cor = map_dbl(data, formula(cor_function)))

  cors_splits
}


#' Piecewise Correlations
#'
#' Independently for each variable of interest `...` and `.target`, calculate
#' the piecewise correlation. Observations included for each correlation
#' calculation are segmented based on local extrema (minima / maxima) of the
#' smoothed model fit.
#'
#' If the smoothing model has no local minima / maxima, correlation is
#' calculated across all observations.
#'
#' @param data Dataframe containing columns of interest.
#' @param .target Column name of target (unquoted). Numeric with
#'   sufficient unique observations.
#' @param ... Column names of variables of interest (unquoted), or a tidy
#'   selection specification (see:
#'   https://dplyr.tidyverse.org/reference/select.html). Separate models are
#'   created for each variable passed into `...`. Numeric with sufficient unique
#'   observations.
#' @param custom_model_spec `NULL` or model specification from the `{parsnip}`
#'   package. Default is `NULL` in which case the model specification will be:
#'
#'   ```
#'   parsnip::gen_additive_mod() %>%
#'     parsnip::set_engine("mgcv", method = "REML") %>%
#'     parsnip::set_mode("regression")
#'   ```
#' (At this point) must be a generalized additive model with `{mgcv}` engine.
#' See for details on parameters:
#' https://parsnip.tidymodels.org/reference/details_gen_additive_mod_mgcv.html
#'
#' @param fit_formula Character string. Pseudo representation of formula passed
#'   to `parsnip::fit()` of model specification. Default is `".target ~
#'   s(...)"`. Note that `.target` and `...` must be included in any custom
#'   `fit_formula`.
#' @param cor_function Character string of a lambda function passed to
#'   `purrr::map*()` to calculate correlation between `.target` and variables of
#'   interest at piecewise segments of observations. Default is:
#'   `"~cor(.x, method = 'spearman')[2,1]"`
#'
#' Should be written in tidyverse friendly shortcut lambda notation (like
#' above), as opposed to traditional lambda function notation: `"function(x)
#' cor(x, method = 'spearman')[2,1]"`
#'
#' Must evaluate to a numeric vector of length one. Could pass in functions
#' specification unrelated to correlation.
#'
#' @return Named list of named list of models `$mods` and named list of
#'   dataframes `$cors` for each variable of interest specified in `...`. Each
#'   dataframe in `.cors` has the columns:
#'
#' `gtoe`: (greater than or equal to)
#' `lt`: (less than)
#' `data`: list column containing dataframes of observations within segment
#' `n_obs`: number of rows in `data` segment
#' `cor`: value returned by `cor_function` for segment of observations
#'
#' @export
#'
#' @examples
#' # See README
piecewise_cors <- function(data,
                           .target,
                           ...,
                           custom_model_spec = NULL,
                           fit_formula = ".target ~ s(...)",
                           cor_function = "~cor(.x, method = 'spearman')[2,1]"){

  # model_spec
  if(is.null(custom_model_spec)){

    model_spec <- parsnip::gen_additive_mod() %>%
      parsnip::set_engine("mgcv", method = "REML") %>%
      parsnip::set_mode("regression")

  } else {
    model_spec <- custom_model_spec
  }

  data_select <- select(data, {{.target}}, ...)

  col_names <- names(data_select)

  fit_formulas <- stringr::str_replace(fit_formula, ".target", col_names[1]) %>%
    stringr::str_replace("\\.\\.\\.", col_names[-1]) %>%
    map(formula)

  mods <- map(
    .x = fit_formulas,
    .f = fit_mod,
    data = data_select,
    model_spec = model_spec
  )

  cors_output <- map2(.x = mods,
                      .y = col_names[-1],
                      .f = find_cuts,
                      data = data_select,
                      cor_function = cor_function) %>%
    purrr::set_names(col_names[-1])

  mods_output <- mods %>%
    purrr::set_names(col_names[-1])

  output <- list(mods = mods_output, cors = cors_output)

  output
}
