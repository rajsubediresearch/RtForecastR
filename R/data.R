#' Weekly measles cases, Ciudad de Mexico (worked example)
#'
#' A two-column weekly case-count series used as the package's worked
#' example. Time is a weekly index, not a calendar date; see
#' `vignette("rtforecastr-walkthrough")` for the recommended
#' `mean_GI`/`var_GI` values to use with this series.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{time}{integer, weekly time index}
#'   \item{cases}{integer, case count for that week}
#' }
#' @source Worked example shipped with RtForecastR; not an official
#'   surveillance release.
"measles_cdmx"

#' Weekly measles cases, Jalisco (larger-outbreak example)
#'
#' A second worked example with a larger peak incidence than
#' [measles_cdmx], useful for exercising the `maxI` prediction-grid
#' argument of [rt_forecast()].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{time}{integer, weekly time index}
#'   \item{cases}{integer, case count for that week}
#' }
#' @source Worked example shipped with RtForecastR; not an official
#'   surveillance release.
"jalisco"
