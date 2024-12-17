## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mixnhreg, .registration = TRUE
## usethis namespace: end
NULL

#' Data set of 2m surface temperature observations and forecasts for station Dillingen/Donau-Fristingen
#'
#' Data set of 2m surface temperature observations, ensemble mean forecasts (\code{*_mean}), control forecast (\code{*_ctrl}) and ensemble standard
#' deviation (\code{*_sd}) of various weather quantities (*) for station Dillingen/Donau-Fristingen between 2016-2020.
#' The ensemble forecasts are bilinearly interpolated from the closest four surrounding grid points to the coordinates of station Dillingen/Donau-Fristingen.
#'
#' @format An object of class \code{data.frame} with 1826 rows and 34 columns. The data set contains the following quantities:
#' \itemize{
#'  \item \code{date}: valid date of observations and forecasts at 1200 UTC
#'  \item \code{doy}: date as day of the year
#'  \item \code{sin1}: sine-transformed day of the year, i.e. sin(2*pi*\code{doy}/365.25)
#'  \item \code{cos1}: cosine-transformed day of the year, i.e. cos(2*pi*\code{doy}/365.25)
#'  \item \code{name}: station name
#'  \item \code{id}: station id
#'  \item \code{lon}: longitude coordinate of the station
#'  \item \code{lat}: latitude coordinate of the station
#'  \item \code{elev}: station elevation
#'  \item \code{obs}: 2m surface temperature observation
#'  \item \code{temp_*}: 2m surface temperature forecasts
#'  \item \code{pres_*}: surface pressure forecasts
#'  \item \code{u_*}: 10m surface u-wind speed component forecasts
#'  \item \code{v_*}: 10m surface v-wind speed component forecasts
#'  \item \code{wspd_*}: 10m surface wind speed forecasts
#'  \item \code{wgust_*}: 10m surface wind gust forecasts
#'  \item \code{tcc_*}: total cloud cover forecasts
#'  \item \code{sh_*}: specific humidity forecasts
#' }
#'
#' @source
#' Observations:
#' \itemize{
#'  \item Data Source: Deutscher Wetterdienst (DWD), Climate Data Center (CDC).
#'  \item Licence: CC BY 4.0
#'  \item URL: https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/
#'  \item Station: Dillingen/Donau-Fristingen (00983),
#'  \item Time Range: Daily observations at 1200 UTC from 2016-01-02 to 2020-12-31.
#' }
#' Ensemble Forecasts:
#' \itemize{
#'  \item Data Source: European Centre for Medium-Range Weather Forecasts (ECMWF).
#'  \item Licence: CC BY 4.0
#'  \item URL: https://www.ecmwf.int
#' }
#'
#'
#' @docType data
#'
#' @usage data(station)
#'
#' @examples
#' data(station)
"station"
