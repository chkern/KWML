#' Example data set
#'
#' A dataset containing a simulated non-probability (n_np = 2000)
#' and probability (n_p = 2000) sample.
#'
#' @format A data frame with 4000 rows and 15 variables:
#' \describe{
#'   \item{w1, ..., w10}{Covariates}
#'   \item{y}{Outcome}
#'   \item{wt}{True weights}
#'   \item{elig_wt}{Observed weights (= 1 for non-probability sample units)}
#'   \item{trt_n}{Sample membership indicator (numeric)}
#'   \item{trt}{Sample membership indicator (factor)}
#' }
"simu_dat"
