#' Calculate optimal 2D velocity
#'
#' For movement such as flight or swimming, we often need to take account of winds or currents. These affect travel not only by speeding or slowing movement (e.g. tailwinds or headwinds), but also push the travel off course (winds blowing from the side). [tc_optimal_2d_velocity()] takes a function that estimates the power required by an animal to move at a given speed relative to its surrounding medium (air or water), along with the speed and direction of the wind or current, and estimates the most efficient speed (minimum energy per unit distance travelled) for the animal to travel at.
#'
#' @param vw scalar: the speed of the wind (or water, or other medium)
#' @param phi scalar: the angle (in radians) of the wind with respect to the desired direction of travel
#' @param va numeric: a vector of airspeeds to evaluate
#' @param powerfun function: a function that takes a single parameter (vector of airspeeds) and returns the power required for the animal to move at each of those speeds
#' @param refine logical: if `TRUE`, use the provided `va` values as a first guess and recursively call [tc_optimal_2d_velocity()] to refine the guess. If `FALSE`, return the optimal result from amongst the provided `va` values
#'
#' @return A data.frame with columns `va_opt`, `theta_opt` (the optimal airspeed and angle relative to the desired direction of travel), `vg` (the resultant ground speed in the desired direction of travel), and `E` the energy required (in J/m)
#'
#' @export
tc_optimal_2d_velocity <- function(vw, phi, va, powerfun, refine = TRUE) {
    if (!is.function(powerfun)) powerfun <- match.fun(powerfun)
    chk <- powerfun(1:10)
    if (length(chk) != 10 || !is.numeric(chk)) warning("powerfun should take a vector of numeric wind (or other) values and return a numeric vector of the same length")
    vw_at <- vw * cos(phi) ## along-track wind
    vw_xt <- vw * sin(phi) ## across-track wind
    ## va is the speed of the bird with respect to the air, at some angle theta to intended direction of motion
    ## the cross-track component va * sin(theta) must balance the cross-track component of wind
    suppressWarnings(theta <- asin(-vw_xt / va)) ## will be NA where not possible, i.e. abs(vw_xt / va) > 1
    theta[abs(vw_xt) < 1e-09] <- 0 ## zero cross-track speed
    ## and then the along-track, ground speed in the intended direction of motion is
    vg <- vw_at + va * cos(theta)
    ## energetic cost of the mechanical effort required to propel an animal at speed va, J/s
    E <- powerfun(va) / vg ## J per metre in the direction of required motion
    E[is.na(E) | E < 0] <- Inf
    idx <- which.min(E)
    if (isTRUE(refine)) {
        ## finer va sequence around our first optimum
        va2 <- seq(va[max(1, idx - 1)], va[min(idx + 1, length(va))], length.out = length(va))
        tc_optimal_2d_velocity(vw, phi, va2, powerfun, refine = FALSE)
    } else {
        data.frame(va_opt = va[idx], theta_opt = theta[idx], vg = vg[idx], E = E[idx])
    }
}
