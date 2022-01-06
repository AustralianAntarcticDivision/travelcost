#' Calculate the travel cost from one or more points to each point on a raster
#'
#' @param x tc_graph: an object of class `tc_graph`, as returned by [tc_build_graph()]
#' @param from numeric or matrix: the points to calculate travel costs from. A two-element numeric vector or a 2-column matrix, or anything else accepted by [raster::cellFromXY()]
#' @param direction string:
#' - "out" - calculate travel cost outwards from the 'from' points to each point on the raster
#' - "in" - calculate travel cost from each point on the raster to each of the 'from' points
#' - "both" the sum of "out" and "in"
#'
#' @return A raster layer (if `from` is one point) or brick (if `from` is more than one point)
#'
#' @seealso [tc_build_graph()]
#'
#' @examples
#' my_raster <- raster::raster(ext = raster::extent(c(153, 155, -56, -54)), res = c(0.5, 0.5),
#'                             crs = "+proj=longlat")
#' g <- tc_build_graph(my_raster)
#'
#' ## set random edge weights
#' g <- tc_set_edge_weights(g, values = runif(sum(tc_adj_matrix(g) > 0)))
#'
#' ## calculate travel cost
#' cx <- c(154, -55)
#' D <- tc_cost(g, from = cx)
#'
#' @export
tc_cost <- function(x, from, direction = "out") {
    direction <- match.arg(direction, c("out", "in", "both"))
    g <- igraph::graph_from_adjacency_matrix(tc_adj_matrix(x), mode = if (x$directed) "directed" else "undirected", weighted = TRUE)

    ## indices of "from" points in our ll grid
    from_idx <- raster::cellFromXY(tc_raster(x), from)
    ## nb for locations on land, we probably want to be able to choose the nearest non-land point instead
    ## e.g. from %>% rowwise %>% dplyr::summarize(idx = which.min((ll[, 1] - lon)^2 + (ll[, 2] - lat)^2)) %>% pull(idx)
    if (any(is.na(from_idx))) stop("at least one \"from\" point is outside of the raster extent")

    ## distance of each "from" point to/from each grid point
    if (direction == "both") {
        D <- igraph::distances(g, v = from_idx, mode = "in", algorithm = "bellman-ford") + igraph::distances(g, v = from_idx, mode = "out", algorithm = "bellman-ford")
    } else {
        D <- igraph::distances(g, v = from_idx, mode = direction, algorithm = "bellman-ford")
    }

    ## into raster form
    out <- raster::brick(lapply(seq_len(nrow(D)), function(i) {
        this <- tc_raster(x)
        values(this) <- D[i, ]
        this
    }))
    if (nrow(D) == 1) out[[1]] else out
}
