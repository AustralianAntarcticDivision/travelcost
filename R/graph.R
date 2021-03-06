#' Build the graph that underlies travel cost calculations on a grid
#'
#' [tc_raster()] is a helper function that returns the raster template used to construct the graph. [tc_adj_matrix()] is a helper function that returns the adjacency matrix that underlies the graph.
#'
#' @param x Raster: a raster layer defining the extent, projection, and resolution of the grid
#' @param directed logical: if `TRUE`, the graph is directed (meaning that the cost of travel from cell x to cell y is not necessarily the same as the cost of travel from cell y to cell x)
#' @param neighbours integer: either 4 (each cell is connected to its cardinal neighbours) or 8 (each cell is connected to all 8 neighbouring cells) or 16 (8 neighbouring cells plus the next 8 that do not lie in cardinal or directly diagonal directions)
#' @param wrap_x logical: if `TRUE`, the graph is wrapped in the x direction (appropriate if the grid represents the full 360 longitude span of the globe)
#'
#' @return An object of class `tc_graph`
#'
#' @seealso [tc_set_edge_weights()], [tc_combine_edge_weights()]
#' @examples
#' my_raster <- raster::raster(ext = raster::extent(c(152.75, 155.25, -56.25, -54.25)),
#'                             res = c(0.5, 0.5), crs = "+proj=longlat")
#' g <- tc_build_graph(my_raster)
#'
#' @export
tc_build_graph <- function(x, directed = TRUE, neighbours = 8, wrap_x = FALSE) {
    if (!inherits(x, "Raster")) stop("x must be a Raster object")
    nx <- ncol(x)
    ny <- nrow(x)
    ngrid <- nx * ny
    if (!neighbours %in% c(4, 8, 16)) stop("neighbours must be 4, 8, or 16")
    wrap_x <- isTRUE(wrap_x)
    ## horrible code to build connected graph
    ## row-columnto index, row-filled
    rc2ind <- function(r, c, sz = c(ny, nx)) c + sz[2] * (r - 1)
    ind2rc <- function(ind, sz = c(ny, nx)) {
        temp <- (ind - 1) / sz[2] + 1
        matrix(c(floor(temp), (temp %% 1) * sz[2] + 1), ncol = 2, byrow = FALSE)
    }
    ind2c <- function(ind, sz = c(ny, nx)) ((ind - 1) %% sz[2]) + 1
    ind2r <- function(ind, sz = c(ny, nx)) floor((ind - 1) / sz[2] + 1)

    ## each cell in our grid gets connected to neighbouring cells
    ## expected_sum is the full number of edges for an undirected graph, or half the number of edges for a directed graph
    ## 4 neighbours
    expected_sum <- ny * if (wrap_x) nx else (nx - 1) ## E neighbours
    expected_sum <- expected_sum + nx * (ny - 1) ## S neighbour
    if (neighbours >= 8) {
        expected_sum <- expected_sum + 2 * ((ny - 1) * if (wrap_x) nx else (nx - 1)) ## *2 because SE, SW neighbour counts are the same
    }
    if (neighbours == 16) {
        ## SSW, SSE neighbours: 2 neighbours but not for last 2 rows
        expected_sum <- expected_sum + 2 * (ny - 2) * if (wrap_x) nx else (nx - 1)
        ## EES, WWS neighbours: 2 neighbours but not for last row
        expected_sum <- expected_sum + 2 * (ny - 1) * if (wrap_x) nx else (nx - 2)
    }

    edge_from <- rep(NA_integer_, expected_sum)
    edge_to <- rep(NA_integer_, expected_sum)
    ptr <- 1

    ## eastern neighbour
    nbfrom <- seq_len(ngrid)
    nbto <- nbfrom + 1
    chk <- ind2r(nbfrom) == ind2r(nbto) ## only those on the same row (same lat)
    edge_from[ptr:(ptr + sum(chk) - 1)] <- nbfrom[chk]
    edge_to[ptr:(ptr + sum(chk) - 1)] <- nbto[chk]
    ptr <- ptr + sum(chk)
    if (wrap_x) {
        ## and wrapped from last col to first col
        nbfrom <- rc2ind(seq_len(ny), nx)
        nbto <- rc2ind(seq_len(ny), 1)
        edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
        edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
        ptr <- ptr + length(nbfrom)
    }

    ## southern neighbour (exclude last row)
    nbfrom <- seq_len(ngrid-nx)
    nbto <- nbfrom+nx
    edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
    edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
    ptr <- ptr+length(nbfrom)

    if (neighbours >= 8) {
        ## SW neighbour (exclude last row and first column)
        nbfrom <- setdiff(seq_len(ngrid - nx), rc2ind(seq_len(ny - 1), 1))
        nbto <- nbfrom + nx - 1
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## first col
            nbfrom <- rc2ind(seq_len(ny - 1), 1)
            nbto <- nbfrom + 2 * nx - 1
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }

        ## SE neighbour (exclude last row and last column)
        nbfrom <- setdiff(seq_len(ngrid - nx), rc2ind(seq_len(ny - 1), nx))
        nbto <- nbfrom + nx + 1
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## last col
            nbfrom <- rc2ind(seq_len(ny - 1), nx)
            nbto <- nbfrom + 1
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }
    }

    if (neighbours == 16) {
        ## SSE neighbour (exclude last 2 rows and last column)
        nbfrom <- setdiff(seq_len(ngrid - 2 * nx), rc2ind(seq_len(ny), nx))
        nbto <- nbfrom + 2 * nx + 1
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## SSE last col
            nbfrom <- rc2ind(seq_len(ny - 2), nx)
            nbto <- nbfrom + nx + 1
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }

        ## SSW neighbour (exclude last 2 rows and first column)
        nbfrom <- setdiff(seq_len(ngrid - 2 * nx), rc2ind(seq_len(ny), 1))
        nbto <- nbfrom + 2 * nx - 1
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## SSW first col
            nbfrom <- rc2ind(seq_len(ny - 2), 1)
            nbto <- nbfrom + 3 * nx - 1
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }

        ## EES neighbour (exclude last row and last 2 columns)
        nbfrom <- setdiff(seq_len(ngrid - nx), rc2ind(seq_len(ny), nx))
        nbfrom <- setdiff(nbfrom, rc2ind(seq_len(ny), nx - 1))
        nbto <- nbfrom + nx + 2
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## EES last 2 cols excluding last row element
            nbfrom <- c(rc2ind(seq_len(ny - 1), nx), rc2ind(seq_len(ny - 1), nx - 1))
            nbto <- nbfrom + 2
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }

        ## WWS neighbour (exclude last row and first 2 columns)
        nbfrom <- setdiff(seq_len(ngrid - nx), rc2ind(seq_len(ny), 1))
        nbfrom <- setdiff(nbfrom, rc2ind(seq_len(ny), 2))
        nbto <- nbfrom + nx - 2
        edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
        edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
        ptr <- ptr + length(nbfrom)
        if (wrap_x) {
            ## WWS first 2 cols excluding last row element
            nbfrom <- c(rc2ind(seq_len(ny - 1), 1), rc2ind(seq_len(ny - 1), 2))
            nbto <- nbfrom + 2 * nx - 2
            edge_from[ptr:(ptr + length(nbfrom) - 1)] <- nbfrom
            edge_to[ptr:(ptr + length(nbfrom) - 1)] <- nbto
            ptr <- ptr + length(nbfrom)
        }
    }

    if (ptr != (expected_sum + 1L)) stop("edge count in A is not as expected")
    A <- sparseMatrix(i = edge_from, j = edge_to, x = 1, dims = c(ngrid, ngrid), repr = "T")
    if (isTRUE(directed)) {
        ## directed graph, so we need each edge to be duplicated in its reverse direction
        A <- A + t(A)
        expected_sum <- expected_sum * 2
    }
    if (sum(A) != expected_sum) stop("edge count in A is not as expected")

    ## if (isTRUE(directed)) {
    ##     ## extract the indices of A, which give us the full edge_from and edge_to vectors
    ##     temp <- which(A, arr.ind = TRUE)
    ##     edge_from <- temp[, 1]
    ##     edge_to <- temp[, 2]
    ## }
    structure(list(A = as(A, "dgTMatrix"), directed = directed, template = raster::raster(x)), class = "tc_graph")
}

#' @method plot tc_graph
#' @export
plot.tc_graph <- function(x, y, ...) {
    rgs <- list(...)
    if (!"edge.curved" %in% names(rgs)) rgs$edge.curved <- x$directed
    if (!"layout" %in% names(rgs)) rgs$layout <- coordinates(tc_raster(x))
    ## some plotting weirdness when weights are very large
    scaled_A <- tc_adj_matrix(x) / max(tc_adj_matrix(x), na.rm = TRUE) * 5
    scaled_A[is.na(scaled_A)] <- 0
    ##g <- igraph::graph_from_adjacency_matrix(as(as(scaled_A, "dgCMatrix"), "dgTMatrix"), mode = if (x$directed) "directed" else "undirected")
    g <- igraph::graph_from_adjacency_matrix(scaled_A, mode = if (x$directed) "directed" else "undirected")
    rgs$x <- g
    do.call(plot, rgs)
}

#' @export
#' @rdname tc_build_graph
tc_raster <- function(x) x$template

#' @export
#' @rdname tc_build_graph
tc_adj_matrix <- function(x) x$A

#' Set the edge weights of a travelcost graph
#'
#' @param x tc_graph: an object of class `tc_graph`, as returned by [tc_build_graph()]
#' @param fun function or string: a function or name of function. `fun` should accept two arguments, each an N x 2 matrix. The first gives the xy-coordinates of starting points, and the second gives the xy-coordinates of ending points. `fun` should return the weight (travel cost) between points. Only one of `fun` or `values` is required
#' @param values numeric: A vector of edge weights. The length of `values` must be the same as the number of non-zero elements of the adjacency matrix. Ignored if `fun` is provided
#'
#' @return A `tc_graph` object with the edge weights set
#'
#' @seealso [tc_build_graph()], [tc_combine_edge_weights()]
#'
#' @examples
#' \dontrun{
#'   my_raster <- raster::raster(ext = raster::extent(c(152.75, 155.25, -56.25, -54.25)),
#'                               res = c(0.5, 0.5), crs = "+proj=longlat")
#'   g <- tc_build_graph(my_raster)
#'   ## set weights according to great-circle distance
#'   g <- tc_set_edge_weights(g, fun = geosphere::distHaversine)
#' }
#' @export
tc_set_edge_weights <- function(x, fun, values) {
    if (!missing(values)) {
        x$A@x <- values
        return(x)
    }
    if (!is.function(fun)) fun <- match.fun(fun)
    ll <- coordinates(tc_raster(x))
    ##temp <- which(tc_adj_matrix(x), arr.ind = TRUE) ## col1 is the "from" node, col2 is the to node
    ##w <- fun(ll[temp[, 1], ], ll[temp[, 2], ])
    ##x$A <- as(sparseMatrix(temp[, 1], temp[, 2], x = w, repr = "T"), "dgTMatrix")
    ## more efficient
    x$A@x <- fun(ll[x$A@i + 1, ], ll[x$A@j + 1, ]) ## @i and @j are zero-based
    x
}

#' Combine the edge weights of two or more travelcost graphs
#'
#' @param ... : two or more `tc_graph` objects, or a list of such objects
#' @param fun function or string: the function (or function name) to be used to combine the edge weights. The function must operate on sparse matrices as returned by [Matrix::sparseMatrix()]
#'
#' @return A `tc_graph` object, which is the first object in `...` but with edge weights modified
#'
#' @seealso [tc_build_graph()], [tc_combine_edge_weights()]
#'
#' @examples
#' my_raster <- raster::raster(ext = raster::extent(c(152.75, 155.25, -56.25, -54.25)),
#'                             res = c(0.5, 0.5), crs = "+proj=longlat")
#' g <- tc_build_graph(my_raster)
#' ## set random edge weights
#' g1 <- tc_set_edge_weights(g, values = runif(sum(tc_adj_matrix(g) > 0)))
#' g2 <- tc_set_edge_weights(g, values = runif(sum(tc_adj_matrix(g) > 0)))
#' ## sum of those two
#' g3 <- tc_combine_edge_weights(g1, g2, fun = "+")
#'
#' @export
tc_combine_edge_weights <- function(..., fun = "*") {
    if (!is.function(fun)) fun <- match.fun(fun)
    tcgs <- list(...)
    if (length(tcgs) == 1 && !inherits(tcgs[[1]], "tc_graph") && is.list(tcgs[[1]])) tcgs <- tcgs[[1]]
    if (!all(vapply(tcgs, inherits, "tc_graph", FUN.VALUE = TRUE))) stop("tc_combine_edge_weights should be called with one or more tc_graph objects, or a list of such")
    for (i in seq_along(tcgs)[-1]) tcgs[[1]]$A <- fun(tcgs[[1]]$A, tcgs[[i]]$A)
    tcgs[[1]]
}
