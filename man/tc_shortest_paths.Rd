% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cost.R
\name{tc_shortest_paths}
\alias{tc_shortest_paths}
\title{Calculate the shortest paths between points on a raster}
\usage{
tc_shortest_paths(x, from, to, direction = "out")
}
\arguments{
\item{x}{tc_graph: an object of class \code{tc_graph}, as returned by \code{\link[=tc_build_graph]{tc_build_graph()}}}

\item{from}{numeric or matrix: the points to calculate travel costs from. A two-element numeric vector or a 2-column matrix, or anything else accepted by \code{\link[raster:cellFrom]{raster::cellFromXY()}}}

\item{to}{numeric or matrix: the points to calculate travel costs to. If missing, paths to all points in the raster will be returned}

\item{direction}{string:
\itemize{
\item "out" - calculate travel cost outwards from the 'from' points to the 'to' points
\item "in" - calculate travel cost from the 'to' points to the 'from' points
}}
}
\value{
A list of paths (two-column matrices giving x and y points)
}
\description{
Calculate the shortest paths between points on a raster
}
\examples{
my_raster <- raster::raster(ext = raster::extent(c(153, 155, -56, -54)), res = c(0.5, 0.5),
                            crs = "+proj=longlat")
g <- tc_build_graph(my_raster)

## set random edge weights
g <- tc_set_edge_weights(g, values = runif(sum(tc_adj_matrix(g) > 0)))

## calculate shortest path
p <- tc_shortest_paths(g, from = c(153, -54), to = c(155, -56))
plot(p[[1]][, 1], p[[1]][, 2], type = "l", xlab = "", ylab = "")

}
\seealso{
\code{\link[=tc_build_graph]{tc_build_graph()}}
}
