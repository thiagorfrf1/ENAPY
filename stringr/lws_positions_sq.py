lws_positions_sq = """
# Ellipsoidal scaling version
lws.positions.sq <- function(enaset) {
  points = as.matrix(enaset$points)
  weights = as.matrix(enaset$line.weights)
  positions = lws_lsq_positions(weights, points, ncol(points));

  node.positions = positions$nodes;
  rownames(node.positions) = enaset$enadata$codes;

  return(list("node.positions" = node.positions, "centroids" = positions$centroids))
}


lws.positions.sq.R6 <- function(enaset) {
  positions = lws_lsq_positions(enaset$line.weights, enaset$points.rotated, enaset$get("dimensions"));

  node.positions = positions$nodes;
  rownames(node.positions) = enaset$enadata$codes;

  return(list("node.positions" = node.positions, "centroids" = positions$centroids))

}
"""