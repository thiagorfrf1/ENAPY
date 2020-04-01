ena_unit_group = """
ena.unit.group = function(set,
  units, method = "mean", scale=T, name = "Group", scaleFactor = 1.0
  # ,keep.dimensions = 1:ncol(set$points.rotated)
) {
  runCIs <- function(pnts) {
    # pntRows = as.matrix(rep(T, nrow(set$points.rotated)))
    # if(length(units.by)>1) {
    #   pntRows = as.data.frame(set$enadata$units[[make.names(units.by[[1]])]]) == x;
    # }
    # pnts = as.matrix(set$points.rotated[pntRows,])
    # dim(pnts) = c(length(which(pntRows)),ncol(set$points.rotated))
    ci = matrix(NA, ncol=ncol(pnts),nrow=2)
    oi = rep(NA, ncol(pnts))
    if(nrow(pnts) > 1) {
      ci = t(apply(pnts,2,function(x) {
        tryCatch(t.test(x, conf.level = 0.95), error = function(e) list(conf.int = rep(x[1],2)))$conf.int
      })) * scaleFactor
      oi = apply(pnts, 2, function(x) { IQR(x) }) * 1.5 * scaleFactor
      # ci = t(matrix(c(
      #   tryCatch(t.test(pnts[, 1], conf.level = 0.95), error = function(e) list(conf.int = c(NA,NA)))$conf.int,
      #   tryCatch(t.test(pnts[, 2], conf.level = 0.95), error = function(e) list(conf.int = c(NA,NA)))$conf.int
      # ), nrow=2)) * scaleFactor;
      # oi = c(IQR(pnts[,1]), IQR(pnts[,2])) * 1.5 * scaleFactor;
    }
    list(
      ci = ci, #[keep.dimensions,],
      oi = oi #[keep.dimensions]
    )
  }

  runMean <- function(x) {
    group = list(
      "names" = as.vector(unique(x)),
      "points" = ena.group(set$points.rotated, x, method=mean),
      "line.weights" = ena.group(set$line.weights, x, method=method)
    )
    group$points = (as.vector(as.matrix(group$points)) * scaleFactor) #[keep.dimensions];

    colnames(group$points) <- NULL;
    if(method == "sum") {
      group$line.weights = group$line.weights * length(which(x == T));
    }
    # cis = runCIs(set$points.rotated[x,]);
    cis = runCIs(matrix(set$points.rotated[x,], ncol = ncol(set$points.rotated)));

    group$conf.ints = cis$ci;
    group$outlier.ints = cis$oi;
    group$line.weights = as.vector(as.matrix(group$line.weights))
    colnames(group$line.weights) = NULL;

    if(scale == T) {
      if(method == "sum") {
        group$line.weights = scales::rescale(group$line.weights, c(0.1,1));
      }
      group$edge.saturation = scales::rescale(group$line.weights, c(0.25,1));
      group$edge.opacity = scales::rescale(group$line.weights, c(0.3,1));
    }

    group
  }

  if(is.list(units)) {
    lapply(names(units), function(x) {
      unitsFnd = set$enadata$unit.names %in% unlist(units[x])
      grp = runMean(unitsFnd);
      grp$names = NULL;
      grp$name = unlist(as.character(x));
      grp$rle = list( lengths = length(which(unitsFnd == T)), values = x );
      grp;
    })
  } else {
    unitsFnd = set$enadata$unit.names %in% units
    grp = runMean(unitsFnd);
    grp$names = NULL;
    grp$name = name
    grp$rle = list( lengths = length(which(unitsFnd == T)), values = name );
    grp;
  }
}
"""