ena_generate = """
# # ##
# # # @title Accumulate and Generate
# # #
# # # @description Accumulate and Generate
# # #
# # # @details [TBD]
# # #
# # # @param file [TBD]
# # # @param window.size.back [TBD]
# # # @param units.by [TBD]
# # # @param conversations.by [TBD]
# # # @param code [TBD]
# # # @param units.used [TBD]
# # # @export
# # # @return list containing the accumulation and set
# # ##
# # ena.generate <- function(
# #   file,
# #   window.size.back,
# #   units.by,
# #   conversations.by,
# #   code,
# #   scale.nodes = T,
# #   units.used = NULL,
# #   ...
# # ) {
# #   args = list(...);
# #   conversations.used = NULL;
# #   weight.by = "binary";
# #   if(!is.null(args$conversations.used)) {
# #     conversations.used = args$conversations.used
# #     file$KEYCOL = merge_columns_c(file,conversations.by)
# #     file = file[file$KEYCOL %in% conversations.used,]
# #   }
# #   if(!is.null(args$weight.by)) {
# #     weight.by = args$weight.by
# #   }
# #
# #   accum = ena.accumulate.data.file(
# #     file = file,
# #     window.size.back = window.size.back,
# #     units.by = make.names(units.by),
# #     units.used = units.used,
# #     model = "EndPoint",
# #     conversations.by = make.names(conversations.by),
# #     codes = make.names(code),
# #     ...
# #   )
# #
# #   rotate.groups = NULL
# #   if(!is.null(args$rotate.by)) {
# #     rotate.meta = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,]
# #     rotate.col = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,][[names(args$rotate.by)[1]]]
# #     rotate.groups = list(
# #       rotate.col == args$rotate.by[[1]][1],
# #       rotate.col == args$rotate.by[[1]][2]
# #     )
# #   }
# #   set = ena.make.set(
# #     enadata = accum,
# #     norm.by = ifelse((is.null(args$sphere.norm) || args$sphere.norm==T),sphere_norm_c,dont_sphere_norm_c),
# #     rotation.by = if(is.null(rotate.groups)) ena.svd else ena.rotate.by.mean, #ifelse(is.null(rotate.groups), NULL, ena.rotate.by.mean),
# #     rotation.params = rotate.groups,
# #     ...
# #   )
# #
# # @param file [TBD]
# # @param window.size.back [TBD]
# # @param units.by [TBD]
# # @param conversations.by [TBD]
# # @param code [TBD]
# # @param units.used [TBD]
# # @export
# # @return list containing the accumulation and set
# ##
# ena.generate <- function(
#   file,
#   window.size.back,
#   units.by,
#   conversations.by,
#   code,
#   scale.nodes = T,
#   units.used = NULL,
#   dimensions = 6,
#   include.meta = F,
#   ...
# ) {
#   startedTime = as.numeric(Sys.time())
#   args = list(...);
#   unit.groups = NULL;
#   conversations.used = NULL;
#   weight.by = "binary";

#   if(!is.null(args$conversations.used)) {
#     conversations.used = args$conversations.used
#     file$KEYCOL = rENA:::merge_columns_c(file, make.names(conversations.by))
#     file = file[file$KEYCOL %in% conversations.used,]
#   }

#   if(!is.null(args$weight.by)) {
#     weight.by = args$weight.by;
#   }
#   if(!is.null(args$unit.groups)){
#     if(is.data.frame(args$unit.groups)) {
#       unit.groups = args$unit.groups$units;
#       names(unit.groups) = args$unit.groups$name;
#     } else {
#       unit.groups = list();
#       group.json = jsonlite::fromJSON(args$unit.groups)
#       for(grp in 1:length(group.json$name)) {
#         unit.groups[group.json$name[grp]] = group.json$units[grp];
#       }
#     }
#   }

#   accum = ena.accumulate.data.file(
#     file = file,
#     window.size.back = window.size.back,
#     units.by = make.names(units.by),
#     units.used = units.used,
#     model = "EndPoint",
#     conversations.by = make.names(conversations.by),
#     codes = make.names(code),
#     include.meta = include.meta,
#     ...
#   )

#   rotate.groups = NULL
#   if(!is.null(args$rotate.by)) {
#     # rotate.meta = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,]
#     # rotate.col = accum$metadata[accum$metadata$ENA_UNIT %in% accum$unit.names,][[make.names(names(args$rotate.by)[1])]]
#     # rotate.groups = list(
#     #   rotate.col == args$rotate.by[[1]][1],
#     #   rotate.col == args$rotate.by[[1]][2]
#     # )
#     rotate.groups = lapply(args$rotate.by, function(x) accum$unit.names %in% x )
#   }

#   use.to.norm = rENA:::fun_skip_sphere_norm;
#   if(is.null(args$sphere.norm) || args$sphere.norm == T) {
#     use.to.norm = rENA:::fun_sphere_norm
#   }
#   rotation.set = NULL
#   if(!is.null(args$rotation.matrix)) {
#     rotation.set = ENARotationSet$new(
#       rotation = args$rotation.matrix$rotation$rotation,
#       node.positions = args$rotation.matrix$rotation$node.positions,
#       codes = args$rotation.matrix$codes
#     );
#     colnames(rotation.set$rotation) = args$rotation.matrix$dimensions;
#   }
#   set = ena.make.set(
#     enadata = accum,
#     norm.by = use.to.norm,
#     rotation.by = if(is.null(rotate.groups)) rENA:::ena.svd else rENA:::ena.rotate.by.mean,
#     rotation.params = rotate.groups,
#     rotation.set = rotation.set,
#     dimensions = dimensions,
#     ...
#   )

#   use.dimensions = 1:2
#   # if(!is.null(args$keep.dimensions)) {
#   #   use.dimensions = which(colnames(set$points.rotated) %in% args$keep.dimensions)
#   # }

#   group.names = NULL;
#   if(length(units.by)>1) {
#     group.names = unique(set$enadata$units[[make.names(units.by)[[1]]]])
#   } else {
#     group.names = units.by
#   }
#   group.cnt = length(group.names);
#   conf.ints = list();
#   outlier.ints = matrix(0, nrow=(group.cnt), ncol=(2));

#   set$points.rotated.scaled = set$points.rotated;
#   scaleFactor = 1.0
#   if(scale.nodes == T) {
#     np.min.x = min(set$node.positions[,use.dimensions[1]])
#     np.min.y = min(set$node.positions[,use.dimensions[2]])
#     rp.min.x = min(set$points.rotated[,use.dimensions[1]])
#     rp.min.y = min(set$points.rotated[,use.dimensions[2]])
#     maxMin = abs(max(np.min.x / rp.min.x, np.min.y / rp.min.y))

#     np.max.x = max(set$node.positions[,use.dimensions[1]])
#     np.max.y = max(set$node.positions[,use.dimensions[2]])
#     rp.max.x = max(set$points.rotated[,use.dimensions[1]])
#     rp.max.y = max(set$points.rotated[,use.dimensions[2]])
#     maxMax = abs(max(np.max.x / rp.max.x, np.max.y / rp.max.y))
#     scaleFactor = min(maxMin, maxMax)
#     # set$points.rotated = set$points.rotated * scaleFactor;
#     set$points.rotated.scaled = set$points.rotated * scaleFactor;
#   }

#   # groups = NULL
#   group.method = "mean"
#   if(!is.null(args$weight.network.by) && (args$weight.network.by %in% c("mean","sum"))) {
#     group.method = args$weight.network.by;
#   }

#   groups = list()
#   if(!is.null(unit.groups) && length(unit.groups) > 0){
#     # for(i in 1:length(names(unit.groups))) {
#     #   groups[[length(groups)+1]] = ena.unit.group(set, set$enadata$unit.names[set$enadata$unit.names %in% unit.groups[[i]]], name = names(unit.groups)[i], method = group.method, scaleFactor = scaleFactor)
#     # }
#     groups = lapply(names(unit.groups), function(nm) {
#       ena.unit.group(
#         set,
#         set$enadata$unit.names[set$enadata$unit.names %in% unit.groups[[nm]]],
#         name = nm,
#         method = group.method,
#         scaleFactor = scaleFactor
#         # ,keep.dimensions = use.dimensions
#       )
#     })
#   }

#   if(
#     !is.null(args$output) && args$output == "save" &&
#     !is.null(args$output.to)
#   ) {
#     setName = tools::file_path_sans_ext(basename(args$output.to))
#     env = environment()
#     assign(x = setName, value = set, envir = env);
#     env[[setName]] = get(x = setName, envir = env)

#     tmp <- tempfile(fileext = ".rdata")
#     on.exit(unlink(tmp))
#     save(list = c(setName), file = tmp, envir = env)
#     bucket <- aws.s3::get_bucketname(args$output.to)
#     object <- aws.s3:::get_objectkey.character(args$output.to)
#     return(aws.s3::put_object(file = tmp, bucket = bucket, object = object));
#   }
#   else {
#     nodes = data.frame(set$node.positions);
#     nodes$weight = rep(0, nrow(nodes))
#     node.rows = rownames(set$node.positions);
#     estimate.over.units = (!(set$enadata$unit.names %in% args$units.exclude))
#     weights = matrix(0, ncol=nrow(set$node.positions), nrow=length(which(estimate.over.units)));

#     colnames(weights) = node.rows
#     network.scaled = set$line.weights[estimate.over.units,];
#     # if(!is.null(scale.weights) && scale.weights == T) {
#     #   network.scaled = network.scaled * (1 / max(abs(network.scaled)));
#     # }

#     mat = set$enadata$adjacency.matrix;
#     # for (x in 1:nrow(network.scaled)) {
#     #   weights[x, ] = sapply(node.rows, function(y) {
#     #     # sum(network.scaled[x,as.logical(colSums(!is.na(apply(mat,2,match, y))))])
#     #     sum(network.scaled[x, as.logical(colSums(mat == y))])
#     #   })
#     #   # network.thickness = network.scaled[x,] #scales::rescale(abs(network.scaled[x,]), thickness);
#     #   # for (i in 1:ncol(mat)) {
#     #   #   weights[x,node.rows==mat[1,i]] = weights[x,node.rows==mat[1,i]] + network.thickness[i];
#     #   # }
#     # }
#     weights = sapply(node.rows, function(x) rowSums(network.scaled[,as.logical(colSums(mat == x) )]))

#     # #weights = t(apply(weights, 1, scales::rescale, c(1,ncol(weights))));
#     weights = scales::rescale(weights, c(1,ncol(weights)));

#     set$line.weights[estimate.over.units,] = set$line.weights[estimate.over.units,];

#     # If not included, remove the weights as to not effect the scaling
#     set$line.weights[!estimate.over.units,] = 0
#     scaleRange = c(min(set$line.weights[estimate.over.units,]) ,1);
#     if(scaleRange[1] < 0.1 && min(set$line.weights)>0) {
#       scaleRange[1] = 0.1;
#     }

#     set$line.weights = scales::rescale(set$line.weights, to=scaleRange, from=range(set$line.weights, na.rm = T, finite = T))
#     # adjRows = triIndices(length(code)) + 1
#     # codedRow1 = code[adjRows[1,]];
#     # codedRow2 = code[adjRows[2,]];

#     tmp = getwd();
#     sess = regexec("temp/(x[^/]*)/workspace", tmp)[[1]]
#     assign("set", set, envir = parent.frame())

#     dimension.names = paste("SVD",1:ncol(set$points.rotated), sep="")
#     if(length(set$function.params$rotation.params) == 2) dimension.names[1] = "MR1"
#     # if(length(args$plotted.nodes) == 2) {
#     #   methods = ena.methods(enaset = set, tool = "webENA", tool.version = "0.1.0", comparison = "parametric", comparison.groups = args$plotted.nodes)
#     # } else {
#     #   methods = ena.methods(enaset = set, tool = "webENA", tool.version = "0.1.0")
#     # }
#     doneTime = as.numeric(Sys.time())

#     ### Limit dimensions
#     # set$points.rotated = set$points.rotated[,use.dimensions]
#     # set$points.rotated.scaled = set$points.rotated.scaled[,use.dimensions]
#     # set$node.positions = set$node.positions[,use.dimensions]
#     # set$rotation.set$rotation = set$rotation.set$rotation[,use.dimensions]
#     return(list(
#       codes = make.names(code),
#       adjacency.matrix = mat, #rbind(codedRow1, codedRow2),
#       set = set,
#       # methods = readChar(methods, file.info(methods)$size),
#       custom.rotation = if(!is.null(rotation.set)) T else F,
#       custom.rotation.set = rotation.set,
#       dimensions = dimension.names, #colnames(set$points.rotated),
#       session = substr(tmp, start=sess[2], stop=sess[2]+attr(sess, "match.length")[2]-1),
#       # groups = groups,
#       groups = groups,
#       scaled = scale.nodes,
#       node.sizes = weights,
#       esitmated.over = args$units.exclude,
#       edge.saturation = scales::rescale(set$line.weights, c(0.25,1)),
#       edge.opacity = scales::rescale(set$line.weights, c(0.3,1)),
#       startedTime = startedTime,
#       doneTime = doneTime,
#       durationTime = doneTime - startedTime
#     ));
#   }
# }
"""