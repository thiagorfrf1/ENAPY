ENAset = """
####
#' ENAset R6class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export

#' @field enadata An \code{\link{ENAdata}} object originally used to create the set
#' @field points.raw A data frame containing accumulated adjacency (co-occurrence) vectors per unit
#' @field points.normed.centered A data frame of centered normed accumulated adjacency (co-occurrence) vectors for each unit
#' @field points.rotated A data frame of point positions for number of dimensions specified in ena.make.set (i.e., the centered, normed, and rotated data)
#' @field line.weights A data frame of connections strengths per unit (Data frame of normed accumu- lated adjacency (co-occurrence) vectors for each unit)
#' @field node.positions - A data frame of positions for each code
#' @field codes - A vector of code names
#' @field rotation.set - An \code{\link{ENARotationSet}} object
#' @field correlation - A data frame of spearman and pearson correlations for each dimension specified
#' @field variance - A vector of variance accounted for by each dimension specified
#' @field centroids - A matrix of the calculated centroid positions
#' @field function.call - The string representation of function called
#' @field function.params - A list of all parameters sent to function call
#'
####

ENAset = R6::R6Class("ENAset",
   public = list(
     #####
     ### Constructor - documented in main class declaration
     #####
     initialize = function(
       enadata,
       dimensions = 2,

       norm.by = fun_sphere_norm,

       rotation.by = ena.svd.R6,
       rotation.params = NULL,
       rotation.set = NULL,

       #center.data = center_data_c,    ### made local to run
       node.position.method = lws.positions.sq.R6,
       endpoints.only = T,
       ...
     ) {
       self$enadata <- enadata;

       private$dimensions <- dimensions;

       self$codes <- enadata$codes;

       self$function.call <- sys.call(-1);

       self$function.params$norm.by <- norm.by;    #was sphere_norm
       #self$function.params$center.data <- center.data;
       self$function.params$node.position.method <- node.position.method;    #was position.method
       self$function.params$rotation.by <- rotation.by;
       self$function.params$rotation.params <- rotation.params;
       self$function.params$rotation.set <- rotation.set;
       self$function.params$endpoints.only <- endpoints.only;

       private$args <- list(...);
     },
     #####
     ### END: Constructor
     #####

     #####
     ## Public Properties
     #####
     rotation_dists = NULL,  #leave for now - to be removed for a temp variable
     enadata = NULL,
     points.raw = NULL,    #was data$raw
     points.normed.centered = NULL,    #was data$centered$normed
     points.rotated = NULL,    #was data$centered$rotated
     points.rotated.scaled = NULL,
     points.rotated.non.zero = NULL,
     line.weights = NULL,   #was data$normed
     line.weights.non.zero = NULL,
     line.weights.unrotated = NULL,
     node.positions = NULL,  #was nodes$positions$scaled
     codes = NULL,
     rotation.set = NULL,   ## new - ENARotation object
     correlations = NULL,   #not formerly listed, comes from optimized node positions in egr.positions
     variance = NULL,     #was self$data$centered$latent
     centroids = NULL,
     function.call = NULL,     #new - string reping function call
     function.params = list(   #list containing parameters function was called with
       norm.by = NULL,
       node.position.method = NULL,
       rotation.by = NULL,
       rotation.params = NULL,
       endpoints.only = NULL
     ),
     #####
     ## END: Public Properties
     #####

     #####
     ## Public Functions
     #####

     ####
     # \code{process()} - Process the ENAset.
     # \preformatted{}
     ####
     process = function() {
       return(private$run())
     },

     ####
     get.data = function(wh = c("normed","centered","rotated"), with.meta = T) {
       wh =  match.arg(wh);
       data = NULL;
       if( wh == "normed" ) {
         data = self$line.weights
       } else if ( wh == "centered" ) {
         data = self$points.normed.centered
       } else if ( wh == "rotated" ) {
         data = self$points.rotated
       }
       df.to.return = NULL;
       if(with.meta == T) {
         data.units = attr(data, opts$UNIT_NAMES);
         df.to.return = merge(
           data.table::data.table(
             data, data.units,
             ENA_UNIT=merge_columns_c(data.units, self$enadata$get("units.by")),
             TRAJ_UNIT=merge_columns_c(data.units, c(self$enadata$get("units.by"), self$enadata$get("trajectory.by")))
           ),
           self$enadata$add.metadata()
         )
       } else {
         df.to.return = data
       }
       df.to.return
     },

     ####
     # \code{get()} - Return a read-only property
     # \preformatted{  Example:
     #     get( x = 'file' )}
     # \preformatted{  Parameters:
     #      x - Property to return. Defaults to 'file', returning the original data}
     ####
     get = function(x = "enadata") {
       return(private[[x]])
     },
     print = function(...) {
       args = list(...);
       fields = NULL;
       to.print = list();

       if (!is.null(args$fields)) {
         fields = args$fields
       } else if(!is.null(private$args$fields)) {
         fields = private$args$fields
       } else {
         #fields = Filter(function(f) { (class(self[[f]]) != "function") }, names(get(class(self))$public_fields))
         fields = Filter(function(f) {
           cls = class(self[[f]]);
           !is(self[[f]], "function") && !is.null(self[[f]])
         }, names(get(class(self))$public_fields))
       }

       for(field in fields) {
         if(grepl("\\$", field)) {
           parts = Filter(function(f) { f!="" }, strsplit(field,"\\$")[[1]])
           to.print[[field]] = Reduce(function(o, i) { o[[i]] }, parts, self)
         } else {
           to.print[[field]] = self[[field]]
         }
       }
       return(to.print);
     }

     #####
     ## END: Public Functions
     #####
   ),

   private = list(
     #####
     ## Private Properties
     #####
     args = NULL,
     data.original = NULL,
     optim = NULL,

     #moved from public
     dimensions = 2,
     #####
     ## END: Private Properties
     #####

     #####
     ## Private Functions
     #####
     run = function() {       
       df = self$enadata$adjacency.vectors;

       ###
       # Backup of ENA data, this is not touched again.
       ###
       #private$data.original = df[,grep("adjacency.code", colnames(df)), with=F];
       private$data.original = df;
       ###
       # Copy of the original data, this is used for all
       # further operations. Unlike, `data.original`, this
       # is likely to be overwritten.
       ###
       self$points.raw = data.table::copy(private$data.original);

       ###
       # Normalize the raw data using self$function.params$norm.by,
       # which defaults to calling rENA::.sphere_norm
       ###
       self$line.weights = self$function.params$norm.by(self$points.raw);

       ###
       # Convert the string vector of code names to their corresponding
       # co-occurence names and set as colnames for the self$line.weights
       ##
       codeNames_tri = svector_to_ut(self$enadata$codes);

       colnames(self$line.weights) = codeNames_tri;
       # set the rownames to that of the original ENAdata file object
       rownames(self$line.weights) = rownames(df);

       attr(self$line.weights, opts$UNIT_NAMES) = attr(df, opts$UNIT_NAMES) #df[, .SD, with=T, .SDcols=self$enadata$get("unitsBy")];
       ###


       ###
       # Center the normed data
       # FIX - store as $data$centered
       ###
       #### ISSUE
       self$points.normed.centered = center_data_c(self$line.weights);

       colnames(self$points.normed.centered) = codeNames_tri;
       rownames(self$points.normed.centered) = rownames(df);
       attr(self$points.normed.centered, opts$UNIT_NAMES) = attr(self$enadata$adjacency.vectors.raw, opts$UNIT_NAMES)

       ###

       ###
       # Generate and Assign the rotation set
       ###
        if(!is.null(self$function.params$rotation.by) && is.null(self$function.params$rotation.set)) {
          self$rotation.set = do.call(self$function.params$rotation.by, list(self, self$function.params$rotation.params));
        } else if (!is.null(self$function.params$rotation.set)) {
          if(is(self$function.params$rotation.set, "ENARotationSet")) {
            print("Using custom rotation.set.")

            self$rotation.set = self$function.params$rotation.set;
          } else {
            stop("Supplied rotation.set is not an instance of ENARotationSet")
          }
        } else {
          stop("Unable to find or create a rotation set")
        }
       ###

       ###
       # Generated the rotated points
       ###
        self$points.rotated = self$points.normed.centered %*% self$rotation.set$rotation;
        private$dimensions = min(private$dimensions, ncol(self$points.rotated))
        attr(self$points.rotated, opts$UNIT_NAMES) = attr(self$points.normed.centered, opts$UNIT_NAMES);
       ###

       ###
       # Calculate node positions
       #  - The supplied methoed is responsible is expected to return a list
       #    with two keys, "node.positions" and "centroids"
       ###
        if(!is.null(self$rotation.set) && is.null(self$function.params$rotation.set)) {
          positions = self$function.params$node.position.method(self);
          if(all(names(positions) %in% c("node.positions","centroids"))) {
            self$node.positions = positions$node.positions
            self$centroids = positions$centroids

            self$rotation.set$node.positions = positions$node.positions
          } else {
            print("The node position method didn't return back the expected objects:")
            print("    Expected: c('node.positions','centroids')");
            print(paste("    Received: ",names(positions),sep=""));
          }
        } else if (!is.null(self$function.params$rotation.set)) {
          self$node.positions = self$function.params$rotation.set$node.positions
        } else {
          stop("Unable to determine the node positions either by calculating
                them using `node.position.method` or using a supplied
                `rotation.set`");
        }
       ###

       ###
       # Variance
       ###
       variance.of.rotated.data = var(self$points.rotated)
       diagonal.of.variance.of.rotated.data = as.vector(diag(variance.of.rotated.data))
       self$variance = diagonal.of.variance.of.rotated.data/sum(diagonal.of.variance.of.rotated.data)

       return(self);
     }
     #####
     ## END: Private Functions
     #####
   )
)
"""