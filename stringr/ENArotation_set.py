ENArotation_set = """
#######
#' ENARotationSet R6class
#
## Node positions
## Rotation matrix
#
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @export
#
# @param windowSize Integer used to select the size of each stanza window within a conversation
# @param binary Logical, whether to convert code values to binary or allow for weigthed values
# @param unitsSelected deprecated
#
# @section Public ENARotationSet methods:
#######
ENARotationSet = R6::R6Class("ENARotationSet",
  public = list(

    #######
    ### Constructor - documented in main class declaration
    #######
    initialize = function(
      rotation,
      codes,
      node.positions,
      eigenvalues = NULL
    ) {
      self$node.positions = node.positions;
      self$rotation = rotation;
      self$codes = codes;
      if(!is.null(codes) && !is.null(self$node.positions)) {
       rownames(self$node.positions) = codes;
      }
      self$eigenvalues = eigenvalues;
    },

    ####
    ## Public Properties
    ####
      rotation = NULL,
      node.positions = NULL,
      codes = NULL,
      eigenvalues = NULL,
    ####
    ## END: Public Properties
    ####,

    ####
    ## Public Functions
    ####
      print = function(...) {
        args = list(...);
        fields = NULL;
        to.print = list();

        if (!is.null(args$fields)) {
          fields = args$fields
        } else if(!is.null(private$args$fields)) {
          fields = private$args$fields
        } else {
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
    ####
    ## END: Public Functions
    ####
  ),
  private = list(
    #####
    ## Private Properties
    #####
      args = NULL
    #####
    ## END: Private Properties
    #####
  )
)
"""