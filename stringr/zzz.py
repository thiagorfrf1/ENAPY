zzz = """
.onLoad <- function(libname, pkgname) {
  globalVariables(c(".","ENA_ROW_IDX","ENA_UNIT","V1","V2","V3","ci.x","ci.y","e","handle","name","unit.groups","V","graph_from_data_frame","%>%","%<>%","X1","X2","dfDT.points","points.raw","lines","KEYCOL","ENA_CONV","..groupCol","..units","..metadata", "..codes", "..conversation"))
#   op <- options()
#   op.rENA <- list(
#     UNIT_NAMES = "ena.unit.names",
#     TRAJ_TYPES = c("accumulated","non-accumulated")
#   );
#
#   toset <- !("rENA" %in% names(op))
#   print(paste("ToSet:", toset));
#
#   if(toset) {
#     options(rENA = op.rENA)
#   }
#
#   invisible()
}
"""