ena_conversations = """
##
#' @title Find conversations by unit
#'
#' @description Find rows of conversations by unit
#'
#' @details [TBD]
#'
#' @param set [TBD]
#' @param units [TBD]
#' @param units.by [TBD]
#' @param codes [TBD]
#' @param conversation.by [TBD]
#' @param window [TBD]
#' @param conversation.exclude [TBD]
#'
#' @examples
#' data(RS.data)
#'
#' codeNames = c('Data','Technical.Constraints','Performance.Parameters',
#'               'Client.and.Consultant.Requests','Design.Reasoning',
#'               'Collaboration');
#'
#' accum = ena.accumulate.data(
#'   units = RS.data[,c("Condition","UserName")],
#'   conversation = RS.data[,c("Condition","GroupName")],
#'   metadata = RS.data[,c("CONFIDENCE.Change","CONFIDENCE.Pre",
#'                         "CONFIDENCE.Post","C.Change")],
#'   codes = RS.data[,codeNames],
#'   model = "EndPoint",
#'   window.size.back = 4
#' );
#' set = ena.make.set(
#'   enadata = accum,
#'   rotation.by = ena.rotate.by.mean,
#'   rotation.params = list(accum$meta.data$Condition=="FirstGame",
#'                          accum$meta.data$Condition=="SecondGame")
#' );
#' ena.conversations(set = RS.data,
#'   units = c("FirstGame.steven z"), units.by=c("Condition","UserName"),
#'   conversation.by = c("Condition","GroupName"),
#'   codes=codeNames, window = 4
#' )
#'
#' @export
#' @return list containing row indices representing conversations
##
ena.conversations = function(set, units, units.by=NULL, codes=NULL, conversation.by = NULL, window = 4, conversation.exclude = c()) {
  # rawData = data.table::copy(set$enadata$raw);
  if(is.null(units.by)) {
    units.by = set$`_function.params`$units.by;
  }
  # conversation.by = set$enadata$function.params$conversations.by;
  # window = set$enadata$function.params$window.size.back;
  # rawAcc = data.table::copy(set$enadata$accumulated.adjacency.vectors);
  if(is(set, "ena.set")) {
    rawAcc2 = set$model$raw.input
  } else {
    rawAcc2 = data.table::data.table(set) #$enadata$raw);
  }

  # rawAcc$KEYCOL = merge_columns_c(rawAcc, conversation.by)
  rawAcc2$KEYCOL = merge_columns_c(rawAcc2, conversation.by)

  # conversationsTable = rawAcc[, paste(.I, collapse = ","), by = c(conversation.by)]
  conversationsTable2 = rawAcc2[, paste(.I, collapse = ","), by = c(conversation.by)]

  # rows = sapply(conversationsTable$V1, function(x) as.numeric(unlist(strsplit(x, split=","))),USE.NAMES = T)
  rows2 = lapply(conversationsTable2$V1, function(x) as.numeric(unlist(strsplit(x, split=","))))
  # browser()
  # names(rows) = merge_columns_c(conversationsTable,conversation.by); #unique(rawAcc[,KEYCOL])
  names(rows2) = merge_columns_c(conversationsTable2,conversation.by); #unique(rawAcc[,KEYCOL])

  # unitRows = merge_columns_c(rawAcc[,c(units.by),with=F], units.by)
  unitRows2 = merge_columns_c(rawAcc2[,c(units.by),with=F], units.by)

  # adjCol = set$enadata$adjacency.matrix[1,] %in%  codes[1] & set$enadata$adjacency.matrix[2,] %in% codes[2]
  # adjColName = paste("adjacency.code.", which(adjCol), sep = "")
  # codedUnitRows = which(unitRows %in% units & rawAcc[[adjColName]] == 1)

  codedRows = rawAcc2[, rowSums(.SD), .SDcols = codes] > 0
  codedUnitRows2 = which(unitRows2 %in% units & codedRows)
  codedUnitRows2 = codedUnitRows2[!(codedUnitRows2 %in% as.vector(unlist(rows2[conversation.exclude])))]
  # codedUnitRowConvs = rawAcc[codedUnitRows,KEYCOL];
  codedUnitRowConvs2 = rawAcc2[codedUnitRows2,KEYCOL];

  # codedUnitRowConvsAll = NULL;
  codedUnitRowConvsAll2 = NULL;
  unitRowsNotCooccurred = c()
  if(length(codedUnitRows2) > 0) {
    codedUnitRowConvsAll = unique(unlist(sapply(X = 1:length(codedUnitRows2), simplify = F, FUN = function(x) {
      thisConvRows = rows2[[codedUnitRowConvs2[x]]]
      thisRowInConv = which(thisConvRows == codedUnitRows2[x])
      winUse = ifelse(is.infinite(window), thisRowInConv, window)
      thisRowAndWindow = rep(thisRowInConv,winUse) - (winUse-1):0
      coOccursFound = all(rawAcc2[thisConvRows[thisRowAndWindow[thisRowAndWindow > 0]], lapply(.SD, sum), .SDcols=codes] > 0)
      if(coOccursFound) {
        thisConvRows[thisRowAndWindow[thisRowAndWindow > 0]]
      } else {
        unitRowsNotCooccurred <<- c(unitRowsNotCooccurred, thisConvRows[thisRowInConv])
        # coOccursFound
        NULL
      }
    })))
  }
  return(list(
    conversations = as.list(rows2),
    unitConvs = unique(rawAcc2[codedUnitRows2,KEYCOL]),
    allRows = codedUnitRowConvsAll,
    unitRows = codedUnitRows2,
    toRemove = unitRowsNotCooccurred
  ));
}

"""