#' BIC method for drc objects 
#' @param object object from which to extract the BIC 
#' @param ... further arguments (none used)
#' @author Tobias Verbeke
#' @importFrom stats4 BIC
#' @export
setMethod("BIC",
    signature(object = "drc"),
    function(object, ...)
    {
        BIC(logLik(object))
    }
)

setOldClass("drc")
