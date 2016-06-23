##' timevaryingPL function
##'
##' A function to 
##'
##' @param formula X 
##' @param t0 X
##' @param dist X 
##' @param data X 
##' @param optimcontrol X 
##' @return ...
##' @export

timevaryingPL <- function(formula,t0,dist,data,optimcontrol=NULL){
    responsename <- as.character(formula[[2]])
    survivaldata <- data[[responsename]]

    ##########
    # This chunk of code borrowed from flexsurvreg    
    ##########  
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0){ 
        stop("A \"formula\" argument is required")
    }
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())

    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    ##########
    # End of borrowed code    
    ##########

    X <- X[, -1, drop = FALSE] 

    #ord <- order(survivaldata[,1])
    #survivaldata <- survivaldata[ord,]
    #X <- X[ord,]

    
}