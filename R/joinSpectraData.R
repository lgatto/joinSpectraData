#' @title Join Spectra Data
#'
#' @description
#' 
#' Individual spectra data variable can be directly added with the
#' `$<-` or `[[<-` syntax. The `joinSpectraData()` function allows to
#' merge a `DataFrame` to the existing spectra data.
#'
#' This function diverges from the `merge()` method in two main ways:
#'
#' - The `by.x` and `by.y` column names must be of length 1. 
#'
#' - If variable names are shared in `x` and `y`, the spectra
#'   variables of `x` are not modified. It's only the `y` variables
#'   that are appended the suffix defined in `suffix.y`. This is to
#'   avoid modifying any core spectra variables that would lead to an
#'   invalid object.
#'
#' @import Spectra
#'
#' @importFrom methods as
#' 
#' @param x A [Spectra()] object.
#' 
#' @param y A `DataFrame` with spectra data to be merged with
#'     `spectraData(x)`.
#' 
#' @param by.x A `character(1)` specifying of the spectra variable
#'     used for merging. Default is `"spectrumId"`.
#' 
#' @param by.y A `character(1)` specifying of the column used for
#'     merging. Set to `by.x` if missing.
#' 
#' @param suffix.y A `character(1)` specifying the suffix to be used
#'     for making the names of columns in the merged spectra variables
#'     unique. This suffix will be used to amend `names(y)`, while
#'     `spectraVariables(x)` will remain unchanged.
#' 
#' @author Laurent Gatto. 
#'
#' @export
#'
#' @examples
#'
#' library("Spectra")
#' library("msdata")
#' library("magrittr")
#' library("PSM")
#'
#' ## Creat a Spectra object
#' ms <- Spectra(msdata::proteomics(pattern = "2014", full.names = TRUE)) %>%
#'     filterMsLevel(2L) %>%
#'     dropNaSpectraVariables()
#' spectraVariables(ms)
#' 
#' ## Additional spectra variables
#' id <- readPSMs(msdata::ident(full.names = TRUE)) %>%
#'     filterPSMs()
#' id$modLocation <- NULL
#' id <- unique(id)
#' names(id)
#'
#' ## Add the new spectra variables
#' ms <- joinSpectraData(ms, id,
#'                       by.x = "spectrumId",
#'                       by.y = "spectrumID")
#' spectraVariables(ms)
joinSpectraData <- function(x, y,
                             by.x = "spectrumId",
                             by.y,
                             suffix.y = ".y") {    
    if (missing(by.y))
        by.y <- by.x
    x_vars <- spectraVariables(x)
    y_vars <- names(y)
    if (length(by.x) != 1 | length(by.y) != 1)
        stop("'by.x' and 'by.y' must be of length 1.")
    if (!is.character(by.x) | !is.character(by.y))
        stop("'by.x' and 'by.y' must be characters.")
    if (any(duplicated(x_vars)))
        stop("Duplicated spectra variables in 'x'.")
    if (any(duplicated(y_vars)))
        stop("Duplicated names in 'y'.")
    if (!by.x %in% x_vars)
        stop("'by.x' not found in spectra variables.")
    if (!by.y %in% y_vars)
        stop("'by.y' not found in 'names(y)'.")    
    ## Keep only rows that (1) have a non-NA by.y and (2) that are in
    ## x[[by.x]] (given that we do a left join here).
    y <- y[!is.na(y[[by.y]]), ]
    y <- y[y[[by.y]] %in% spectraData(x)[[by.x]], ]
    k <- match(y[[by.y]], spectraData(x)[[by.x]])
    n <- length(x)
    ## Don't need by.y anymore
    y_vars <- y_vars[-grep(by.y, y_vars)]
    y <- y[, y_vars]
    ## Check if there are any shared x_vars and y_vars. If so, the
    ## y_vars are appended suffix.y.
    if (length(xy_vars <- intersect(x_vars, y_vars))) {
        y_vars[y_vars %in% xy_vars] <- paste0(y_vars[y_vars %in% xy_vars], suffix.y[1])
        names(y) <- y_vars
    }    
    for (y_var in y_vars) {
        ## Initialise the new column as a list, vector or List (in
        ## that order!)
        if (is.list(y[[y_var]])) {
            x[[y_var]] <- vector("list", length = n)
        } else if (is.vector(y[[y_var]])) {
            x[[y_var]] <- rep(NA, n)
        } else if (inherits(y[[y_var]], "List")) {
            x[[y_var]] <- as(vector("list", length = n), class(y[[y_var]]))
        }
        ## Populate the new column
        x[[y_var]][k] <- y[[y_var]]            
    }
    x
}
