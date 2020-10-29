##' This is the implementation for `merge,DFrame,DFrame`. It works
##' exactly as `merge,data.frame,data.frame` with the additional check
##' that the variables used for merging must be of the same class.
##'
##' @noRd
##' 
##' @author Laurent Gatto, based on
##'     `base::merge,data.frame,data.frame`.
mergeDFrame <- function (x, y,
                         by = intersect(names(x), names(y)),
                         by.x = by, 
                         by.y = by,
                         all = FALSE,
                         all.x = all,
                         all.y = all,
                         suffixes = c(".x", ".y"),
                         no.dups = TRUE,
                         incomparables = NULL,
                         ...) {
    fix.by <- function(by, df) {
        ## fix up 'by' to be a valid set of cols by number: 0 is row.names
        if (is.null(by)) by <- numeric()
        by <- as.vector(by)
        nc <- ncol(df)
        if (is.character(by)) {
            poss <- c("row.names", names(df))
            ## names(df) are not necessarily unique, so check for multiple matches.
            if (any(bad <- !charmatch(by, poss, 0L)))
                stop(ngettext(sum(bad),
                              "'by' must specify a uniquely valid column",
                              "'by' must specify uniquely valid columns"),
                     domain = NA)
            by <- match(by, poss) - 1L
        } else if (is.numeric(by)) {
            if (any(by < 0L) || any(by > nc))
                stop("'by' must match numbers of columns")
        } else if (is.logical(by)) {
            if (length(by) != nc) stop("'by' must match number of columns")
            by <- seq_along(by)[by]
        } else stop("'by' must specify one or more columns as numbers, names or logical")
        if (any(bad <- is.na(by)))
            stop(ngettext(sum(bad),
                          "'by' must specify a uniquely valid column",
                          "'by' must specify uniquely valid columns"),
                 domain = NA)
        unique(by)
    }
    nx <- nrow(x)
    ny <- nrow(y)
    if (nx >= 2^31 || ny >= 2^31) stop("long vectors are not supported")
    by.x <- fix.by(by.x, x)
    by.y <- fix.by(by.y, y)
    if ((l.b <- length(by.x)) != length(by.y))
        stop("'by.x' and 'by.y' specify different numbers of columns")
    if (l.b == 0L) { 
        ## return the cartesian product of x and y, fixing up common names
        nm <- nm.x <- names(x)
        nm.y <- names(y)
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if (has.common.nms) {
            names(x)[cnm] <- paste0(nm.x[cnm], suffixes[1L])
            cnm <- nm.y %in% nm
            names(y)[cnm] <- paste0(nm.y[cnm], suffixes[2L])
        }
        if (nx == 0L || ny == 0L) {
            res <- cbind(x[FALSE, ], y[FALSE, ])
        } else {
            ij <- expand.grid(seq_len(nx), seq_len(ny))
            res <- cbind(x[ij[, 1L], , drop = FALSE],
                         y[ij[, 2L], , drop = FALSE])
        }
    }
    else {
        if (any(by.x == 0L)) {
            x <- cbind(Row.names = I(row.names(x)), x)
            by.x <- by.x + 1L
        }
        if (any(by.y == 0L)) {
            y <- cbind(Row.names = I(row.names(y)), y)
            by.y <- by.y + 1L
        }
        row.names(x) <- NULL
        row.names(y) <- NULL
        ## create keys from 'by' columns:
        if (l.b == 1L) {                  # (be faster)
            bx <- x[, by.x]; if (is.factor(bx)) bx <- as.character(bx)
            by <- y[, by.y]; if (is.factor(by)) by <- as.character(by)
        } else {
            if (!is.null(incomparables))
                stop("'incomparables' is supported only for merging on a single column")
            ## Do these together for consistency in as.character.
            ## Use same set of names.
            bx <- x[, by.x, drop=FALSE]; by <- y[, by.y, drop=FALSE]
            names(bx) <- names(by) <- paste0("V", seq_len(ncol(bx)))
            ## Check that bx and by columns are of same class
            cl.bx <- vapply(bx, class, character(1))
            cl.by <- vapply(by, class, character(1))
            if (!all(cl.bx == cl.by))
                stop("The marching columns in 'x' and 'y' used for merging must ",
                     "be of identical classes.")
            ## Convert list|List columns to chars for subsequent matching
            bx_list_cols <- vapply(bx, function(.x) is.list(.x) | inherits(.x, "List"), TRUE)
            if (any(bx_list_cols))
                for (i in which(bx_list_cols)) 
                    bx[[i]] <- sapply(bx[[i]], paste, collapse = " ")
            by_list_cols <- vapply(by, function(.x) is.list(.x) | inherits(.x, "List"), TRUE)
            if (any(by_list_cols))
                for (i in which(by_list_cols)) 
                    by[[i]] <- sapply(by[[i]], paste, collapse = " ")

            bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
            bx <- bz[seq_len(nx)]
            by <- bz[nx + seq_len(ny)]
        }
        comm <- match(bx, by, 0L)
        bxy <- bx[comm > 0L]             # the keys which are in both
        xinds <- match(bx, bxy, 0L, incomparables)
        yinds <- match(by, bxy, 0L, incomparables)
        if (nx > 0L && ny > 0L)
            m <- .Internal(merge(xinds, yinds, all.x, all.y))
        else
            m <- list(xi = integer(), yi = integer(),
                      x.alone = seq_len(nx), y.alone = seq_len(ny))
        nm <- nm.x <- names(x)[-by.x]
        nm.by <- names(x)[by.x]
        nm.y <- names(y)[-by.y]
        ncx <- ncol(x)
        if (all.x) all.x <- (nxx <- length(m$x.alone)) > 0L
        if (all.y) all.y <- (nyy <- length(m$y.alone)) > 0L
        lxy <- length(m$xi)             # == length(m$yi)
        ## x = [ by | x ] :
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if (has.common.nms && nzchar(suffixes[1L]))
            nm.x[cnm] <- paste0(nm.x[cnm], suffixes[1L])
        x <- x[c(m$xi, if(all.x) m$x.alone),
               c(by.x, seq_len(ncx)[-by.x]), drop=FALSE]
        names(x) <- c(nm.by, nm.x)
        if (all.y) { ## add the 'y.alone' rows to x[]
            ## need to have factor levels extended as well -> using [cr]bind
            ya <- y[m$y.alone, by.y, drop = FALSE]
            names(ya) <- nm.by
            ## this used to use a logical matrix, but that was not good
            ## enough as x could be zero-row.
            ## workaround possibly duplicated names: PR#15618
            xa <- x[rep.int(NA_integer_, nyy), nm.x, drop=FALSE ]
            names(xa) <- nm.x
            x <- rbind(x, cbind(ya, xa))
        }
        ## y (w/o 'by'):
        if (has.common.nms && nzchar(suffixes[2L])) {
            cnm <- nm.y %in% nm
            nm.y[cnm] <- paste0(nm.y[cnm], suffixes[2L])
        }
        y <- y[c(m$yi, if(all.x) rep.int(1L, nxx), if(all.y) m$y.alone),
               -by.y, drop = FALSE]
        if (all.x) {
            zap <- (lxy+1L):(lxy+nxx)
            for(i in seq_along(y)) {
                ## do it this way to invoke methods for e.g. factor
                if (is.matrix(y[[i]])) y[[i]][zap, ] <- NA
                else is.na(y[[i]]) <- zap
            }
        }

        if (has.common.nms) names(y) <- nm.y
        ## If by.x %in% names(y) , duplicate column names may still arise.
        ## Apply suffixes to just y - keeping back compatibility when referring to by.x:
        if (no.dups && any((mi <- match(nm.by, names(y), 0L)) > 0L) && nzchar(suffixes[2L]))
            names(y)[mi] <- paste0(names(y)[mi], suffixes[2L])
        ## if(!no.dups) ## <<- speedup if we were sure we'd caught all dups above 
        nm <- c(names(x), names(y))
        if (any(d <- duplicated(nm)))
            if (sum(d) > 1L)
                warning("column names ",
                        paste(sQuote(nm[d]), collapse = ", "),
                        " are duplicated in the result", domain = NA)
            else
                warning("column name ", sQuote(nm[d]),
                        " is duplicated in the result", domain = NA)
        res <- cbind(x, y)
    }
    res
}



#' @title Join Spectra Data
#'
#' @description
#' 
#' Individual spectra data variable can be directly added with the
#' `$<-` syntax (see examples below).
#'
#' The `joinSpectraData()` function allows to merge a `DataFrame` to
#' the existing spectra data. It works like the `base::merge()`
#' function: the `by.x` (default is `"spectrumId"`) and `by.y`
#' arguments define the variable names in the current spectra data and
#' the news `DataFrame`. In addition, the `key` argument (default is
#' also `"spectrumId"`) is used the match the merged spectra data to
#' the original spectra.
#'
#' For details and background:
#' - https://github.com/rformassspectrometry/Spectra/issues/149
#' - https://github.com/Bioconductor/S4Vectors/issues/73
#'
#' @import Spectra
#' 
#' @export joinSpectraData
#' 
#' @param object A [Spectra()] object.
#' 
#' @param y A [DataFrame()] with spectra data to be merged with
#'     `spectraData(object)`.
#' 
#' @param by.x by.x Specification of the columns used for merging. See
#'     `base::merge()` for details. Default is `"spectrumId"`.
#' 
#' @param by.y by.y Specification of the columns used for merging. See
#'     `base::merge()` for details.
#' 
#' @param key A `character(1)` used to match the merged spectra data
#'     back to the spectra. Default is `"spectrumId"`.
#'
#' @author Laurent Gatto. Code of `joinSpectraData:::mergeDFrame()`
#'     based on `base::merge,data.frame,data.frame`.
joinSpectraData <- function(object, y,
                            by.x = "spectrumId",
                            by.y = NULL,
                            key = "spectrumId") {
    stopifnot(inherits(object, "Spectra"))
    if (is.null(by.y)) by.y <- by.x
    message("Merging by '", by.x, "' and '",
            by.y, "'.")
    x <- spectraData(object)
    xy <- mergeDFrame(x, y, all.x = TRUE, by.x = by.x, by.y = by.y)
    stopifnot(identical(nrow(xy), length(object)))
    message("Matching by '", key, "'.")
    o <- match(xy[, key], x[, key])
    xy <- xy[o, ]
    suppressWarnings(spectraData(object) <- xy)
    object
}