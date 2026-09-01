# Encoding-robust delimited-text reader
#
# These helpers provide deterministic and heuristic decoding of UTF-8,
# UTF-16, and Windows-1252 text files without additional dependencies.


#' Heuristically detect the encoding of a text file
#'
#' Detects UTF-16 byte-order marks and byte patterns, then distinguishes
#' UTF-8 from Windows-1252 by attempting to read the file as UTF-8.
#'
#' @param path Character path to the text file.
#'
#' @return A character encoding name suitable for the `encoding` argument of
#'   [base::file()].
#'
#' @keywords internal
#' @noRd
.detect_text_encoding <- function(path) {
  rb <- readBin(
    path,
    what = "raw",
    n = file.info(path)$size
  )

  n <- length(rb)

  if (!n) {
    return("UTF-8")
  }

  # Detect byte-order marks first.
  if (n >= 2L) {
    b2 <- as.integer(rb[1:2])

    if (identical(b2, c(0xFF, 0xFE))) {
      return("UTF-16LE")
    }

    if (identical(b2, c(0xFE, 0xFF))) {
      return("UTF-16BE")
    }
  }

  if (n >= 3L) {
    b3 <- as.integer(rb[1:3])

    if (identical(b3, c(0xEF, 0xBB, 0xBF))) {
      return("UTF-8")
    }
  }

  # Detect likely UTF-16 without a byte-order mark. Text in the ASCII range
  # produces many NUL bytes in one parity of a UTF-16 byte sequence.
  zeros <- which(rb == as.raw(0x00))
  p_zero <- length(zeros) / n

  if (p_zero > 0.20) {
    even_zeros <- sum(zeros %% 2L == 0L)
    odd_zeros <- sum(zeros %% 2L == 1L)

    if (even_zeros >= odd_zeros) {
      return("UTF-16LE")
    }

    return("UTF-16BE")
  }

  read_with <- function(enc) {
    con <- file(
      path,
      open = "r",
      encoding = enc
    )

    on.exit(close(con), add = TRUE)

    had_warn <- FALSE
    warn_msg <- NULL

    txt <- withCallingHandlers(
      readLines(con, warn = FALSE),
      warning = function(w) {
        had_warn <<- TRUE
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )

    list(
      txt = txt,
      had_warn = had_warn,
      warn_msg = warn_msg
    )
  }

  res <- read_with("UTF-8")

  has_replacement <- length(res$txt) &&
    any(grepl("\ufffd", res$txt))

  if (
    isTRUE(res$had_warn) ||
    isTRUE(has_replacement) ||
    !length(res$txt)
  ) {
    return("windows-1252")
  }

  "UTF-8"
}


#' Clean a decoded table column name
#'
#' Removes a leading byte-order mark, leading replacement characters, and
#' non-breaking spaces, then trims surrounding whitespace.
#'
#' @param x Character vector of names.
#'
#' @return A cleaned character vector.
#'
#' @keywords internal
#' @noRd
.clean_name <- function(x) {
  x <- sub("^\ufeff", "", x)
  x <- sub("^\ufffd+", "", x)
  x <- gsub("\u00a0", " ", x, useBytes = FALSE)

  trimws(x)
}


#' Read a delimited text file normalized to UTF-8
#'
#' Reads a delimited text file after detecting likely UTF-8, UTF-16, or
#' Windows-1252 encoding. Decoded text is normalized to UTF-8 before being
#' parsed with [utils::read.table()].
#'
#' Detection proceeds in the following order:
#'
#' * UTF-16 little-endian or big-endian byte-order mark
#' * UTF-8 byte-order mark
#' * UTF-16 without a byte-order mark, inferred from NUL-byte patterns
#' * UTF-8
#' * Windows-1252 as a fallback
#'
#' Column names are cleaned of byte-order marks, replacement characters,
#' non-breaking spaces, and surrounding whitespace. Character columns are
#' returned in UTF-8.
#'
#' @param path Character path to the input file.
#' @param sep Character field separator. The default is a tab character.
#' @param header Logical. Whether the file contains a header line.
#' @param quote Character string specifying quoting characters. An empty string
#'   disables quote processing.
#' @param dec Character decimal separator.
#' @param fix_html_entities Logical. If `TRUE`, replace the HTML entity
#'   `&amp;` with `&` in character columns.
#' @param ... Additional arguments passed to [utils::read.table()].
#'
#' @return A data frame whose character columns are encoded in UTF-8.
#'
#' @examples
#' \dontrun{
#' plan <- read_table_utf8(
#'   "grafting-plan.tsv",
#'   sep = "\t",
#'   header = TRUE,
#'   quote = ""
#' )
#' }
#'
#' @export
read_table_utf8 <- function(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    dec = ".",
    fix_html_entities = TRUE,
    ...
) {
  if (
    length(path) != 1L ||
    !is.character(path) ||
    is.na(path) ||
    !nzchar(path)
  ) {
    stop(
      "`path` must be a single non-empty character value.",
      call. = FALSE
    )
  }

  if (!file.exists(path)) {
    stop(
      "File not found: ",
      path,
      call. = FALSE
    )
  }

  enc <- .detect_text_encoding(path)

  con <- file(
    path,
    open = "r",
    encoding = enc
  )

  on.exit(close(con), add = TRUE)

  txt <- suppressWarnings(
    readLines(con, warn = FALSE)
  )

  # If the UTF-8 attempt unexpectedly produces no text, try Windows-1252.
  if (!length(txt) && identical(enc, "UTF-8")) {
    con2 <- file(
      path,
      open = "r",
      encoding = "windows-1252"
    )

    on.exit(close(con2), add = TRUE)

    txt <- suppressWarnings(
      readLines(con2, warn = FALSE)
    )
  }

  if (!length(txt)) {
    stop(
      "No readable text lines found in file: ",
      path,
      call. = FALSE
    )
  }

  # Strip an optional byte-order mark and normalize text to UTF-8.
  txt[1] <- sub("^\ufeff", "", txt[1])
  txt <- enc2utf8(txt)
  txt <- gsub("\u00a0", " ", txt, useBytes = FALSE)

  tc <- textConnection(txt)
  on.exit(close(tc), add = TRUE)

  df <- utils::read.table(
    file = tc,
    sep = sep,
    header = header,
    quote = quote,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    dec = dec,
    comment.char = "",
    fill = TRUE,
    ...
  )

  names(df) <- vapply(
    names(df),
    .clean_name,
    character(1)
  )

  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x <- enc2utf8(x)

      if (isTRUE(fix_html_entities)) {
        x <- gsub(
          "&amp;",
          "&",
          x,
          fixed = TRUE
        )
      }
    }

    x
  })

  df
}