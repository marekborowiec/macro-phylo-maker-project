# Lightweight logging utilities (internal)
#
# A minimal logger used across the package to:
# - create a per-run log file,
# - mirror log messages to the console,
# - write titled section headers for readability.
#
# Typical workflow:
# - Call `.make_logger()` once at the start of a task.
# - Use `log_msg()` for one-line messages.
# - Use `log_section()` to emit a boxed section title.
# - Call `.close_logger()` when finished (typically via `on.exit()`).
#
# Logging behavior:
# - `enable` controls whether messages are written to a log file.
# - `console` controls whether messages are mirrored to the console.
# - When `enable = FALSE`, a logger object is still created and a
#   logfile path is reserved, but no file is opened or written.
#
# Box-drawing characters use UTF-8 by default. To force ASCII output, set
# either:
# - `options(MacroPhyloMaker.ascii = TRUE)` or
# - `options(AntPhyloMaker.ascii = TRUE)` (backward-compatible alias).

# -------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------

#' Resolve the package ASCII logger option (with backward-compatible alias)
#' @keywords internal
#' @noRd
.ascii_logging_enabled <- function() {
  getOption(
    "MacroPhyloMaker.ascii",
    getOption("AntPhyloMaker.ascii", FALSE)
  )
}

#' Should this logger mirror to console?
#' @keywords internal
#' @noRd
.logger_wants_console <- function(logger, .console = NULL) {
  if (!is.null(.console)) {
    return(isTRUE(.console))
  }

  if (inherits(logger, "smart_logger")) {
    return(isTRUE(logger$console))
  }

  # Non-logger input (e.g. NULL) defaults to console output
  TRUE
}

#' Should this logger write to file?
#' @keywords internal
#' @noRd
.logger_wants_file <- function(logger) {
  inherits(logger, "smart_logger") &&
    isTRUE(logger$enabled) &&
    inherits(logger$con, "connection") &&
    isOpen(logger$con)
}

# -------------------------------------------------------------------
# Box-drawing character set
# -------------------------------------------------------------------

#' Box-drawing characters with ASCII fallback
#' @keywords internal
#' @noRd
.box_chars <- function(ascii = .ascii_logging_enabled()) {
  if (isTRUE(ascii)) {
    list(
      h = "-",
      tl = "+",
      tr = "+",
      bl = "+",
      br = "+",
      v = "|"
    )
  } else {
    list(
      h = "\u2500",
      tl = "\u250c",
      tr = "\u2510",
      bl = "\u2514",
      br = "\u2518",
      v = "\u2502"
    )
  }
}

# -------------------------------------------------------------------
# Create/open a logger
# -------------------------------------------------------------------

#' Create a file-backed smart logger
#'
#' @details
#' Creates a per-run log file path in \code{outdir} and returns a
#' \code{"smart_logger"} object containing:
#' \itemize{
#'   \item the planned logfile path,
#'   \item an open connection when \code{enable = TRUE},
#'   \item the selected \code{mode},
#'   \item the logger name and header,
#'   \item console/file logging flags.
#' }
#'
#' If \code{enable = TRUE}, the logfile is opened immediately and the session
#' header is written once to the file. If \code{console = TRUE}, the same
#' header is also printed to the console.
#'
#' If \code{enable = FALSE}, no file is opened or written, but the log object
#' still records the intended logfile path. This is useful when downstream code
#' wants a stable log-derived file stem without actually writing a log.
#'
#' The log filename stem is either \code{file_prefix} or a mode-derived default,
#' followed by a timestamp in \code{YYYYMMDD-HHMMSS} format.
#'
#' @keywords internal
#' @noRd
.make_logger <- function(outdir,
                         genus,
                         enable = TRUE,
                         console = TRUE,
                         mode = c("generic", "extract", "graft"),
                         file_prefix = NULL) {
  mode <- match.arg(mode)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")

  stem <- if (!is.null(file_prefix) && nzchar(file_prefix)) {
    file_prefix
  } else {
    switch(mode,
      extract = sprintf("%s_extract", genus),
      graft   = sprintf("%s_graft", genus),
      generic = sprintf("%s_log", genus)
    )
  }

  logfile <- file.path(outdir, sprintf("%s_%s.log", stem, ts))
  con <- NULL

  title <- switch(mode,
    extract = "Genus extraction",
    graft   = "Genus grafting",
    generic = "MacroPhyloMaker"
  )

  header <- sprintf(
    "=== %s log: %s | %s ===",
    title, genus, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  if (isTRUE(enable)) {
    con <- file(logfile, open = "a", encoding = "UTF-8")
    writeLines(header, con)
    flush(con)
  }

  if (isTRUE(console)) {
    cat(header, "\n")
  }

  structure(
    list(
      enabled = isTRUE(enable),
      console = isTRUE(console),
      file    = logfile,
      con     = con,
      mode    = mode,
      name    = genus,
      header  = header
    ),
    class = "smart_logger"
  )
}

# -------------------------------------------------------------------
# Close a logger safely
# -------------------------------------------------------------------

#' Close a smart logger safely
#'
#' @details
#' Closes the logger's file connection if it exists and is open. Safe to call
#' repeatedly or on loggers with no active file connection.
#'
#' @keywords internal
#' @noRd
.close_logger <- function(logger) {
  if (inherits(logger, "smart_logger") &&
      inherits(logger$con, "connection")) {
    try(
      {
        if (isOpen(logger$con)) {
          close(logger$con)
        }
      },
      silent = TRUE
    )
  }
  invisible(NULL)
}

# -------------------------------------------------------------------
# Write a one-line message
# -------------------------------------------------------------------

#' Write a one-line log message
#'
#' @details
#' Concatenates all supplied components into a single string. If
#' \code{.timestamp = TRUE}, the line is prefixed with the current time in
#' \code{HH:MM:SS} format.
#'
#' Messages are written to the console when console mirroring is enabled, and
#' to the logfile when file logging is enabled. For non-logger inputs
#' (e.g. \code{NULL}), output defaults to the console only.
#'
#' @keywords internal
#' @noRd
log_msg <- function(logger, ..., .timestamp = FALSE, .console = NULL) {
  line <- paste0(..., collapse = "")

  if (isTRUE(.timestamp)) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), line)
  }

  if (.logger_wants_console(logger, .console = .console)) {
    cat(line, "\n")
  }

  if (.logger_wants_file(logger)) {
    writeLines(line, logger$con)
    flush(logger$con)
  }

  invisible(NULL)
}

# -------------------------------------------------------------------
# Draw a titled section box
# -------------------------------------------------------------------

#' Write a boxed section title
#'
#' @details
#' Produces a three-line boxed section header using either UTF-8 box-drawing
#' characters or ASCII fallbacks from \code{.box_chars()}. Section output is
#' routed through \code{log_msg()}, so it respects the current console/file
#' logging settings.
#'
#' @keywords internal
#' @noRd
log_section <- function(logger, title) {
  chars <- .box_chars()
  bar <- paste(rep(chars$h, max(10, nchar(title) + 2L)), collapse = "")

  log_msg(logger, "")
  log_msg(logger, chars$tl, bar, chars$tr)
  log_msg(logger, chars$v, " ", title, " ", chars$v)
  log_msg(logger, chars$bl, bar, chars$br)

  invisible(NULL)
}