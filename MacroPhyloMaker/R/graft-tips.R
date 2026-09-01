# ---- small utilities --------------------------------------------------------

#' Fallback-if-null operator
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Normalize a function label from the plan
#' @keywords internal
#' @noRd
normalize_fun <- function(x) {
  tolower(gsub("[^A-Za-z]", "", x))
}

#' Parse Sister column into character vector of tip labels
#' @keywords internal
#' @noRd
parse_sister <- function(x) {
  if (is.null(x) || !nzchar(as.character(x))) {
    return(character(0L))
  }
  s <- gsub("[\r\n]", "", as.character(x))
  s <- unlist(strsplit(s, "[;,|]"))
  s <- trimws(s)
  unique(s[nzchar(s)])
}

#' Get numeric value from data.frame row with default
#' @keywords internal
#' @noRd
get_num_or <- function(row, col, default) {
  if (col %in% names(row)) {
    v <- suppressWarnings(as.numeric(row[[col]]))
    if (!is.na(v)) {
      return(v)
    }
  }
  default
}

#' Get logical value from data.frame row with default
#' @keywords internal
#' @noRd
get_logi_or <- function(row, col, default) {
  if (col %in% names(row)) {
    v <- row[[col]]
    if (is.logical(v)) {
      return(v)
    }
    if (is.numeric(v)) {
      return(v != 0)
    }
    if (is.character(v)) {
      vv <- tolower(trimws(v))
      if (vv %in% c("t", "true", "1", "yes", "y")) {
        return(TRUE)
      }
      if (vv %in% c("f", "false", "0", "no", "n")) {
        return(FALSE)
      }
    }
  }
  default
}

#' Set global seed according to a mode (fixed/random/integer)
#' @keywords internal
#' @noRd
set_global_seed <- function(seed_mode = "fixed", out_dir = getwd()) {
  if (is.numeric(seed_mode) && length(seed_mode) == 1L && is.finite(seed_mode)) {
    s <- as.integer(seed_mode)
  } else {
    sm <- tolower(as.character(seed_mode))
    if (identical(sm, "fixed")) {
      p <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
      # Convert path bytes to integers before hashing
      bytes <- as.integer(charToRaw(p))
      if (!length(bytes)) bytes <- 0L
      # Simple stable hash from bytes
      h <- sum((seq_along(bytes) %% 1024L + 1L) * bytes)
      # Map hash into a safe integer seed range
      rng_max <- .Machine$integer.max - 1e6
      s <- as.integer(1e6 + (abs(h) %% rng_max))
    } else if (identical(sm, "random")) {
      s <- as.integer(sample.int(.Machine$integer.max, 1L))
    } else {
      s <- as.integer(12345)
    }
  }
  set.seed(s)
  s
}

# ---- I/O: read trees --------------------------------------------------------

#' Read one or more phylogenies from a Newick file or text
#'
#' @param path Path to file with Newick content; can also be a single Newick string.
#' @return A list of `phylo` (length 1 for single-tree files; >1 for multi-tree files)
#' @keywords internal
#' @noRd
read_trees_any <- function(path) {
  # If path looks like Newick text (has '(' and ')' and ends with ';'), read directly
  if (is.character(path) && length(path) == 1L && grepl("[()]", path) && grepl(";\n?$", path)) {
    return(list(ape::read.tree(text = path)))
  }
  obj <- tryCatch(ape::read.tree(path), error = function(e) e)
  if (inherits(obj, "phylo")) {
    return(list(obj))
  }
  if (inherits(obj, "multiPhylo")) {
    return(unclass(obj))
  } # list of phylo
  if (is.list(obj) && all(vapply(obj, inherits, logical(1L), what = "phylo"))) {
    return(obj)
  }
  if (is.character(obj)) {
    s <- obj[nzchar(obj)]
    if (length(s) >= 1L) {
      return(list(ape::read.tree(text = s[1])))
    }
  }
  if (inherits(obj, "error")) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) >= 1L) {
      return(list(ape::read.tree(text = lines[1])))
    }
    stop("Could not read a phylogeny from file or text: ", conditionMessage(obj))
  }
  stop("Backbone file did not yield a usable phylogeny.")
}

# ---- pruning / ingroup extraction ------------------------------------------

#' Drop only the specified tip labels (ignore non-existent)
#' @keywords internal
#' @noRd
safe_drop_tips <- function(tree, tips) {
  drop <- intersect(tree$tip.label, unique(as.character(tips)))
  if (!length(drop)) {
    return(tree)
  }
  ape::drop.tip(tree, drop)
}

#' Extract ingroup defined by MRCA of available anchors (if ≥ 2 present)
#' @keywords internal
#' @noRd
extract_ingroup_by_anchors <- function(tree, anchors) {
  if (is.null(anchors) || !length(anchors)) {
    return(list(tree = tree, msg = "anchors: none requested"))
  }
  present <- intersect(anchors, tree$tip.label)
  if (length(present) < 2L) {
    return(list(
      tree = tree,
      msg = sprintf(
        "anchors: %d/%d present (%s) -> no extraction",
        length(present), length(anchors),
        if (length(present)) paste(present, collapse = ",") else "none"
      )
    ))
  }
  node <- ape::getMRCA(tree, present)
  if (is.null(node) || length(node) == 0L) {
    return(list(tree = tree, msg = "anchors: MRCA not found -> no extraction"))
  }
  ext <- ape::extract.clade(tree, node)
  list(
    tree = ext,
    msg = sprintf(
      "anchors: %d/%d present (%s) -> extracted %d tips",
      length(present), length(anchors), paste(present, collapse = ","),
      ape::Ntip(ext)
    )
  )
}

# ---- core: apply plan -------------------------------------------------------

#' Apply grafting plan to a tree from a parameter table
#'
#' \strong{Grafting modes}
#'
#' The \code{Function} column in the plan controls how taxa are grafted:
#'
#' \itemize{
#'   \item \code{"sister_to_tip"}:
#'     Attach the new taxon as a sister to a specified tip.
#'
#'   \item \code{"sister_to_clade"}:
#'     Attach the new taxon as sister to a clade defined by one or more tips.
#'
#'   \item \code{"within_clade_random"}:
#'     Insert the new taxon at a random position within a clade,
#'     optionally constrained by Beta-distributed branch position priors.
#' }
#'
#' \strong{Branch position priors}
#'
#' When placing a taxon along a branch, its position is sampled from a
#' Beta distribution defined by \code{shape1} and \code{shape2},
#' and constrained to lie between \code{min_frac} and \code{max_frac}
#' of branch length.
#'
#' @param tree A `phylo` object.
#' @param plan_df Data frame with at least columns: `Function`, `GraftedTip`, `Sister`.
#' @param default_shape1,default_shape2 Beta-shape parameters for randomization (if used by graft funcs).
#' @param default_min_frac,default_max_frac Range on edge where to attach.
#' @param default_include_terminal,default_include_stem Flags for within-clade grafting.
#' @param seed Optional integer seed.
#' @param strict If `TRUE`, stop on first error; otherwise continue with warnings.
#' @param return_log If `TRUE`, return list with `tree` and `log`; else just the tree.
#' @param verbose Logical. If `TRUE`, print progress and diagnostic messages
#'   while applying grafting operations.
#' @return `list(tree=phylo, log=data.frame)` or `phylo`.
#' @export
apply_grafts_from_table <- function(tree,
                                    plan_df,
                                    default_shape1 = 1,
                                    default_shape2 = 1,
                                    default_min_frac = 0.4,
                                    default_max_frac = 0.6,
                                    default_include_terminal = TRUE,
                                    default_include_stem = FALSE,
                                    seed = NULL,
                                    strict = TRUE,
                                    return_log = TRUE,
                                    verbose = FALSE) {
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo.")
  if (!is.null(seed)) set.seed(as.integer(seed))
  plan <- as.data.frame(plan_df, stringsAsFactors = FALSE, check.names = FALSE)
  log_list <- vector("list", nrow(plan))

  # Ensure dependent graft functions exist; if not, instruct the user
  required_funs <- c("graft_sister_to_clade", "graft_sister_to_tip", "graft_within_clade_random")
  missing <- required_funs[!vapply(required_funs, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Missing grafting function(s): ", paste(missing, collapse = ", "),
      "\nPlease make sure they are defined in your package or imported."
    )
  }

  current <- tree
  n_total <- nrow(plan)

  for (i in seq_len(nrow(plan))) {
    row <- plan[i, , drop = FALSE]
    fn_raw <- as.character(row$Function)
    if (is.null(fn_raw) || !nzchar(fn_raw)) stop("Row ", i, ": 'Function' is empty.")
    fn <- normalize_fun(fn_raw)

    new_label <- as.character(row$GraftedTip)
    if (!nzchar(new_label)) stop("Row ", i, ": 'GraftedTip' is empty.")

    sister <- parse_sister(row$Sister)
    shape1 <- get_num_or(row, "shape1", default_shape1)
    shape2 <- get_num_or(row, "shape2", default_shape2)
    min_f <- get_num_or(row, "min_frac", default_min_frac)
    max_f <- get_num_or(row, "max_frac", default_max_frac)
    inc_t <- get_logi_or(row, "include_terminal", default_include_terminal)
    inc_s <- get_logi_or(row, "include_stem", default_include_stem)

    # Pre-drop an existing tip named as GraftedTip (avoid label collision)
    pre_dropped <- FALSE
    if (new_label %in% current$tip.label) {
      current <- ape::drop.tip(current, new_label)
      pre_dropped <- TRUE
    }

    res <- try(
      {
        if (fn %in% c("graftsistertoclade")) {
          out <- graft_sister_to_clade(
            tr = current, new_label = new_label, clade_tips = sister,
            shape1 = shape1, shape2 = shape2, min_frac = min_f, max_frac = max_f,
            .capture = TRUE
          )
          current <- out$tree
          meta <- out$log
          target_txt <- paste0("(", paste(sister, collapse = ","), ")")
          mode_txt <- "graft sister to clade"
          edge_type <- NA_character_
        } else if (fn %in% c("graftsistertotip")) {
          if (length(sister) != 1L) {
            stop("Row ", i, ": 'graft sister to tip' requires exactly one Sister label.")
          }
          out <- graft_sister_to_tip(
            tr = current, new_label = new_label, sister_to = sister,
            shape1 = shape1, shape2 = shape2, min_frac = min_f, max_frac = max_f,
            .capture = TRUE
          )
          current <- out$tree
          meta <- out$log
          target_txt <- sister[1]
          mode_txt <- "graft sister to tip"
          edge_type <- NA_character_
        } else if (fn %in% c("graftwithincladerandom")) {
          out <- graft_within_clade_random(
            tr = current, new_label = new_label, clade_tips = sister,
            include_terminal = inc_t, include_stem = inc_s,
            min_frac = min_f, max_frac = max_f, return_log = TRUE
          )
          current <- out$tree
          meta <- list(
            edge_index = out$log$chosen_edge_index,
            edge_length = out$log$edge_length,
            chosen_frac = out$log$chosen_frac,
            position_from_child = out$log$position_from_child,
            attach_age = out$log$attach_age,
            where_node = out$log$where_node,
            new_tip_length = out$log$new_tip_length,
            tree_depth = out$log$tree_depth,
            edge_type = out$log$edge_type,
            is_stem_included = out$log$is_stem_included,
            is_terminal_included = out$log$is_terminal_included
          )
          target_txt <- paste0("(", paste(sister, collapse = ","), ")")
          mode_txt <- "graft within clade random"
          edge_type <- meta$edge_type %||% NA_character_
        } else {
          stop("Row ", i, ": Unrecognized Function value: ", fn_raw)
        }
        list(ok = TRUE, meta = meta, target_txt = target_txt, mode_txt = mode_txt, edge_type = edge_type)
      },
      silent = TRUE
    )

    if (inherits(res, "try-error") || !isTRUE(res$ok)) {
      msg <- paste0(
        "Row ", i, " (", fn_raw, ", new_label=", new_label, "): ",
        as.character(if (inherits(res, "try-error")) res else "unknown error")
      )
      if (strict) stop(msg)
      warning(msg)
    }

    meta <- if (!inherits(res, "try-error")) res$meta else list()
    pos_str <- sprintf("%.6g", meta$position_from_child %||% NA_real_)
    age_str <- sprintf("%.6g", meta$attach_age %||% NA_real_)
    drop_str <- if (pre_dropped) " (pre-dropped existing tip)" else ""
    etype <- if (!is.null(res$edge_type) && !is.na(res$edge_type)) paste0(" edge=", res$edge_type) else ""
    if (isTRUE(verbose)) {
      cat(sprintf(
        "[%d/%d] %s: %s  ⟂  %s  at pos=%s (attach_age=%s)%s%s\n",
        i, n_total, res$mode_txt, new_label, res$target_txt, pos_str, age_str, etype, drop_str
      ))
    }

    log_list[[i]] <- data.frame(
      step = i,
      function_raw = fn_raw,
      new_label = new_label,
      pre_dropped = pre_dropped,
      sister_input = as.character(row$Sister),
      shape1 = shape1,
      shape2 = shape2,
      min_frac = min_f,
      max_frac = max_f,
      include_terminal = if (identical(fn, "graftwithincladerandom")) inc_t else NA,
      include_stem = if (identical(fn, "graftwithincladerandom")) inc_s else NA,
      edge_index = meta$edge_index %||% NA_integer_,
      edge_length = meta$edge_length %||% NA_real_,
      chosen_frac = meta$chosen_frac %||% NA_real_,
      position_from_child = meta$position_from_child %||% NA_real_,
      attach_age = meta$attach_age %||% NA_real_,
      where_node = meta$where_node %||% NA_integer_,
      new_tip_length = meta$new_tip_length %||% NA_real_,
      tree_depth = meta$tree_depth %||% NA_real_,
      edge_type = meta$edge_type %||% NA_character_,
      reference = if ("Reference" %in% names(row)) row$Reference else NA,
      doi = if ("DOI" %in% names(row)) row$DOI else NA,
      stringsAsFactors = FALSE
    )
  }

  current <- ape::ladderize(current)
  if (return_log) {
    log_df <- do.call(rbind, log_list)
    return(list(tree = current, log = log_df))
  } else {
    return(current)
  }
}

# ---- post-process & outputs -------------------------------------------------

#' Ensure ultrametric and ladderize
#' @keywords internal
#' @noRd
pipeline_postprocess <- function(tree) {
  if (!ape::is.ultrametric(tree)) {
    tree <- phytools::force.ultrametric(tree, method = "extend")
  }
  stopifnot(ape::is.ultrametric(tree))
  ape::ladderize(tree)
}

#' Write outputs (tree, optional plot PDF, graft log, tips list)
#' Note: plotting is delegated to plot_tree_autosize() in the plotting module.
#' @keywords internal
#' @noRd
pipeline_write_outputs <- function(
  tree, log_df, out_prefix = NULL,
  plot_pdf = TRUE, pdf_width = NULL, pdf_height = NULL,
  pdf_auto = TRUE, plot_cex = 0.4, plot_fn = NULL
) {
  if (is.null(out_prefix)) {
    tree_path <- "grafted_tree.tre"
    pdf_path <- "grafted_tree.pdf"
    log_path <- "graft_log.tsv"
    tips_path <- "grafted_tree_tips.txt"
  } else {
    tree_path <- paste0(out_prefix, ".tre")
    pdf_path <- paste0(out_prefix, ".pdf")
    log_path <- paste0(out_prefix, "_graft_log.tsv")
    tips_path <- paste0(out_prefix, "_tips.txt")
    out_dir <- dirname(normalizePath(tree_path, mustWork = FALSE))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ape::write.tree(tree, file = tree_path)
  if (isTRUE(plot_pdf)) {
    w <- if (is.null(pdf_width)) 6 else as.numeric(pdf_width)
    h <- if (is.null(pdf_height)) 30 else as.numeric(pdf_height)
    if (isTRUE(pdf_auto) && (is.null(pdf_width) || is.null(pdf_height))) {
      dims <- suggest_pdf_size(tree, pdf_width = w, cex = plot_cex)
      if (is.null(pdf_height)) h <- dims$height
      if (is.null(pdf_width)) w <- dims$width
    }
    grDevices::pdf(pdf_path, width = w, height = h)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    if (is.function(plot_fn)) plot_fn(tree) else graphics::plot(ape::ladderize(tree), cex = plot_cex)
    grDevices::dev.off()
  }

  utils::write.table(log_df, file = log_path, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(data.frame(Tip = tree$tip.label),
    file = tips_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  invisible(list(
    tree = tree_path,
    pdf = if (isTRUE(plot_pdf)) pdf_path else NULL,
    log = log_path, tips = tips_path
  ))
}

# ---- top-level API ----------------------------------------------------------

#' Apply grafting plan to backbones using grafting rules
#'
#' Reads one or more backbone trees, optionally extracts an ingroup defined by anchors,
#' optionally prunes tips, then applies a TSV "plan" that specifies grafting actions.
#' Produces a ladderized ultrametric result per backbone (optionally saved).
#'
#' @param backbone_path Path to a Newick file containing a `phylo` or `multiPhylo`.
#' @param plan_path TSV with columns like `Function`, `GraftedTip`, `Sister`, etc.
#' @param out_prefix Optional prefix for outputs (per backbone).
#' @param seed_mode "fixed", "random", or an integer.
#' @param ingroup_anchors Character vector of at least two tips present in each backbone to define an MRCA.
#' @param prefer_tree_index Ignored (kept for compatibility).
#' @param plot_pdf,pdf_width,pdf_height,pdf_auto,plot_cex,plot_fn Plotting controls.
#' @return Invisibly returns `TRUE` (side effects: writes tree(s), logs, optional PDF).
#' @examples
#' \dontrun{
#' run_tip_grafting(
#'   backbone_path = system.file("extdata", "example_backbone.tre", package = "AntPhyloMaker"),
#'   plan_path = system.file("extdata", "example_graft_plan.tsv", package = "AntPhyloMaker"),
#'   out_prefix = tempfile("grafted_"),
#'   seed_mode = "fixed",
#'   ingroup_anchors = c("Martialis", "Camponotus"),
#'   plot_pdf = FALSE
#' )
#' }
#' @export
run_tip_grafting <- function(backbone_path = NULL,
                             plan_path = NULL,
                             out_prefix = NULL,
                             seed_mode = "fixed",
                             ingroup_anchors = NULL,
                             prefer_tree_index = 2L, # kept for compatibility, ignored
                             drop_tips = NULL,
                             plot_pdf = TRUE,
                             pdf_width = NULL,
                             pdf_height = NULL,
                             pdf_auto = TRUE,
                             plot_cex = 0.5,
                             plot_fn = NULL) {
  # RNG / out dir prep
  out_dir <- if (is.null(out_prefix)) {
    "."
  } else {
    p <- dirname(normalizePath(paste0(out_prefix, ".tre"), mustWork = FALSE))
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    p
  }
  run_seed <- set_global_seed(seed_mode, out_dir = out_dir)

  # Load plan (use package's UTF-8 normalizer if available)
  plan <- tryCatch(
    {
      if (exists("read_table_utf8", mode = "function")) {
        read_table_utf8(plan_path, sep = "\t", header = TRUE, quote = "")
      } else {
        utils::read.delim(plan_path, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
      }
    },
    error = function(e) stop("Failed to read plan TSV: ", conditionMessage(e))
  )

  # Load 1..N backbone trees
  backbones <- read_trees_any(backbone_path)
  nB <- length(backbones)

  for (b in seq_len(nB)) {
    cat(sprintf("\n=== Processing backbone %d/%d ===\n", b, nB))

    # (0) take backbone b
    tree1 <- backbones[[b]]

    # (1) optional ingroup extraction by anchors
    ext <- extract_ingroup_by_anchors(tree1, ingroup_anchors)
    tree1 <- ext$tree
    cat(sprintf(" - %s\n", ext$msg))

    # (2) optional drop-tips
    if (!is.null(drop_tips)) {
      drop_now <- intersect(tree1$tip.label, unique(as.character(drop_tips)))
      tree1 <- if (length(drop_now)) ape::drop.tip(tree1, drop_now) else tree1
      cat(sprintf(" - drop-tips: removed %d\n", length(drop_now)))
    }

    # (3) apply plan on THIS backbone
    res <- apply_grafts_from_table(
      tree = tree1, plan_df = plan,
      default_shape1 = 1, default_shape2 = 1,
      default_min_frac = 0.2, default_max_frac = 0.8,
      default_include_terminal = TRUE, default_include_stem = TRUE,
      seed = run_seed, strict = TRUE, return_log = TRUE
    )
    tree_final <- pipeline_postprocess(res$tree)

    # (4) outputs
    suffix <- if (nB > 1) sprintf("__t%02d", b) else ""
    op <- if (is.null(out_prefix)) paste0("grafted", suffix) else paste0(out_prefix, suffix)

    pipeline_write_outputs(tree_final, res$log,
      out_prefix = op,
      plot_pdf = plot_pdf,
      pdf_width = pdf_width, pdf_height = pdf_height,
      pdf_auto = pdf_auto, plot_cex = plot_cex, plot_fn = plot_fn
    )
  }
  invisible(TRUE)
}
