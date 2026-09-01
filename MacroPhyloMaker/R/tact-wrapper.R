# tact_wrapper.R
# MacroPhyloMaker TACT integration wrapper
#
# This file prepares inputs for TACT and calls either a Docker-based or
# system-installed TACT executable. It does not distribute or vendor TACT.

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------

#' Create a TACT workflow logger
#'
#' Creates a small timestamped log file next to the requested output prefix.
#'
#' @param out_prefix Character. Output prefix used by the TACT wrapper.
#' @param file_prefix Character. Prefix for the log filename.
#'
#' @return A logger object with a `file` element.
#' @keywords internal
.tact_make_logger <- function(out_prefix, file_prefix = "tact_grafting") {
  log_dir <- file.path(dirname(normalizePath(out_prefix, mustWork = FALSE)), "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(
    log_dir,
    paste0(file_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  cat("TACT grafting log\n", file = logfile)
  cat("Started: ", format(Sys.time()), "\n\n", file = logfile, append = TRUE, sep = "")
  structure(list(file = logfile), class = "tact_logger")
}

#' Write a message to console and logger
#'
#' Uses a package-level `log_msg()` if present, otherwise prints to console and
#' appends to the logger file.
#'
#' @param ... Components of the message.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns the message string.
#' @keywords internal
.tact_msg <- function(..., logger = NULL) {
  msg <- paste0(...)
  if (exists("log_msg", mode = "function", inherits = TRUE)) {
    get("log_msg", mode = "function", inherits = TRUE)(logger, msg)
  } else {
    cat(msg, "\n", sep = "")
    if (!is.null(logger) && !is.null(logger$file)) {
      cat(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n",
        file = logger$file, append = TRUE, sep = ""
      )
    }
  }
  invisible(msg)
}

#' Write a section heading to the TACT log
#'
#' @param title Character section title.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns `title`.
#' @keywords internal
.tact_section <- function(title, logger = NULL) {
  if (exists("log_section", mode = "function", inherits = TRUE)) {
    get("log_section", mode = "function", inherits = TRUE)(logger, title)
  } else {
    line <- paste(rep("-", nchar(title) + 4), collapse = "")
    .tact_msg(line, logger = logger)
    .tact_msg("| ", title, " |", logger = logger)
    .tact_msg(line, logger = logger)
  }
  invisible(title)
}

#' Require R packages used by the TACT wrapper
#'
#' @param pkgs Character vector of package names.
#'
#' @return Invisibly returns `TRUE` if all packages are available.
#' @keywords internal
.tact_require <- function(pkgs = c("ape", "stringr", "phytools")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Required R package(s) not installed: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Taxon-name helpers
# -----------------------------------------------------------------------------

#' Extract genus-like prefixes from tip labels
#'
#' Assumes tree labels use `Genus_epithet` style. Labels without underscores are
#' returned unchanged.
#'
#' @param labels Character vector of tip labels.
#'
#' @return Character vector of genus prefixes.
#' @keywords internal
.tact_get_genus <- function(labels) {
  has_us <- grepl("_", labels, fixed = TRUE)
  out <- labels
  out[has_us] <- sub("_.*$", "", labels[has_us])
  out
}

#' Extract epithets from tip labels
#'
#' @param labels Character vector of labels.
#'
#' @return Character vector of everything after the first underscore, or `NA`.
#' @keywords internal
.tact_get_epithet <- function(labels) {
  ifelse(grepl("_", labels, fixed = TRUE), sub("^[^_]+_", "", labels), NA_character_)
}

#' Detect genus-only labels
#'
#' @param labels Character vector of labels.
#'
#' @return Logical vector.
#' @keywords internal
.tact_is_genus_only <- function(labels) {
  !grepl("_", labels, fixed = TRUE)
}

#' Detect code-like species epithets
#'
#' Detects labels such as `Eburopone_MG04` or `Recurvidris_TH01` where the
#' epithet is uppercase letters followed by optional digits.
#'
#' @param labels Character vector of labels.
#'
#' @return Logical vector.
#' @keywords internal
.tact_is_code_epithet <- function(labels) {
  ep <- .tact_get_epithet(labels)
  !is.na(ep) & grepl("^[A-Z]{1,5}[0-9]{0,4}$", ep)
}

#' Convert underscores to spaces
#'
#' @param x Character vector.
#'
#' @return Character vector with all underscores converted to spaces.
#' @keywords internal
.tact_to_space <- function(x) gsub("_", " ", x, fixed = TRUE)

#' Convert spaces to underscores
#'
#' @param x Character vector.
#' @return Character vector with spaces converted to underscores.
#' @keywords internal
.tact_to_underscore <- function(x) gsub(" ", "_", x, fixed = TRUE)

#' Make a label unique among used labels
#'
#' @param x Proposed label.
#' @param used Character vector of already-used labels.
#'
#' @return Unique label.
#' @keywords internal
.tact_make_unique <- function(x, used) {
  if (!(x %in% used)) return(x)
  i <- 1L
  cand <- paste0(x, "_dup", i)
  while (cand %in% used) {
    i <- i + 1L
    cand <- paste0(x, "_dup", i)
  }
  cand
}

# -----------------------------------------------------------------------------
# Taxonomy readers/builders
# -----------------------------------------------------------------------------

#' Read and normalize taxonomy for TACT
#'
#' Reads AntWiki-style, simple, or TACT-style CSV taxonomy into a standard table
#' used by the wrapper. The returned table includes both a space-separated
#' species label and an underscore-separated species label.
#'
#' @param taxonomy Character path to a taxonomy file, or a data frame.
#' @param taxonomy_format One of `"antwiki"`, `"simple"`, or `"tact_csv"`.
#'   `"antwiki"` expects columns `TaxonName`, `Genus`, and `Species`.
#'   `"simple"` accepts columns `genus` and `species`, or `genus.species`.
#'   `"tact_csv"` expects columns `Family`, `genus`, and `genus.species`.
#' @param family Character family name to assign when absent.
#'
#' @return A data frame with columns `Family`, `genus`, `species`, and
#'   `species_underscore`.
#' @export
read_tact_taxonomy <- function(taxonomy,
                               taxonomy_format = c("antwiki", "simple", "tact_csv"),
                               family = "Formicidae") {
  taxonomy_format <- match.arg(taxonomy_format)

  if (is.data.frame(taxonomy)) {
    tx <- taxonomy
  } else if (taxonomy_format == "antwiki") {
    tx <- utils::read.table(
      taxonomy,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  } else {
    tx <- utils::read.table(
      taxonomy,
      header = TRUE,
      sep = if (grepl("\\.csv$", taxonomy, ignore.case = TRUE)) "," else "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }

  if (taxonomy_format == "antwiki") {
    req <- c("TaxonName", "Genus", "Species")
    miss <- setdiff(req, names(tx))
    if (length(miss)) {
      stop(
        "AntWiki taxonomy is missing required column(s): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }

    # AntWiki export columns are not fully consistent for species parsing:
    # `TaxonName` is the authoritative binomial/trinomial string, whereas some
    # `Species` cells contain a full binomial rather than an epithet. Trim all
    # character columns before any tokenization to avoid trailing-space labels.
    char_cols <- vapply(tx, is.character, logical(1))
    tx[char_cols] <- lapply(tx[char_cols], trimws)

    taxon_tokens <- strsplit(tx$TaxonName, "\\s+")
    is_species_rank <- vapply(taxon_tokens, length, integer(1)) == 2L

    sp <- tx[is_species_rank, , drop = FALSE]
    toks <- strsplit(sp$TaxonName, "\\s+")

    genus_from_taxon <- vapply(toks, `[`, character(1), 1)
    epithet_from_taxon <- vapply(toks, `[`, character(1), 2)

    genus_clean <- trimws(sp$Genus)
    empty_genus <- is.na(genus_clean) | genus_clean == ""
    genus_clean[empty_genus] <- genus_from_taxon[empty_genus]

    disagree <- genus_clean != genus_from_taxon
    if (any(disagree, na.rm = TRUE)) {
      warning(
        "Omitting AntWiki rows where Genus does not match the TaxonName genus: ",
        paste(utils::head(sp$TaxonName[disagree], 20), collapse = ", "),
        if (sum(disagree, na.rm = TRUE) > 20) " ..." else "",
        call. = FALSE
      )
    }

    keep_rows <- !disagree

    out <- data.frame(
      Family = family,
      genus = genus_clean[keep_rows],
      species = paste(genus_clean[keep_rows], epithet_from_taxon[keep_rows]),
      species_underscore = paste(genus_clean[keep_rows], epithet_from_taxon[keep_rows], sep = "_"),
      stringsAsFactors = FALSE
    )
  } else if (taxonomy_format == "tact_csv") {
    req <- c("Family", "genus", "genus.species")
    miss <- setdiff(req, names(tx))
    if (length(miss)) {
      stop(
        "TACT CSV taxonomy is missing required column(s): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }

    out <- data.frame(
      Family = tx$Family,
      genus = tx$genus,
      species = tx[["genus.species"]],
      species_underscore = .tact_to_underscore(tx[["genus.species"]]),
      stringsAsFactors = FALSE
    )
  } else {
    if (!"Family" %in% names(tx)) tx$Family <- family
    if (!"genus" %in% names(tx) && "Genus" %in% names(tx)) tx$genus <- tx$Genus

    if (!"species" %in% names(tx)) {
      if ("genus.species" %in% names(tx)) {
        tx$species <- tx[["genus.species"]]
      } else if (all(c("genus", "Species") %in% names(tx))) {
        tx$species <- paste(tx$genus, tx$Species)
      } else {
        stop("Simple taxonomy needs columns genus and species, or genus.species.", call. = FALSE)
      }
    }

    out <- data.frame(
      Family = tx$Family,
      genus = tx$genus,
      species = tx$species,
      species_underscore = .tact_to_underscore(tx$species),
      stringsAsFactors = FALSE
    )
  }

  out <- out[!is.na(out$genus) & !is.na(out$species_underscore), , drop = FALSE]
  out <- out[nchar(out$genus) > 0 & nchar(out$species_underscore) > 0, , drop = FALSE]
  out <- out[!duplicated(out$species_underscore), , drop = FALSE]
  out <- out[order(out$genus, out$species_underscore), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Write a TACT rank CSV
#'
#' Writes taxonomy as a rank-ordered CSV for `tact_build_taxonomic_tree`. The
#' species column uses spaces rather than underscores because DendroPy treats
#' unquoted Newick underscores as spaces when reading the backbone and taxonomy
#' tree. For example, `TACTTMPProcryptocerus01_carbonarius` is written as
#' `TACTTMPProcryptocerus01 carbonarius`.
#'
#' @param tax Data frame produced by [read_tact_taxonomy()] or by the wrapper's
#'   preprocessing steps.
#' @param file Character output CSV path.
#' @param family Character family name. Present for API consistency; the value
#'   in `tax$Family` is used.
#'
#' @return Invisibly returns `file`.
#' @export
write_tact_taxonomy_csv <- function(tax, file, family = "Formicidae") {
  species_label <- if ("species" %in% names(tax)) {
    tax$species
  } else {
    tax$species_underscore
  }

  species_label <- .tact_to_space(species_label)

  out <- data.frame(
    Family = tax$Family,
    genus = tax$genus,
    species = species_label,
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$Family) & !is.na(out$genus) & !is.na(out$species),
    ,
    drop = FALSE
  ]
  out <- out[
    nchar(out$Family) > 0 & nchar(out$genus) > 0 & nchar(out$species) > 0,
    ,
    drop = FALSE
  ]
  out <- out[!duplicated(out$species), , drop = FALSE]
  out <- out[order(out$Family, out$genus, out$species), , drop = FALSE]

  utils::write.table(
    out,
    file = file,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )

  invisible(file)
}

#' Write a simple taxonomy Newick tree
#'
#' This function is retained for debugging and fallback use. The recommended
#' production workflow is to write a CSV with [write_tact_taxonomy_csv()] and
#' build the taxonomy tree using TACT's `tact_build_taxonomic_tree` command.
#'
#' @param tax Data frame with taxonomy.
#' @param file Character output Newick path.
#' @param family Character root label.
#'
#' @return Invisibly returns `file`.
#' @keywords internal
write_tact_taxonomy_newick <- function(tax, file, family = "Formicidae") {
  tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  genus_groups <- split(tax$species_underscore, tax$genus)
  genus_strings <- vapply(names(genus_groups), function(g) {
    tips <- genus_groups[[g]]
    if (length(tips) == 1L) {
      paste0("(", tips, ")", g)
    } else {
      paste0("(", paste(tips, collapse = ","), ")", g)
    }
  }, character(1))
  nwk <- paste0("(", paste(genus_strings, collapse = ","), ")", family, ";")
  writeLines(nwk, file)
  invisible(file)
}

# -----------------------------------------------------------------------------
# Exclusion utilities
# -----------------------------------------------------------------------------

#' Resolve taxa and MRCA clades to exclude from TACT grafting
#'
#' Exclusions are used in two ways: tips present in the backbone are protected
#' from receiving missing taxa during non-monophyletic genus splitting; taxa not
#' present in the backbone can be removed from the TACT taxonomy by
#' `.remove_excluded_taxa_not_in_backbone()`.
#'
#' @param tree An object of class `phylo`.
#' @param exclude_taxa Optional character vector of exact tip labels to protect.
#' @param exclude_mrca Optional character vector or list of character vectors.
#'   Each vector is resolved to an MRCA and all descendant tips are protected.
#' @param exclude_missing One of `"warn"`, `"error"`, or `"ignore"`, controlling
#'   behavior when an excluded tip is absent from the backbone.
#'
#' @return Character vector of backbone tip labels to protect.
#' @export
resolve_tact_exclusions <- function(tree,
                                    exclude_taxa = NULL,
                                    exclude_mrca = NULL,
                                    exclude_missing = c("warn", "error", "ignore")) {
  exclude_missing <- match.arg(exclude_missing)
  out <- character(0)

  if (!is.null(exclude_taxa) && length(exclude_taxa)) {
    exclude_taxa <- unique(.tact_to_underscore(exclude_taxa))
    missing <- setdiff(exclude_taxa, tree$tip.label)
    if (length(missing)) {
      msg <- paste("Excluded tip(s) not found in backbone:", paste(missing, collapse = ", "))
      if (exclude_missing == "error") stop(msg, call. = FALSE)
      if (exclude_missing == "warn") warning(msg, call. = FALSE)
    }
    out <- union(out, intersect(exclude_taxa, tree$tip.label))
  }

  if (!is.null(exclude_mrca) && length(exclude_mrca)) {
    sets <- if (is.list(exclude_mrca)) exclude_mrca else list(exclude_mrca)
    for (anchors in sets) {
      anchors <- .tact_to_underscore(anchors)
      missing <- setdiff(anchors, tree$tip.label)
      if (length(missing)) {
        stop("MRCA exclusion anchor(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
      }

      node <- if (length(anchors) == 1L) {
        match(anchors, tree$tip.label)
      } else {
        ape::getMRCA(tree, anchors)
      }

      if (is.null(node) || is.na(node)) {
        stop("Could not resolve MRCA exclusion for: ", paste(anchors, collapse = ", "), call. = FALSE)
      }

      desc <- if (node <= ape::Ntip(tree)) node else phytools::getDescendants(tree, node)
      desc <- desc[desc <= ape::Ntip(tree)]
      out <- union(out, tree$tip.label[desc])
    }
  }

  sort(unique(out))
}

#' Remove excluded taxa that are absent from the backbone from taxonomy
#'
#' @param tax TACT taxonomy data frame.
#' @param exclude_taxa Character vector of taxa requested for exclusion.
#' @param backbone_tips Character vector of current backbone tip labels.
#' @param logger Optional logger object.
#'
#' @return A list with `taxonomy` and `removed` elements.
#' @keywords internal
.remove_excluded_taxa_not_in_backbone <- function(tax, exclude_taxa, backbone_tips, logger = NULL) {
  empty <- data.frame(
    species_underscore = character(),
    species = character(),
    genus = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )

  if (is.null(exclude_taxa) || !length(exclude_taxa)) {
    return(list(taxonomy = tax, removed = empty))
  }

  exclude_taxa <- unique(.tact_to_underscore(exclude_taxa))
  absent_from_backbone <- setdiff(exclude_taxa, backbone_tips)
  if (!length(absent_from_backbone)) {
    return(list(taxonomy = tax, removed = empty))
  }

  hit <- tax$species_underscore %in% absent_from_backbone
  removed <- tax[hit, , drop = FALSE]
  removed_log <- data.frame(
    species_underscore = removed$species_underscore,
    species = removed$species,
    genus = removed$genus,
    reason = "excluded_taxon_absent_from_backbone_removed_from_tact_taxonomy",
    stringsAsFactors = FALSE
  )

  if (nrow(removed_log)) {
    .tact_msg(
      "Excluded taxa removed from TACT taxonomy because absent from backbone: ",
      paste(removed_log$species_underscore, collapse = ", "),
      logger = logger
    )
  }

  list(taxonomy = tax[!hit, , drop = FALSE], removed = removed_log)
}

# -----------------------------------------------------------------------------
# Monophyly diagnostic utilities
# -----------------------------------------------------------------------------

#' Detect non-monophyletic genera for TACT preprocessing
#'
#' Uses Eddie's `nonmono_genera()` helper if it is available. Otherwise it falls
#' back to `MonoPhy::AssessMonophyly()`.
#'
#' @param tree An object of class `phylo`.
#' @param exclude_tips Optional tips to drop before testing. For TACT splitting,
#'   this should usually be `character(0)` so excluded tips can still trigger
#'   pseudo-genus creation.
#' @param genus_sep Character genus separator passed to Eddie's helper.
#'
#' @return Character vector of non-monophyletic genera.
#' @export
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#' nonmono <- detect_nonmono_genera_for_tact(tree)
#' }
detect_nonmono_genera_for_tact <- function(tree, exclude_tips = character(0), genus_sep = "_") {
  tr <- tree
  if (length(exclude_tips)) {
    tr <- ape::drop.tip(tr, intersect(exclude_tips, tr$tip.label))
  }

  if (exists("nonmono_genera", mode = "function", inherits = TRUE)) {
    return(get("nonmono_genera", mode = "function", inherits = TRUE)(tr, genus_sep = genus_sep))
  }

  if (!requireNamespace("MonoPhy", quietly = TRUE)) {
    stop("Package MonoPhy is required unless Eddie's nonmono_genera() is available.", call. = FALSE)
  }

  genera <- .tact_get_genus(tr$tip.label)
  tab <- table(genera)
  multi <- names(tab)[tab >= 2]
  tax <- data.frame(species = tr$tip.label, genus = genera, stringsAsFactors = FALSE)
  ass <- MonoPhy::AssessMonophyly(tr, tax)
  res <- MonoPhy::GetResultMonophyly(ass)
  if (!"Taxa" %in% names(res)) return(character(0))

  status_col <- intersect(c("Monophyly", "Status", "monophyly"), names(res))[1]
  if (is.na(status_col)) return(character(0))

  sort(intersect(res$Taxa[!tolower(res[[status_col]]) %in% "monophyletic"], multi))
}

#' Return descendant tip labels from a node
#'
#' @param tree An object of class `phylo`.
#' @param node Numeric node number.
#'
#' @return Character vector of descendant tip labels.
#' @keywords internal
.tact_desc_tips <- function(tree, node) {
  if (node <= ape::Ntip(tree)) return(tree$tip.label[node])
  dd <- phytools::getDescendants(tree, node)
  tree$tip.label[dd[dd <= ape::Ntip(tree)]]
}

#' Partition a non-monophyletic genus into TACT-safe singleton pseudo-genera
#'
#' The earlier MRCA-child partitioning strategy can create pseudo-genera whose
#' local placement is too broad in large trees. For TACT completion, each
#' existing graftable backbone representative is therefore treated as its own
#' temporary pseudo-genus. Missing species can then be allocated across these
#' single-tip anchors without allowing placement on unrelated high-level branches.
#'
#' @param tree An object of class `phylo`.
#' @param genus Character genus name.
#' @param exclude_tips Character vector of protected tips that should not receive
#'   missing species.
#'
#' @return A list of character vectors. Each element contains one graftable tip.
#' @keywords internal
.partition_genus_by_mrca_children <- function(tree, genus, exclude_tips = character(0)) {
  all_tips <- grep(paste0("^", genus, "_"), tree$tip.label, value = TRUE)
  graftable <- setdiff(all_tips, exclude_tips)
  if (!length(graftable)) return(list())
  as.list(graftable)
}

# -----------------------------------------------------------------------------
# Backbone preparation
# -----------------------------------------------------------------------------

#' Prepare backbone labels for TACT
#'
#' Handles genus-only tips, code-like epithets, and backbone tips that are absent
#' from the taxonomic table. Code-like tips are kept by default and added to the
#' taxonomy as known terminals so TACT does not treat them as missing.
#'
#' @param tree An object of class `phylo`.
#' @param tax Taxonomy data frame from [read_tact_taxonomy()].
#' @param genus_only One of `"replace_random_species"`, `"temporary_species"`,
#'   or `"keep"`.
#' @param species_code One of `"keep"`, `"temporary_species"`, or `"drop"`.
#' @param sample_with_replacement Logical. Whether genus-only tips can sample an
#'   already-used species if no unused valid species remain.
#' @param fallback_prefix Character prefix for temporary species labels.
#' @param seed Integer random seed.
#' @param logger Optional logger object.
#'
#' @return A list with `tree`, `taxonomy`, `label_map`, and `dropped` elements.
#' @export
prepare_tact_backbone_labels <- function(tree,
                                         tax,
                                         genus_only = c("replace_random_species", "temporary_species", "keep"),
                                         species_code = c("keep", "temporary_species", "drop"),
                                         sample_with_replacement = TRUE,
                                         fallback_prefix = "sp",
                                         seed = 1,
                                         logger = NULL) {
  genus_only <- match.arg(genus_only)
  species_code <- match.arg(species_code)
  set.seed(seed)

  labels <- tree$tip.label
  used <- labels
  map <- data.frame(old = character(), new = character(), reason = character(), stringsAsFactors = FALSE)
  dropped <- character(0)
  species_by_genus <- split(tax$species_underscore, tax$genus)

  idx_go <- which(.tact_is_genus_only(labels))
  if (length(idx_go)) {
    .tact_msg("Genus-only backbone tips detected: ", length(idx_go), logger = logger)
    for (i in idx_go) {
      g <- labels[i]
      if (genus_only == "keep") next

      if (genus_only == "replace_random_species") {
        pool <- setdiff(species_by_genus[[g]], used)
        if (length(pool) == 0L && isTRUE(sample_with_replacement)) pool <- species_by_genus[[g]]
        if (is.null(pool) || length(pool) == 0L) {
          proposal <- paste0(g, "_", fallback_prefix, sprintf("%03d", i))
          reason <- "genus_only_fallback"
        } else {
          proposal <- sample(pool, 1L)
          reason <- "genus_only_random_valid_species"
        }
      } else {
        proposal <- paste0(g, "_", fallback_prefix, sprintf("%03d", i))
        reason <- "genus_only_temporary_species"
      }

      new <- .tact_make_unique(proposal, used)
      map <- rbind(map, data.frame(old = labels[i], new = new, reason = reason, stringsAsFactors = FALSE))
      labels[i] <- new
      used <- c(used, new)
    }
  }

  idx_code <- which(.tact_is_code_epithet(labels))
  if (length(idx_code)) {
    .tact_msg("Species-code-like backbone tips detected: ", length(idx_code), logger = logger)
  }

  if (length(idx_code) && species_code == "drop") {
    dropped <- labels[idx_code]
    labels <- labels[-idx_code]
    tree <- ape::drop.tip(tree, dropped)
  } else if (length(idx_code) && species_code == "temporary_species") {
    for (i in idx_code) {
      g <- .tact_get_genus(labels[i])
      ep <- .tact_get_epithet(labels[i])
      new <- paste0(g, "_", fallback_prefix, ep)
      new <- .tact_make_unique(new, setdiff(labels, labels[i]))
      map <- rbind(map, data.frame(
        old = labels[i], new = new,
        reason = "species_code_temporary_species",
        stringsAsFactors = FALSE
      ))
      labels[i] <- new
    }
  }

  tree$tip.label <- labels

  missing_from_tax <- setdiff(tree$tip.label, tax$species_underscore)
  if (length(missing_from_tax)) {
    add <- data.frame(
      Family = tax$Family[1],
      genus = .tact_get_genus(missing_from_tax),
      species = .tact_to_space(missing_from_tax),
      species_underscore = missing_from_tax,
      stringsAsFactors = FALSE
    )
    if ("realm" %in% names(tax)) {
      add$realm <- "Unknown"
    }
    tax <- rbind(tax, add)
    tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  }

  list(tree = tree, taxonomy = tax, label_map = map, dropped = dropped)
}


# -----------------------------------------------------------------------------
# Biogeography helpers
# -----------------------------------------------------------------------------

#' Sanitize realm labels for tree tip names
#'
#' @param x Character vector.
#' @return Character vector safe to append to Newick tip labels.
#' @keywords internal
.tact_sanitize_realm <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(x), "Unknown", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "Unknown")
}

#' Read country-to-realm map
#'
#' @param country_realm_map A path or data frame with country and realm columns.
#' @param logger Optional logger.
#' @return Named character vector mapping country names to realms.
#' @keywords internal
.tact_read_country_realm_map <- function(country_realm_map, logger = NULL) {
  if (is.null(country_realm_map)) {
    return(setNames(character(0), character(0)))
  }

  rm <- if (is.data.frame(country_realm_map)) {
    country_realm_map
  } else {
    utils::read.table(
      country_realm_map,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }

  names(rm) <- trimws(names(rm))
  country_col <- intersect(c("Country", "country"), names(rm))[1]
  realm_col <- intersect(c("UdvardyRealm", "Realm", "realm", "BiogeographicRealm"), names(rm))[1]

  if (is.na(country_col) || is.na(realm_col)) {
    stop(
      "Country-realm map must contain columns `Country` and `UdvardyRealm` or `Realm`.",
      call. = FALSE
    )
  }

  country <- trimws(rm[[country_col]])
  realm <- trimws(rm[[realm_col]])
  keep <- nzchar(country) & nzchar(realm)
  country <- country[keep]
  realm <- realm[keep]

  dup <- duplicated(country)
  if (any(dup)) {
    warning(
      "Duplicate country entries in realm map; using the first occurrence for: ",
      paste(unique(country[dup]), collapse = ", "),
      call. = FALSE
    )
    country <- country[!dup]
    realm <- realm[!dup]
  }

  out <- setNames(realm, country)
  .tact_msg("Country-realm map entries loaded: ", length(out), logger = logger)
  out
}

#' Assign biogeographic realms to normalized taxonomy
#'
#' Uses the AntWiki `TypeLocalityCountry` column and an external country-realm
#' map. The function mirrors the AntWiki binomial filtering in
#' `read_tact_taxonomy()`, including omission of rows where `Genus` disagrees
#' with the genus parsed from `TaxonName`.
#'
#' @param tax Normalized taxonomy data frame.
#' @param taxonomy Original taxonomy input path or data frame.
#' @param taxonomy_format Taxonomy format.
#' @param country_realm_map Country-to-realm map path or data frame.
#' @param unknown_realm Label for missing/unmapped countries.
#' @param logger Optional logger.
#' @return Taxonomy data frame with a `realm` column.
#' @keywords internal
.tact_add_biogeo_realms <- function(tax,
                                    taxonomy,
                                    taxonomy_format,
                                    country_realm_map,
                                    unknown_realm = "Unknown",
                                    logger = NULL) {
  if ("realm" %in% names(tax)) return(tax)

  tax$realm <- unknown_realm

  if (!identical(taxonomy_format, "antwiki")) {
    .tact_msg(
      "Biogeography enabled, but taxonomy_format is not 'antwiki'; assigning all taxa to `",
      unknown_realm, "`.",
      logger = logger
    )
    return(tax)
  }

  realm_map <- .tact_read_country_realm_map(country_realm_map, logger = logger)
  if (!length(realm_map)) {
    warning("Biogeography enabled but no usable country-realm map was supplied; all realms set to Unknown.", call. = FALSE)
    return(tax)
  }

  tx <- if (is.data.frame(taxonomy)) {
    taxonomy
  } else {
    utils::read.table(
      taxonomy,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }

  req <- c("TaxonName", "Genus")
  miss <- setdiff(req, names(tx))
  if (length(miss)) {
    warning(
      "Cannot assign biogeographic realms; AntWiki taxonomy is missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
    return(tax)
  }

  char_cols <- vapply(tx, is.character, logical(1))
  tx[char_cols] <- lapply(tx[char_cols], trimws)

  # Some AntWiki exports have fewer fields per row than the full header. In
  # those files the type-locality country values can be shifted into another
  # column name, often `Author`. Choose the candidate column with the largest
  # exact overlap with the supplied country-realm map.
  country_candidates <- intersect(
    c("TypeLocalityCountry", "Author", "Year", "ChangedComb", "TrinomialAuthority"),
    names(tx)
  )
  if (!length(country_candidates)) {
    warning("Cannot assign biogeographic realms; no plausible country column found.", call. = FALSE)
    return(tax)
  }
  country_overlap <- vapply(country_candidates, function(cc) {
    sum(unique(tx[[cc]]) %in% names(realm_map))
  }, integer(1))
  country_col <- country_candidates[which.max(country_overlap)]
  if (!length(country_col) || country_overlap[[country_col]] == 0L) {
    warning("No AntWiki country column values matched the country-realm map; all realms set to Unknown.", call. = FALSE)
    return(tax)
  }
  .tact_msg("Using AntWiki column for type-locality country realm matching: ", country_col, logger = logger)

  tokens <- strsplit(tx$TaxonName, "\\s+")
  is_species_rank <- vapply(tokens, length, integer(1)) == 2L
  sp <- tx[is_species_rank, , drop = FALSE]
  toks <- strsplit(sp$TaxonName, "\\s+")
  genus_from_taxon <- vapply(toks, `[`, character(1), 1)
  epithet_from_taxon <- vapply(toks, `[`, character(1), 2)
  genus_clean <- trimws(sp$Genus)
  disagree <- genus_clean != genus_from_taxon
  sp <- sp[!disagree, , drop = FALSE]
  genus_from_taxon <- genus_from_taxon[!disagree]
  epithet_from_taxon <- epithet_from_taxon[!disagree]

  species_underscore <- paste(genus_from_taxon, epithet_from_taxon, sep = "_")
  country <- trimws(sp[[country_col]])
  realm <- unname(realm_map[country])
  realm[is.na(realm) | !nzchar(realm)] <- unknown_realm

  realm_lookup <- setNames(realm, species_underscore)
  idx <- match(tax$species_underscore, names(realm_lookup))
  hit <- !is.na(idx)
  tax$realm[hit] <- realm_lookup[idx[hit]]

  unmapped_country <- sort(unique(country[(is.na(unname(realm_map[country])) | !nzchar(unname(realm_map[country]))) & nzchar(country)]))
  if (length(unmapped_country)) {
    warning(
      "Type-locality countries not found in realm map; assigning `", unknown_realm, "`: ",
      paste(utils::head(unmapped_country, 30), collapse = ", "),
      if (length(unmapped_country) > 30) " ..." else "",
      call. = FALSE
    )
  }

  .tact_msg(
    "Biogeographic realms assigned to taxonomy rows: ",
    sum(tax$realm != unknown_realm), " / ", nrow(tax),
    logger = logger
  )
  tax
}

#' Append biogeographic realm names to tree tip labels
#'
#' @param tree A `phylo` or `multiPhylo` object.
#' @param tax Taxonomy data frame with `species_underscore` and `realm` columns.
#' @param file Output tree path.
#' @param unknown_realm Realm label used for unknown or scaffold taxa.
#' @return Invisibly returns `file`.
#' @keywords internal
.tact_write_realm_labelled_tree <- function(tree, tax, file, unknown_realm = "Unknown") {
  realm_lookup <- if ("realm" %in% names(tax)) {
    setNames(tax$realm, tax$species_underscore)
  } else {
    setNames(rep(unknown_realm, nrow(tax)), tax$species_underscore)
  }

  append_one <- function(tr) {
    lab <- tr$tip.label
    realm <- unname(realm_lookup[lab])
    realm[is.na(realm) | !nzchar(realm)] <- unknown_realm
    tr$tip.label <- paste0(lab, "__", .tact_sanitize_realm(realm))
    tr
  }

  out_tree <- if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, append_one)
    class(out) <- "multiPhylo"
    out
  } else {
    append_one(tree)
  }

  ape::write.tree(phy = out_tree, file = file)
  invisible(file)
}

# -----------------------------------------------------------------------------
# Nonmonophyletic genus splitting
# -----------------------------------------------------------------------------

#' Split non-monophyletic genera into temporary pseudo-genera for TACT
#'
#' Non-monophyletic genera cannot be safely represented as a single taxonomic
#' unit for TACT grafting. This helper temporarily renames separate backbone
#' occurrences of a non-monophyletic genus as pseudo-genera such as
#' `TACTTMPNeoponera01`. Protected excluded tips are renamed to `TACTEXCL...`
#' pseudo-genera and receive no missing species. After TACT, temporary labels are
#' removed by [restore_tact_temp_names()].
#'
#' @param tree An object of class `phylo`.
#' @param tax Taxonomy data frame.
#' @param nonmono One of `"split"`, `"skip"`, or `"error"`.
#' @param nonmono_allocation Missing-species allocation among pseudo-genera:
#'   `"proportional"`, `"equal"`, or `"random"`.
#' @param exclude_tips Character vector of backbone tips to protect from grafting.
#' @param seed Integer random seed.
#' @param logger Optional logger object.
#' @param biogeo Logical. If `TRUE`, preferentially allocate missing species to
#'   temporary anchors assigned to the same biogeographic realm.
#' @param realm_col Character. Name of the taxonomy column containing
#'   biogeographic realm assignments. If absent, the column is created and
#'   populated with `biogeo_unknown_realm`.
#' @param biogeo_unknown_realm Character. Label used for missing or unknown
#'   realm assignments.
#' @param biogeo_apply_to_all_genera Logical. If `TRUE` and `biogeo = TRUE`,
#'   create realm-aware singleton anchors for all genera represented in the
#'   backbone that also have missing species in the taxonomy, rather than only
#'   for genera detected as non-monophyletic.
#'
#' @return A list with `tree`, `taxonomy`, `map`, `nonmono`, and `skipped`.
#' @export
split_nonmono_genera_for_tact <- function(tree,
                                          tax,
                                          nonmono = c("split", "skip", "error"),
                                          nonmono_allocation = c("proportional", "equal", "random"),
                                          exclude_tips = character(0),
                                          seed = 1,
                                          logger = NULL,
                                          biogeo = FALSE,
                                          realm_col = "realm",
                                          biogeo_unknown_realm = "Unknown",
                                          biogeo_apply_to_all_genera = TRUE) {
  nonmono <- match.arg(nonmono)
  nonmono_allocation <- match.arg(nonmono_allocation)
  set.seed(seed)

  if (!(realm_col %in% names(tax))) {
    tax[[realm_col]] <- biogeo_unknown_realm
  }
  tax[[realm_col]][is.na(tax[[realm_col]]) | !nzchar(tax[[realm_col]])] <- biogeo_unknown_realm

  nonmono_g <- detect_nonmono_genera_for_tact(tree, exclude_tips = character(0))
  excluded_genera <- unique(.tact_get_genus(exclude_tips))
  nonmono_g <- sort(unique(c(
    nonmono_g,
    intersect(excluded_genera, unique(.tact_get_genus(tree$tip.label)))
  )))

  if (isTRUE(biogeo) && isTRUE(biogeo_apply_to_all_genera)) {
    tree_genera <- unique(.tact_get_genus(tree$tip.label))
    tax_genera <- unique(tax$genus)
    candidate_genera <- intersect(tree_genera, tax_genera)
    candidate_genera <- candidate_genera[vapply(candidate_genera, function(g) {
      species_in_tax <- tax$species_underscore[tax$genus == g]
      species_in_tree <- grep(paste0("^", g, "_"), tree$tip.label, value = TRUE)
      length(species_in_tree) > 0 && length(setdiff(species_in_tax, species_in_tree)) > 0
    }, logical(1))]
    nonmono_g <- sort(unique(c(nonmono_g, candidate_genera)))
  }

  if (!length(nonmono_g)) {
    return(list(tree = tree, taxonomy = tax, map = data.frame(), nonmono = data.frame(), skipped = data.frame()))
  }

  .tact_msg(
    "Non-monophyletic/biogeographic genera detected for TACT handling: ",
    paste(nonmono_g, collapse = ", "),
    logger = logger
  )
  .tact_msg("Using singleton pseudo-genus anchors for handled genera.", logger = logger)
  if (isTRUE(biogeo)) {
    .tact_msg("Biogeographic realm-aware allocation enabled.", logger = logger)
  }

  if (nonmono == "error") {
    stop("Non-monophyletic genera present: ", paste(nonmono_g, collapse = ", "), call. = FALSE)
  }

  map <- data.frame(
    original_genus = character(), temp_genus = character(), old_label = character(),
    new_label = character(), cluster = integer(), protected = logical(),
    realm = character(), stringsAsFactors = FALSE
  )
  skipped <- data.frame(
    genus = character(), species_underscore = character(), reason = character(),
    stringsAsFactors = FALSE
  )
  labels <- tree$tip.label
  realm_lookup <- setNames(tax[[realm_col]], tax$species_underscore)

  for (g in nonmono_g) {
    species_in_tax <- tax$species_underscore[tax$genus == g]
    species_in_tree <- grep(paste0("^", g, "_"), labels, value = TRUE)
    species_missing <- setdiff(species_in_tax, species_in_tree)
    protected <- intersect(species_in_tree, exclude_tips)

    if (nonmono == "skip") {
      if (length(species_missing)) {
        skipped <- rbind(skipped, data.frame(
          genus = g,
          species_underscore = species_missing,
          reason = "handled_genus_skipped",
          stringsAsFactors = FALSE
        ))
      }
      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      next
    }

    groups <- .partition_genus_by_mrca_children(tree, g, exclude_tips = exclude_tips)
    if (!length(groups)) {
      if (length(species_missing)) {
        skipped <- rbind(skipped, data.frame(
          genus = g,
          species_underscore = species_missing,
          reason = "only_excluded_branches_available",
          stringsAsFactors = FALSE
        ))
      }
      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      next
    }

    anchor_realm <- vapply(groups, function(gr) {
      rr <- realm_lookup[gr]
      rr <- rr[!is.na(rr) & nzchar(rr)]
      if (!length(rr)) return(biogeo_unknown_realm)
      names(sort(table(rr), decreasing = TRUE))[1]
    }, character(1))
    anchor_realm[is.na(anchor_realm) | !nzchar(anchor_realm)] <- biogeo_unknown_realm

    temp_genera <- paste0(
      "TACTTMP", g, "_", .tact_sanitize_realm(anchor_realm),
      sprintf("%02d", seq_along(groups))
    )

    for (k in seq_along(groups)) {
      old <- groups[[k]]
      new <- sub(paste0("^", g, "_"), paste0(temp_genera[k], "_"), old)
      labels[match(old, labels)] <- new
      map <- rbind(map, data.frame(
        original_genus = g,
        temp_genus = temp_genera[k],
        old_label = old,
        new_label = new,
        cluster = k,
        protected = FALSE,
        realm = anchor_realm[k],
        stringsAsFactors = FALSE
      ))
    }

    if (length(protected)) {
      prot_realm <- realm_lookup[protected]
      prot_realm[is.na(prot_realm) | !nzchar(prot_realm)] <- biogeo_unknown_realm
      temp_ex <- paste0(
        "TACTEXCL", g, "_", .tact_sanitize_realm(prot_realm),
        sprintf("%02d", seq_along(protected))
      )
      for (j in seq_along(protected)) {
        old <- protected[j]
        new <- sub(paste0("^", g, "_"), paste0(temp_ex[j], "_"), old)
        labels[match(old, labels)] <- new
        map <- rbind(map, data.frame(
          original_genus = g,
          temp_genus = temp_ex[j],
          old_label = old,
          new_label = new,
          cluster = NA_integer_,
          protected = TRUE,
          realm = prot_realm[j],
          stringsAsFactors = FALSE
        ))
      }
    }

    tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_in_tree), , drop = FALSE]
    g_map <- map[map$original_genus == g, , drop = FALSE]
    if (nrow(g_map)) {
      add_existing <- data.frame(
        Family = tax$Family[1],
        genus = g_map$temp_genus,
        species = .tact_to_space(g_map$new_label),
        species_underscore = g_map$new_label,
        stringsAsFactors = FALSE
      )
      add_existing[[realm_col]] <- g_map$realm
      tax <- rbind(tax, add_existing)
    }

    if (length(species_missing)) {
      missing_realm <- realm_lookup[species_missing]
      missing_realm[is.na(missing_realm) | !nzchar(missing_realm)] <- biogeo_unknown_realm
      weights <- vapply(groups, length, integer(1))
      if (nonmono_allocation == "equal") weights <- rep(1, length(groups))

      assign <- integer(length(species_missing))
      for (m in seq_along(species_missing)) {
        candidates <- seq_along(groups)
        if (isTRUE(biogeo) && !identical(missing_realm[m], biogeo_unknown_realm)) {
          same_realm <- which(anchor_realm == missing_realm[m])
          if (length(same_realm)) candidates <- same_realm
        }

        # Important: `sample()` treats a single numeric `x` as `seq_len(x)`.
        # If `candidates` has length 1 and equals, for example, 5, then
        # `sample(candidates, 1, prob = weights[candidates])` incorrectly
        # expects five probabilities. Sample over candidate positions instead.
        if (length(candidates) == 1L) {
          assign[m] <- candidates
        } else if (nonmono_allocation == "random") {
          assign[m] <- candidates[sample.int(length(candidates), 1L)]
        } else {
          assign[m] <- candidates[sample.int(
            length(candidates),
            1L,
            prob = weights[candidates]
          )]
        }
      }

      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      missing_epithet <- sub(paste0("^", g, "_"), "", species_missing)
      missing_temp_labels <- paste0(temp_genera[assign], "_", missing_epithet)

      add_missing <- data.frame(
        Family = tax$Family[1],
        genus = temp_genera[assign],
        species = .tact_to_space(missing_temp_labels),
        species_underscore = missing_temp_labels,
        stringsAsFactors = FALSE
      )
      add_missing[[realm_col]] <- missing_realm
      tax <- rbind(tax, add_missing)
    }
  }

  tree$tip.label <- labels
  tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  nonmono_log <- data.frame(
    genus = nonmono_g,
    action = if (isTRUE(biogeo)) "split_biogeo" else nonmono,
    stringsAsFactors = FALSE
  )

  list(tree = tree, taxonomy = tax, map = map, nonmono = nonmono_log, skipped = skipped)
}

#' Restore temporary TACT pseudo-genus names
#'
#' Converts temporary labels such as `TACTTMPNeoponera01_apicalis` and
#' `TACTEXCLNeoponera01_bucki` back to `Neoponera_apicalis` and
#' `Neoponera_bucki`.
#'
#' @param tree An object of class `phylo` or `multiPhylo`.
#' @param map Mapping data frame returned by [split_nonmono_genera_for_tact()].
#'
#' @return A tree object of the same class with temporary prefixes removed.
#' @export
restore_tact_temp_names <- function(tree, map) {
  if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, restore_tact_temp_names, map = map)
    class(out) <- "multiPhylo"
    return(out)
  }

  if (is.null(map) || !nrow(map)) return(tree)
  labels <- tree$tip.label
  for (i in seq_len(nrow(map))) {
    tg <- map$temp_genus[i]
    og <- map$original_genus[i]
    labels <- sub(paste0("^", tg, "_"), paste0(og, "_"), labels)
  }
  tree$tip.label <- labels
  tree
}

#' Normalize TACT output labels from spaces to underscores
#'
#' TACT/DendroPy may emit labels with spaces because underscores in unquoted
#' Newick labels are interpreted as spaces. This helper supports both `phylo`
#' and `multiPhylo` outputs.
#'
#' @param tree A `phylo` or `multiPhylo` object.
#'
#' @return A tree object of the same class with spaces converted to underscores
#'   in tip labels.
#' @keywords internal
.tact_normalize_output_labels <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, .tact_normalize_output_labels)
    class(out) <- "multiPhylo"
    return(out)
  }
  tree$tip.label <- .tact_to_underscore(tree$tip.label)
  tree
}

# -----------------------------------------------------------------------------
# TACT execution
# -----------------------------------------------------------------------------

#' Run `tact_build_taxonomic_tree`
#'
#' @param work_dir Working directory mounted into Docker or used locally.
#' @param taxonomy_csv TACT rank CSV path.
#' @param taxonomy_tree Output taxonomy tree path.
#' @param tact_runner One of `"docker"` or `"system"`.
#' @param docker_image Docker image name.
#' @param tact_build_bin Local executable for taxonomy-tree construction.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list with stdout, stderr, status, and taxonomy path.
#' @export
run_tact_build_taxonomy_external <- function(work_dir,
                                             taxonomy_csv,
                                             taxonomy_tree,
                                             tact_runner = c("docker", "system"),
                                             docker_image = "jonchang/tact",
                                             tact_build_bin = "tact_build_taxonomic_tree",
                                             logger = NULL) {
  tact_runner <- match.arg(tact_runner)

  stdout_file <- file.path(work_dir, "tact_build_taxonomy_stdout.log")
  stderr_file <- file.path(work_dir, "tact_build_taxonomy_stderr.log")

  if (tact_runner == "docker") {
    args <- c(
      "run", "--rm",
      "-v", paste0(normalizePath(work_dir), ":/workdir"),
      "-w", "/workdir",
      docker_image,
      "tact_build_taxonomic_tree",
      basename(taxonomy_csv),
      "--output", basename(taxonomy_tree)
    )

    .tact_msg("TACT taxonomy-build Docker command: docker ", paste(args, collapse = " "), logger = logger)
    status <- system2("docker", args = args, stdout = stdout_file, stderr = stderr_file)
  } else {
    args <- c(taxonomy_csv, "--output", taxonomy_tree)
    .tact_msg("TACT taxonomy-build command: ", tact_build_bin, " ", paste(args, collapse = " "), logger = logger)
    status <- system2(tact_build_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  }

  if (!identical(status, 0L)) {
    err <- if (file.exists(stderr_file)) {
      paste(tail(readLines(stderr_file, warn = FALSE), 80), collapse = "\n")
    } else {
      "<no stderr>"
    }
    stop(
      "TACT taxonomy-tree construction failed with exit status ", status,
      ".\nStderr tail:\n", err,
      "\nFull stderr: ", stderr_file,
      call. = FALSE
    )
  }

  invisible(list(stdout = stdout_file, stderr = stderr_file, status = status, taxonomy_tree = taxonomy_tree))
}

#' Run `tact_add_taxa`
#'
#' Calls TACT using either Docker or a system-installed executable.
#'
#' @param work_dir Working directory mounted into Docker or used locally.
#' @param backbone_file TACT-ready backbone tree.
#' @param taxonomy_file TACT taxonomy tree produced by `tact_build_taxonomic_tree`.
#' @param output_prefix TACT output prefix.
#' @param outgroups Optional outgroup labels to pass to TACT.
#' @param tact_runner One of `"docker"` or `"system"`.
#' @param docker_image Docker image name.
#' @param tact_bin Local `tact_add_taxa` executable.
#' @param extra_args Additional command-line arguments for TACT.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list with stdout, stderr, and status.
#' @export
run_tact_external <- function(work_dir,
                              backbone_file,
                              taxonomy_file,
                              output_prefix,
                              outgroups = NULL,
                              tact_runner = c("docker", "system"),
                              docker_image = "jonchang/tact",
                              tact_bin = "tact_add_taxa",
                              extra_args = character(0),
                              logger = NULL) {
  tact_runner <- match.arg(tact_runner)
  stdout_file <- file.path(work_dir, "tact_stdout.log")
  stderr_file <- file.path(work_dir, "tact_stderr.log")

  if (tact_runner == "docker") {
    args <- c(
      "run", "--rm",
      "-v", paste0(normalizePath(work_dir), ":/workdir"),
      "-w", "/workdir",
      docker_image,
      "tact_add_taxa",
      "--backbone", basename(backbone_file),
      "--taxonomy", basename(taxonomy_file),
      "--output", basename(output_prefix)
    )
    if (!is.null(outgroups) && length(outgroups) && all(nzchar(outgroups))) {
      args <- c(args, "--outgroups", paste(outgroups, collapse = ","))
    }
    args <- c(args, extra_args)
    .tact_msg("TACT Docker command: docker ", paste(args, collapse = " "), logger = logger)
    status <- system2("docker", args = args, stdout = stdout_file, stderr = stderr_file)
  } else {
    args <- c(
      "--backbone", backbone_file,
      "--taxonomy", taxonomy_file,
      "--output", output_prefix
    )
    if (!is.null(outgroups) && length(outgroups) && all(nzchar(outgroups))) {
      args <- c(args, "--outgroups", paste(outgroups, collapse = ","))
    }
    args <- c(args, extra_args)
    .tact_msg("TACT command: ", tact_bin, " ", paste(args, collapse = " "), logger = logger)
    status <- system2(tact_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  }

  if (!identical(status, 0L)) {
    err <- if (file.exists(stderr_file)) {
      paste(tail(readLines(stderr_file, warn = FALSE), 80), collapse = "\n")
    } else {
      "<no stderr>"
    }
    stop(
      "TACT failed with exit status ", status,
      ".\nStderr tail:\n", err,
      "\nFull stderr: ", stderr_file,
      call. = FALSE
    )
  }

  invisible(list(stdout = stdout_file, stderr = stderr_file, status = status))
}

#' Find a TACT output tree in a work directory
#'
#' TACT commonly writes `{output}.newick.tre` and `{output}.nexus.tre`.
#'
#' @param work_dir Character work directory.
#' @param output_prefix Character output prefix passed to TACT.
#'
#' @return Character path or `NA_character_`.
#' @keywords internal
.find_tact_output_tree <- function(work_dir, output_prefix) {
  stem <- basename(output_prefix)

  # Prefer Newick over Nexus. TACT writes both, and Nexus is often modified last,
  # so choosing by mtime can select the less convenient output.
  candidates <- c(
    file.path(work_dir, paste0(stem, ".newick.tre")),
    file.path(work_dir, paste0(stem, ".tre")),
    file.path(work_dir, paste0(stem, ".nwk")),
    file.path(work_dir, paste0(stem, ".newick")),
    file.path(work_dir, paste0(stem, ".nexus.tre")),
    file.path(work_dir, paste0(stem, ".nex"))
  )

  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  candidates[1]
}


#' Get tip counts for phylo or multiPhylo objects
#'
#' @param tree A `phylo`, `multiPhylo`, or other object.
#'
#' @return Integer vector of tip counts, or `NA_integer_`.
#' @keywords internal
.tact_n_tips <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    return(vapply(tree, ape::Ntip, integer(1)))
  }
  if (inherits(tree, "phylo")) {
    return(ape::Ntip(tree))
  }
  NA_integer_
}

#' Detect temporary TACT labels in phylo or multiPhylo objects
#'
#' @param tree A `phylo`, `multiPhylo`, or other object.
#'
#' @return Integer count of labels matching `TACTTMP` or `TACTEXCL`.
#' @keywords internal
.tact_n_temp_labels <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    return(vapply(tree, .tact_n_temp_labels, integer(1)))
  }
  if (!inherits(tree, "phylo")) {
    return(NA_integer_)
  }
  sum(grepl("TACTTMP|TACTEXCL", tree$tip.label))
}

#' Detect scaffold taxa that should not remain in final tree
#'
#' Detects code-like species and placeholder taxa such as `Uwari_sp`.
#'
#' @param labels Character vector of tip labels.
#'
#' @return Logical vector.
#' @keywords internal
.tact_is_scaffold_taxon <- function(labels) {
  ep <- .tact_get_epithet(labels)

  code_like <- !is.na(ep) & grepl("^[A-Z]{1,5}[0-9]{0,4}$", ep)

  placeholder <- !is.na(ep) & tolower(ep) %in% c(
    "sp", "spp", "species", "nr", "cf", "aff"
  )

  trailing_empty <- grepl("_$", labels)

  code_like | placeholder | trailing_empty
}

#' Drop scaffold taxa from TACT output
#'
#' Removes specimen-code labels and placeholder labels such as `Uwari_sp`.
#'
#' @param tree A `phylo` or `multiPhylo` object.
#'
#' @return A list with `tree` and `dropped` elements.
#' @keywords internal
.tact_drop_code_species <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    pieces <- lapply(tree, .tact_drop_code_species)
    out <- lapply(pieces, `[[`, "tree")
    class(out) <- "multiPhylo"
    dropped <- sort(unique(unlist(lapply(pieces, `[[`, "dropped"), use.names = FALSE)))
    return(list(tree = out, dropped = dropped))
  }

  if (!inherits(tree, "phylo")) {
    return(list(tree = tree, dropped = character(0)))
  }

  dropped <- tree$tip.label[.tact_is_scaffold_taxon(tree$tip.label)]
  if (length(dropped)) {
    tree <- ape::drop.tip(tree, dropped)
  }

  list(tree = tree, dropped = dropped)
}

# -----------------------------------------------------------------------------
# Main wrapper
# -----------------------------------------------------------------------------

#' Run TACT taxonomy-based grafting for a MacroPhyloMaker backbone
#'
#' Prepares a backbone tree and taxonomy for TACT, optionally handles
#' non-monophyletic genera by splitting them into temporary pseudo-genera,
#' protects excluded taxa or clades from receiving grafted species, runs TACT via
#' Docker or a system executable, and restores temporary names in the final tree.
#'
#' TACT itself is not distributed with MacroPhyloMaker. This function only
#' prepares inputs and calls external TACT commands.
#'
#' @param backbone_tree Character path to a Newick backbone tree, or a `phylo`
#'   object.
#' @param taxonomy Character path to a taxonomy file, or a data frame.
#' @param out_prefix Character output prefix. The wrapper writes audit files and
#'   a TACT work directory using this prefix.
#' @param family Character family/root label for taxonomy rows.
#' @param taxonomy_format One of `"antwiki"`, `"simple"`, or `"tact_csv"`.
#' @param taxonomy_output Currently retained for compatibility. The recommended
#'   path always writes CSV and lets TACT build the taxonomy tree.
#' @param tact_runner One of `"docker"`, `"system"`, or `"none"`. Use `"none"`
#'   to prepare files without running TACT.
#' @param docker_image Docker image name for TACT.
#' @param tact_bin Local executable name or path for `tact_add_taxa`.
#' @param outgroups Optional character vector of outgroup taxa to pass to TACT.
#' @param seed Integer random seed for stochastic preprocessing decisions.
#' @param nonmono One of `"split"`, `"skip"`, or `"error"`, controlling treatment
#'   of non-monophyletic genera.
#' @param nonmono_allocation One of `"proportional"`, `"equal"`, or `"random"`,
#'   controlling allocation of missing species among pseudo-genera.
#' @param genus_only One of `"replace_random_species"`, `"temporary_species"`, or
#'   `"keep"`, controlling genus-only backbone tips.
#' @param species_code One of `"keep"`, `"temporary_species"`, or `"drop"`,
#'   controlling code-like labels such as `Eburopone_MG04`.
#' @param drop_code_species_after_tact Logical. If `TRUE`, remove code-like and `_sp` scaffold terminals from the final TACT tree after grafting is complete.
#' @param enforce_taxonomy_tip_count Logical. If `TRUE`, stop on a final
#'   tree/taxonomy tip-count mismatch. If `FALSE` (default), issue a warning,
#'   write mismatch reports, and return the full tree.
#' @param biogeo Logical. If `TRUE`, use type-locality country mapped to
#'   biogeographic realm to allocate missing species preferentially to same-realm
#'   temporary genus anchors.
#' @param country_realm_map Optional path or data frame mapping `Country` to
#'   `UdvardyRealm` or `Realm`.
#' @param biogeo_unknown_realm Character label for missing or unmapped realm.
#' @param biogeo_apply_to_all_genera Logical. If `TRUE`, realm-aware singleton
#'   anchors are created for all genera with backbone representatives and missing
#'   species, not only genera detected as non-monophyletic.
#' @param write_biogeo_labelled_trees Logical. If `TRUE`, write realm-appended
#'   input/backbone and final cleaned trees for visual inspection.
#' @param exclude_taxa Optional character vector of exact taxa to protect from
#'   grafting. If present in the backbone, the corresponding branch is renamed to
#'   a protected `TACTEXCL...` pseudo-genus. If absent from the backbone but
#'   present in taxonomy, it is removed from TACT taxonomy.
#' @param exclude_mrca Optional clade exclusions. Supply a character vector of
#'   tip labels or a list of such vectors; each vector is resolved to an MRCA and
#'   all descendant tips are protected.
#' @param exclude_missing One of `"warn"`, `"error"`, or `"ignore"`, controlling
#'   behavior when `exclude_taxa` are not found in the backbone.
#' @param sample_with_replacement Logical. Whether genus-only labels can sample
#'   an already-used species if no unused species remain.
#' @param fallback_prefix Character prefix for temporary placeholder species.
#' @param extra_tact_args Character vector of additional arguments passed to
#'   `tact_add_taxa`.
#' @param keep_temp Logical. Whether to keep the TACT work directory after a
#'   successful run.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list of class `tact_grafting_result` containing
#'   the final tree, prepared tree and taxonomy, audit tables, paths, and
#'   parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- run_tact_grafting(
#'   backbone_tree = "project/results/grafted/final_tree.tre",
#'   taxonomy = "project/tables/antwiki-valid-species-8Mar2026.txt",
#'   out_prefix = "project/results/tact/Formicidae_complete_tact",
#'   taxonomy_format = "antwiki",
#'   tact_runner = "docker",
#'   docker_image = "jonchang/tact",
#'   outgroups = NULL,
#'   nonmono = "split",
#'   exclude_taxa = "Neoponera_bucki"
#' )
#' }
run_tact_grafting <- function(backbone_tree,
                              taxonomy,
                              out_prefix,
                              family = "Formicidae",
                              taxonomy_format = c("antwiki", "simple", "tact_csv"),
                              taxonomy_output = c("newick", "csv"),
                              tact_runner = c("docker", "system", "none"),
                              docker_image = "jonchang/tact",
                              tact_bin = "tact_add_taxa",
                              outgroups = NULL,
                              seed = 42,
                              nonmono = c("split", "skip", "error"),
                              nonmono_allocation = c("proportional", "equal", "random"),
                              genus_only = c("replace_random_species", "temporary_species", "keep"),
                              species_code = c("keep", "temporary_species", "drop"),
                              drop_code_species_after_tact = TRUE,
                              enforce_taxonomy_tip_count = FALSE,
                              biogeo = FALSE,
                              country_realm_map = NULL,
                              biogeo_unknown_realm = "Unknown",
                              biogeo_apply_to_all_genera = TRUE,
                              write_biogeo_labelled_trees = TRUE,
                              exclude_taxa = NULL,
                              exclude_mrca = NULL,
                              exclude_missing = c("warn", "error", "ignore"),
                              sample_with_replacement = TRUE,
                              fallback_prefix = "sp",
                              extra_tact_args = character(0),
                              keep_temp = TRUE,
                              logger = NULL) {
  .tact_require(c("ape", "stringr", "phytools"))

  taxonomy_format <- match.arg(taxonomy_format)
  taxonomy_output <- match.arg(taxonomy_output)
  tact_runner <- match.arg(tact_runner)
  nonmono <- match.arg(nonmono)
  nonmono_allocation <- match.arg(nonmono_allocation)
  genus_only <- match.arg(genus_only)
  species_code <- match.arg(species_code)
  exclude_missing <- match.arg(exclude_missing)

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  if (is.null(logger)) logger <- .tact_make_logger(out_prefix)

  .tact_section("TACT grafting: input", logger = logger)
  .tact_msg("Backbone tree: ", if (inherits(backbone_tree, "phylo")) "<phylo object>" else backbone_tree, logger = logger)
  .tact_msg("Taxonomy format: ", taxonomy_format, logger = logger)
  .tact_msg("TACT runner: ", tact_runner, logger = logger)

  tree <- if (inherits(backbone_tree, "phylo")) backbone_tree else ape::read.tree(backbone_tree)
  tax <- read_tact_taxonomy(taxonomy, taxonomy_format = taxonomy_format, family = family)
  if (isTRUE(biogeo)) {
    tax <- .tact_add_biogeo_realms(
      tax = tax,
      taxonomy = taxonomy,
      taxonomy_format = taxonomy_format,
      country_realm_map = country_realm_map,
      unknown_realm = biogeo_unknown_realm,
      logger = logger
    )
  } else if (!("realm" %in% names(tax))) {
    tax$realm <- biogeo_unknown_realm
  }
  target_taxonomy_species_count <- nrow(tax)

  .tact_msg("Backbone tips: ", ape::Ntip(tree), logger = logger)
  .tact_msg("Taxonomy species rows: ", nrow(tax), logger = logger)

  .tact_section("TACT grafting: label preparation", logger = logger)
  prep <- prepare_tact_backbone_labels(
    tree = tree,
    tax = tax,
    genus_only = genus_only,
    species_code = species_code,
    sample_with_replacement = sample_with_replacement,
    fallback_prefix = fallback_prefix,
    seed = seed,
    logger = logger
  )
  tree <- prep$tree
  tax <- prep$taxonomy

  exclude_tips <- resolve_tact_exclusions(
    tree,
    exclude_taxa = exclude_taxa,
    exclude_mrca = exclude_mrca,
    exclude_missing = exclude_missing
  )
  if (length(exclude_tips)) {
    .tact_msg("TACT-excluded backbone tips: ", length(exclude_tips), logger = logger)
  }

  tax_excl <- .remove_excluded_taxa_not_in_backbone(
    tax = tax,
    exclude_taxa = exclude_taxa,
    backbone_tips = tree$tip.label,
    logger = logger
  )
  tax <- tax_excl$taxonomy
  excluded_taxonomy <- tax_excl$removed
  target_taxonomy_species_count <- target_taxonomy_species_count - nrow(excluded_taxonomy)
  target_taxonomy_tips <- sort(unique(tax$species_underscore))

  backbone_realm_labelled_path <- paste0(out_prefix, "_backbone_biogeo_labels.tre")
  final_realm_labelled_path <- paste0(out_prefix, "_tacted_cleaned_biogeo_labels.tre")
  if (isTRUE(write_biogeo_labelled_trees)) {
    .tact_write_realm_labelled_tree(
      tree = tree,
      tax = tax,
      file = backbone_realm_labelled_path,
      unknown_realm = biogeo_unknown_realm
    )
    .tact_msg("Backbone tree with realm labels written: ", backbone_realm_labelled_path, logger = logger)
  }

  .tact_section("TACT grafting: non-monophyletic genera", logger = logger)
  spl <- split_nonmono_genera_for_tact(
    tree = tree,
    tax = tax,
    nonmono = nonmono,
    nonmono_allocation = nonmono_allocation,
    exclude_tips = exclude_tips,
    seed = seed,
    logger = logger,
    biogeo = biogeo,
    realm_col = "realm",
    biogeo_unknown_realm = biogeo_unknown_realm,
    biogeo_apply_to_all_genera = biogeo_apply_to_all_genera
  )
  tree_tact <- spl$tree
  tax_tact <- spl$taxonomy

  outgroups_tact <- outgroups
  if (!is.null(outgroups_tact) && length(outgroups_tact)) {
    all_map <- rbind(
      if (nrow(prep$label_map)) {
        data.frame(old_label = prep$label_map$old, new_label = prep$label_map$new)
      } else {
        data.frame(old_label = character(), new_label = character())
      },
      if (nrow(spl$map)) {
        spl$map[, c("old_label", "new_label")]
      } else {
        data.frame(old_label = character(), new_label = character())
      }
    )
    if (nrow(all_map)) {
      idx <- match(outgroups_tact, all_map$old_label)
      outgroups_tact[!is.na(idx)] <- all_map$new_label[idx[!is.na(idx)]]
    }
  }

  work_dir <- paste0(out_prefix, "_tact_work")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  backbone_file <- file.path(work_dir, "backbone_tact_ready.tre")
  taxonomy_csv <- file.path(work_dir, "taxonomy_tact_ready.csv")
  taxonomy_file <- file.path(work_dir, "taxonomy_tact_ready.tre")
  raw_output_prefix <- file.path(work_dir, "tact_output")

  ape::write.tree(tree_tact, backbone_file)
  write_tact_taxonomy_csv(tax_tact, taxonomy_csv, family = family)

  if (tact_runner != "none") {
    run_tact_build_taxonomy_external(
      work_dir = work_dir,
      taxonomy_csv = taxonomy_csv,
      taxonomy_tree = taxonomy_file,
      tact_runner = tact_runner,
      docker_image = docker_image,
      logger = logger
    )
  } else {
    .tact_msg(
      "TACT runner set to 'none'; taxonomy CSV written but taxonomy tree not built by TACT.",
      logger = logger
    )
  }

  utils::write.table(
    prep$label_map,
    paste0(out_prefix, "_label_replacements.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$map,
    paste0(out_prefix, "_nonmono_temp_name_map.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$nonmono,
    paste0(out_prefix, "_nonmono_genera.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$skipped,
    paste0(out_prefix, "_skipped_taxa.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    excluded_taxonomy,
    paste0(out_prefix, "_excluded_taxonomy.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  writeLines(
    if (length(exclude_tips)) exclude_tips else character(0),
    paste0(out_prefix, "_excluded_tips.txt")
  )

  .tact_section("TACT grafting: run external TACT", logger = logger)
  .tact_msg("TACT-ready backbone: ", backbone_file, logger = logger)
  .tact_msg("TACT-ready taxonomy CSV: ", taxonomy_csv, logger = logger)
  .tact_msg("TACT-ready taxonomy tree: ", taxonomy_file, logger = logger)

  if (tact_runner != "none") {
    run_tact_external(
      work_dir = work_dir,
      backbone_file = backbone_file,
      taxonomy_file = taxonomy_file,
      output_prefix = raw_output_prefix,
      outgroups = outgroups_tact,
      tact_runner = tact_runner,
      docker_image = docker_image,
      tact_bin = tact_bin,
      extra_args = extra_tact_args,
      logger = logger
    )
  } else {
    .tact_msg("TACT runner set to 'none'; prepared files only.", logger = logger)
  }

  raw_tree_path <- NA_character_
  final_tree <- NA
  final_path <- NA_character_
  dropped_code_species <- character(0)

  if (tact_runner != "none") {
    raw_tree_path <- .find_tact_output_tree(work_dir, raw_output_prefix)
    if (!is.na(raw_tree_path) && file.exists(raw_tree_path)) {
      raw_tree <- ape::read.tree(raw_tree_path)
      raw_tree <- .tact_normalize_output_labels(raw_tree)
      final_tree <- restore_tact_temp_names(raw_tree, spl$map)

      if (isTRUE(drop_code_species_after_tact)) {
        dropped_res <- .tact_drop_code_species(final_tree)
        final_tree <- dropped_res$tree
        dropped_code_species <- dropped_res$dropped
        if (length(dropped_code_species)) {
          .tact_msg(
            "Dropped code-like/scaffold taxa after TACT: ",
            paste(dropped_code_species, collapse = ", "),
            logger = logger
          )
        }
      }

      final_path <- paste0(out_prefix, "_tacted_cleaned.tre")
      ape::write.tree(phy = final_tree, file = final_path)
      if (isTRUE(write_biogeo_labelled_trees)) {
        .tact_write_realm_labelled_tree(
          tree = final_tree,
          tax = tax,
          file = final_realm_labelled_path,
          unknown_realm = biogeo_unknown_realm
        )
        .tact_msg("Cleaned TACT tree with realm labels written: ", final_realm_labelled_path, logger = logger)
      }

      dropped_code_species_path <- paste0(out_prefix, "_dropped_code_species.tsv")
      utils::write.table(
        data.frame(tip = dropped_code_species, reason = "code_like_or_placeholder_scaffold_taxon", stringsAsFactors = FALSE),
        dropped_code_species_path,
        sep = " ",
        row.names = FALSE,
        quote = FALSE
      )

      validation_path <- paste0(out_prefix, "_validation.tsv")
      final_tip_counts <- .tact_n_tips(final_tree)
      taxonomy_count_match <- final_tip_counts == target_taxonomy_species_count
      validation <- data.frame(
        output_tree = seq_along(final_tip_counts),
        raw_output_path = raw_tree_path,
        raw_n_tips = .tact_n_tips(raw_tree),
        cleaned_n_tips = final_tip_counts,
        remaining_temp_labels = .tact_n_temp_labels(final_tree),
        dropped_code_species = length(dropped_code_species),
        expected_taxonomy_species = target_taxonomy_species_count,
        taxonomy_tip_count_match = taxonomy_count_match,
        taxonomy_species_rows_after_tact_preprocessing = nrow(tax_tact),
        backbone_tips = ape::Ntip(tree),
        stringsAsFactors = FALSE
      )
      utils::write.table(
        validation,
        validation_path,
        sep = " ",
        row.names = FALSE,
        quote = FALSE
      )

      .tact_msg("Raw TACT output tree detected: ", raw_tree_path, logger = logger)
      .tact_msg("Cleaned TACT tree written: ", final_path, logger = logger)
      .tact_msg("Cleaned TACT tip count(s): ", paste(final_tip_counts, collapse = ", "), logger = logger)
      tree_tips_final <- if (inherits(final_tree, "multiPhylo")) {
        sort(unique(unlist(lapply(final_tree, `[[`, "tip.label"), use.names = FALSE)))
      } else {
        sort(final_tree$tip.label)
      }
      tree_not_taxonomy <- sort(setdiff(tree_tips_final, target_taxonomy_tips))
      taxonomy_not_tree <- sort(setdiff(target_taxonomy_tips, tree_tips_final))

      tree_not_taxonomy_path <- paste0(out_prefix, "_tree_not_taxonomy.tsv")
      taxonomy_not_tree_path <- paste0(out_prefix, "_taxonomy_not_tree.tsv")
      mismatch_summary_path <- paste0(out_prefix, "_taxonomy_mismatch_summary.tsv")

      utils::write.table(
        data.frame(species = tree_not_taxonomy, stringsAsFactors = FALSE),
        tree_not_taxonomy_path,
        sep = " ",
        row.names = FALSE,
        quote = FALSE
      )
      utils::write.table(
        data.frame(species = taxonomy_not_tree, stringsAsFactors = FALSE),
        taxonomy_not_tree_path,
        sep = " ",
        row.names = FALSE,
        quote = FALSE
      )
      utils::write.table(
        data.frame(
          expected_taxonomy_species = target_taxonomy_species_count,
          observed_cleaned_tree_tips = paste(final_tip_counts, collapse = ","),
          tree_not_taxonomy = length(tree_not_taxonomy),
          taxonomy_not_tree = length(taxonomy_not_tree),
          taxonomy_tip_count_match = all(taxonomy_count_match),
          stringsAsFactors = FALSE
        ),
        mismatch_summary_path,
        sep = " ",
        row.names = FALSE,
        quote = FALSE
      )

      .tact_msg("TACT validation table written: ", validation_path, logger = logger)
      .tact_msg("TACT mismatch summary written: ", mismatch_summary_path, logger = logger)
      .tact_msg("Tree-only species report written: ", tree_not_taxonomy_path, logger = logger)
      .tact_msg("Taxonomy-only species report written: ", taxonomy_not_tree_path, logger = logger)

      if (any(!taxonomy_count_match)) {
        msg <- paste0(
          "Final TACT tree tip count does not match taxonomy species count. ",
          "Expected ", target_taxonomy_species_count,
          "; observed ", paste(final_tip_counts, collapse = ", "),
          ". Validation table: ", validation_path,
          ". Mismatch summary: ", mismatch_summary_path
        )
        if (isTRUE(enforce_taxonomy_tip_count)) {
          stop(msg, call. = FALSE)
        } else {
          warning(msg, call. = FALSE)
        }
      }
    } else {
      .tact_msg(
        "TACT completed, but no output tree was detected. Check work directory: ",
        work_dir,
        logger = logger
      )
    }
  } else {
    .tact_msg(
      "TACT runner set to 'none'; no TACT output tree will be read or cleaned.",
      logger = logger
    )
  }

  if (!keep_temp && tact_runner != "none") {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }

  paths <- list(
    tact_ready_backbone = backbone_file,
    tact_ready_taxonomy_csv = taxonomy_csv,
    tact_ready_taxonomy = taxonomy_file,
    raw_tact_output_tree = raw_tree_path,
    cleaned_tree = final_path,
    label_replacements = paste0(out_prefix, "_label_replacements.tsv"),
    nonmono_map = paste0(out_prefix, "_nonmono_temp_name_map.tsv"),
    nonmono_genera = paste0(out_prefix, "_nonmono_genera.tsv"),
    skipped_taxa = paste0(out_prefix, "_skipped_taxa.tsv"),
    excluded_taxonomy = paste0(out_prefix, "_excluded_taxonomy.tsv"),
    dropped_code_species = paste0(out_prefix, "_dropped_code_species.tsv"),
    validation = paste0(out_prefix, "_validation.tsv"),
    biogeo_taxonomy_realms = paste0(out_prefix, "_biogeo_taxonomy_realms.tsv"),
    backbone_biogeo_labels = backbone_realm_labelled_path,
    cleaned_tree_biogeo_labels = final_realm_labelled_path,
    taxonomy_mismatch_summary = paste0(out_prefix, "_taxonomy_mismatch_summary.tsv"),
    tree_not_taxonomy = paste0(out_prefix, "_tree_not_taxonomy.tsv"),
    taxonomy_not_tree = paste0(out_prefix, "_taxonomy_not_tree.tsv"),
    excluded_tips = paste0(out_prefix, "_excluded_tips.txt"),
    log = logger$file,
    work_dir = work_dir
  )

  res <- list(
    tree = final_tree,
    tact_ready_tree = tree_tact,
    tact_ready_taxonomy = tax_tact,
    label_map = prep$label_map,
    nonmono_map = spl$map,
    nonmono_genera = spl$nonmono,
    skipped_taxa = spl$skipped,
    excluded_taxonomy = excluded_taxonomy,
    dropped_code_species = dropped_code_species,
    excluded_tips = exclude_tips,
    paths = paths,
    parameters = list(
      seed = seed,
      nonmono = nonmono,
      nonmono_allocation = nonmono_allocation,
      genus_only = genus_only,
      species_code = species_code,
      drop_code_species_after_tact = drop_code_species_after_tact,
      enforce_taxonomy_tip_count = enforce_taxonomy_tip_count,
      biogeo = biogeo,
      country_realm_map = country_realm_map,
      biogeo_unknown_realm = biogeo_unknown_realm,
      biogeo_apply_to_all_genera = biogeo_apply_to_all_genera,
      write_biogeo_labelled_trees = write_biogeo_labelled_trees,
      tact_runner = tact_runner,
      docker_image = docker_image,
      taxonomy_output = taxonomy_output
    )
  )
  class(res) <- c("tact_grafting_result", class(res))
  invisible(res)
}
