#' Grafting utilities for clade‑ and tip‑level integration into a backbone chronogram
#'
#' This module provides utilities to graft donor phylogenies onto an ultrametric
#' backbone chronogram using either **clade (MRCA‑based)** placement or
#' **terminal‑tip** replacement. Grafting behavior is fully specified by a
#' tabular plan that associates each donor phylogeny with a placement target
#' and is designed to support reproducible megaphylogeny assembly from
#' heterogeneous published sources.
#'
#' Placement targets are interpreted as follows:
#' \itemize{
#'   \item A placement token containing multiple labels (comma‑separated)
#'         specifies a **clade graft**, identifying the MRCA of those tips
#'         in the backbone.
#'   \item A placement token containing a single label specifies a
#'         **terminal‑tip graft**, replacing the branch leading to that tip
#'         in the backbone.
#' }
#'
#' Donor trees are processed according to an explicit plan column named
#' \code{stem_mode}, which controls whether donor stem geometry is used.
#'
#' \describe{
#' \item{\code{stem_mode = "outgroup"}}{Donor has a meaningful stem above the ingroup
#' crown; stem fraction \eqn{r} is inferred and used (default, backward‑compatible).}
#' \item{\code{stem_mode = "crown"}}{Donor is ingroup‑only or stem depth is not trusted;
#' grafting is performed using crown‑only geometry.}
#' }
#'
#' Missing, empty, or \code{NA} values default to \code{"outgroup"}.
#'
#' The donor crown subtree is then rescaled and
#' grafted according to the placement mode:
#'
#' \describe{
#'   \item{Clade (MRCA) grafting}{
#'     The donor crown is inserted along the incoming edge to the target MRCA
#'     such that the backbone stem above the splice equals
#'     \eqn{r \times (L_{in} + C_{backbone})}, and the donor crown height equals
#'     \eqn{(1 - r) \times (L_{in} + C_{backbone})}, where \eqn{L_{in}} is the
#'     length of the MRCA’s incoming edge and \eqn{C_{backbone}} is the backbone
#'     crown depth of the clade. Existing descendant tips of the MRCA are removed
#'     after grafting.
#'   }
#'   \item{Terminal‑tip grafting}{
#'     The donor crown replaces the terminal branch leading to the specified tip.
#'     The donor is rescaled so that its crown height equals
#'     \eqn{(1 - r) \times L}, where \eqn{L} is the terminal branch length, and is
#'     grafted onto that edge followed by removal of the original tip. This
#'     operation is implemented as terminal‑edge replacement rather than node
#'     insertion, ensuring stable branch‑length geometry.
#'   }
#' }
#'
#' Chronogram conversion of donor trees is performed using
#' \pkg{ape}::\code{chronos()} when enabled. Converged chronos fits are preferred,
#' but if no model converges while still producing finite likelihoods, the best
#' non‑converged likelihood solution is retained with an explicit warning.
#' Non‑parametric rate smoothing (NNLS) is used only as a last resort when all
#' likelihood‑based fits fail.
#'
#' @section Dependencies:
#' These functions rely on \pkg{ape} and \pkg{phytools}. The surrounding workflow
#' assumes the availability of project‑level utilities such as
#' \code{read_trees_any()}, \code{extract_ingroup_by_anchors()},
#' \code{safe_drop_tips()}, \code{set_global_seed()},
#' \code{plot_tree_autosize()}, \code{log_msg()}, and \code{log_section()}.
#'
#' @keywords phylogenetics grafting chronogram ultrametric
#' @name graft_clades
NULL

# -----------------------------------------------------------------------------
# Utils (internal) -------------------------------------------------------------
# -----------------------------------------------------------------------------

#' Graft a donor crown onto a terminal backbone edge
#'
#' Internal helper to graft a donor subtree onto the terminal branch leading
#' to a single backbone tip. The donor is interpreted as a \emph{crown-only}
#' subtree (scaled to unit height), and is rescaled and attached such that its
#' crown depth equals \eqn{(1 - r) \times L}, where \eqn{L} is the length of the
#' terminal backbone edge and \eqn{r} is the donor stem fraction inferred from
#' the source phylogeny.
#'
#' Conceptually, this operation \emph{replaces the terminal edge} rather than
#' inserting a clade into the backbone: the original tip is dropped after
#' grafting, and the donor crown becomes the sole descendant of the former
#' parent node. This mirrors the behavior used in genus-level tip grafting and
#' avoids edge-collapse issues common to node-based insertion on terminal edges.
#'
#' Placement geometry:
#' \itemize{
#'   \item Terminal edge length: \eqn{L}
#'   \item Donor crown height: \eqn{(1 - r) \times L}
#'   \item Graft position: \eqn{(1 - r) \times L} measured from the tip towards
#'         the parent node (as required by \code{ape::bind.tree()} semantics)
#' }
#'
#' A small epsilon is used to clamp graft positions away from the exact endpoints
#' of the edge to prevent zero-length branches.
#'
#' @param X A backbone \code{phylo} object (assumed ultrametric).
#' @param tip Character scalar giving the name of the terminal tip to be replaced.
#' @param donor_unit A \code{phylo} donor crown subtree scaled to unit height
#'   (rooted at the donor crown; no outgroup).
#' @param r Numeric donor stem fraction in \eqn{(0, 1)}.
#' @param eps Numeric epsilon used to clamp graft positions (default: \code{1e-10}).
#'
#' @details
#' This function assumes that donor stem length has already been inferred using
#' an explicit outgroup during donor template preparation. As such, \code{r}
#' must reflect a meaningful stem-to-total ratio for the donor clade.
#'
#' The helper intentionally avoids explicit edge splitting and instead relies on
#' terminal-edge replacement via \code{bind.tree()} followed by tip removal.
#' This design preserves intended branch-length geometry under ape's tree
#' normalization rules.
#'
#' @return A \code{phylo} tree with the specified tip replaced by the donor crown.
#'
#' @keywords internal grafting chronogram terminal-edge
#' @noRd
.graft_on_terminal_edge <- function(X, tip, donor_unit, r, eps = 1e-10) {
  H <- phytools::nodeHeights(X)
  tip_idx <- match(tip, X$tip.label)
  ei <- which(X$edge[, 2] == tip_idx)

  L <- H[ei, 2] - H[ei, 1]
  pos <- (1 - r) * L

  pos_eps <- max(eps, 1e-6 * L)
  pos <- min(max(pos, pos_eps), L - pos_eps)

  donor_scaled <- .rescale_to_height(donor_unit, pos)

  X2 <- ape::bind.tree(
    X,
    donor_scaled,
    where = tip_idx,
    position = pos
  )

  # Always drop the backbone tip by *index* (never by label)
  out <- ape::drop.tip(X2, tip_idx)

  .ensure_positive_edges(out)
}
#' Sample an attachment fraction from a Beta distribution (internal)
#' @keywords internal
#' @noRd
.sample_beta_frac <- function(shape1, shape2, min_frac, max_frac) {
  shape1 <- max(1e-8, as.numeric(shape1))
  shape2 <- max(1e-8, as.numeric(shape2))
  min_frac <- max(0, min(1, as.numeric(min_frac)))
  max_frac <- max(0, min(1, as.numeric(max_frac)))
  if (max_frac < min_frac) {
    tmp <- min_frac
    min_frac <- max_frac
    max_frac <- tmp
  }
  u <- stats::rbeta(1L, shape1 = shape1, shape2 = shape2)
  min_frac + u * (max_frac - min_frac)
}

#' Graft donor crown onto a terminal backbone edge at a Beta-random position
#' (internal; used when donor has no meaningful stem fraction)
#' @keywords internal
#' @noRd
.graft_on_terminal_edge_beta <- function(
  X, tip, donor_unit,
  shape1 = 1, shape2 = 1,
  min_frac = 0.1, max_frac = 0.9,
  eps = 1e-10
) {
  H <- phytools::nodeHeights(X)
  tip_idx <- match(tip, X$tip.label)
  ei <- which(X$edge[, 2] == tip_idx)
  if (!length(ei)) stop("Tip not found in edge matrix: ", tip)

  L <- H[ei, 2] - H[ei, 1]
  if (!is.finite(L) || L <= 0) stop("Non-positive terminal edge length for tip: ", tip)

  frac <- .sample_beta_frac(shape1, shape2, min_frac, max_frac)
  pos <- frac * L

  # clamp away from endpoints
  pos_eps <- max(eps, 1e-6 * L)
  pos <- min(max(pos, pos_eps), L - pos_eps)

  donor_scaled <- .rescale_to_height(donor_unit, pos)
  X2 <- ape::bind.tree(
    X,
    donor_scaled,
    where = tip_idx,
    position = pos
  )

  # Always drop the backbone tip by *index* (never by label)
  out <- ape::drop.tip(X2, tip_idx)

  .ensure_positive_edges(out)
}

#' Crown-depth clade replacement: drop old clade tips and bind rescaled donor at backbone MRCA node depth
#' @keywords internal
#' @noRd
# Crown-mode MRCA replacement via terminal-edge grafting at crown depth
# @keywords internal
# @noRd
.graft_replace_clade_at_crown_depth <- function(
  X, mrca, donor_unit, eps = 1e-10
) {
  ## Measure backbone crown depth
  H <- phytools::nodeHeights(X)
  root_age <- max(H[, 2])
  nh <- phytools::nodeheight(X, mrca)
  crown_depth <- root_age - nh

  if (!is.finite(crown_depth) || crown_depth <= 0) {
    stop("Invalid crown depth at MRCA.")
  }

  ## Collapse MRCA clade to exactly one real tip
  collapsed <- .collapse_mrca_to_existing_tip(X, mrca)
  Xc <- collapsed$tree
  keep_tip <- collapsed$tip

  ## Get terminal edge length of retained tip
  Hc <- phytools::nodeHeights(Xc)
  tip_idx <- match(keep_tip, Xc$tip.label)
  ei <- which(Xc$edge[, 2] == tip_idx)
  L <- Hc[ei, 2] - Hc[ei, 1]

  if (crown_depth >= L) {
    stop("Crown depth exceeds retained terminal edge length.")
  }

  ## Graft position: exactly crown_depth from the tip
  pos <- crown_depth
  pos_eps <- max(eps, 1e-6 * L)
  pos <- min(max(pos, pos_eps), L - pos_eps)

  ## Scale donor crown to match backbone crown depth
  donor_scaled <- .rescale_to_height(donor_unit, crown_depth)

  ## Bind donor and remove retained backbone tip
  X2 <- ape::bind.tree(
    Xc,
    donor_scaled,
    where = tip_idx,
    position = pos
  )
  out <- ape::drop.tip(X2, keep_tip)
  out <- .ensure_positive_edges(out)

  list(
    tree = out,
    retained_tip = keep_tip
  )
}
#' Ensure strictly positive edge lengths (internal)
#'
#' Replaces non-finite or non-positive edge lengths with a small epsilon
#' computed from the minimum positive edge.
#'
#' @param tr A \code{phylo} object.
#' @param base_eps Numeric; minimum epsilon used if no positive edges exist.
#' @return The \code{phylo} with non-positive edges nudged to epsilon.
#' @keywords internal
#' @noRd
.ensure_positive_edges <- function(tr, base_eps = 1e-12) {
  if (is.null(tr$edge.length)) {
    return(tr)
  }
  el <- tr$edge.length
  pos <- el[is.finite(el) & el > 0]
  eps <- if (length(pos)) max(base_eps, 1e-8 * min(pos)) else base_eps
  el[!is.finite(el) | el <= 0] <- eps
  tr$edge.length <- el
  tr
}
#' Collapse an MRCA clade into a single terminal tip (internal)
#'
#' Drops all descendant tips of a specified MRCA node except for one
#' terminal. This is used to convert clade grafting into terminal-edge
#' grafting.
#'
#' @param X A \code{phylo} backbone tree (assumed ultrametric).
#' @param mrca Integer node number of the MRCA to collapse.
#' @param label Character label for the synthetic tip.
#' @return A modified \code{phylo} with the clade collapsed to one tip.
#' @keywords internal
#' @noRd
.collapse_mrca_to_existing_tip <- function(X, mrca) {
  desc <- phytools::getDescendants(X, mrca)
  tip_desc <- desc[desc <= ape::Ntip(X)]
  labels <- X$tip.label[tip_desc]

  keep <- sort(labels)[1] # deterministic choice
  drop <- setdiff(labels, keep)

  X2 <- ape::drop.tip(X, drop)

  list(tree = X2, tip = keep)
}
#' Try multiple methods to force ultrametricity (internal)
#'
#' Attempts \code{phytools::force.ultrametric()} with methods in order until
#' ultrametric or applies a minimal nudge for tips.
#'
#' @param tr A \code{phylo} object.
#' @param method_order Character vector of methods (default: c("nnls","extend")).
#' @return An ultrametric \code{phylo}.
#' @keywords internal
#' @noRd
.safe_force_ultrametric <- function(tr, method_order = c("nnls", "extend")) {
  if (ape::is.ultrametric(tr)) {
    return(tr)
  }
  for (m in method_order) {
    tt <- try(phytools::force.ultrametric(tr, method = m), silent = TRUE)
    if (!inherits(tt, "try-error") && ape::is.ultrametric(tt)) {
      return(tt)
    }
  }
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
  eps <- max(1e-8, 1e-6 * max(tr$edge.length, na.rm = TRUE))
  te <- which(tr$edge[, 2] <= ape::Ntip(tr))
  tr$edge.length[te] <- tr$edge.length[te] + eps
  tr
}

#' Rescale a tree to a specific height (internal)
#'
#' @param tr A \code{phylo} object.
#' @param target_height Positive numeric target height.
#' @return The \code{phylo} rescaled so its height equals \code{target_height}.
#' @keywords internal
#' @noRd
.rescale_to_height <- function(tr, target_height) {
  H <- suppressWarnings(max(phytools::nodeHeights(tr)))
  if (!is.finite(H) || H <= 0) {
    if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
    H <- sum(tr$edge.length)
  }
  if (!is.finite(H) || H <= 0) H <- 1
  tr$edge.length <- tr$edge.length / H * max(1e-8, target_height)
  tr
}
#' Compute stem/crown metrics for an ingroup (internal)
#'
#' Given a rooted tree and an ingroup tip set, returns the MRCA crown node,
#' root age, stem length (root->crown), crown depth (crown->tips), total depth,
#' and the clamped stem fraction r = stem/total.
#'
#' @param TREE A \code{phylo}.
#' @param ig Character vector of ingroup tip labels.
#' @return A list with: ok, reason, crown_node, root_age, nh_crown, stem_len, crown_depth, total, r
#' @keywords internal
#' @noRd
.compute_crown_metrics <- function(TREE, ig) {
  cr <- ape::getMRCA(TREE, ig)
  if (is.na(cr)) {
    return(list(ok = FALSE, reason = "No ingroup MRCA"))
  }

  H <- phytools::nodeHeights(TREE)
  root_age <- max(H[, 2])
  nh <- phytools::nodeheight(TREE, cr)

  stem_len <- nh
  crown_depth <- root_age - nh
  total <- stem_len + crown_depth

  if (!is.finite(total) || total <= 0) {
    return(list(ok = FALSE, reason = "Non-positive donor total depth"))
  }

  r <- min(max(stem_len / total, 1e-8), 1 - 1e-8)

  list(
    ok = TRUE,
    crown_node = cr,
    root_age = root_age,
    nh_crown = nh,
    stem_len = stem_len,
    crown_depth = crown_depth,
    total = total,
    r = r
  )
}
#' Infer a singleton outgroup if present (internal)
#'
#' Tries to find a root child that is a single tip (common in donor trees with
#' one explicit outgroup).
#'
#' @param tr A \code{phylo}.
#' @return A tip label (character) or \code{NA_character_}.
#' @keywords internal
#' @noRd
.infer_single_outgroup <- function(tr) {
  N <- ape::Ntip(tr)
  roots <- setdiff(unique(tr$edge[, 1]), tr$edge[, 2])
  root <- roots[1]
  kids <- tr$edge[tr$edge[, 1] == root, 2]
  if (length(kids) == 2) {
    size <- function(node) if (node <= N) 1L else length(ape::extract.clade(tr, node)$tip.label)
    n1 <- size(kids[1])
    n2 <- size(kids[2])
    if (xor(n1 == 1L, n2 == 1L)) {
      return(if (n1 == 1L) tr$tip.label[kids[1]] else tr$tip.label[kids[2]])
    }
  }
  NA_character_
}
#' Fallback: choose exactly one outgroup tip if the root-singleton heuristic fails (internal)
#'
#' Strategy:
#'   (a) If a root child is a singleton tip, return it (same as primary heuristic).
#'   (b) Otherwise, choose the tip with the largest mean cophenetic distance to all other tips
#'       (a common “farthest-from-centroid” proxy for a single outgroup).
#'
#' @param tr A phylo tree (assumed rooted).
#' @return A single tip label or NA_character_ if something is wrong.
#' @keywords internal
#' @noRd
.fallback_single_outgroup <- function(tr) {
  N <- ape::Ntip(tr)
  roots <- setdiff(unique(tr$edge[, 1]), tr$edge[, 2])
  if (!length(roots)) {
    return(NA_character_)
  }
  root <- roots[1]
  kids <- tr$edge[tr$edge[, 1] == root, 2]

  # (a) root child singleton (mirror .infer_single_outgroup)
  if (length(kids) == 2) {
    size <- function(node) if (node <= N) 1L else length(ape::extract.clade(tr, node)$tip.label)
    n1 <- size(kids[1])
    n2 <- size(kids[2])
    if (xor(n1 == 1L, n2 == 1L)) {
      return(if (n1 == 1L) tr$tip.label[kids[1]] else tr$tip.label[kids[2]])
    }
  }

  # (b) farthest-from-centroid tip as a single outgroup proxy
  D <- ape::cophenetic.phylo(tr)
  if (!is.matrix(D) || !nrow(D)) {
    return(NA_character_)
  }
  diag(D) <- NA_real_
  means <- rowMeans(D, na.rm = TRUE)
  out <- names(which.max(means))
  if (length(out)) out[1] else NA_character_
}
#' Resolve authority data (Genus+Species or Binomial)
#'
#' Accepts a data.frame, a TSV path, or a character vector of binomials and
#' returns a unique character vector of valid \code{Genus_species} names.
#'
#' @param authority Data frame with columns \code{Genus} and \code{Species} or
#'   column \code{Binomial}; OR a single TSV path to such a table; OR a character
#'   vector of already-resolved binomials; OR \code{NULL}.
#' @return A character vector of unique binomials, or \code{NULL} if input is
#'   \code{NULL}.
#' @export
resolve_authority_binomials <- function(authority) {
  if (is.null(authority)) {
    return(NULL)
  }
  if (is.character(authority) && length(authority) > 1L) {
    return(unique(authority[nzchar(authority)]))
  }
  if (is.data.frame(authority)) {
    cn <- names(authority)
    if (all(c("Genus", "Species") %in% cn)) {
      out <- paste(authority$Genus, authority$Species, sep = "_")
      return(unique(out[nzchar(out)]))
    }
    if ("Binomial" %in% cn) {
      out <- as.character(authority$Binomial)
      return(unique(out[nzchar(out)]))
    }
    stop("Authority data.frame must have columns Genus+Species or Binomial.")
  }
  if (is.character(authority) && length(authority) == 1L && nzchar(authority)) {
    path <- normalizePath(authority, mustWork = TRUE)
    df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", fileEncoding = "UTF-8")
    return(resolve_authority_binomials(df))
  }
  stop("`authority` must be a data.frame, a readable file path, or a character vector of binomials.")
}
.normalize_stem_mode <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    return("outgroup")
  }

  v <- tolower(trimws(x))

  if (v %in% c("outgroup", "stem", "with_outgroup")) {
    return("outgroup")
  }

  if (v %in% c("crown", "crown_only", "no_outgroup", "ingroup_only")) {
    return("crown")
  }

  warning("Unknown stem_mode '", x, "'; defaulting to 'outgroup'")
  "outgroup"
}

# -----------------------------------------------------------------------------
# Donor template preparation ---------------------------------------------------
# -----------------------------------------------------------------------------

#' Prepare a reusable donor template for grafting
#'
#' Constructs a reusable donor template suitable for clade‑ or tip‑level grafting
#' into a backbone chronogram. The function standardizes donor phylogenies by
#' resolving polytomies (optionally), cleaning tip labels, applying authority
#' filters, converting donors to time trees when needed, and computing the donor
#' stem fraction \eqn{r} used for geometric graft placement.
#'
#' Donor trees may be provided either as phylograms or as chronograms:
#' \itemize{
#'   \item If a donor tree is already ultrametric (within tolerance), it is treated
#'         as a chronogram and used directly after validation and rescaling.
#'   \item If the donor tree is not ultrametric, a chronogram is estimated using
#'         \pkg{ape}::\code{chronos()} according to \code{chronos_select}.
#' }
#'
#' Likelihood‑based chronos fits are preferred whenever possible. If no model
#' converges but one or more fits produce finite likelihoods, the best
#' non‑converged solution is retained with an explicit warning. Non‑parametric
#' smoothing (NNLS) is used only if all likelihood‑based fits fail.
#'
#' When \code{expand_ingroup_to_full_donor = TRUE}, the supplied
#' \code{ingroup_tips} are treated as a hint only: a single outgroup is inferred
#' and the ingroup is expanded to all remaining donor tips. This is required for
#' meaningful stem‑length estimation.
#'
#' @param donor_tree A rooted \code{phylo} donor tree (phylogram or chronogram).
#' @param ingroup_tips Character vector of donor tips belonging to the ingroup.
#'   These may be expanded depending on \code{expand_ingroup_to_full_donor}.
#' @param template_path Optional character path for the prepared donor-tree
#'   template.
#' @param write_template Logical. If `TRUE`, write the prepared template to
#'   `template_path`.
#' @param reuse_template Logical. If `TRUE` and an existing template is
#'   available at `template_path`, reuse it rather than rebuilding it.
#' @param stem_mode Character. Controls how donor-tree stem information is
#'   represented when preparing the grafting template.
#' @param authority Optional species authority data accepted by
#'   \code{resolve_authority_binomials()}.
#' @param chronos_select One of \code{"auto"}, \code{"fixed"}, or \code{"off"}.
#'   In \code{"auto"} mode, chronos is only applied if the donor is not already
#'   ultrametric.
#' @param chronos_model Chronos model used if \code{chronos_select = "fixed"}.
#' @param lambda,lambda_grid,models_grid,rate_cats,calib_df Parameters controlling
#'   chronos model fitting and selection.
#' @param resolve_polytomies One of \code{"none"}, \code{"donor"},
#'   \code{"backbone"}, or \code{"both"}.
#' @param resolve_random Logical; if TRUE, polytomies are resolved randomly.
#' @param poly_eps Numeric epsilon used to nudge non‑positive edges.
#' @param seed Numeric seed used upstream for reproducibility logging.
#' @param logger Optional logger object.
#' @param expand_ingroup_to_full_donor Logical; whether to infer an outgroup and
#'   expand the ingroup to the full donor tree.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{ok}: logical indicating success
#'     \item \code{donor_unit}: donor crown subtree scaled to unit height
#'     \item \code{stem_fraction}: inferred donor stem fraction \eqn{r}
#'     \item \code{template_log}: data.frame summarizing donor preprocessing
#'     \item \code{logfile}: path to the logfile, if logging is enabled
#'   }
#'
#' @keywords phylogenetics grafting chronogram donor‑processing
#' @export
prepare_clade_template <- function(
  donor_tree, ingroup_tips,
  template_path = NULL,
  write_template = TRUE,
  reuse_template = TRUE,
  authority = NULL,
  stem_mode = NULL,
  chronos_select = c("auto", "fixed", "off"),
  chronos_model = "correlated",
  lambda = 1,
  lambda_grid = c(1e-1, 1, 10),
  models_grid = c("correlated", "relaxed", "discrete"),
  rate_cats = c(1, 2, 5),
  resolve_polytomies = c("none", "donor", "backbone", "both"),
  resolve_random = TRUE,
  poly_eps = 1e-8,
  calib_df = NULL,
  seed = 42,
  logger = NULL,
  expand_ingroup_to_full_donor = TRUE
) {
  stem_mode <- .normalize_stem_mode(stem_mode)
  if (!is.null(template_path) &&
    reuse_template &&
    file.exists(template_path)) {
    tpl <- readRDS(template_path)

    if (isTRUE(tpl$ok)) {
      log_msg(logger, "Reusing cached donor template: ", template_path)
      return(tpl)
    }
  }
  chronos_select <- match.arg(chronos_select)
  resolve_polytomies <- match.arg(resolve_polytomies)
  tr <- donor_tree

  log_section(logger, "Template: input summary")

  ## Resolve donor polytomies if requested
  if (resolve_polytomies %in% c("donor", "both")) {
    tr <- ape::multi2di(tr, random = isTRUE(resolve_random))
    if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
    tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- poly_eps
  }

  ## -------------------------------------------------------------------------
  ## Clean labels
  ## -------------------------------------------------------------------------

  clean_species_label <- function(x) {
    ## 1. normalize dots in qualifiers (cf., nr., aff., sp.)
    x <- gsub("(cf|nr|aff|sp)\\.", "\\1", x)

    ## 2. extract Genus_species safely, drop suffixes/IDs
    sub(
      "^([A-Z][a-z]+)_(?:cf|nr|aff|sp_nr)?_?([a-z]{3,}).*",
      "\\1_\\2",
      x
    )
  }

  ## Store ORIGINAL donor labels so we can protect the outgroup if needed
  donor_labels_raw <- tr$tip.label

  ## Apply cleaning to tree tips
  tr$tip.label <- vapply(tr$tip.label, clean_species_label, character(1))

  ## Clean ingroup hints in the same namespace
  ingroup_tips <- intersect(
    vapply(ingroup_tips, clean_species_label, character(1)),
    tr$tip.label
  )

  ## Fallback: if ingroup collapsed entirely, treat all as ingroup for now
  if (!length(ingroup_tips)) ingroup_tips <- tr$tip.label


  ## -------------------------------------------------------------------------
  ## Outgroup handling
  ## -------------------------------------------------------------------------

  out_tip <- NA_character_
  out_tip_raw <- NA_character_

  if (isTRUE(expand_ingroup_to_full_donor) ||
    length(ingroup_tips) == length(tr$tip.label)) {
    og <- .infer_single_outgroup(tr)

    if (!is.na(og) && og %in% tr$tip.label) {
      out_tip <- og
      out_tip_raw <- og

      ingroup_tips <- setdiff(tr$tip.label, out_tip)
      log_msg(logger, "Outgroup inferred as root singleton: ", out_tip)
    }
  }


  ## -------------------------------------------------------------------------
  ## Authority filtering
  ## -------------------------------------------------------------------------

  valid <- if (is.null(authority)) NULL else resolve_authority_binomials(authority)

  if (!is.null(valid)) {
    keep_ig <- intersect(ingroup_tips, valid)

    log_section(logger, "Template: authority filter")
    log_msg(
      logger,
      "Valid ingroup species kept: ", length(keep_ig),
      "; dropped absent: ", length(setdiff(ingroup_tips, valid))
    )

    if (!length(keep_ig)) {
      return(list(ok = FALSE, reason = "No ingroup retained after authority filter"))
    }

    ## ---- PROTECT OUTGROUP FROM AUTHORITY FILTER ----
    keep <- keep_ig
    if (!is.na(out_tip_raw)) {
      keep <- union(keep, out_tip_raw)
    }

    tr <- ape::keep.tip(tr, keep)
    ingroup_tips <- keep_ig
  }

  ## -------------------------------------------------------------------------
  ## Donor ultrametric conversion
  ## -------------------------------------------------------------------------
  log_section(logger, "Template: donor ultrametric conversion")

  G <- tr

  tree_type <- attr(donor_tree, "input_tree_type", exact = TRUE)
  if (is.null(tree_type)) tree_type <- "phylogram"

  tree_type <- tolower(tree_type)
  is_ultra <- ape::is.ultrametric(G, tol = 1e-6)

  log_msg(
    logger,
    "Donor tree type: ", tree_type,
    "; ultrametric: ", is_ultra
  )

  safe_num1 <- function(x) {
    if (inherits(x, "try-error") || is.null(x) || !length(x)) {
      return(NA_real_)
    }
    x <- suppressWarnings(as.numeric(x[1]))
    if (is.finite(x)) x else NA_real_
  }

  extract_chronos_metrics <- function(fit) {
    phi_info <- attr(fit, "PHIIC", exact = TRUE)

    phiic <- NA_real_
    loglik <- NA_real_

    if (is.list(phi_info)) {
      phiic <- safe_num1(phi_info$PHIIC)
      loglik <- safe_num1(phi_info$logLik)
    }

    list(
      PHIIC = phiic,
      logLik = loglik,
      has_phiic_attr = !is.null(phi_info),
      phiic_attr_class = if (is.null(phi_info)) {
        NA_character_
      } else {
        paste(class(phi_info), collapse = ",")
      }
    )
  }

  fmt_num <- function(x) {
    if (length(x) == 1L && is.finite(x)) {
      sprintf("%.6f", x)
    } else {
      "NA"
    }
  }

empty_chronos_meta <- function(status = "not_run", reason = NA_character_) {
  list(
    status = status,
    reason = reason,
    selected = data.frame(),
    trials = data.frame()
  )
}

chronos_meta <- empty_chronos_meta()

  ## -------------------------------------------------------------------------
  ## Case 1: Declared chronogram — NEVER run chronos
  ## -------------------------------------------------------------------------
  if (identical(tree_type, "chronogram")) {
    log_msg(logger, "Input tree declared as chronogram; skipping chronos.")

    if (!is_ultra) {
      log_msg(
        logger,
        "Chronogram not ultrametric (likely numerical jitter); ",
        "enforcing ultrametricity using 'extend'."
      )
      G_ultra <- phytools::force.ultrametric(G, method = "extend")
      chronos_meta <- empty_chronos_meta(
        status = "not_run",
        reason = "Input tree declared as chronogram."
      )
    } else {
      G_ultra <- G
      chronos_meta <- empty_chronos_meta(
        status = "not_run",
        reason = "Input tree declared as chronogram."
      )
    }
  } else {
    ## -----------------------------------------------------------------------
    ## Case 2: Phylogram — allow chronos logic
    ## -----------------------------------------------------------------------
    if (is_ultra && chronos_select != "fixed") {
      log_msg(
        logger,
        "Phylogram already ultrametric; skipping chronos and using as chronogram."
      )
      G_ultra <- G
      chronos_meta <- empty_chronos_meta(
        status = "not_run",
        reason = "Phylogram already ultrametric; chronos skipped."
      )
    } else {
      log_msg(
        logger,
        "Phylogram detected; estimating chronogram with chronos."
      )

      cal <- if (!is.null(calib_df)) calib_df else ape::makeChronosCalib(G)

      run_chronos_safely <- function(expr) {
        warnings <- character(0)

        old_warn <- getOption("warn")
        on.exit(options(warn = old_warn), add = TRUE)

        # Prevent user/test-session options(warn = 2) from converting chronos
        # warnings into errors.
        options(warn = 0)

        fit <- withCallingHandlers(
          tryCatch(
            expr,
            error = function(e) e
          ),
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

        list(
          fit = fit,
          warnings = unique(warnings),
          error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_
        )
      }

      pick_best <- function(G, logger, ig) {
        candidates <- list()
        trials <- list()

        add_trial <- function(
          mdl, lam, kcat, fit = NULL, status = "ok",
          error_message = NA_character_,
          warning_message = NA_character_
        ) {
          if (!is.null(fit) && !inherits(fit, "try-error")) {
            info <- extract_chronos_metrics(fit)
            conv <- isTRUE(attr(fit, "convergence"))

            crown_metrics <- try(.compute_crown_metrics(fit, ig), silent = TRUE)
            crown_depth <- if (!inherits(crown_metrics, "try-error") &&
                               isTRUE(crown_metrics$ok)) {
              crown_metrics$crown_depth
            } else {
              NA_real_
            }

            is_ultra_fit <- try(ape::is.ultrametric(fit, tol = 1e-6), silent = TRUE)
            is_ultra_fit <- if (inherits(is_ultra_fit, "try-error")) {
              NA
            } else {
              isTRUE(is_ultra_fit)
            }
            edge_lengths <- fit$edge.length
            min_edge <- if (is.null(edge_lengths)) {
              NA_real_
            } else {
              suppressWarnings(min(edge_lengths, na.rm = TRUE))
            }
            max_edge <- if (is.null(edge_lengths)) {
              NA_real_
            } else {
              suppressWarnings(max(edge_lengths, na.rm = TRUE))
            }
            n_nonpositive <- if (is.null(edge_lengths)) {
              NA_integer_
            } else {
              sum(!is.finite(edge_lengths) | edge_lengths <= 0)
            }
            row <- data.frame(
              model = mdl,
              lambda = lam,
              nb_rate_cat = if (is.na(kcat)) NA_integer_ else as.integer(kcat),
              status = status,
              converged = conv,
              PHIIC = info$PHIIC,
              logLik = info$logLik,
              has_phiic_attr = info$has_phiic_attr,
              phiic_attr_class = info$phiic_attr_class,
              crown_depth = crown_depth,
              is_ultrametric = is_ultra_fit,
              n_tip = ape::Ntip(fit),
              n_node = fit$Nnode,
              min_edge = min_edge,
              max_edge = max_edge,
              n_nonpositive_edges = n_nonpositive,
              error_message = error_message,
              warning_message = warning_message,
              stringsAsFactors = FALSE
            )

            log_msg(
              logger,
              "chronos trial: ",
              "model=", mdl,
              " lambda=", lam,
              " rate_cats=", if (is.na(kcat)) "NA" else kcat,
              " converged=", conv,
              " PHIIC=", fmt_num(info$PHIIC),
              " logLik=", fmt_num(info$logLik),
              " crown_depth=", fmt_num(crown_depth),
              " ultrametric=", is_ultra_fit
            )

            trials[[length(trials) + 1L]] <<- row

            ## Keep any finite-likelihood fit, even if PHIIC is NA/NaN.
            if (is.finite(info$logLik)) {
              candidates[[length(candidates) + 1L]] <<- list(
                tree = fit,
                model = mdl,
                lambda = lam,
                nb_rate_cat = if (is.na(kcat)) NA_integer_ else as.integer(kcat),
                PHIIC = info$PHIIC,
                logLik = info$logLik,
                converged = conv,
                crown = crown_depth,
                is_ultrametric = is_ultra_fit
              )
            }
          } else {
            row <- data.frame(
              model = mdl,
              lambda = lam,
              nb_rate_cat = if (is.na(kcat)) NA_integer_ else as.integer(kcat),
              status = status,
              converged = FALSE,
              PHIIC = NA_real_,
              logLik = NA_real_,
              has_phiic_attr = FALSE,
              phiic_attr_class = NA_character_,
              crown_depth = NA_real_,
              is_ultrametric = NA,
              n_tip = ape::Ntip(G),
              n_node = G$Nnode,
              min_edge = if (is.null(G$edge.length)) {
                NA_real_
              } else {
                suppressWarnings(min(G$edge.length, na.rm = TRUE))
              },
              max_edge = if (is.null(G$edge.length)) {
                NA_real_
              } else {
                suppressWarnings(max(G$edge.length, na.rm = TRUE))
              },
              n_nonpositive_edges = if (is.null(G$edge.length)) {
                NA_integer_
              } else {
                sum(!is.finite(G$edge.length) | G$edge.length <= 0)
              },
              error_message = error_message,
              warning_message = warning_message,
              stringsAsFactors = FALSE
            )

            log_msg(
              logger,
              "chronos trial FAILED: ",
              "model=", mdl,
              " lambda=", lam,
              " rate_cats=", if (is.na(kcat)) "NA" else kcat,
              if (!is.na(warning_message)) paste0(" warning=", warning_message) else "",
              if (!is.na(error_message)) paste0(" error=", error_message) else ""
            )

            trials[[length(trials) + 1L]] <<- row
          }
        }

        for (mdl in models_grid) {
          if (mdl == "discrete") {
            for (kcat in rate_cats) {
              for (lam in lambda_grid) {

                ctrl <- ape::chronos.control()
                ctrl$nb.rate.cat <- kcat

                chr <- run_chronos_safely(
                  ape::chronos(
                    G,
                    lambda = lam,
                    model = mdl,
                    quiet = FALSE,
                    calibration = cal,
                    control = ctrl
                  )
                )

                fit <- chr$fit

                warning_message <- if (length(chr$warnings)) {
                  paste(chr$warnings, collapse = " | ")
                } else {
                  NA_character_
                }

                if (inherits(fit, "error")) {
                  add_trial(
                    mdl = mdl,
                    lam = lam,
                    kcat = kcat,
                    fit = NULL,
                    status = "error",
                    error_message = chr$error,
                    warning_message = warning_message
                  )
                } else {
                  add_trial(
                    mdl = mdl,
                    lam = lam,
                    kcat = kcat,
                    fit = fit,
                    status = "ok",
                    warning_message = warning_message
                  )
                }
              }
            }
          } else {
            for (lam in lambda_grid) {
              chr <- run_chronos_safely(
                ape::chronos(
                  G,
                  lambda = lam,
                  model = mdl,
                  quiet = FALSE,
                  calibration = cal
                )
              )

              fit <- chr$fit

              warning_message <- if (length(chr$warnings)) {
                paste(chr$warnings, collapse = " | ")
              } else {
                NA_character_
              }

              if (inherits(fit, "error")) {
                add_trial(
                  mdl = mdl,
                  lam = lam,
                  kcat = NA_integer_,
                  fit = NULL,
                  status = "error",
                  error_message = chr$error,
                  warning_message = warning_message
                )
              } else {
                add_trial(
                  mdl = mdl,
                  lam = lam,
                  kcat = NA_integer_,
                  fit = fit,
                  status = "ok",
                  warning_message = warning_message
                )
              }
            }
          }
        }

        trials_df <- if (length(trials)) {
          do.call(rbind, trials)
        } else {
          data.frame()
        }

        if (!length(candidates)) {
          log_msg(
            logger,
            "ERROR: No chronos fit produced a finite logLik; ",
            "falling back to NNLS."
          )

          selected_df <- data.frame(
            model = "nnls",
            lambda = NA_real_,
            nb_rate_cat = NA_integer_,
            PHIIC = NA_real_,
            logLik = NA_real_,
            converged = FALSE,
            ranking_criterion = "none",
            note = "No finite chronos logLik; fallback to NNLS.",
            stringsAsFactors = FALSE
          )

          return(list(
            tree = phytools::force.ultrametric(G, method = "nnls"),
            model = "nnls",
            lambda = NA_real_,
            nb_rate_cat = NA_integer_,
            PHIIC = NA_real_,
            logLik = NA_real_,
            converged = FALSE,
            chronos_status = "nnls_fallback",
            chronos_trials = trials_df,
            chronos_selected = selected_df
          ))
        }

        conv_flag <- vapply(candidates, `[[`, logical(1), "converged")

        if (any(conv_flag)) {
          pool <- candidates[conv_flag]
          note <- "Selected best CONVERGED chronos model."
        } else {
          pool <- candidates
          note <- "WARNING: No chronos model converged; using best NON-converged model."
        }

        phiics <- vapply(pool, function(x) x$PHIIC, numeric(1))
        logliks <- vapply(pool, function(x) x$logLik, numeric(1))

        if (any(is.finite(phiics))) {
          rank_score <- ifelse(is.finite(phiics), phiics, Inf)
          j <- which.min(rank_score)
          ranking_criterion <- "PHIIC"
        } else {
          rank_score <- ifelse(is.finite(logliks), logliks, -Inf)
          j <- which.max(rank_score)
          ranking_criterion <- "logLik"
          note <- paste(note, "PHIIC unavailable/NaN; ranked by logLik.")
        }

        best <- pool[[j]]

        selected_df <- data.frame(
          model = best$model,
          lambda = best$lambda,
          nb_rate_cat = best$nb_rate_cat,
          PHIIC = best$PHIIC,
          logLik = best$logLik,
          converged = best$converged,
          ranking_criterion = ranking_criterion,
          note = note,
          stringsAsFactors = FALSE
        )

        log_msg(
          logger,
          "Selected chronos:",
          "\n  model=", best$model,
          "\n  lambda=", best$lambda,
          "\n  rate_cats=", best$nb_rate_cat,
          "\n  PHIIC=", fmt_num(best$PHIIC),
          "\n  logLik=", fmt_num(best$logLik),
          "\n  converged=", best$converged,
          "\n  ranking_criterion=", ranking_criterion,
          "\n  note=", note
        )

        best$chronos_status <- if (best$converged) {
          "selected_converged"
        } else {
          "selected_nonconverged"
        }
        best$chronos_trials <- trials_df
        best$chronos_selected <- selected_df

        best
      }

      sel <- switch(chronos_select,
        "off" = list(
          tree = phytools::force.ultrametric(G, method = "nnls"),
          model = "nnls",
          lambda = NA_real_,
          nb_rate_cat = NA_integer_,
          PHIIC = NA_real_,
          logLik = NA_real_,
          converged = FALSE,
          chronos_status = "chronos_off_nnls",
          chronos_trials = data.frame(),
          chronos_selected = data.frame(
            model = "nnls",
            lambda = NA_real_,
            nb_rate_cat = NA_integer_,
            PHIIC = NA_real_,
            logLik = NA_real_,
            converged = FALSE,
            ranking_criterion = "none",
            note = "chronos_select='off'; used NNLS.",
            stringsAsFactors = FALSE
          )
        ),
        "fixed" = {
          mdl <- tolower(chronos_model)
          ctrl <- ape::chronos.control()
          chr <- run_chronos_safely(
            ape::chronos(
              G,
              lambda = lambda,
              model = mdl,
              quiet = FALSE,
              calibration = cal,
              control = ctrl
            )
          )

          fit <- chr$fit

          warning_message <- if (length(chr$warnings)) {
            paste(chr$warnings, collapse = " | ")
          } else {
            NA_character_
          }

          if (!inherits(fit, "error") && isTRUE(attr(fit, "convergence"))) {
            info <- extract_chronos_metrics(fit)

            log_msg(
              logger,
              "Fixed chronos:",
              "\n  model=", mdl,
              "\n  lambda=", lambda,
              "\n  PHIIC=", fmt_num(info$PHIIC),
              "\n  logLik=", fmt_num(info$logLik),
              "\n  converged=", isTRUE(attr(fit, "convergence"))
            )

            trial_df <- data.frame(
              model = mdl,
              lambda = lambda,
              nb_rate_cat = ctrl$nb.rate.cat,
              status = "ok",
              converged = isTRUE(attr(fit, "convergence")),
              PHIIC = info$PHIIC,
              logLik = info$logLik,
              has_phiic_attr = info$has_phiic_attr,
              phiic_attr_class = info$phiic_attr_class,
              crown_depth = .compute_crown_metrics(fit, ingroup_tips)$crown_depth,
              is_ultrametric = ape::is.ultrametric(fit, tol = 1e-6),
              n_tip = ape::Ntip(fit),
              n_node = fit$Nnode,
              min_edge = suppressWarnings(min(fit$edge.length, na.rm = TRUE)),
              max_edge = suppressWarnings(max(fit$edge.length, na.rm = TRUE)),
              n_nonpositive_edges = sum(!is.finite(fit$edge.length) | fit$edge.length <= 0),
              error_message = NA_character_,
              warning_message = warning_message,
              stringsAsFactors = FALSE
            )

            selected_df <- data.frame(
              model = mdl,
              lambda = lambda,
              nb_rate_cat = ctrl$nb.rate.cat,
              PHIIC = info$PHIIC,
              logLik = info$logLik,
              converged = isTRUE(attr(fit, "convergence")),
              ranking_criterion = "fixed",
              note = "Fixed chronos model requested.",
              stringsAsFactors = FALSE
            )

            list(
              tree = fit,
              model = mdl,
              lambda = lambda,
              nb_rate_cat = ctrl$nb.rate.cat,
              PHIIC = info$PHIIC,
              logLik = info$logLik,
              converged = isTRUE(attr(fit, "convergence")),
              chronos_status = "fixed",
              chronos_trials = trial_df,
              chronos_selected = selected_df
            )
          } else {
            log_msg(
              logger,
              "Fixed chronos failed; using NNLS.",
              if (!is.na(warning_message)) paste0(" warning=", warning_message) else "",
              if (inherits(fit, "error")) paste0(" error=", chr$error) else ""
            )

            selected_df <- data.frame(
              model = "nnls",
              lambda = NA_real_,
              nb_rate_cat = NA_integer_,
              PHIIC = NA_real_,
              logLik = NA_real_,
              converged = FALSE,
              ranking_criterion = "none",
              note = "Fixed chronos failed; fallback to NNLS.",
              stringsAsFactors = FALSE
            )

            list(
              tree = phytools::force.ultrametric(G, method = "nnls"),
              model = "nnls",
              lambda = NA_real_,
              nb_rate_cat = NA_integer_,
              PHIIC = NA_real_,
              logLik = NA_real_,
              converged = FALSE,
              chronos_status = "fixed_failed_nnls",
              chronos_trials = data.frame(),
              chronos_selected = selected_df
            )
          }
        },
        "auto" = pick_best(G, logger, ingroup_tips)
      )

      if (is.null(sel$chronos_trials)) {
        sel$chronos_trials <- data.frame()
      }

      if (is.null(sel$chronos_selected)) {
        sel$chronos_selected <- data.frame(
          model = sel$model %||% NA_character_,
          lambda = sel$lambda %||% NA_real_,
          nb_rate_cat = sel$nb_rate_cat %||% NA_integer_,
          PHIIC = sel$PHIIC %||% NA_real_,
          logLik = sel$logLik %||% NA_real_,
          converged = sel$converged %||% NA,
          ranking_criterion = NA_character_,
          note = sel$chronos_status %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }

      chronos_meta <- list(
        status = sel$chronos_status %||% "selected",
        reason = NA_character_,
        selected = sel$chronos_selected,
        trials = sel$chronos_trials
      )

      G_ultra <- sel$tree
    }
  }
  ## -------------------------------------------------------------------------
  ## Stem / crown metrics and donor crown extraction
  ## -------------------------------------------------------------------------
  metrics <- .compute_crown_metrics(G_ultra, ingroup_tips)
  if (!isTRUE(metrics$ok)) {
    return(list(ok = FALSE, reason = metrics$reason))
  }

  if (stem_mode == "outgroup") {
    r <- metrics$r
    stem_len <- metrics$stem_len
  } else {
    r <- NA_real_
    stem_len <- NA_real_
  }
  crown_node <- metrics$crown_node

  donor_crown <- ape::extract.clade(G_ultra, crown_node)
  if (!is.na(out_tip) && out_tip %in% donor_crown$tip.label) {
    donor_crown <- ape::drop.tip(donor_crown, out_tip)
  }

  donor_crown <- .safe_force_ultrametric(donor_crown)

  unique_species <- unique(donor_crown$tip.label)
  needs_suffix <- length(unique_species) == 1L
  orig_label <- unique_species
  donor_label <- paste0(orig_label, "_donor")

  if (needs_suffix) {
    donor_crown$tip.label[] <- donor_label
  }

  donor_unit <- .rescale_to_height(donor_crown, 1.0)
  donor_unit$root.edge <- NULL

  label_map <- if (needs_suffix) {
    data.frame(from = donor_label, to = orig_label, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  log_section(logger, "Template: output")
  log_msg(
    logger, "Donor stem_fraction=", sprintf("%.6f", r),
    " crown_fraction=", sprintf("%.6f", 1 - r)
  )

  selected_chronos <- chronos_meta$selected

  chronos_model <- if (nrow(selected_chronos)) {
    as.character(selected_chronos$model[1])
  } else {
    NA_character_
  }

  chronos_lambda <- if (nrow(selected_chronos)) {
    selected_chronos$lambda[1]
  } else {
    NA_real_
  }

  chronos_rate_cat <- if (nrow(selected_chronos)) {
    selected_chronos$nb_rate_cat[1]
  } else {
    NA_integer_
  }

  chronos_phiic <- if (nrow(selected_chronos)) {
    selected_chronos$PHIIC[1]
  } else {
    NA_real_
  }

  chronos_loglik <- if (nrow(selected_chronos)) {
    selected_chronos$logLik[1]
  } else {
    NA_real_
  }

  chronos_converged <- if (nrow(selected_chronos)) {
    selected_chronos$converged[1]
  } else {
    NA
  }

  chronos_rank <- if (nrow(selected_chronos)) {
    as.character(selected_chronos$ranking_criterion[1])
  } else {
    NA_character_
  }

  tpl_row <- data.frame(
    n_ingroup_kept = length(ingroup_tips),
    outgroup_used = if (!is.na(out_tip)) out_tip else NA_character_,
    stem_len_inferred = stem_len,
    stem_fraction = r,
    crown_fraction = 1 - r,
    chronos_status = chronos_meta$status,
    chronos_reason = chronos_meta$reason,
    chronos_model = chronos_model,
    chronos_lambda = chronos_lambda,
    chronos_rate_cat = chronos_rate_cat,
    chronos_PHIIC = chronos_phiic,
    chronos_logLik = chronos_loglik,
    chronos_converged = chronos_converged,
    chronos_ranking_criterion = chronos_rank,
    n_chronos_trials = nrow(chronos_meta$trials),
    stringsAsFactors = FALSE
  )

  out <- list(
    ok = TRUE,
    donor_unit = donor_unit,
    stem_fraction = r,
    stem_len = stem_len,
    label_map = label_map,
    template_log = tpl_row,
    chronos_meta = chronos_meta,
    logfile = if (inherits(logger, "smart_logger")) logger$file else NA_character_
  )
  
  if (!is.null(template_path) && write_template) {
    dir.create(
      dirname(normalizePath(template_path, mustWork = FALSE)),
      recursive = TRUE,
      showWarnings = FALSE
    )

    tmp_path <- paste0(template_path, ".tmp")

    saveRDS(out, tmp_path)

    if (!file.rename(tmp_path, template_path)) {
      unlink(tmp_path)
      stop("Failed to move temporary template into place: ", template_path)
    }

    log_msg(logger, "Saved donor template to: ", template_path)
  }

  return(out)

}

# -----------------------------------------------------------------------------
# Placement (tip or clade) -----------------------------------------------------
# -----------------------------------------------------------------------------

#' Graft a donor template onto an ultrametric backbone
#'
#' Grafts a donor crown subtree onto an ultrametric backbone tree using
#' **terminal-edge replacement**. Two placement modes are supported:
#'
#' \describe{
#'   \item{Tip replacement}{
#'     The donor is grafted directly onto the terminal branch leading to a
#'     specified backbone tip. The donor stem fraction \eqn{r} determines the
#'     relative depth of grafting within that branch.
#'   }
#'   \item{Clade (MRCA) replacement}{
#'     All descendant tips of the target MRCA clade are dropped except for
#'     a single retained backbone tip (chosen deterministically). The donor
#'     is then grafted onto the terminal branch of that retained real tip using
#'     the same logic as tip replacement. No synthetic tips or internal-edge
#'     insertion is used.
#'   }
#' }
#'
#' This design ensures that clade grafting preserves:
#' \itemize{
#'   \item valid donor stem geometry (via \eqn{r}),
#'   \item true ancestral context for grafting,
#'   \item strict ultrametricity without post-hoc adjustment.
#' }
#'
#' @param backbone_ultra An ultrametric \code{phylo} backbone tree.
#' @param placement A list describing placement:
#'   \itemize{
#'     \item \code{list(type = "tip", tip = <label>)}
#'     \item \code{list(type = "clade", anchors = <character vector>)}
#'   }
#' @param template A donor template produced by \code{prepare_clade_template()}.
#' @param eps Numeric epsilon used for positional clamping.
#' @param logger Optional logger.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{tree}: the grafted \code{phylo}
#'     \item \code{log}: a data.frame describing the graft
#'   }
#'
#' @export
graft_template_onto_backbone <- function(
  backbone_ultra,
  placement,
  template,
  eps = 1e-10,
  logger = NULL
) {
  stem_mode <- .normalize_stem_mode(placement$stem_mode)
  stopifnot(isTRUE(template$ok))
  X <- ape::reorder.phylo(backbone_ultra, "postorder")

  donor_unit <- template$donor_unit
  r <- template$stem_fraction

  ## -------------------------------------------------------------------------
  ## TIP REPLACEMENT
  ## -------------------------------------------------------------------------
  if (identical(placement$type, "tip")) {
    if (stem_mode == "crown") {
      out <- .graft_on_terminal_edge_beta(
        X,
        placement$tip,
        donor_unit,
        shape1 = placement$shape1 %||% 6,
        shape2 = placement$shape2 %||% 2,
        min_frac = placement$min_frac %||% 0.1,
        max_frac = placement$max_frac %||% 0.9,
        eps = eps
      )
    } else {
      out <- .graft_on_terminal_edge(
        X,
        placement$tip,
        donor_unit,
        r,
        eps
      )
    }

    log_row <- transform(
      template$template_log,
      stem_mode = stem_mode,
      target_type = "tip",
      target_label = placement$tip
    )

    return(list(tree = out, log = log_row))
  }

  ## -------------------------------------------------------------------------
  ## CLADE (MRCA) REPLACEMENT
  ## -------------------------------------------------------------------------
  if (identical(placement$type, "clade")) {
    anchors <- as.character(placement$anchors)
    if (!all(anchors %in% X$tip.label)) {
      stop(
        "Anchors not found in backbone: ",
        paste(setdiff(anchors, X$tip.label), collapse = ", ")
      )
    }

    mrca <- ape::getMRCA(X, anchors)
    if (is.na(mrca)) {
      stop("Could not determine MRCA for anchors.")
    }

    ## identify all descendant tips
    desc <- unique(phytools::getDescendants(X, mrca))
    tip_desc <- desc[desc <= ape::Ntip(X)]
    clade_tips <- X$tip.label[tip_desc]

    ## deterministically retain one real backbone tip
    keep_tip <- sort(clade_tips)[1]
    drop_tips <- setdiff(clade_tips, keep_tip)

    log_msg(
      logger,
      "Collapsing MRCA(", paste(anchors, collapse = ","),
      ") by retaining backbone tip ", keep_tip,
      " and grafting donor as terminal replacement."
    )

    if (stem_mode == "crown") {
      crown_res <- .graft_replace_clade_at_crown_depth(
        X,
        mrca,
        donor_unit,
        eps
      )

      out <- crown_res$tree
      retained_tip <- crown_res$retained_tip

      log_row <- template$template_log
      log_row$stem_mode <- stem_mode
      log_row$target_type <- "clade"
      log_row$target_label <- paste(anchors, collapse = ",")
      log_row$retained_tip <- retained_tip

      return(list(tree = out, log = log_row))
    }

    ## prune backbone clade
    X_pruned <- ape::drop.tip(X, drop_tips)

    ## graft donor exactly as a tip replacement
    out <- .graft_on_terminal_edge(
      X = X_pruned,
      tip = keep_tip,
      donor_unit = donor_unit,
      r = r,
      eps = eps
    )
    log_row <- template$template_log
    log_row$stem_mode <- stem_mode
    log_row$target_type <- "clade"
    log_row$target_label <- paste(anchors, collapse = ",")
    log_row$retained_tip <- keep_tip

    return(list(tree = out, log = log_row))
  }

  stop("placement$type must be 'tip' or 'clade'.")
}

# -----------------------------------------------------------------------------
# Batch driver & wrapper -------------------------------------------------------
# -----------------------------------------------------------------------------

#' Batch grafting across one or many backbones
#'
#' Prepares donor templates (once per donor) and grafts them into each backbone
#' tree. Optionally enforces ultrametricity at the end (default: \code{"extend"}).
#'
#' @param backbones A \code{phylo} or \code{multiPhylo} backbone chronogram.
#' @param donors Named list of donor \code{phylo} trees.
#' @param placements Named list of placement descriptors (one per donor name).
#' @param authority Resolved authority vector or object accepted by
#'   [resolve_authority_binomials()].
#' @param chronos_select One of \code{"auto"}, \code{"fixed"}, \code{"off"}.
#' @param chronos_model Chronos model if \code{chronos_select="fixed"}.
#' @param lambda,lambda_grid,models_grid,rate_cats,calib_df Chronos controls.
#' @param seed Numeric seed used upstream.
#' @param eps Small epsilon for positional clamping.
#' @param ultrametric_final One of \code{"extend"}, \code{"nnls"}, \code{"none"}.
#' @param resolve_polytomies \code{"none"|"donor"|"backbone"|"both"}; threaded
#'   into donor preparation.
#' @param resolve_random Logical; if TRUE, \code{multi2di()} resolves randomly.
#' @param poly_eps Numeric epsilon used to nudge edges after \code{multi2di()}.
#' @param out_tree Output Newick path.
#' @param out_log Output TSV log path.
#' @param logger Optional logger.
#' @return A list with \code{tree} (phylo|multiPhylo) and \code{log} (data.frame).
#' @export
graft_many_clades <- function(
  backbones, donors, placements,
  authority = NULL,
  chronos_select = c("auto", "fixed", "off"),
  chronos_model = "correlated", lambda = 1,
  lambda_grid = c(1e-1, 1, 10), models_grid = c("correlated", "relaxed", "discrete"),
  rate_cats = c(1, 2, 5), calib_df = NULL,
  seed = 42, eps = 1e-10,
  ultrametric_final = c("extend", "nnls", "none"),
  resolve_polytomies = c("none", "donor", "backbone", "both"),
  resolve_random = TRUE,
  poly_eps = 1e-8,
  out_tree = "clade_grafted.tre",
  out_log = "clade_grafting_log.tsv",
  logger = NULL
) {
  chronos_select <- match.arg(chronos_select)
  ultrametric_final <- match.arg(ultrametric_final)
  resolve_polytomies <- match.arg(resolve_polytomies)

  is_multi <- inherits(backbones, "multiPhylo")
  Blist <- if (is_multi) backbones else list(backbones)
  if (!all(vapply(Blist, ape::is.ultrametric, logical(1)))) {
    stop("All backbone trees must be ultrametric.")
  }

  valid <- resolve_authority_binomials(authority)

  # Prepare donor templates once
  templates <- list()
  for (nm in names(donors)) {
    dtr <- donors[[nm]]
    pl <- placements[[nm]]
    ## Require ingroup_tips ONLY for raw donor trees
    if (!is.list(dtr)) {
      if (is.null(attr(dtr, "ingroup_tips"))) {
        stop("Donor '", nm, "' must carry attr(, 'ingroup_tips') from extractor.")
      }
      ig <- attr(dtr, "ingroup_tips")
    } else {
      ig <- NULL # templates do not need ingroup_tips
    }
    if (is.list(dtr) && isTRUE(dtr$ok) && !is.null(dtr$donor_unit)) {
      # Already a prepared template
      tpl <- dtr
    } else {
      template_path <- NULL
      if (is.character(attr(dtr, "source_path"))) {
        template_path <- sub("\\.[^.]+$", ".template.rds", attr(dtr, "source_path"))
      }

      # Raw donor tree → prepare template
      tpl <- prepare_clade_template(
        donor_tree = dtr, ingroup_tips = ig, template_path = template_path, authority = valid,
        chronos_select = chronos_select, chronos_model = chronos_model,
        lambda = lambda, lambda_grid = lambda_grid, models_grid = models_grid,
        rate_cats = rate_cats, calib_df = calib_df, seed = seed, logger = logger,
        resolve_polytomies = resolve_polytomies, resolve_random = resolve_random,
        poly_eps = poly_eps, expand_ingroup_to_full_donor = pl$type %in% c("clade", "tip")
      )
    }

    if (!isTRUE(tpl$ok)) {
      log_msg(
        logger, "[SKIP] ", nm, ": template could not be prepared (",
        tpl$reason, ")."
      )
    } else {
      templates[[nm]] <- tpl
    }
  }
  if (!length(templates)) stop("No donor templates could be prepared; aborting.")

  logs <- list()
  out_trees <- vector("list", length(Blist))
  for (i in seq_along(Blist)) {
    backbone_i <- Blist[[i]]
    for (nm in names(templates)) {
      tpl <- templates[[nm]]
      pl <- placements[[nm]]
      ok_target <- if (pl$type == "tip") {
        pl$tip %in% backbone_i$tip.label
      } else {
        length(intersect(pl$anchors, backbone_i$tip.label)) >= 2
      }
      if (!ok_target) {
        log_msg(logger, "[SKIP tree ", i, "] ", nm, ": target missing in this backbone.")
        next
      }
      res <- graft_template_onto_backbone(backbone_i,
        placement = pl, template = tpl,
        eps = eps, logger = logger
      )
      backbone_i <- res$tree
      res$log$backbone_index <- i
      logs[[length(logs) + 1]] <- res$log
    }

    if (ultrametric_final != "none" && !ape::is.ultrametric(backbone_i)) {
      backbone_i <- phytools::force.ultrametric(
        backbone_i,
        method = ultrametric_final
      )
    }

    if (ultrametric_final != "none") {
      stopifnot(ape::is.ultrametric(backbone_i))
    }

    backbone_i$node.label <- NULL
    out_trees[[i]] <- backbone_i
  }

  # Name-aware log binding ----------------------------------------------------
  if (length(logs)) {
    all_cols <- unique(unlist(lapply(logs, names)))
    logs_std <- lapply(logs, function(df) {
      miss <- setdiff(all_cols, names(df))
      if (length(miss)) df[miss] <- NA
      df[, all_cols, drop = FALSE]
    })
    log_df <- do.call(rbind, logs_std)
  } else {
    log_df <- data.frame()
  }
  if (ncol(log_df) > 0) {
    cn <- colnames(log_df)
    cn[is.na(cn) | cn == ""] <- "X"
    colnames(log_df) <- make.names(cn, unique = TRUE)
  }

  dir.create(dirname(normalizePath(out_log, mustWork = FALSE)),
    recursive = TRUE,
    showWarnings = FALSE
  )
  utils::write.table(log_df,
    file = out_log, sep = "\t", row.names = FALSE,
    quote = FALSE
  )
  dir.create(dirname(normalizePath(out_tree, mustWork = FALSE)),
    recursive = TRUE,
    showWarnings = FALSE
  )
  if (is_multi) {
    names(out_trees) <- NULL
    attr(out_trees, "order") <- NULL
    class(out_trees) <- "multiPhylo"
    ape::write.tree(out_trees, file = out_tree)
  } else {
    ape::write.tree(out_trees[[1]], file = out_tree)
  }

  invisible(list(
    tree = if (is_multi) out_trees else out_trees[[1]],
    log = log_df
  ))
}

#' Run clade grafting from a plan (MRCA or TIP in one column)
#'
#' Reads a backbone (\code{phylo} or \code{multiPhylo}), parses a placement
#' plan with either \code{MRCA}+\code{Phylogeny_file_path} (new schema) or
#' \code{Genus}+\code{Phylogeny_file_path} (legacy), prepares donor templates
#' and performs grafting. Placement is chosen automatically: cells in \code{MRCA}
#' containing commas become clade anchors; otherwise a single TIP replacement is
#' performed using the cell verbatim as the target tip label (including raw tip
#' names or a bare Genus).
#'
#' @param backbone_path Path to a backbone Newick (single tree or multiPhylo).
#' @param plan_path Path to a TSV plan: either with columns
#'   \code{MRCA}+\code{Phylogeny_file_path} or \code{Genus}+\code{Phylogeny_file_path}.
#' @param authority Species authority object accepted by
#'   [resolve_authority_binomials()].
#' @param out_prefix Output prefix (no extension).
#' @param seed_mode "fixed", "random" or a single integer; passed to your
#'   project-level seeding utility.
#' @param ingroup_anchors Optional anchors to extract an ingroup from the
#'   backbone before grafting.
#' @param drop_tips Optional character vector of tips to drop from backbone.
#' @param ensure_ultrametric One of \code{"nnls"}, \code{"extend"}, \code{"none"}
#'   (applied to the backbone before grafting).
#' @param chronos_select One of \code{"auto"}, \code{"fixed"}, \code{"off"}.
#' @param chronos_model Chronos model when \code{chronos_select="fixed"}.
#' @param lambda,lambda_grid,models_grid,rate_cats,calib_df Chronos controls.
#' @param ultrametric_final One of \code{"extend"}, \code{"nnls"}, \code{"none"}
#'   (applied after grafting).
#' @param resolve_polytomies \code{"none"|"donor"|"backbone"|"both"}; if
#'   includes \code{"donor"} or \code{"both"}, resolves donor polytomies before
#'   chronos; if includes \code{"backbone"} or \code{"both"}, resolves backbone
#'   polytomies before grafting.
#' @param resolve_random Logical; if TRUE, \code{multi2di()} resolves randomly.
#' @param poly_eps Numeric epsilon to nudge non-positive edges created by
#'   \code{multi2di()}.
#' @param plot_pdf,pdf_width,pdf_height,pdf_auto,plot_cex Plotting controls.
#' @param logger Optional logger.
#' @return A list with \code{tree} (phylo|multiPhylo) and \code{log} (data.frame).
#' @examples
#' \dontrun{
#' res <- run_clade_grafting(
#'   backbone_path = "./backbone.tre",
#'   plan_path = "./graft-plan.tsv", # MRCA + Phylogeny_file_path
#'   authority = NULL,
#'   out_prefix = "./out/grafted",
#'   seed_mode = 42,
#'   chronos_select = "off",
#'   ultrametric_final = "extend",
#'   resolve_polytomies = "none"
#' )
#' }
#' @export
run_clade_grafting <- function(
  backbone_path,
  plan_path,
  authority = NULL,
  out_prefix = "ant_clade_chrono_grafted",
  seed_mode = "fixed",
  ingroup_anchors = NULL,
  drop_tips = NULL,
  ensure_ultrametric = c("nnls", "extend", "none"),
  chronos_select = c("auto", "fixed", "off"),
  chronos_model = "correlated",
  lambda = 1,
  lambda_grid = c(1e-1, 1, 10),
  models_grid = c("correlated", "relaxed", "discrete"),
  rate_cats = c(1, 2, 5),
  ultrametric_final = c("extend", "nnls", "none"),
  resolve_polytomies = c("none", "donor", "backbone", "both"),
  resolve_random = TRUE,
  poly_eps = 1e-8,
  calib_df = NULL,
  plot_pdf = TRUE,
  pdf_width = NULL,
  pdf_height = NULL,
  pdf_auto = TRUE,
  plot_cex = 0.35,
  logger = NULL
) {
  ensure_ultrametric <- match.arg(ensure_ultrametric)
  chronos_select <- match.arg(chronos_select)
  ultrametric_final <- match.arg(ultrametric_final)
  resolve_polytomies <- match.arg(resolve_polytomies)

  out_dir <- dirname(normalizePath(paste0(out_prefix, ".tre"), mustWork = FALSE))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_seed <- set_global_seed(seed_mode, out_dir = out_dir)

  if (is.null(logger) && exists(".make_logger", mode = "function")) {
    logger <- .make_logger(
      outdir = file.path(out_dir, "logs"),
      genus = "CladeGraft", mode = "graft",
      file_prefix = "clade_grafting"
    )
    on.exit(try(.close_logger(logger), silent = TRUE), add = TRUE)
  }

  log_section(logger, "load plan")
  plan_df <- tryCatch(
    {
      if (exists("read_table_utf8", mode = "function")) {
        read_table_utf8(plan_path, sep = "\t", header = TRUE, quote = "")
      } else {
        utils::read.delim(plan_path,
          sep = "\t", header = TRUE, stringsAsFactors = FALSE,
          quote = "", check.names = FALSE
        )
      }
    },
    error = function(e) stop("Failed to read plan TSV: ", conditionMessage(e))
  )
  log_msg(logger, "Plan rows: ", nrow(plan_df))

  has_col <- function(nm) tolower(nm) %in% tolower(names(plan_df))
  get_col <- function(nm) plan_df[[names(plan_df)[which(tolower(names(plan_df)) == tolower(nm))[1]]]]

  has_mrca <- has_col("MRCA")
  has_genus <- has_col("Genus")
  has_path <- has_col("Phylogeny_file_path")
  has_type <- has_col("Input_tree_type")

  if (!(has_path && (has_mrca || has_genus))) {
    stop("Plan must contain either (MRCA + Phylogeny_file_path) or (Genus + Phylogeny_file_path).")
  }

  log_section(logger, "assemble donors & placements")
  donors <- list()
  placements <- list()
  plan_dir <- dirname(normalizePath(plan_path, mustWork = FALSE))

  parse_csv <- function(x) unique(trimws(unlist(strsplit(x, ","))))

  n_ok <- 0
  n_skip <- 0
  for (i in seq_len(nrow(plan_df))) {
    if (has_mrca) {
      token <- trimws(as.character(get_col("MRCA")[i]))
      path <- trimws(as.character(get_col("Phylogeny_file_path")[i]))
      tree_type <- if (has_type) {
        x <- tolower(trimws(as.character(get_col("Input_tree_type")[i])))
        if (!nzchar(x) || is.na(x)) "phylogram" else x
      } else {
        "phylogram"
      }
      log_msg(logger, "Row ", i, ": Input_tree_type = ", tree_type)
      if (!tree_type %in% c("phylogram", "chronogram")) {
        stop(
          "Invalid Input_tree_type in row ", i,
          " (got '", tree_type,
          "'); must be 'phylogram' or 'chronogram'."
        )
      }
      if (!nzchar(token) || !nzchar(path)) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] Missing MRCA or file path")
        next
      }
      path <- trimws(path)
      if (!grepl("^(/|[A-Za-z]:)", path)) {
        if (grepl("^project/", path)) {
          path <- here::here(path)
        } else {
          path <- file.path(plan_dir, path)
        }
      }
      if (!file.exists(path)) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] file not found: ", path)
        next
      }
      if (grepl(",", token)) {
        anchors <- parse_csv(token)
        placement <- list(type = "clade", anchors = anchors)

        stem_mode <- if (has_col("Stem_mode")) {
          .normalize_stem_mode(get_col("Stem_mode")[i])
        } else {
          "outgroup"
        }
        placement$stem_mode <- stem_mode

        label <- paste0("MRCA{", paste(anchors, collapse = ","), "}")
      } else {
        placement <- list(type = "tip", tip = token)

        stem_mode <- if (has_col("Stem_mode")) {
          .normalize_stem_mode(get_col("Stem_mode")[i])
        } else {
          "outgroup"
        }
        placement$stem_mode <- stem_mode

        label <- paste0("TIP{", token, "}")
      }

      if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        tpl <- readRDS(path)

        if (!is.list(tpl) || !isTRUE(tpl$ok)) {
          stop("Template RDS is invalid: ", path)
        }

        donors[[label]] <- tpl
        placements[[label]] <- placement
        n_ok <- n_ok + 1L
        next
      }

      dtr <- tryCatch(ape::read.tree(path), error = function(e) e)
      attr(dtr, "source_path") <- path
      attr(dtr, "input_tree_type") <- tree_type
      if (inherits(dtr, "error")) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] failed to read donor tree: ", dtr$message)
        next
      }

      # If the MRCA cell is a single bare Genus (=> TIP mode), restrict ingroup to that genus only.
      # Otherwise (clade anchors or arbitrary label), fall back to all Genus_species-like tips.
      if (!grepl(",", token) && grepl("^[A-Z][a-z]+$", token)) {
        ig_hint <- grep(paste0("^", token, "_"), dtr$tip.label, value = TRUE)
      } else {
        ig_hint <- grep("^[A-Z][a-z]+_[a-z]+", dtr$tip.label, value = TRUE)
      }
      if (!length(ig_hint)) ig_hint <- dtr$tip.label
      attr(dtr, "ingroup_tips") <- unique(ig_hint)

      donors[[label]] <- dtr
      placements[[label]] <- placement
      n_ok <- n_ok + 1L
    } else {
      genus <- trimws(as.character(get_col("Genus")[i]))
      path <- trimws(as.character(get_col("Phylogeny_file_path")[i]))

      if (!nzchar(genus) || !nzchar(path)) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] Missing Genus or file path")
        next
      }

      if (!grepl("^(/|[A-Za-z]:)", path)) {
        if (grepl("^project/", path)) {
          path <- here::here(path)
        } else {
          path <- file.path(plan_dir, path)
        }
      }

      if (!file.exists(path)) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] file not found: ", path)
        next
      }

      ## ✅ placement ALWAYS defined first
      placement <- list(type = "tip", tip = genus)
      label <- genus

      ## ✅ safe early exit for templates
      if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        tpl <- readRDS(path)
        if (!is.list(tpl) || !isTRUE(tpl$ok)) {
          stop("Template RDS is invalid: ", path)
        }
        donors[[label]] <- tpl
        placements[[label]] <- placement
        n_ok <- n_ok + 1L
        next
      }

      ## raw donor tree
      dtr <- tryCatch(ape::read.tree(path), error = function(e) e)
      if (inherits(dtr, "error")) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] failed to read donor tree: ", dtr$message)
        next
      }

      tip_labels <- dtr$tip.label
      ig <- tip_labels[grepl(paste0("^", genus, "_"), tip_labels)]
      if (!length(ig)) {
        n_skip <- n_skip + 1L
        log_msg(logger, "[skip row ", i, "] no donor ingroup detected for Genus_ prefix: ", genus)
        next
      }

      attr(dtr, "ingroup_tips") <- unique(ig)

      donors[[label]] <- dtr
      placements[[label]] <- placement
      n_ok <- n_ok + 1L
    }
  }
  if (!length(donors)) stop("No usable donor rows found in plan (check file paths).")
  log_msg(logger, "Donors ready: ", n_ok, " (skipped: ", n_skip, ")")

  # Backbones -----------------------------------------------------------------
  log_section(logger, "load backbones")
  backbones <- read_trees_any(backbone_path)
  log_msg(logger, "Backbones: ", length(backbones))

  for (b in seq_along(backbones)) {
    tr <- backbones[[b]]

    # Optional: resolve backbone polytomies BEFORE grafting
    if (resolve_polytomies %in% c("backbone", "both")) {
      tr <- ape::multi2di(tr, random = isTRUE(resolve_random))
      if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
      tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- poly_eps
    }

    if (!is.null(ingroup_anchors) && length(ingroup_anchors)) {
      ext <- extract_ingroup_by_anchors(tr, ingroup_anchors)
      log_msg(logger, sprintf("Backbone %d: %s", b, ext$msg))
      tr <- ext$tree
    }

    if (!is.null(drop_tips) && length(drop_tips)) {
      before <- ape::Ntip(tr)
      tr <- safe_drop_tips(tr, drop_tips)
      log_msg(logger, sprintf("Backbone %d: drop-tips removed %d", b, before - ape::Ntip(tr)))
    }

    if (ensure_ultrametric != "none" && !ape::is.ultrametric(tr)) {
      tr <- phytools::force.ultrametric(tr, method = ensure_ultrametric)
      stopifnot(ape::is.ultrametric(tr))
      log_msg(logger, sprintf("Backbone %d: ultrametricized via %s", b, ensure_ultrametric))
    }
    backbones[[b]] <- tr
  }

  if (length(backbones) > 1) {
    class(backbones) <- "multiPhylo"
    attr(backbones, "order") <- NULL
  } else {
    backbones <- backbones[[1]]
  }

  # Grafting ------------------------------------------------------------------
  log_section(logger, "graft many clades")
  res <- graft_many_clades(
    backbones = backbones,
    donors = donors,
    placements = placements,
    authority = authority,
    chronos_select = chronos_select,
    chronos_model = chronos_model,
    lambda = lambda,
    lambda_grid = lambda_grid,
    models_grid = models_grid,
    rate_cats = rate_cats,
    calib_df = calib_df,
    seed = run_seed,
    eps = 1e-10,
    ultrametric_final = ultrametric_final,
    resolve_polytomies = resolve_polytomies,
    resolve_random = resolve_random,
    poly_eps = poly_eps,
    out_tree = paste0(out_prefix, ".tre"),
    out_log = paste0(out_prefix, "_graft_log.tsv"),
    logger = logger
  )

  # Outputs -------------------------------------------------------------------
  log_section(logger, "write plots & tip lists")
  write_one <- function(tr, idx = NULL) {
    suffix <- if (is.null(idx)) "" else sprintf("__t%02d", idx)
    base <- paste0(out_prefix, suffix)
    if (isTRUE(plot_pdf)) {
      plot_tree_autosize(
        tree = tr,
        pdf_path = paste0(base, ".pdf"),
        cex = plot_cex,
        pdf_width = pdf_width,
        pdf_height = pdf_height,
        pdf_auto = pdf_auto
      )
    }
    utils::write.table(
      data.frame(Tip = tr$tip.label),
      file = paste0(base, "_tips.txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
    )
  }

  if (inherits(res$tree, "multiPhylo")) {
    for (i in seq_along(res$tree)) write_one(res$tree[[i]], idx = i)
  } else {
    write_one(res$tree, idx = NULL)
  }

  invisible(res)
}
