# ---- internal helpers -------------------------------------------------------

#' Rescale phylogenetic branch lengths by a constant factor
#'
#' @param tree A `phylo` object or path to a Newick tree file.
#' @param factor Numeric scaling factor applied to all branch lengths.
#' @param out_path Optional output path for the rescaled tree.
#'
#' @return A rescaled `phylo` object.
#' @export
rescale_tree_time_units <- function(tree, factor = 100, out_path = NULL) {
  if (is.character(tree)) {
    tree <- ape::read.tree(tree)
  }

  tree$edge.length <- tree$edge.length * factor

  if (!is.null(tree$root.edge)) {
    tree$root.edge <- tree$root.edge * factor
  }

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ape::write.tree(phy = tree, file = out_path)
  }

  tree
}

#' Draw a fraction along an edge using a Beta(shape1,shape2) scaled to \[min,max\]
#' @keywords internal
#' @noRd
draw_depth_fraction <- function(shape1 = 1, shape2 = 1, min_frac = 0, max_frac = 1) {
  if (!is.finite(min_frac) || !is.finite(max_frac) || min_frac > max_frac) {
    min_frac <- 0
    max_frac <- 1
  }
  min_frac <- max(0, min(1, min_frac))
  max_frac <- max(0, min(1, max_frac))
  if (shape1 <= 0 || shape2 <= 0) {
    # fallback uniform
    u <- stats::runif(1L)
  } else {
    u <- stats::rbeta(1L, shape1 = shape1, shape2 = shape2)
  }
  min_frac + u * (max_frac - min_frac)
}

#' Total tree depth (max node height)
#' @keywords internal
#' @noRd
tree_depth <- function(tr) {
  max(phytools::nodeHeights(tr))
}

#' Bind a new tip at a given node and position; edge length computed from attach age
#' @keywords internal
#' @noRd
bind_tip_at <- function(tr, new_label, where_node, position, attach_age) {
  Htot <- tree_depth(tr)
  newlen <- max(0, Htot - attach_age)
  pos <- max(0, as.numeric(position))
  phytools::bind.tip(tr, tip.label = new_label, where = where_node, position = pos, edge.length = newlen)
}

#' Find the incoming edge of a node
#'
#' @param tr An object of class `phylo`.
#' @param node Integer node number.
#'
#' @return The row index of the incoming edge, or `integer(0)` for the root.
#' @keywords internal
incoming_edge_index <- function(tr, node) {
  which(tr$edge[, 2] == node)
}

# ---- exported primitives ----------------------------------------------------

#' Graft a new tip as sister to an existing tip
#'
#' Attaches a new tip labeled `new_label` as sister to the tip `sister_to` by splitting
#' the incoming edge of that tip at a random fraction (Beta-distributed, optionally
#' constrained to the interval defined by `min_frac` and `max_frac`). If the target tip is
#' an immediate child of the root (no incoming edge), the function attaches the new tip as
#' a polytomy at that node.
#'
#' @param tr A `phylo` tree.
#' @param new_label New tip label to add.
#' @param sister_to Existing tip label to which the new tip will be sister.
#' @param shape1,shape2 Beta distribution shape parameters.
#' @param min_frac,max_frac Min/max fraction along the edge from child toward parent.
#' @param .capture If `TRUE`, returns list with `tree` and detailed `log` (positions, ages).
#' @details
#' The new tip is attached by splitting the edge leading to \code{sister_to}.
#' A fraction \eqn{f} is drawn from a Beta(\code{shape1}, \code{shape2})
#' distribution and linearly scaled into the interval
#' \code{[min_frac, max_frac]}. The attachment position is then
#' \eqn{f \times L}, where \eqn{L} is the length of the target edge,
#' measured from the child tip toward its parent.
#'
#' If the target tip has no detectable incoming edge (a pathological or
#' degenerate case), the function falls back to attaching the new tip as a
#' polytomy at that node.
#'
#' When \code{.capture = TRUE}, the returned log reports the chosen edge,
#' sampled fraction, attachment age, insertion position, and the implied
#' new tip length.
#' @return Either a `phylo` tree or a list with `tree` and `log` if `.capture=TRUE`.
#' @examples
#' \dontrun{
#' tr <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
#' tr2 <- graft_sister_to_tip(tr, new_label = "X", sister_to = "A")
#' }
#' @export
graft_sister_to_tip <- function(tr, new_label, sister_to,
                                shape1 = 1, shape2 = 1,
                                min_frac = 0, max_frac = 1,
                                .capture = FALSE) {
  if (!sister_to %in% tr$tip.label) {
    stop(sprintf("Tip '%s' not found.", sister_to))
  }
  where <- which(tr$tip.label == sister_to)
  ei <- incoming_edge_index(tr, where)

  if (length(ei) == 0) {
    # attach as polytomy at the node
    Htot <- tree_depth(tr)
    tr2 <- phytools::bind.tip(tr, tip.label = new_label, where = where, position = 0, edge.length = Htot)
    if (.capture) {
      return(list(tree = tr2, log = list(
        edge_index = NA_integer_, edge_length = NA_real_,
        chosen_frac = NA_real_, position_from_child = 0,
        attach_age = Htot, where_node = where,
        new_tip_length = Htot, tree_depth = Htot
      )))
    } else {
      return(tr2)
    }
  }

  H <- phytools::nodeHeights(tr)
  L <- H[ei, 2] - H[ei, 1]
  frac <- draw_depth_fraction(shape1, shape2, min_frac, max_frac)
  pos <- if (L <= 0) 0 else frac * L
  aage <- H[ei, 2] - pos
  tr2 <- bind_tip_at(tr, new_label, where_node = where, position = pos, attach_age = aage)

  if (.capture) {
    Htot <- max(H)
    newlen <- Htot - aage
    return(list(tree = tr2, log = list(
      edge_index = ei, edge_length = as.numeric(L), chosen_frac = as.numeric(frac),
      position_from_child = as.numeric(pos), attach_age = as.numeric(aage),
      where_node = as.integer(where), new_tip_length = as.numeric(newlen),
      tree_depth = as.numeric(Htot)
    )))
  } else {
    return(tr2)
  }
}

#' Graft a new tip as sister to an existing clade (at the MRCA)
#'
#' Attaches a new tip labeled `new_label` as sister to the MRCA of `clade_tips` by splitting
#' the incoming edge of that MRCA at a random fraction. If the MRCA is the root (no incoming
#' edge), attaches the new tip as a polytomy at the root.
#'
#' @inheritParams graft_sister_to_tip
#' @param clade_tips Character vector of ≥2 tips defining the clade.
#' @details
#' The target clade is defined by the MRCA of \code{clade_tips}. The new tip
#' is attached by splitting the edge entering that MRCA. The split position is
#' drawn exactly as in \code{graft_sister_to_tip()}, i.e. from a Beta-distributed
#' fraction scaled to \code{[min_frac, max_frac]} and measured from the child
#' node (the MRCA) toward its parent.
#'
#' If the MRCA is the root and therefore has no incoming edge, the new tip is
#' attached as a root polytomy.
#'
#' When \code{.capture = TRUE}, the returned log reports the chosen edge,
#' sampled fraction, attachment age, insertion position, and the implied
#' new tip length.

#' @export
graft_sister_to_clade <- function(tr, new_label, clade_tips,
                                  shape1 = 1, shape2 = 1,
                                  min_frac = 0, max_frac = 1,
                                  .capture = FALSE) {
  clade_tips <- intersect(clade_tips, tr$tip.label)
  if (length(clade_tips) < 2) {
    stop("clade_tips must include at least two existing tips from the tree.")
  }
  where <- ape::getMRCA(tr, clade_tips)
  if (is.null(where) || length(where) == 0) {
    stop("Could not find MRCA for clade_tips.")
  }

  ei <- incoming_edge_index(tr, where)
  if (length(ei) == 0) {
    # MRCA is the root -> polytomy at root
    Htot <- tree_depth(tr)
    tr2 <- phytools::bind.tip(tr, tip.label = new_label, where = where, position = 0, edge.length = Htot)
    if (.capture) {
      return(list(tree = tr2, log = list(
        edge_index = NA_integer_, edge_length = NA_real_, chosen_frac = NA_real_,
        position_from_child = 0, attach_age = Htot, where_node = where,
        new_tip_length = 0, tree_depth = Htot
      )))
    } else {
      return(tr2)
    }
  }

  H <- phytools::nodeHeights(tr)
  L <- H[ei, 2] - H[ei, 1]
  frac <- draw_depth_fraction(shape1, shape2, min_frac, max_frac)
  pos <- if (L <= 0) 0 else frac * L
  aage <- H[ei, 2] - pos
  tr2 <- bind_tip_at(tr, new_label, where_node = where, position = pos, attach_age = aage)

  if (.capture) {
    Htot <- max(H)
    newlen <- Htot - aage
    return(list(tree = tr2, log = list(
      edge_index = ei, edge_length = as.numeric(L), chosen_frac = as.numeric(frac),
      position_from_child = as.numeric(pos), attach_age = as.numeric(aage),
      where_node = as.integer(where), new_tip_length = as.numeric(newlen),
      tree_depth = as.numeric(Htot)
    )))
  } else {
    return(tr2)
  }
}

#' Graft a new tip at a random edge within a clade
#'
#' Samples an edge inside the MRCA-defined clade (optionally including terminal and/or the stem
#' edge) with probability proportional to the edge length, then attaches the new tip at a random
#' position along that edge.
#'
#' @param include_terminal Logical; if `FALSE`, excludes terminal edges.
#' @param include_stem Logical; if `TRUE`, allows grafting on the edge entering the MRCA (stem).
#' @param return_log If `TRUE`, returns `list(tree, log=list(...))` with details.
#' @details
#' Randomness in this function is controlled by the global RNG state.
#' To obtain reproducible placements, call \code{set.seed()} (or a package-level
#' workflow seed helper) before invoking this function or any higher-level wrapper
#' that uses it.
#'
#' Candidate edges are sampled from the clade descending from the MRCA of
#' \code{clade_tips}. By default, all descendant edges are eligible, including
#' terminal edges and the edge entering the MRCA (the stem edge).
#'
#' If \code{include_terminal = FALSE}, terminal edges are excluded. If
#' \code{include_stem = TRUE}, the edge entering the MRCA is added to the
#' candidate set even though it lies immediately above the crown of the clade.
#'
#' Candidate edges are sampled with probability proportional to their branch
#' lengths. Once an edge is selected, an attachment fraction is drawn from a
#' Beta(\code{shape1}, \code{shape2}) distribution scaled to the interval
#' \code{[min_frac, max_frac]}, and the new tip is inserted at that position
#' measured from the child toward the parent.
#'
#' When \code{return_log = TRUE}, the function returns a list containing the
#' grafted tree and a detailed log, including the chosen edge index, edge type
#' (\code{"stem"}, \code{"internal"}, or \code{"terminal"}), sampled fraction,
#' attachment age, insertion node, and implied new tip length.
#' @examples
#' \dontrun{
#' set.seed(42)
#' tr <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
#' tr2 <- graft_within_clade_random(tr, new_label = "X", clade_tips = c("A", "B"))
#' }
#' @inheritParams graft_sister_to_tip
#' @export

graft_within_clade_random <- function(tr, new_label, clade_tips,
                                      include_terminal = TRUE,
                                      include_stem = TRUE,
                                      shape1 = 1, shape2 = 1,
                                      min_frac = 0, max_frac = 1,
                                      return_log = FALSE) {
  if (new_label %in% tr$tip.label) {
    stop(sprintf("Label '%s' already exists in the tree.", new_label))
  }
  clade_tips <- intersect(clade_tips, tr$tip.label)
  if (length(clade_tips) < 2) stop("clade_tips must include at least two existing tips from the tree.")
  mrca <- ape::getMRCA(tr, clade_tips)
  if (is.null(mrca) || length(mrca) == 0) stop("Could not find MRCA for clade_tips.")

  H <- phytools::nodeHeights(tr)
  S <- unique(stats::na.omit(c(mrca, phytools::getDescendants(tr, mrca))))
  e_inside <- which(tr$edge[, 1] %in% S & tr$edge[, 2] %in% S)

  if (!include_terminal) {
    tip_nodes <- seq_len(ape::Ntip(tr))
    e_inside <- e_inside[!(tr$edge[e_inside, 2] %in% tip_nodes)]
    if (length(e_inside) == 0 && !include_stem) {
      stop("No internal edges in this clade (consider include_terminal=TRUE or include_stem=TRUE).")
    }
  }
  stem_ei <- integer(0)
  if (include_stem) stem_ei <- which(tr$edge[, 2] == mrca)
  candidates <- unique(c(e_inside, stem_ei))
  if (length(candidates) == 0) stop("No candidate edges to sample.")

  Ls <- H[candidates, 2] - H[candidates, 1]
  if (any(Ls <= 0)) stop("Non-positive edge length among candidates.")
  probs <- Ls / sum(Ls)
  i <- sample(candidates, size = 1, prob = probs)
  Li <- H[i, 2] - H[i, 1]
  frac <- draw_depth_fraction(shape1, shape2, min_frac, max_frac)
  pos <- if (Li <= 0) 0 else frac * Li
  aage <- H[i, 2] - pos
  where <- tr$edge[i, 2]
  Htot <- max(H)
  newlen <- Htot - aage
  tr2 <- phytools::bind.tip(tr, tip.label = new_label, where = where, position = pos, edge.length = newlen)

  if (return_log) {
    edge_type <- if (length(stem_ei) == 1L && i == stem_ei) "stem" else if (tr$edge[i, 2] <= ape::Ntip(tr)) "terminal" else "internal"
    out <- list(
      tree = tr2,
      log = list(
        chosen_edge_index = i,
        edge_type = edge_type,
        is_stem_included = include_stem,
        is_terminal_included = include_terminal,
        edge_length = Li,
        chosen_frac = frac,
        position_from_child = pos,
        attach_age = aage,
        where_node = where,
        tree_depth = Htot,
        new_tip_length = newlen,
        mrca_node = mrca,
        stem_edge_index = if (length(stem_ei)) stem_ei else NA_integer_
      )
    )
    return(out)
  } else {
    return(tr2)
  }
}

#' @importFrom ape getMRCA drop.tip Ntip write.tree ladderize read.tree extract.clade is.ultrametric
#' @importFrom phytools nodeHeights getDescendants force.ultrametric bind.tip
#' @importFrom stats rbeta runif
NULL
