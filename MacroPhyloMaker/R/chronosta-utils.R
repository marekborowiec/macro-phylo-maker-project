# =============================================================================
# Namespace aliases for standalone sourcing and package-internal use
# =============================================================================
# The original ChronoSTA helper functions were written as project scripts and used
# unqualified calls such as str_split_fixed(), Ntip(), nodeHeights(), etc. When the
# functions are moved into a package, packages listed in DESCRIPTION are installed
# but not necessarily attached. These aliases keep the original code working
# without requiring library(stringr), library(ape), etc. in the user's session.
str_detect       <- stringr::str_detect
str_split_fixed  <- stringr::str_split_fixed
str_split_1      <- stringr::str_split_1
str_remove_all   <- stringr::str_remove_all
str_replace_all  <- stringr::str_replace_all

Ntip             <- ape::Ntip
Nnode            <- ape::Nnode
read.tree        <- ape::read.tree
write.tree       <- ape::write.tree
drop.tip         <- ape::drop.tip
extract.clade    <- ape::extract.clade
getMRCA          <- ape::getMRCA
nodepath         <- ape::nodepath
bind.tree        <- ape::bind.tree
is.ultrametric   <- ape::is.ultrametric
chronos          <- ape::chronos
chronos.control  <- ape::chronos.control

nodeHeights      <- phytools::nodeHeights
force.ultrametric <- phytools::force.ultrametric
getDescendants   <- phytools::getDescendants

AssessMonophyly  <- MonoPhy::AssessMonophyly
GetResultMonophyly <- MonoPhy::GetResultMonophyly
PlotMonophyly    <- MonoPhy::PlotMonophyly

graph_from_data_frame <- igraph::graph_from_data_frame
bipartite_projection  <- igraph::bipartite_projection
components            <- igraph::components
V                     <- igraph::V


check_and_recalibrate <- function(tree_list,
                                  tree_ref,
                                  tolerance_pct   = 0.05,  # 5% relative tolerance for the check
                                  tolerance_min   = 1.0,   # min Ma window passed to recalibration
                                  calib_ratio     = 0.3,
                                  calib_min_nodes = 3,
                                  calib_max_nodes = 50,
                                  lambda          = 1) {
  
  if (inherits(tree_list, "phylo")) tree_list <- list(tree_list)
  
  tree_names <- names(tree_list)
  if (is.null(tree_names)) tree_names <- paste0("tree_", seq_along(tree_list))
  
  recalibrated <- vector("list", length(tree_list))
  names(recalibrated) <- tree_names
  
  recal_names <- character(0)   # track which trees were recalibrated
  
  for (i in seq_along(tree_list)) {
    
    nm  <- tree_names[i]
    tre <- tree_list[[i]]
    
    # --- Identify shared nodes with the reference tree -----------------------
    shared <- suppressMessages(get_shared_nodes(
      tree_ref        = tree_ref,
      tree_target     = tre,
      calib_ratio     = calib_ratio,
      calib_min_nodes = calib_min_nodes,
      calib_max_nodes = calib_max_nodes
    ))
    
    if (is.null(shared) || nrow(shared) == 0) {
      recalibrated[[i]] <- tre
      next
    }
    
    # --- Compare ages: relative deviation per shared node --------------------
    rel_dev <- ifelse(shared$age_ref > 0,
                      abs(shared$age_target - shared$age_ref) / shared$age_ref,
                      0)
    
    n_off   <- sum(rel_dev > tolerance_pct)
    max_dev <- max(rel_dev, na.rm = TRUE)
    
    if (n_off == 0) {
      recalibrated[[i]] <- tre
      next
    }
    
    # --- Recalibration needed: one info line per recalibrated tree -----------
    cat("  Recalibrating '", nm, "' (", i, "/", length(tree_list), "): ",
        n_off, "/", nrow(shared), " node(s) off, max deviation ",
        round(max_dev * 100, 1), "%\n", sep = "")
    
    recalibrated[[i]] <- suppressMessages(suppressWarnings(
      recalibrate_with_chronos(
        tree_ref        = tree_ref,
        tree_target     = tre,
        shared_nodes    = shared,
        tolerance_pct   = tolerance_pct,
        tolerance_min   = tolerance_min,
        calib_min_nodes = calib_min_nodes,
        lambda          = lambda
      )
    ))
    recal_names <- c(recal_names, nm)
  }
  
  # --- Final summary line ---------------------------------------------------
  if (length(recal_names) == 0) {
    cat("No recalibration needed: all trees within ",
        round(tolerance_pct * 100), "% tolerance.\n", sep = "")
  } else {
    cat(length(recal_names), "/", length(tree_list),
        " tree(s) recalibrated.\n", sep = "")
  }
  
  return(recalibrated)
}

####################################### Filter consistent calibrations (no incompatibilities) #######################################


filter_consistent_calibrations <- function(tree, calib, min_margin = 0.5) {
  
  if (nrow(calib) <= 1) return(calib)
  
  node_to_row <- setNames(seq_len(nrow(calib)), calib$node)
  
  get_parent <- function(nd) {
    row <- which(tree$edge[, 2] == nd)
    if (length(row) == 0) return(NA_integer_)  # nd is root
    tree$edge[row, 1]
  }
  
  # Nearest calibrated ancestor of a node (climbs until it hits a calib node)
  nearest_calib_ancestor <- function(nd, calib_nodes) {
    p <- get_parent(nd)
    while (!is.na(p)) {
      if (p %in% calib_nodes) return(p)
      p <- get_parent(p)
    }
    NA_integer_
  }
  
  # Node ages in target tree (for discordance scoring + midpoints)
  calib$node_age_target <- sapply(calib$node, function(nd) node_age(tree, nd))
  calib$mid             <- (calib$age.min + calib$age.max) / 2
  calib$discordance     <- abs(calib$node_age_target - calib$mid)
  
  removed_hard <- character(0)
  removed_tight <- character(0)
  
  # ---------------------------------------------------------------------------
  # PASS 1 — remove HARD conflicts (child age.min > parent age.max)
  # (unchanged behaviour from the original function)
  # ---------------------------------------------------------------------------
  max_iter <- nrow(calib)
  for (i in seq_len(max_iter)) {
    if (nrow(calib) <= 1) break
    
    conflicts <- lapply(seq_len(nrow(calib)), function(r) {
      nd     <- calib$node[r]
      parent <- get_parent(nd)
      if (is.na(parent)) return(NULL)
      p_row <- which(calib$node == parent)
      if (length(p_row) == 0) return(NULL)
      if (calib$age.min[r] > calib$age.max[p_row]) {
        return(data.frame(r_child = r, r_parent = p_row))
      }
      NULL
    })
    conflicts <- do.call(rbind, Filter(Negate(is.null), conflicts))
    if (is.null(conflicts) || nrow(conflicts) == 0) break
    
    conflict_nodes <- unique(c(conflicts$r_child, conflicts$r_parent))
    worst <- conflict_nodes[which.max(calib$discordance[conflict_nodes])]
    removed_hard <- c(removed_hard, as.character(calib$node[worst]))
    calib <- calib[-worst, ]
  }
  
  # ---------------------------------------------------------------------------
  # PASS 2 — enforce minimum ordering margin between nested calibrations
  # Remove the CHILD of any parent-child calibrated pair whose midpoints are
  # closer than min_margin. Iterate until all nested pairs have enough margin.
  # ---------------------------------------------------------------------------
  repeat {
    if (nrow(calib) <= 1) break
    calib_nodes <- calib$node
    
    tight <- lapply(seq_len(nrow(calib)), function(r) {
      nd  <- calib$node[r]
      anc <- nearest_calib_ancestor(nd, calib_nodes)
      if (is.na(anc)) return(NULL)
      a_row <- which(calib$node == anc)
      margin <- calib$mid[a_row] - calib$mid[r]   # parent_mid - child_mid
      # Too tight if the ordering margin is below the threshold
      # (covers both small positive margins and inverted midpoints)
      if (margin < min_margin) {
        return(data.frame(r_child = r, margin = margin))
      }
      NULL
    })
    tight <- do.call(rbind, Filter(Negate(is.null), tight))
    if (is.null(tight) || nrow(tight) == 0) break
    
    # Remove the child of the tightest pair (smallest margin).
    # The child is shallower (smaller clade) and less reliable than its
    # calibrated ancestor, so dropping it preserves the deeper anchor.
    worst_child <- tight$r_child[which.min(tight$margin)]
    removed_tight <- c(removed_tight, as.character(calib$node[worst_child]))
    calib <- calib[-worst_child, ]
  }
  
  # ---------------------------------------------------------------------------
  # Reporting
  # ---------------------------------------------------------------------------
  if (length(removed_hard) > 0) {
    message("  consistency filter removed ", length(removed_hard),
            " hard-conflict node(s): ", paste(removed_hard, collapse = ", "))
  }
  if (length(removed_tight) > 0) {
    message("  margin filter removed ", length(removed_tight),
            " tightly-nested node(s) (< ", min_margin, " Ma apart): ",
            paste(removed_tight, collapse = ", "))
  }
  
  # Drop helper columns
  calib$node_age_target <- NULL
  calib$mid             <- NULL
  calib$discordance     <- NULL
  
  return(calib)
}

#####################################################################################################################################
######################################## Get the shared nodes of trees ##############################################################
get_shared_nodes <- function(tree_ref, tree_target,
                             calib_ratio, 
                             calib_min_nodes,   
                             calib_max_nodes) 
{
  
  tips_ref    <- tree_ref$tip.label
  tips_target <- tree_target$tip.label
  shared_tips <- intersect(tips_ref, tips_target)
  
  if (length(shared_tips) < 2) {
    message("  fewer than 2 shared taxa — recalibration skipped.")
    return(NULL)
  }
  
  n_tips_target  <- length(tree_target$tip.label)
  internal_nodes <- (n_tips_target + 1):max(tree_target$edge)
  
  results <- lapply(internal_nodes, function(nd) {
    
    # Tip descendants of this node in the target tree
    desc_target <- tryCatch(
      tree_target$tip.label[
        phangorn::Descendants(tree_target, nd, type = "tips")[[1]]
      ],
      error = function(e) character(0)
    )
    
    desc_shared <- intersect(desc_target, shared_tips)
    if (length(desc_shared) < 2) return(NULL)
    
    # Corresponding MRCA in the reference tree
    mrca_ref <- getMRCA(tree_ref, desc_shared)
    if (is.null(mrca_ref)) return(NULL)
    
    data.frame(
      node_target   = nd,
      node_ref      = mrca_ref,
      age_target    = node_age(tree_target, nd),
      age_ref       = node_age(tree_ref, mrca_ref),
      n_desc_shared = length(desc_shared)
    )
  })
  
  results <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(results) || nrow(results) == 0) return(NULL)
  
  # Deduplicate: if two target nodes map to the same reference MRCA,
  # keep the one with the most shared descendants (maximum information)
  # In case of exact duplication, keep the youngest node as it mean that it is a pair of species 
  results <- results[order(results$node_ref, -results$n_desc_shared, results$age_target), ]
  
  results <- results[!duplicated(results$node_ref), ]
  
  # Compute adaptive cap from tree size and per-level ratio, then clamp
  n_internal  <- tree_target$Nnode
  ratio       <- calib_ratio
  n_cap       <- max(calib_min_nodes,
                     min(calib_max_nodes, round(n_internal * ratio)))
  
  n_candidates <- nrow(results)
  
  if (n_candidates > n_cap) {
    # Retain the n_cap most informative nodes (most shared descendants)
    results <- results[order(-results$n_desc_shared), ][seq_len(n_cap), ]
    message("  calibration nodes: ", n_cap, " retained out of ", n_candidates,
            " candidates ( ratio = ", ratio,
            ", n_internal = ", n_internal, ")")
  } else {
    message("  calibration nodes: ", n_candidates, " (all retained — below cap of ",
            n_cap)
  }
  
  return(results)
}
#####################################################################################################################################
# ---------------------------------------------------------------------------
# Helper: extract genus from a tip label (prefix before `genus_sep`)
# ---------------------------------------------------------------------------

get_genus <- function(labels, genus_sep = "_") {
  sub(paste0("\\Q", genus_sep, "\\E", ".*$"), "", labels, perl = TRUE)
}

make_taxonomy <- function(tree, genus_sep = "_") {
  data.frame(
    tip   = tree$tip.label,
    genus = get_genus(tree$tip.label, genus_sep),
    stringsAsFactors = FALSE
  )
}
# ---------------------------------------------------------------------------
# Helper: which genera are non-monophyletic in a single tree?
#   - Only genera with >= 2 tips can be non-monophyletic (monotypic excluded).
#   - Returns a character vector of genus names.
# ---------------------------------------------------------------------------
nonmono_genera <- function(tree, genus_sep = "_") {
  if (is.null(tree) || Ntip(tree) < 3) return(character(0))
  tax <- make_taxonomy(tree, genus_sep)
  # genera represented by >= 2 tips in THIS tree
  multi <- names(which(table(tax$genus) >= 2))
  if (length(multi) == 0) return(character(0))
  sol <- suppressWarnings(
    MonoPhy::AssessMonophyly(tree, taxonomy = tax, verbosity = 0)
  )
  res <- MonoPhy::GetResultMonophyly(sol)[[1]]
  # Result table: rownames = genera, column "Monophyly" in {"Yes","No","Monotypic"}
  status <- res[, "Monophyly"]
  names(status) <- rownames(res)
  bad <- names(status)[status == "No"]
  intersect(bad, multi)
}

# ---------------------------------------------------------------------------
# Build the protected baseline:
#   union of non-monophyletic genera across reference + all source trees.
# Labels must already be harmonised across trees (flagged, not automated).
# ---------------------------------------------------------------------------
compute_baseline_nonmono <- function(ref_tree,
                                      source_trees = list(),
                                      genus_sep = "_") {
  baseline <- nonmono_genera(ref_tree, genus_sep)
  for (st in source_trees) {
    baseline <- union(baseline, nonmono_genera(st, genus_sep))
  }
  sort(unique(baseline))
}


# ---------------------------------------------------------------------------
# Core optimiser: minimal-removal plan for ONE non-monophyletic genus.
#
# Strategy (exact over candidate cores):
#   The retained "core" clade of genus G must be the subtree rooted at some
#   internal node that is an ancestor of >=1 tip of G. For each such node we
#   compute:
#       drop = (G tips OUTSIDE the clade)  U  (non-G tips INSIDE the clade)
#   subject to:
#       - core keeps >= min_tips_per_genus tips of G
#       - removing the intruders never empties their own genus elsewhere
#   We pick the node minimising length(drop); ties broken by the SMALLEST
#   resulting clade (most specific), then by node id (determinism).
#
# Returns list(core_node, to_remove, n_removed, kept_G) or NULL if no plan
# satisfies the no-empty-genus constraints.
# ---------------------------------------------------------------------------
min_removal_plan <- function(tree, genus, genus_sep = "_",
                             min_tips_per_genus = 1) {
  tip_genus <- get_genus(tree$tip.label, genus_sep)
  G_tips    <- tree$tip.label[tip_genus == genus]
  if (length(G_tips) < 2) return(NULL)
  
  G_idx <- which(tree$tip.label %in% G_tips)
  
  # genus sizes (to protect intruders' own genera from being emptied)
  genus_counts <- table(tip_genus)
  
  # candidate core nodes = all internal nodes that are ancestors of any G tip,
  # plus the G tips' direct context. Practically: all nodes on the paths from
  # each G tip to the root. Dedup.
  root <- Ntip(tree) + 1
  anc_nodes <- unique(unlist(lapply(G_idx, function(ti) {
    nodepath(tree, from = ti, to = root)
  })))
  # keep only internal nodes (a single-tip "clade" can't host >=2 G tips
  # unless min_tips==1; we still allow tip-level cores when min==1)
  anc_nodes <- sort(unique(anc_nodes))
  
  best <- NULL
  for (nd in anc_nodes) {
    if (nd <= Ntip(tree)) {
      clade_tips <- tree$tip.label[nd]
    } else {
      desc <- phytools::getDescendants(tree, nd)
      clade_tips <- tree$tip.label[desc[desc <= Ntip(tree)]]
    }
    G_in   <- intersect(clade_tips, G_tips)
    if (length(G_in) < min_tips_per_genus) next
    G_out  <- setdiff(G_tips, G_in)            # outliers of G to drop
    nonG_in <- setdiff(clade_tips, G_in)       # intruders to drop
    
    # constraint: dropping intruders must not empty their own genus
    intr_genus <- get_genus(nonG_in, genus_sep)
    drop_per_genus <- table(intr_genus)
    empties <- vapply(names(drop_per_genus), function(h) {
      as.integer(drop_per_genus[[h]]) >= as.integer(genus_counts[[h]])
    }, logical(1))
    if (any(empties)) next   # this core would wipe out another genus -> skip
    
    to_remove <- c(G_out, nonG_in)
    cost <- length(to_remove)
    clade_size <- length(clade_tips)
    
    cand <- list(core_node = nd, to_remove = to_remove,
                 n_removed = cost, kept_G = length(G_in),
                 clade_size = clade_size)
    if (is.null(best) ||
        cost < best$n_removed ||
        (cost == best$n_removed && clade_size < best$clade_size) ||
        (cost == best$n_removed && clade_size == best$clade_size &&
         nd < best$core_node)) {
      best <- cand
    }
  }
  best
}


####################################### Get the age of a node ###########################################################
node_age <- function(tr, node) {
  nh       <- nodeHeights(tr)
  root_h   <- max(nh)
  node_row <- which(tr$edge[, 2] == node)
  if (length(node_row) == 0) return(root_h)  # node is the root
  return(root_h - nh[node_row[1], 2])
}
#####################################################################################################################################
####################################### Recalibrate a tree with penalized likelihood (chronos) ######################################
recalibrate_with_chronos <- function(tree_ref,
                                     tree_target,
                                     shared_nodes,
                                     tolerance_pct = tolerance_pct,   
                                     tolerance_min = tolerance_min,
                                     calib_min_nodes = calib_min_nodes,
                                     lambda = lambda,
                                     margin_nod = 0.5)    
{
  
  shared <- shared_nodes
  
  if (is.null(shared) || nrow(shared) == 0) {
    warning("  no shared nodes found — tree returned without recalibration.")
    return(tree_target)
  }
  
  # Build calibration constraints with adaptive tolerance window
  half_window <- pmax(shared$age_ref * tolerance_pct, tolerance_min)
  
  calib <- data.frame(
    node        = shared$node_target,
    age.min     = pmax(shared$age_ref - half_window, 0),
    age.max     = shared$age_ref + half_window,
    soft.bounds = FALSE
  )
  calib <- calib[!duplicated(calib$node), ]
  
  # Remove parent-child age conflicts before any chronos call
  calib <- filter_consistent_calibrations(tree_target, calib, min_margin = margin_nod)
  
  if (nrow(calib) < 2) {
    warning("  fewer than 2 consistent calibration nodes remain after ",
            "conflict filtering — tree returned without recalibration.")
    return(tree_target)
  }
  
  # Helper: attempt one chronos() call, return NULL on failure
  try_chronos <- function(tr, cal, lam) {
    tryCatch(
      chronos(tr, lambda = lambda, model = "corr", calibration = cal,
              control = chronos.control(nb.rate.cat = 10,
                                        dual.iter.max = 100)),
      error = function(e) NULL
    )
  }
  
  # --- Fallback 1: hard constraints, fixed lambda ----------------------------
  best_lambda <- lambda
  result      <- try_chronos(tree_target, calib, best_lambda)
  if (!is.null(result)) return(result)
  message("  fallback 1 failed — trying soft bounds...")
  
  # --- Fallback 2: soft bounds, same constraints and lambda -------------------
  calib_soft              <- calib
  calib_soft$soft.bounds  <- TRUE
  result <- try_chronos(tree_target, calib_soft, best_lambda)
  if (!is.null(result)) {
    message("  fallback 2 succeeded (soft bounds).")
    return(result)
  }
  message("  fallback 2 failed — trying deepest nodes only...")
  
  # --- Fallback 3: keep only the calib_min_nodes deepest calibration nodes ---
  # Deep nodes span the largest clades and are least likely to conflict.
  # Use age_ref (reference age) as proxy for depth — older = deeper.
  calib_deep <- calib[order(-shared$age_ref[match(calib$node,
                                                  shared$node_target)]), ]
  calib_deep <- calib_deep[seq_len(min(calib_min_nodes, nrow(calib_deep))), ]
  calib_deep$soft.bounds <- FALSE
  
  result <- try_chronos(tree_target, calib_deep, lambda)
  if (!is.null(result)) {
    message("  fallback 3 succeeded (", nrow(calib_deep), " deepest nodes, ",
            "hard constraints, fallback lambda = ", lambda, ").")
    return(result)
  }
  
  # --- Fallback 4: give up, return uncalibrated tree -------------------------
  warning("  all chronos() fallbacks failed for this tree.\n",
          "  Returning uncalibrated tree — CHECK THIS TREE MANUALLY.\n",
          "  n_calib = ", nrow(calib), ", lambda tried = ", best_lambda)
  return(tree_target)
}

#####################################################################################################################################

# =============================================================================
# restore_monophyly.R
# -----------------------------------------------------------------------------
# Restore genus-level monophyly on a grafted supertree (post Chrono-STA),
# under four constraints:
#   1. "Legitimate" non-monophylies (present in reference OR any source tree)
#      are PRESERVED untouched (both the offending tips and the affected clade).
#   2. NEW non-monophylies introduced by the merge are removed.
#   3. No genus is ever fully eliminated (>= min_tips_per_genus tips kept).
#   4. Removal is MINIMAL: for each genus we drop the fewest tips that restore
#      monophyly (exact search over candidate "core" clades), and we re-assess
#      after each correction so nested cases collapse naturally.
#
# Deterministic: no random component. Genera processed in a fixed order.


restore_monophyly <- function(final_tree,
                                 ref_tree,
                                 source_trees = list_tre,
                                 genus_sep = "_",
                                 min_tips_per_genus = 1,
                                 max_passes = 1000,
                                 verbose = TRUE,
                                 exc = NA) {
  
  stopifnot(inherits(final_tree, "phylo"), inherits(ref_tree, "phylo"))
  
  # 1. protected baseline (legitimate non-monophylies)
  protected <- compute_baseline_nonmono(ref_tree, source_trees, genus_sep)
  
  
  #Exceptions
  if(length(exc) >0){
    protected <- protected[!protected %in% exc]
  }
  
  if (verbose) message("Protected (legitimate) non-mono genera: ",
                       if (length(protected)) paste(protected, collapse = ", ")
                       else "<none>")
  
  tree <- final_tree
  log_rows <- list()
  frozen   <- character(0)   # genera skipped (uncorrectable w/o genus loss)
  pass <- 0L
  
  repeat {
    pass <- pass + 1L
    if (pass > max_passes) {
      warning("max_passes reached; stopping with remaining violations.")
      break
    }
    
    current_bad <- nonmono_genera(tree, genus_sep)
    # exclude protected baseline and frozen genera
    new_bad <- setdiff(current_bad, union(protected, frozen))
    if (length(new_bad) == 0) break
    
    # process smallest violation first (helps nested cases resolve cleanly).
    # "smallest" = fewest tips of that genus; ties -> alphabetical (determinism)
    tip_genus <- get_genus(tree$tip.label, genus_sep)
    sizes <- vapply(new_bad, function(g) sum(tip_genus == g), integer(1))
    ord   <- order(sizes, new_bad)
    g     <- new_bad[ord[1]]
    
    plan <- min_removal_plan(tree, g, genus_sep, min_tips_per_genus)
    
    if (is.null(plan) || plan$n_removed == 0) {
      frozen <- c(frozen, g)
      log_rows[[length(log_rows) + 1L]] <- data.frame(
        pass = pass, genus = g, action = "skipped",
        reason = "no plan without emptying a genus",
        tip = NA_character_, tip_genus = NA_character_,
        core_node = NA_integer_, stringsAsFactors = FALSE)
      if (verbose) message("[pass ", pass, "] ", g,
                           " -> skipped (uncorrectable without genus loss)")
      next
    }
    
    rem_genus <- get_genus(plan$to_remove, genus_sep)
    log_rows[[length(log_rows) + 1L]] <- data.frame(
      pass = pass, genus = g, action = "drop_tip",
      reason = sprintf("minimal removal (kept %d of %s; core node %d)",
                       plan$kept_G, g, plan$core_node),
      tip = plan$to_remove, tip_genus = rem_genus,
      core_node = plan$core_node, stringsAsFactors = FALSE)
    
    tree <- ape::drop.tip(tree, plan$to_remove)
    if (verbose) message("[pass ", pass, "] ", g, " -> removed ",
                         plan$n_removed, " tip(s): ",
                         paste(plan$to_remove, collapse = ", "))
  }
  
  removed_df <- if (length(log_rows)) do.call(rbind, log_rows) else
    data.frame(pass = integer(0), genus = character(0), action = character(0),
               reason = character(0), tip = character(0),
               tip_genus = character(0), core_node = integer(0))
  
  # sanity: confirm no genus was eliminated
  before_genera <- unique(get_genus(final_tree$tip.label, genus_sep))
  after_genera  <- unique(get_genus(tree$tip.label, genus_sep))
  lost <- setdiff(before_genera, after_genera)
  if (length(lost)) warning("Genera lost (should be none): ",
                            paste(lost, collapse = ", "))
  
  list(
    tree       = tree,
    removed    = removed_df,
    protected  = protected,
    frozen     = unique(frozen),
    lost_genera = lost,
    n_tips_before = Ntip(final_tree),
    n_tips_after  = Ntip(tree)
  )
}

#####################################################################################################################################
################################################## Chrono-STA merging ###############################################################
patch_chronosta_readonly_array <- function(path, logger = NULL) {
  if (!file.exists(path)) {
    stop("Cannot patch ChronoSTA script; file not found: ", path, call. = FALSE)
  }

  x <- readLines(path, warn = FALSE)

  # Already patched
  if (any(grepl("m\\.values\\.copy\\(\\)", x)) ||
      any(grepl("m\\.to_numpy\\(copy=True\\)", x))) {
    .csta_msg("ChronoSTA read-only-array patch already present: ", path, logger = logger)
    return(invisible(TRUE))
  }

  old <- "np.fill_diagonal(m.values, 0.0)"
  hit <- grep(old, x, fixed = TRUE)

  if (!length(hit)) {
    warning(
      "Could not find expected ChronoSTA line to patch: ",
      old,
      "\nFile: ",
      path
    )
    return(invisible(FALSE))
  }

  indent <- sub("^(\\s*).*", "\\1", x[hit[1]])

  replacement <- c(
    paste0(indent, "mat = m.values.copy()"),
    paste0(indent, "np.fill_diagonal(mat, 0.0)"),
    paste0(indent, "m.iloc[:, :] = mat")
  )

  x <- append(x[-hit], values = replacement, after = hit[1] - 1L)
  writeLines(x, path)

  .csta_msg(
    "Applied ChronoSTA read-only-array compatibility patch to temporary script: ",
    path,
    logger = logger
  )

  invisible(TRUE)
}

# Run Chrono-STA
# The reference tree is duplicated ref_weight times to increase its weightin the averaged supermatrix.

run_chrono_sta <- function(tree_list,
                           ref_name = "NONE",
                           newname,
                           ref_weight = 1,
                           sta_way = NA,
                           csta_folder = NA, 
                           outp_way = NA,
                           outp_mes = TRUE) {
  
  sta_way <- paste0(sta_way, "temp/")
  
  #Should output be displayed
  if(outp_mes == TRUE){
    outp_mes <- ""
  }
  
  #Create a new dir to store the phylogeny
  if(!dir.exists(paste0(sta_way)))
  {
    dir.create(paste0(sta_way))
  }
  if(!dir.exists(paste0(sta_way,newname,"/")))
  {
    dir.create(paste0(sta_way,newname,"/"))
  }
  if(!dir.exists(paste0(outp_way,"/")))
  {
    dir.create(paste0(outp_way,"/"))
  }
  
  run_chronosta_script <- file.path(run_dir, "chronosta.py")

  file.copy(
    chronosta_script,
    run_chronosta_script,
    overwrite = TRUE
  )

  patch_chronosta_readonly_array(
    run_chronosta_script,
    logger = logger
  )
  
  #Copy the tree files in the directory
  for(i in 1:length(tree_list)){
    
    #In case of weighting backbone
    tre <- tree_list[[i]]  
    
    if(names(tree_list)[i] == ref_name)
    {
      for(j in 1:ref_weight)
      {
        write.tree(tre,
                   paste0(sta_way,newname,"/",ref_name,j,".nwk"))
      }
    } else { #Other files
      write.tree(tre,
                 paste0(sta_way,newname,"/",names(tree_list)[i],".nwk"))
    }
  }
  
  # Directory for bash, in case there is not the full way on the argument
  if(str_detect(sta_way, "^[.]"))
  {bsta_way <- paste0(getwd(), str_replace_all(sta_way, "^[.]*", ""))} else {
    bsta_way <- sta_way
  }
  
  bsta_way <- str_replace_all(bsta_way, "\\/", "\\\\")
  she_way <- paste0(bsta_way,newname)
  
  # Launch the python code
  system2(
    "cmd",
    args = c("/c", paste0(
      "cd /d \"", she_way, "\" && python -u chronosta.py"
    )),
    stdout = outp_mes,
    stderr = FALSE #can be switched to "" in case of error
  )
  
  #Rename the file
  file.rename(from = paste0(sta_way,newname,"/chronosta_supertree.newick"),
              to = paste0(sta_way,newname,"/",newname,".nwk"))
  
  #Copy it to the nwk folder 
  file.copy(paste0(sta_way,newname,"/",newname,".nwk"),
            paste0(outp_way,"/",newname,".nwk"),
            overwrite = T)
  
  supertree <- read.tree(paste0(outp_way,"/",newname,".nwk"))
  
  #Remove temp directory
  unlink(paste0(sta_way,"*"), recursive = T)
  
  return(supertree)
}

#####################################################################################################################################
#Split trees in genus level subtrees, except for genus with no overlapping specie sin the reference tree

split_gen_trees <- function(tre, ref, nam_tre){
  
  all_gen <- unique(str_split_fixed(tre$tip.label, "_", n = 2)[,1])
  dat_gen_ovl <- data.frame(tre = all_gen , ovl = FALSE, grp = NA)
  
  for(i in 1:nrow(dat_gen_ovl))
  {
    gen <- dat_gen_ovl$tre[i]
    sp_tre <- tre$tip.label[str_detect(tre$tip.label, paste0(gen,"_"))]
    sp_ref <- ref$tip.label[str_detect(ref$tip.label, paste0(gen,"_"))]
    
    if(length(intersect(sp_tre, sp_ref))==0)
    {
      dat_gen_ovl$ovl[i] <- FALSE
    } else {
      dat_gen_ovl$ovl[i] <- TRUE
    }
  }
  #Case with no overlap in any genus, return the full tree
  if(sum(dat_gen_ovl$ovl == FALSE) == nrow(dat_gen_ovl))
  {
    cat(paste0("No overlap with reference tree for ",nam_tre,". Returning the fulltree.\n"))
    list_tre <- list()
    list_tre[[1]] <- tre
    names(list_tre) <- nam_tre
    return(list_tre)
    
  } else {
    dat_gen_ovl$grp[dat_gen_ovl$ovl == TRUE] <- dat_gen_ovl$tre[dat_gen_ovl$ovl == TRUE]
    
    #Find the closest genus for overlap
    no_ovlp_gen <- dat_gen_ovl$tre[dat_gen_ovl$ovl == FALSE]
    
    for(i in no_ovlp_gen)
    {
      #List all species
      sp_tre  <- tre$tip.label[str_detect(tre$tip.label, paste0(i,"_"))]
      
      if(length(sp_tre) == 1)
      {
       mrca <-  which(tre$tip.label == sp_tre)
      } else {
      #Get the mrca
      mrca <- getMRCA(tre, sp_tre)
      }
      
      #Get the nodepath
      wayt <- ape::nodepath(tre, from = mrca, to = Ntip(tre)+1)
      
      is_ovl <- FALSE
      j <- 2
      while(is_ovl == FALSE)
      {
        #Get descendants
        desc <- extract.clade(tre, wayt[j])$tip.label
        
        #Get genus of descendants
        gen_des <- unique(str_split_fixed(desc, n =2, "_")[,1])
        
        if(sum(gen_des %in% dat_gen_ovl$tre[dat_gen_ovl$ovl == TRUE]) > 0)
        {
          dat_gen_ovl$grp[dat_gen_ovl$tre == i] <- gen_des[gen_des %in% dat_gen_ovl$tre[dat_gen_ovl$ovl == TRUE]][1]
          is_ovl <- TRUE
        } else {
          j <- j+1
        }
        
      }
      
    }
      
      #Split all the trees 
      all_grp <- unique(dat_gen_ovl$grp)
      list_tre <- list()
      
      for(i in 1:length(all_grp))
      {
       grp <- all_grp[i]
        
       all_sp <- tre$tip.label[str_detect(tre$tip.label, paste(paste0(dat_gen_ovl$tre[dat_gen_ovl$grp == grp],"_"), collapse = "|"))] 
       
       #Case with one species, means it is already in the ref, so we can drop it
       if(length(all_sp) == 1){
         next
       }
       
       #Extract clade
       sbt <- drop.tip(tre, tre$tip.label[!tre$tip.label %in% all_sp])
       
       list_tre[[paste0(nam_tre,"_",grp)]] <- sbt
      }
      
      return(list_tre)
  }
}
####################################### Split trees in subtrees of missing clades, with two outgroups in each #######################

split_tree_in_subtrees <- function(tree_graft,
                                   tree_ref,
                                   tre_name,
                                   n.iter.spl = 100,
                                   seed = 1
)
{
  trl <- list()
  trl[[tre_name]] <- tree_graft
  tree_graft <- trl
  
  sp_ref <- tree_ref$tip.label
  sp_gra <- lapply(tree_graft ,function(x){x$tip.label})
  sp_mis <- lapply(sp_gra, function(x){setdiff(x, sp_ref)})
  
  if(do.call(sum, lapply(sp_mis, function(x){ length(x)==0})) >0)
  {
    cat(paste0("No missing species in ", tre_name,"\n"))
    return()
  }
  
  prog <- progress::progress_bar$new(total = n.iter.spl,  format = " progress [:bar] :percent running for: :elapsed remaining: :eta")
  
  cnt_sim <- 0
  set.seed(seed)
  tab_sbt <- data.frame(tree = character(),focal_sp = character(), nb_mis_sp = numeric(), nb_ovl_sp = numeric(), iter =numeric())  
  
  for(i in 1:n.iter.spl){
    prog$tick()
    
    sp_ref <- tree_ref$tip.label
    sp_gra <- lapply(tree_graft ,function(x){x$tip.label})
    
    #Matching species
    sp_mat <- lapply(sp_gra, function(x){intersect(sp_ref, x)})
    #Missing species
    sp_mis <- lapply(sp_gra, function(x){setdiff(x, sp_ref)})
    
    #Table to store overlap information
    sbt_list <- list()
    
    #Loop for each tree in the list
    
    for(j in 1:length(tree_graft))
    {
      mis <- sp_mis[[j]]
      matc <- sp_mat[[j]]
      tre_g <- tree_graft[[j]]
      
      #To skip in case it has already been reduced to maximum
      if(i > 1 & length(mis) == 1){
        nam_sbt <- names(tree_graft)[j]
        sbt_list[[nam_sbt]] <-  tre_g
        
        next
      }
      
      #Shuffle the order of species
      mis <- sample(mis, length(mis), replace = F)
      
      nbsbt <- 0 #Counter for the subtree names
      
      while(length(mis)!= 0){
        
        nbsbt <- nbsbt+1
        #Target species
        trg_sp <- mis[1]  
        #Way until the root
        wayt <- ape::nodepath(tre_g, from = which(tre_g$tip.label == trg_sp), to = Ntip(tre_g)+1)
        
        cnt_otg <- 0 #Counting outgroups
        cnt_chk <- 2
        nd_chk <- wayt[cnt_chk] #first node to check
        
        while (cnt_otg < 2 & !is.na(nd_chk)) {
          
          clad <- extract.clade(tre_g, nd_chk)
          
          #Check the species in the other clade
          des_nod <- getDescendants(tre_g, nd_chk)[1:2][getDescendants(tre_g, nd_chk)[1:2] != wayt[cnt_chk-1]]
          
          if(length(getDescendants(tre_g, des_nod)) == 1){
            xclad <- tre_g$tip.label[getDescendants(tre_g, des_nod)]
          } else {
            xclad <- extract.clade(tre_g, des_nod)$tip.label
          } 
          
          #Check if there is at least one species that is not in the reference tree
          if(sum(xclad %in% mis) > 0){
            cnt_otg <- 0
          } else{
            cnt_otg <- cnt_otg + 1
          }
          
          if(nd_chk == wayt[length(wayt)]){
            nd_chk <- NA #Stop if we reach the root
          } else {
            cnt_chk <- cnt_chk+1 #Else, iterate to the next node
            nd_chk <- wayt[cnt_chk]
          }
        }
        
        #Work on the subtree to keep 1 species for clades already in tree ref
        nods_cla <- (Ntip(clad)+1):(Nnode(clad)+Ntip(clad))
        sp_to_rem <- c()
        
        for(k in nods_cla){
          cla_chk <- extract.clade(clad, k)$tip.label
          if(sum(cla_chk %in% mis) > 0) {next} else {
            sp_to_rem <- c(sp_to_rem, cla_chk[2:length(cla_chk)])
          }
        }
        
        
        clad <- drop.tip(clad, sp_to_rem)
        
        if(!identical(clad, tre_g))
        {
          #Get the iteration of subtree for the name
          num_m <- as.numeric(str_split_fixed(c(names(tree_graft), names(sbt_list)), "_s", n =2)[,2])
          if(length(num_m) == 1 & sum(is.na(num_m)) == length(num_m)){
            max_sbt <- 0 #First iteration case
          }   else{
          max_sbt <- max(num_m, na.rm = T)
          }

          #Name
          nam_sbt <- paste0(tre_name, "_s",max_sbt+1)
          
          tab_l <- data.frame(
            tree = nam_sbt,
            focal_sp = trg_sp,
            nb_mis_sp = length(setdiff(clad$tip.label, tree_ref$tip.label)),
            nb_ovl_sp = length(intersect(clad$tip.label, tree_ref$tip.label)),
            iter = i
          )
          
          tab_sbt <- rbind(tab_sbt, tab_l) #Add this to the tab
          sbt_list[[nam_sbt]] <-  clad
          
        } else {
          #Case for the first tree
          if(i == 1)
          {
            nam_sbt <- paste0(tre_name,"_s0")
            tab_l <- data.frame(
              tree = nam_sbt,
              focal_sp = trg_sp,
              nb_mis_sp = length(setdiff(clad$tip.label, tree_ref$tip.label)),
              nb_ovl_sp = length(intersect(clad$tip.label, tree_ref$tip.label)),
              iter = i
            )
            
            tab_sbt <- rbind(tab_sbt, tab_l) #Add this to the tab
            sbt_list[[nam_sbt]] <-  clad
          } else {
            
            #Get the iteration of subtree for the name
            nam_sbt <- names(tree_graft)[j]
            
            tab_l <- data.frame(
              tree = nam_sbt,
              focal_sp = trg_sp,
              nb_mis_sp = length(setdiff(clad$tip.label, tree_ref$tip.label)),
              nb_ovl_sp = length(intersect(clad$tip.label, tree_ref$tip.label)),
              iter = i
            )
            
            tab_sbt <- rbind(tab_sbt, tab_l) #Add this to the tab
            sbt_list[[nam_sbt]] <-  clad
          }
          
        }
        
        #Once we have the result, remove the species from tre_g and mis
        tre_g <- drop.tip(tre_g, clad$tip.label[clad$tip.label %in% mis])
        mis <- mis[!mis %in% clad$tip.label]
        
      }
      
      
    }
    
    # if(identical(tree_graft, sbt_list)){
    #   cnt_sim <- cnt_sim +1
    # } else {
    #   cnt_sim <- 0
    # }
    # 
    # if(cnt_sim > 9)
    # {
    #   warning(paste0("Optimum reached after ", i, " iterations" ))
    #   stop()
    # }
    tree_graft <- sbt_list
    
    
  }
  tab_sbt <- tab_sbt[tab_sbt$tree %in% names(tree_graft),]
  tab_sbt <- tab_sbt[!duplicated(tab_sbt[,c("tree", "nb_mis_sp", "nb_ovl_sp")]), ]
  sbt_list <- list(tab_sbt, tree_graft)
  names(sbt_list) <- c("data", "subtrees")
  return(sbt_list)
  
}

#####################################################################################################################################

# =============================================================================
# MacroPhyloMaker ChronoSTA integration helpers (minimal adaptation)
# =============================================================================
# These helpers are appended after the original utilities so they override only
# the pieces that need to be package-safe: ChronoSTA download/discovery, Python
# checking, and platform-independent execution.

.csta_require_packages <- function(pkgs = c("ape", "phytools", "phangorn", "stringr", "progress", "MonoPhy")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Required R package(s) not installed: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.csta_msg <- function(..., logger = NULL) {
  msg <- paste0(...)
  if (exists("log_msg", mode = "function", inherits = TRUE)) {
    get("log_msg", mode = "function", inherits = TRUE)(logger, msg)
  } else {
    cat(msg, "\n", sep = "")
    if (!is.null(logger) && !is.null(logger$file)) {
      cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", file = logger$file, append = TRUE, sep = "")
    }
  }
  invisible(msg)
}

.csta_section <- function(title, logger = NULL) {
  if (exists("log_section", mode = "function", inherits = TRUE)) {
    get("log_section", mode = "function", inherits = TRUE)(logger, title)
  } else {
    line <- paste(rep("-", nchar(title) + 4), collapse = "")
    .csta_msg(line, logger = logger)
    .csta_msg("| ", title, " |", logger = logger)
    .csta_msg(line, logger = logger)
  }
  invisible(title)
}

.csta_make_logger <- function(out_prefix = NULL, output_folder = NULL, file_prefix = "chronosta_grafting") {
  if (is.null(output_folder)) {
    if (!is.null(out_prefix)) {
      output_folder <- dirname(normalizePath(out_prefix, mustWork = FALSE))
    } else {
      output_folder <- getwd()
    }
  }
  log_dir <- file.path(output_folder, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(log_dir, paste0(file_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  cat("ChronoSTA grafting log\n", file = logfile)
  cat("Started: ", format(Sys.time()), "\n\n", file = logfile, append = TRUE, sep = "")
  structure(list(file = logfile), class = "csta_logger")
}

chronosta_default_dir <- function() {
  if (requireNamespace("rappdirs", quietly = TRUE)) {
    rappdirs::user_data_dir("MacroPhyloMaker", "chronosta")
  } else {
    file.path(path.expand("~"), ".macro_phylo_maker", "chronosta")
  }
}

#' Download the ChronoSTA Python script
#'
#' Downloads chronosta.py from the public Chrono-STA GitHub repository into a
#' user-writable cache directory or a user-specified directory. This avoids
#' requiring users to manually place chronosta.py in the working directory.
#'
#' @param dest_dir Directory where chronosta.py will be stored.
#' @param url Raw GitHub URL for chronosta.py.
#' @param overwrite Logical; overwrite existing file?
#' @return Normalized path to chronosta.py.
download_chronosta <- function(
  dest_dir = chronosta_default_dir(),
  url = "https://raw.githubusercontent.com/josebarbamontoya/chrono-sta/main/code/chronosta.py",
  overwrite = FALSE
) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dest_dir, "chronosta.py")
  if (file.exists(dest) && !isTRUE(overwrite)) return(normalizePath(dest))
  utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    stop("Failed to download chronosta.py to ", dest, call. = FALSE)
  }
  normalizePath(dest)
}

find_chronosta_script <- function(chronosta_script = NULL, csta_folder = NULL, download = TRUE) {
  if (!is.null(chronosta_script) && !is.na(chronosta_script) && nzchar(chronosta_script)) {
    if (!file.exists(chronosta_script)) stop("ChronoSTA script not found: ", chronosta_script, call. = FALSE)
    return(normalizePath(chronosta_script))
  }

  if (!is.null(csta_folder) && !is.na(csta_folder) && nzchar(csta_folder)) {
    cand <- file.path(csta_folder, "chronosta.py")
    if (file.exists(cand)) return(normalizePath(cand))
    if (isTRUE(download)) return(download_chronosta(dest_dir = csta_folder))
  }

  bundled <- system.file("python", "chronosta.py", package = "MacroPhyloMaker")
  if (nzchar(bundled) && file.exists(bundled)) return(normalizePath(bundled))

  if (isTRUE(download)) return(download_chronosta())

  stop("Could not find chronosta.py. Provide `chronosta_script`, `csta_folder`, or set `download = TRUE`.", call. = FALSE)
}

find_python <- function(python = NULL) {
  if (!is.null(python) && nzchar(python)) {
    p <- Sys.which(python)
    if (nzchar(p)) return(unname(p))
    if (file.exists(python)) return(normalizePath(python))
    stop("Python executable not found: ", python, call. = FALSE)
  }
  p <- Sys.which("python3")
  if (!nzchar(p)) p <- Sys.which("python")
  if (!nzchar(p)) stop("No Python executable found. Install Python 3 or provide `python`.", call. = FALSE)
  unname(p)
}

#' Check the Python environment for ChronoSTA
#'
#' Verifies that a selected Python executable can locate the Python packages
#' required by the ChronoSTA-enabled grafting workflow. By default, the function
#' checks for Biopython, pandas, NumPy, SciPy, and Matplotlib.
#'
#' @param python Character. Optional path or command name for a Python
#'   executable. If `NULL`, `find_python()` searches for `python3` and then
#'   `python` on the system path.
#' @param packages Character vector of Python import names to check. The default
#'   values correspond to Biopython (`Bio`), pandas, NumPy, SciPy, and
#'   Matplotlib.
#'
#' @return Invisibly returns the resolved path to the Python executable when all
#'   requested packages are available. Stops with an informative error if the
#'   executable cannot be found or one or more packages are unavailable.
#'
#' @details
#' This function only checks an existing Python environment. It does not create
#' a virtual environment, install packages, activate an environment, or change
#' the Python executable used by later calls. Pass the returned executable, or
#' the same value supplied to `python`, to [run_chronosta_grafting()].
#'
#' Each package is checked in a separate temporary Python script using
#' `importlib.util.find_spec()`.
#'
#' @examples
#' \dontrun{
#' # Automatically locate Python.
#' check_chronosta_python()
#'
#' # Check a project-specific virtual environment on Unix-like systems.
#' check_chronosta_python(".venv-chronosta/bin/python")
#'
#' # Native Windows example.
#' check_chronosta_python(".venv-chronosta/Scripts/python.exe")
#' }
#'
#' @export
check_chronosta_python <- function(
  python = NULL,
  packages = c("Bio", "pandas", "numpy", "scipy", "matplotlib")
) {
  py <- find_python(python)

  missing <- character(0)

  for (pkg in packages) {
    tf <- tempfile(fileext = ".py")

    writeLines(
      c(
        "import importlib.util",
        "import sys",
        sprintf("pkg = %s", encodeString(pkg, quote = "\"")),
        "sys.exit(0 if importlib.util.find_spec(pkg) else 1)"
      ),
      con = tf
    )

    status <- system2(
      py,
      args = tf,
      stdout = TRUE,
      stderr = TRUE
    )

    exit_status <- attr(status, "status")
    if (is.null(exit_status)) {
      exit_status <- 0L
    }

    unlink(tf)

    if (!identical(as.integer(exit_status), 0L)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing)) {
    stop(
      "Python executable found, but missing ChronoSTA dependency/dependencies: ",
      paste(missing, collapse = ", "),
      ". Install them with e.g. `",
      py,
      " -m pip install biopython pandas numpy scipy matplotlib`.",
      call. = FALSE
    )
  }

  invisible(py)
}

#' Set up the Python environment and script for ChronoSTA
#'
#' Prepares the external Python components required by the ChronoSTA-enabled
#' grafting workflow. The function can install the required Python packages,
#' download `chronosta.py`, and verify that the selected Python executable can
#' locate all required dependencies.
#'
#' @param python Character. Optional path or command name for a Python
#'   executable. If `NULL`, `find_python()` searches for `python3` and then
#'   `python` on the system path.
#' @param install Logical. If `TRUE`, install Biopython, pandas, NumPy, SciPy,
#'   and Matplotlib into the environment associated with `python` using
#'   `python -m pip install`.
#' @param download Logical. If `TRUE`, download `chronosta.py` into `dest_dir`.
#'   If `FALSE`, use `find_chronosta_script()` to locate an existing script
#'   without downloading it.
#' @param dest_dir Character. Directory in which `chronosta.py` is stored when
#'   `download = TRUE`. Defaults to the user-writable ChronoSTA directory
#'   returned by `chronosta_default_dir()`.
#'
#' @return Invisibly returns a list with two elements: `python`, the resolved
#'   Python executable path, and `chronosta_script`, the resolved path to
#'   `chronosta.py`.
#'
#' @details
#' The function operates on the Python environment associated with the selected
#' executable. It does not activate that environment or permanently configure R.
#' Use the same Python executable for setup, verification with
#' [check_chronosta_python()], and [run_chronosta_grafting()].
#'
#' Set `install = FALSE` when the Python dependencies are already installed. Set
#' `download = FALSE` only when an existing `chronosta.py` can be located by
#' `find_chronosta_script()`.
#'
#' @examples
#' \dontrun{
#' # Automatically locate Python, install dependencies, and download ChronoSTA.
#' env <- setup_chronosta_env()
#'
#' # Use a project-specific virtual environment on Unix-like systems.
#' env <- setup_chronosta_env(
#'   python = ".venv-chronosta/bin/python"
#' )
#'
#' # Native Windows example.
#' env <- setup_chronosta_env(
#'   python = ".venv-chronosta/Scripts/python.exe"
#' )
#'
#' # Verify an environment without reinstalling its Python packages.
#' env <- setup_chronosta_env(
#'   python = ".venv-chronosta/bin/python",
#'   install = FALSE
#' )
#' }
#'
#' @export
setup_chronosta_env <- function(
  python = NULL,
  install = TRUE,
  download = TRUE,
  dest_dir = chronosta_default_dir()
) {
  py <- find_python(python)

  message("Using Python: ", py)

  if (isTRUE(install)) {
    status <- system2(
      py,
      c(
        "-m", "pip", "install",
        "biopython", "pandas", "numpy", "scipy", "matplotlib"
      )
    )

    if (!identical(status, 0L)) {
      stop(
        "Python dependency installation failed for: ", py,
        "\nTry running manually:\n",
        py, " -m pip install biopython pandas numpy scipy matplotlib",
        call. = FALSE
      )
    }
  }

  script <- if (isTRUE(download)) {
    download_chronosta(dest_dir = dest_dir)
  } else {
    find_chronosta_script(download = FALSE)
  }

  check_chronosta_python(py)

  message("ChronoSTA setup complete.")
  message("Python: ", py)
  message("ChronoSTA script: ", script)

  invisible(list(
    python = py,
    chronosta_script = script
  ))
}

# Cross-platform replacement for the original run_chrono_sta(). This preserves
# the original argument names where possible, but adds python/chronosta_script and
# writes stdout/stderr logs for each ChronoSTA call.
run_chrono_sta <- function(tree_list,
                           ref_name = "NONE",
                           newname,
                           ref_weight = 1,
                           sta_way = tempdir(),
                           csta_folder = NULL,
                           outp_way = NULL,
                           outp_mes = TRUE,
                           chronosta_script = NULL,
                           python = NULL,
                           logger = NULL,
                           keep_temp = FALSE,
                           download = TRUE) {
  .csta_require_packages()
  if (missing(newname) || is.null(newname) || !nzchar(newname)) stop("`newname` is required.", call. = FALSE)
  if (is.null(outp_way) || is.na(outp_way) || !nzchar(outp_way)) outp_way <- sta_way
  dir.create(sta_way, recursive = TRUE, showWarnings = FALSE)
  dir.create(outp_way, recursive = TRUE, showWarnings = FALSE)

  chronosta_script <- find_chronosta_script(chronosta_script = chronosta_script, csta_folder = csta_folder, download = download)
  py <- find_python(python)

  .csta_msg("ChronoSTA source script: ", chronosta_script, logger = logger)
  .csta_msg("ChronoSTA Python executable: ", py, logger = logger)

  # Fail early with a clear message if Python dependencies are unavailable.
  check_chronosta_python(py)

  run_dir <- file.path(sta_way, "temp", newname)
  if (dir.exists(run_dir)) unlink(run_dir, recursive = TRUE, force = TRUE)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_chronosta_script <- file.path(run_dir, "chronosta.py")

  file.copy(
    chronosta_script,
    run_chronosta_script,
    overwrite = TRUE
  )

  .csta_msg("ChronoSTA run script copied to: ", run_chronosta_script, logger = logger)

  patch_chronosta_readonly_array(
    run_chronosta_script,
    logger = logger
  )

  nms <- names(tree_list)
  if (is.null(nms)) nms <- paste0("tree_", seq_along(tree_list))
  nms[!nzchar(nms)] <- paste0("tree_", which(!nzchar(nms)))
  names(tree_list) <- nms

  for (i in seq_along(tree_list)) {
    tre <- tree_list[[i]]
    nm <- names(tree_list)[i]
    safe_nm <- gsub("[^A-Za-z0-9_.-]", "_", nm)
    if (identical(nm, ref_name)) {
      for (j in seq_len(max(1L, as.integer(ref_weight)))) {
        ape::write.tree(tre, file.path(run_dir, paste0(safe_nm, j, ".nwk")))
      }
    } else {
      ape::write.tree(tre, file.path(run_dir, paste0(safe_nm, ".nwk")))
    }
  }

  stdout_file <- file.path(run_dir, "chronosta_stdout.log")
  stderr_file <- file.path(run_dir, "chronosta_stderr.log")
  .csta_msg("Running ChronoSTA: ", newname, " in ", run_dir, logger = logger)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run_dir)

  status <- system2(
    command = py,
    args = c("-u", "chronosta.py"),
    stdout = if (isTRUE(outp_mes)) "" else stdout_file,
    stderr = stderr_file,
    wait = TRUE
  )

  if (!identical(status, 0L)) {
    stop("ChronoSTA failed for ", newname, " with exit status ", status,
         ". See stderr log: ", stderr_file, call. = FALSE)
  }

  raw_out <- file.path(run_dir, "chronosta_supertree.newick")
  if (!file.exists(raw_out)) {
    stop("ChronoSTA completed for ", newname, " but did not create chronosta_supertree.newick in ", run_dir, call. = FALSE)
  }

  final_out <- file.path(run_dir, paste0(newname, ".nwk"))
  file.rename(raw_out, final_out)
  out_copy <- file.path(outp_way, paste0(newname, ".nwk"))
  file.copy(final_out, out_copy, overwrite = TRUE)
  st <- ape::read.tree(out_copy)

  if (!isTRUE(keep_temp)) unlink(file.path(sta_way, "temp"), recursive = TRUE, force = TRUE)
  .csta_msg("ChronoSTA output written: ", out_copy, logger = logger)
  st
}
