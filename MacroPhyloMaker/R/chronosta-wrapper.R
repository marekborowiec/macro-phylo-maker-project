# =============================================================================
# ChronoSTA-enabled species grafting wrapper for MacroPhyloMaker
# Minimal package adaptation of sta_sp_grafting()
# =============================================================================
# For package use, place this file under R/ along with chronosta_utils_adapted.R.
# For standalone use, source chronosta_utils_adapted.R before this file.

#' Run time-informed species grafting using ChronoSTA
#'
#' This is a non-interactive, package-oriented wrapper around the original
#' ChronoSTA species-grafting workflow. It assumes the reference tree has already
#' been produced by the MacroPhyloMaker tip/clade grafting steps, identifies
#' species present in donor trees but absent from the reference, optionally
#' recalibrates donor trees to the reference timescale, splits trees to reduce
#' topology disturbance, runs ChronoSTA, and optionally restores genus-level
#' monophyly introduced by the final merge.
#'
#' @param tree_ref Legacy alias for `reference_tree`. Use `reference_tree` in
#'   new code.
#' @param tree_folder Legacy alias for `donor_tree_dir`. Use `donor_tree_dir`
#'   in new code.
#' @param recalibrate Logical. If `TRUE`, recalibrate eligible donor trees to
#'   the timescale of the reference tree before Chrono-STA merging.
#' @param split_gen Logical. If `TRUE`, split donor trees into genus-level
#'   subtrees before identifying missing-species subtrees.
#' @param split_sbt Logical. If `TRUE`, extract minimal donor subtrees that
#'   contain missing species and sufficient overlap with the reference tree.
#' @param monoph_restore Logical. If `TRUE`, apply genus-level monophyly
#'   restoration after Chrono-STA grafting.
#' @param plot_monoph Logical. If `TRUE`, produce diagnostic monophyly plots
#'   when monophyly restoration is requested.
#' @param prefuse Logical. If `TRUE`, merge compatible or overlapping donor
#'   subtrees with Chrono-STA before merging them with the reference tree.
#' @param ultr_tol Numeric tolerance used when testing whether donor trees are
#'   ultrametric.
#' @param calib_ratio Numeric proportion used when selecting the number of
#'   shared calibration nodes for donor-tree recalibration.
#' @param calib_min_nodes Integer minimum number of shared calibration nodes
#'   targeted during donor-tree recalibration.
#' @param calib_max_nodes Integer maximum number of shared calibration nodes
#'   retained during donor-tree recalibration.
#' @param margin_nod Numeric minimum temporal separation required between
#'   nested calibration constraints.
#' @param tolerance_pct Numeric proportional tolerance applied around
#'   reference-tree calibration ages.
#' @param tolerance_min Numeric minimum absolute tolerance applied to
#'   calibration ages.
#' @param lambda Numeric smoothing parameter passed to [ape::chronos()] during
#'   donor-tree recalibration.
#' @param split_seed Integer random seed used during donor-subtree splitting
#'   and other stochastic preprocessing.
#' @param ref_weight Integer weighting assigned to the reference tree during
#'   Chrono-STA merging.
#' @param paraph_exc Optional character vector of genera that should not be
#'   protected as baseline non-monophyletic genera during post-grafting
#'   monophyly restoration.
#' @param download_chronosta Logical. If `TRUE`, download `chronosta.py` when
#'   no local script can be found.
#' @param keep_temp Logical. If `TRUE`, retain temporary Chrono-STA run
#'   directories and intermediate files.
#' @param reference_tree A phylo object or path to a Newick tree. Alias: tree_ref.
#' @param donor_tree_dir Folder containing donor .nwk/.tre trees. Alias: tree_folder.
#' @param out_prefix Output file prefix. If omitted, output_folder is used.
#' @param output_folder Output folder for legacy compatibility.
#' @param chronosta_script Optional path to chronosta.py.
#' @param csta_folder Optional folder containing chronosta.py. If missing and
#'   download_chronosta = TRUE, chronosta.py is downloaded there.
#' @param python Optional Python executable.
#' @param include_trees Optional character vector of donor tree names to include.
#'   Names are file basenames without extension. NULL means all eligible trees.
#' @param logger Optional MacroPhyloMaker logger. If NULL, a simple logfile is made.
#' @return A list with final tree, logs, and output paths. The final usable tree is
#'   always available as result$tree.
#' @export
run_chronosta_grafting <- function(
  reference_tree = NULL,
  donor_tree_dir = NULL,
  out_prefix = NULL,
  output_folder = NULL,
  tree_ref = NULL,
  tree_folder = NULL,
  csta_folder = NULL,
  chronosta_script = NULL,
  python = NULL,
  include_trees = NULL,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  monoph_restore = TRUE,
  plot_monoph = TRUE,
  prefuse = TRUE,
  ultr_tol = 0.05,
  calib_ratio = 0.3,
  calib_min_nodes = 3,
  calib_max_nodes = 70,
  margin_nod = 0.5,
  tolerance_pct = 0.05,
  tolerance_min = 0.2,
  lambda = 1,
  split_seed = 1,
  ref_weight = 1,
  paraph_exc = NULL,
  logger = NULL,
  download_chronosta = TRUE,
  keep_temp = FALSE
) {
  .csta_require_packages()

  if (is.null(reference_tree)) reference_tree <- tree_ref
  if (is.null(donor_tree_dir)) donor_tree_dir <- tree_folder

  if (is.null(reference_tree) || length(reference_tree) == 0 || all(is.na(reference_tree))) {
    stop("Please provide `reference_tree`/`tree_ref` as a phylo object or Newick path.", call. = FALSE)
  }
  if (is.null(donor_tree_dir) || length(donor_tree_dir) != 1 || is.na(donor_tree_dir) || !nzchar(donor_tree_dir)) {
    stop("Please provide `donor_tree_dir`/`tree_folder`.", call. = FALSE)
  }

  if (inherits(reference_tree, "phylo")) {
    tree_ref <- reference_tree
  } else if (is.character(reference_tree) && length(reference_tree) == 1 && file.exists(reference_tree)) {
    tree_ref <- ape::read.tree(reference_tree)
  } else {
    stop("`reference_tree` must be a phylo object or an existing Newick path.", call. = FALSE)
  }

  if (is.null(out_prefix)) {
    if (is.null(output_folder) || is.na(output_folder) || !nzchar(output_folder)) {
      output_folder <- getwd()
    }
    out_prefix <- file.path(output_folder, "chronosta_gapfilled")
  }
  output_folder <- dirname(normalizePath(out_prefix, mustWork = FALSE))
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  if (is.null(logger)) logger <- .csta_make_logger(out_prefix = out_prefix, file_prefix = "chronosta_grafting")

  chronosta_script <- find_chronosta_script(
    chronosta_script = chronosta_script,
    csta_folder = csta_folder,
    download = download_chronosta
  )

  .csta_section("ChronoSTA grafting: input", logger = logger)
  .csta_msg("Reference tips: ", ape::Ntip(tree_ref), logger = logger)
  .csta_msg("Donor tree folder: ", donor_tree_dir, logger = logger)
  .csta_msg("ChronoSTA script: ", chronosta_script, logger = logger)

  tree_files <- list.files(donor_tree_dir, pattern = "\\.(nwk|tre)$", full.names = TRUE, ignore.case = TRUE)
  tree_names <- tools::file_path_sans_ext(basename(tree_files))
  if (!length(tree_files)) stop("No .nwk/.tre donor trees found in `donor_tree_dir`: ", donor_tree_dir, call. = FALSE)

  list_tre <- lapply(tree_files, ape::read.tree)
  names(list_tre) <- tree_names

  if (!is.null(include_trees)) {
    missing_requested <- setdiff(include_trees, names(list_tre))
    if (length(missing_requested)) {
      stop("Requested donor tree(s) not found: ", paste(missing_requested, collapse = ", "), call. = FALSE)
    }
    list_tre <- list_tre[include_trees]
  }

  # Check separators. Allow genus-only labels on the reference, then add _sp below.
  sep_ref <- sum(stringr::str_detect(tree_ref$tip.label, "_")) == 0
  sep_trs <- unlist(lapply(list_tre, function(x) sum(stringr::str_detect(x$tip.label, "_")) == 0))
  if (sum(sep_trs) > 0) {
    stop("Donor tree tip labels should use `Genus_species` format with `_` separators.", call. = FALSE)
  }
  if (sep_ref) {
    .csta_msg("Reference labels appear genus-only; appending `_sp` to genus-only reference tips.", logger = logger)
  }
  tree_ref$tip.label[!stringr::str_detect(tree_ref$tip.label, "_")] <- paste0(tree_ref$tip.label[!stringr::str_detect(tree_ref$tip.label, "_")], "_sp")

  .csta_section("ChronoSTA grafting: overlap", logger = logger)
  data_overlap <- data.frame(tre = names(list_tre), total = NA_integer_, overlap = NA_integer_, missing = NA_integer_)
  for (i in names(list_tre)) {
    tre <- list_tre[[i]]
    data_overlap$total[data_overlap$tre == i] <- ape::Ntip(tre)
    data_overlap$overlap[data_overlap$tre == i] <- length(intersect(tre$tip.label, tree_ref$tip.label))
    data_overlap$missing[data_overlap$tre == i] <- length(setdiff(tre$tip.label, tree_ref$tip.label))
  }
  utils::write.table(data_overlap, paste0(out_prefix, "_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  .csta_msg("Donor trees read: ", length(list_tre), logger = logger)

  data_overlap <- data_overlap[data_overlap$missing != 0, , drop = FALSE]
  list_tre <- list_tre[names(list_tre) %in% data_overlap$tre]
  if (nrow(data_overlap) == 0 || !length(list_tre)) {
    .csta_msg("All donor species are already in the reference tree. No ChronoSTA grafting performed.", logger = logger)
    ape::write.tree(tree_ref, paste0(out_prefix, "_final_tree.nwk"))
    return(list(tree = tree_ref, raw_tree = tree_ref, overlap = data_overlap, paths = list(final_tree = paste0(out_prefix, "_final_tree.nwk"), log = logger$file)))
  }
  .csta_msg("Trees with missing species retained: ", length(list_tre), logger = logger)

  tree_to_add <- list_tre

  if (isTRUE(recalibrate)) {
    .csta_section("ChronoSTA grafting: recalibration", logger = logger)
    cm_nd_tab <- lapply(tree_to_add, function(x) {
      get_shared_nodes(
        tree_ref = tree_ref,
        tree_target = x,
        calib_ratio = calib_ratio,
        calib_min_nodes = calib_min_nodes,
        calib_max_nodes = calib_max_nodes
      )
    })

    recal_trees <- vector("list", length(tree_to_add))
    names(recal_trees) <- names(tree_to_add)
    for (nm in names(tree_to_add)) {
      .csta_msg("Starting recalibration for ", nm, " (", match(nm, names(tree_to_add)), "/", length(tree_to_add), ").", logger = logger)
      recal_trees[[nm]] <- recalibrate_with_chronos(
        tree_ref = tree_ref,
        tree_target = tree_to_add[[nm]],
        shared_nodes = cm_nd_tab[[nm]],
        tolerance_pct = tolerance_pct,
        tolerance_min = tolerance_min,
        calib_min_nodes = calib_min_nodes,
        lambda = lambda,
        margin_nod = margin_nod
      )
    }
  } else {
    .csta_msg("No recalibration performed (recalibrate = FALSE).", logger = logger)
    recal_trees <- tree_to_add
  }

  test_chrono <- vapply(recal_trees, function(x) !ape::is.ultrametric(x, tol = ultr_tol), logical(1))
  if (sum(test_chrono) > 0) {
    .csta_msg("Removing ", sum(test_chrono), " trees due to failing ultrametricity test.", logger = logger)
    recal_trees <- recal_trees[!test_chrono]
  }
  if (!length(recal_trees)) stop("No recalibrated/chronogram donor trees remain for ChronoSTA grafting.", call. = FALSE)

  .csta_section("ChronoSTA grafting: split donor trees", logger = logger)
  if (isTRUE(split_gen)) {
    gen_sub_trees_nested <- lapply(names(recal_trees), function(nm) {
      .csta_msg("Genus subtree splitting for ", nm, ".", logger = logger)
      split_gen_trees(tre = recal_trees[[nm]], ref = tree_ref, nam_tre = nm)
    })
    gen_sub_trees <- unlist(gen_sub_trees_nested, recursive = FALSE, use.names = TRUE)
  } else {
    .csta_msg("No genus splitting performed (split_gen = FALSE).", logger = logger)
    gen_sub_trees <- recal_trees
  }

  if (!length(gen_sub_trees)) stop("No genus subtrees remain after splitting.", call. = FALSE)
  test_dups <- duplicated(lapply(gen_sub_trees, function(x) x$tip.label))
  if (any(test_dups)) gen_sub_trees <- gen_sub_trees[!test_dups]

  missing_spe <- unique(unlist(lapply(gen_sub_trees, function(x) x$tip.label)))
  missing_spe <- missing_spe[!missing_spe %in% tree_ref$tip.label]
  test_missing <- vapply(gen_sub_trees, function(x) sum(x$tip.label %in% missing_spe) > 0, logical(1))
  gen_sub_trees <- gen_sub_trees[test_missing]
  if (!length(gen_sub_trees)) stop("No genus subtrees contain missing species.", call. = FALSE)

  if (isTRUE(split_sbt)) {
    .csta_msg("Starting extraction of minimal missing-species subtrees.", logger = logger)
    sub_trees_list <- lapply(names(gen_sub_trees), function(nm) {
      split_tree_in_subtrees(tree_graft = gen_sub_trees[[nm]], tree_ref = tree_ref, tre_name = nm, seed = split_seed)
    })
    sub_trees_list <- Filter(Negate(is.null), sub_trees_list)
    splits_data <- if (length(sub_trees_list)) do.call(rbind, lapply(sub_trees_list, function(x) x$data)) else data.frame()
    sub_trees <- if (length(sub_trees_list)) unlist(lapply(sub_trees_list, function(x) x$subtrees), recursive = FALSE, use.names = TRUE) else list()
    .csta_msg("Subtree splits completed: ", length(sub_trees), " subtree(s).", logger = logger)
  } else {
    .csta_msg("No missing-species subtree splitting performed (split_sbt = FALSE).", logger = logger)
    sub_trees <- gen_sub_trees
    splits_data <- data.frame()
  }
  utils::write.table(splits_data, paste0(out_prefix, "_splits.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!length(sub_trees)) stop("No subtrees remain for ChronoSTA merging.", call. = FALSE)

  .csta_section("ChronoSTA grafting: prefusion", logger = logger)
  if (isTRUE(prefuse) && length(sub_trees) > 1) {
    ovl_tab <- data.frame(t(utils::combn(names(sub_trees), 2)), stringsAsFactors = FALSE)
    colnames(ovl_tab) <- c("tA", "tB")
    ovl_tab$nb_ovl <- apply(ovl_tab, 1, function(x) {
      tA <- sub_trees[[unlist(x["tA"])]][["tip.label"]]
      tB <- sub_trees[[unlist(x["tB"])]][["tip.label"]]
      length(intersect(tA, tB))
    })
    utils::write.table(ovl_tab, paste0(out_prefix, "_prefusion_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

    if (sum(ovl_tab$nb_ovl) == 0) {
      .csta_msg("0 species overlap among subtrees; prefusion skipped.", logger = logger)
      tree_list_fusion <- sub_trees
    } else {
      merg_df <- data.frame(tree = names(sub_trees), stringsAsFactors = FALSE)
      df_long <- utils::stack(lapply(sub_trees, function(x) x$tip.label))
      colnames(df_long) <- c("species", "tree")

      g_biparti <- igraph::graph_from_data_frame(df_long, directed = FALSE)
      g_elements <- igraph::bipartite_projection(
        g_biparti,
        types = igraph::V(g_biparti)$name %in% df_long$species
      )$proj1
      compo <- igraph::components(g_elements)
      df_resultat_groupes <- data.frame(
        tree = names(compo$membership),
        group = paste0("Groupe_", compo$membership),
        stringsAsFactors = FALSE
      )
      df_groups <- merge(merg_df, df_resultat_groupes, by = "tree", all.x = TRUE)
      df_groups$group[is.na(df_groups$group)] <- paste0("Groupe_al_", seq_len(sum(is.na(df_groups$group))))
      utils::write.table(df_groups, paste0(out_prefix, "_prefusion_groups.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

      .csta_msg("Prefusion groups found: ", length(unique(df_groups$group)), logger = logger)
      prog_pref <- progress::progress_bar$new(total = length(unique(df_groups$group)), format = " ChronoSTA prefusion [:bar] :percent elapsed: :elapsed eta: :eta")
      fused_trees <- list()
      for (m in unique(df_groups$group)) {
        prog_pref$tick()
        if (stringr::str_detect(m, "_al_") || sum(df_groups$group == m) == 1) {
          fused_trees[[m]] <- sub_trees[[df_groups$tree[df_groups$group == m][1]]]
        } else {
          tre_in_gr <- sub_trees[names(sub_trees) %in% df_groups$tree[df_groups$group == m]]
          tre_merg <- run_chrono_sta(
            tree_list = tre_in_gr,
            ref_name = "NONE",
            newname = m,
            ref_weight = 1,
            sta_way = output_folder,
            outp_way = file.path(output_folder, "prefusions"),
            csta_folder = csta_folder,
            chronosta_script = chronosta_script,
            python = python,
            outp_mes = FALSE,
            logger = logger,
            keep_temp = keep_temp,
            download = download_chronosta
          )
          fused_trees[[m]] <- tre_merg
        }
      }
      tree_list_fusion <- fused_trees
      .csta_msg("Prefusion complete. Trees entering final merge: ", length(tree_list_fusion), logger = logger)
    }
  } else {
    .csta_msg("No prefusion performed.", logger = logger)
    tree_list_fusion <- sub_trees
  }

  .csta_section("ChronoSTA grafting: final merge", logger = logger)

  # One-species genus placeholders in the reference (Genus_sp) are renamed to a
  # real species from the same genus when available in the ChronoSTA subtrees.
  gen_list <- unique(stringr::str_split_fixed(tree_ref$tip.label, "_", n = 2)[, 1])
  sp_per_gen <- data.frame(
    gen = gen_list,
    nbsp = vapply(gen_list, function(x) sum(stringr::str_detect(tree_ref$tip.label, paste0("^", x, "_"))), integer(1)),
    stringsAsFactors = FALSE
  )
  sp_per_gen <- sp_per_gen[sp_per_gen$nbsp == 1, , drop = FALSE]
  if (nrow(sp_per_gen) > 0) {
    mnt_und_sp <- tree_ref$tip.label[
      stringr::str_detect(tree_ref$tip.label, paste(sp_per_gen$gen, collapse = "|")) &
        stringr::str_detect(tree_ref$tip.label, "_sp$")
    ]
    for (u in mnt_und_sp) {
      gen <- stringr::str_split_fixed(u, "_", n = 2)[, 1]
      list_sp_gen <- unlist(lapply(tree_list_fusion, function(x) x$tip.label[stringr::str_detect(x$tip.label, paste0("^", gen, "_"))]))
      list_sp_gen <- list_sp_gen[!stringr::str_detect(list_sp_gen, "_sp$")]
      if (length(list_sp_gen) > 0) tree_ref$tip.label[tree_ref$tip.label == u] <- list_sp_gen[1]
    }
  }

  all_taxa <- unique(c(unlist(lapply(tree_list_fusion, function(x) x$tip.label)), tree_ref$tip.label))
  extr_valid <- all_taxa[!stringr::str_detect(all_taxa, "_sp$")]
  gen_val <- unique(stringr::str_split_fixed(extr_valid, "_", n = 2)[, 1])
  gen_tree <- stringr::str_split_fixed(all_taxa, "_", n = 2)[, 1]
  gen_novalid_sp <- !gen_tree %in% gen_val
  valid_taxa <- c(extr_valid, all_taxa[gen_novalid_sp])
  tree_list_fusion <- lapply(tree_list_fusion, function(x) ape::drop.tip(x, x$tip.label[!x$tip.label %in% valid_taxa]))
  tree_list_fusion <- tree_list_fusion[vapply(tree_list_fusion, function(x) length(x$tip.label) > 1, logical(1))]

  ovl_ref <- vapply(tree_list_fusion, function(x) sum(x$tip.label %in% tree_ref$tip.label), integer(1))
  tree_list_fusion <- tree_list_fusion[ovl_ref > 0]
  if (!length(tree_list_fusion)) stop("No subtrees overlap with the reference tree.", call. = FALSE)

  prog <- progress::progress_bar$new(total = length(tree_list_fusion), format = " ChronoSTA final merge [:bar] :percent elapsed: :elapsed eta: :eta")
  tree_bkb <- tree_ref

  for (i in seq_along(tree_list_fusion)) {
    prog$tick()
    tre <- tree_list_fusion[i]
    ref <- tree_bkb
    shared <- tre[[1]]$tip.label[tre[[1]]$tip.label %in% ref$tip.label]
    if (length(shared) == 1) {
      nd <- which(ref$tip.label == shared)
    } else {
      nd <- ape::getMRCA(ref, shared)
    }
    tax_ref <- ape::nodepath(tree_bkb, to = ape::Ntip(tree_bkb) + 1, from = nd)
    nd_ref <- if (length(tax_ref) > 1) tax_ref[2] else tax_ref[1]
    
    gen_ref <- unique(str_split_fixed(extract.clade(ref, nd_ref)$tip.label, "_", n = 2)[,1])
    gen_bra <- unique(str_split_fixed(tre[[1]]$tip.label, "_", n = 2)[,1])
    
    if(sum(!gen_bra %in% gen_ref) > 0){
      nd_ref <- getMRCA(ref, ref$tip.label[str_detect(ref$tip.label, paste(paste0("^", gen_bra), collapse = "|"))])
    }
    
    ref <- ape::extract.clade(tree_bkb, nd_ref)

    nam_tr <- paste0("fsub_", names(tree_list_fusion)[i])
    tre[["ref"]] <- ref

    sta_tree <- run_chrono_sta(
      tree_list = tre,
      ref_name = "ref",
      newname = nam_tr,
      ref_weight = ref_weight,
      sta_way = output_folder,
      outp_way = file.path(output_folder, "prefusions"),
      csta_folder = csta_folder,
      chronosta_script = chronosta_script,
      python = python,
      outp_mes = FALSE,
      logger = logger,
      keep_temp = keep_temp,
      download = download_chronosta
    )

    tree_length <- phytools::nodeHeights(ref)
    tree_length <- tree_length[nrow(tree_length), 2]
    sta_length <- phytools::nodeHeights(sta_tree)
    sta_length <- sta_length[nrow(sta_length), 2]
    if (is.finite(tree_length) && is.finite(sta_length) && sta_length > 0) {
      sta_tree$edge.length <- sta_tree$edge.length * (tree_length / sta_length)
    }

    tree_bkb$tip.label[tree_bkb$tip.label %in% sta_tree$tip.label] <- paste0(tree_bkb$tip.label[tree_bkb$tip.label %in% sta_tree$tip.label], "_DELETE")
    tree_bkb <- ape::bind.tree(tree_bkb, sta_tree, where = nd_ref)
    tree_bkb <- ape::drop.tip(tree_bkb, tree_bkb$tip.label[stringr::str_detect(tree_bkb$tip.label, "_DELETE$")])
  }

  final_tree <- phytools::force.ultrametric(tree_bkb, message = FALSE)
  final_tree_path <- paste0(out_prefix, "_final_tree.nwk")
  ape::write.tree(final_tree, final_tree_path)
  .csta_msg("Final ChronoSTA-grafted tree written: ", final_tree_path, logger = logger)

  monoph_tree <- NULL
  monophyletic_tree_path <- NA_character_
  if (isTRUE(monoph_restore)) {
    .csta_section("ChronoSTA grafting: monophyly cleanup", logger = logger)
    taxo <- data.frame(species = final_tree$tip.label, genus = stringr::str_split_fixed(final_tree$tip.label, "_", n = 2)[, 1])
    monophy <- MonoPhy::AssessMonophyly(final_tree, taxo)
    monoph_tree <- restore_monophyly(
      final_tree = final_tree,
      ref_tree = tree_ref,
      source_trees = list_tre,
      genus_sep = "_",
      min_tips_per_genus = 1,
      max_passes = 1000,
      exc = paraph_exc
    )
    monophyletic_tree_path <- paste0(out_prefix, "_final_tree_monophyletic.nwk")
    ape::write.tree(monoph_tree$tree, monophyletic_tree_path)
    utils::write.table(monoph_tree$removed, paste0(out_prefix, "_monophyly_removed.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    .csta_msg("Monophyly-cleaned tree written: ", monophyletic_tree_path, logger = logger)

    if (isTRUE(plot_monoph)) {
      mono2 <- MonoPhy::AssessMonophyly(monoph_tree$tree, taxo[taxo$species %in% monoph_tree$tree$tip.label, ])
      MonoPhy::PlotMonophyly(monophy, final_tree, taxlevels = 1, plot.type = "monophyly", monocoll = TRUE,
                             cex = 0.4, label.offset = 2, ladderize = TRUE, PDF = TRUE,
                             PDF_filename = paste0(out_prefix, "_monoph_plot_final_tree.pdf"))
      MonoPhy::PlotMonophyly(mono2, monoph_tree$tree, taxlevels = 1, plot.type = "monophyly", monocoll = TRUE,
                             cex = 0.4, label.offset = 2, ladderize = TRUE, PDF = TRUE,
                             PDF_filename = paste0(out_prefix, "_monoph_plot_final_tree_monoph.pdf"))
    }
  }

  result_tree <- if (isTRUE(monoph_restore)) monoph_tree$tree else final_tree
  result <- list(
    tree = result_tree,
    raw_tree = final_tree,
    overlap = data_overlap,
    split_log = splits_data,
    monophyly = monoph_tree,
    paths = list(
      final_tree = final_tree_path,
      monophyletic_tree = monophyletic_tree_path,
      overlap = paste0(out_prefix, "_overlap.tsv"),
      splits = paste0(out_prefix, "_splits.tsv"),
      log = logger$file
    ),
    parameters = list(
      recalibrate = recalibrate, split_gen = split_gen, split_sbt = split_sbt,
      monoph_restore = monoph_restore, prefuse = prefuse, ultr_tol = ultr_tol,
      calib_ratio = calib_ratio, calib_min_nodes = calib_min_nodes,
      calib_max_nodes = calib_max_nodes, margin_nod = margin_nod,
      tolerance_pct = tolerance_pct, tolerance_min = tolerance_min,
      lambda = lambda, split_seed = split_seed, ref_weight = ref_weight
    )
  )
  class(result) <- c("chronosta_grafting_result", class(result))
  .csta_msg("ChronoSTA grafting complete. Final tips: ", ape::Ntip(result_tree), logger = logger)
  result
}

# Backward-compatible alias. Existing scripts can still call sta_sp_grafting(),
# but the function is now non-interactive and creates logs/structured output.
sta_sp_grafting <- function(tree_ref = NULL,
                            tree_folder = NULL,
                            csta_folder = NULL,
                            output_folder = NULL,
                            ...) {
  run_chronosta_grafting(
    tree_ref = tree_ref,
    tree_folder = tree_folder,
    csta_folder = csta_folder,
    output_folder = output_folder,
    ...
  )
}
