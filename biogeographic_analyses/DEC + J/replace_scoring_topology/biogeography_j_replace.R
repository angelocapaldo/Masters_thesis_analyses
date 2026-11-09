# =========================================================================
# BioGeoBEARS DEC+J — Eusuchian biogeography (tree: replace)
#
# Data: croco_occ.txt (PHYLIP), biogeography_ages.txt (8 matrices),
#       calibrated_newmix_last_tree_acutus.tre
# =========================================================================

library(ape)
library(BioGeoBEARS)

# ---- PATHS ----
work_dir <- "your_directory_here"

setwd(work_dir)

# =========================================================================
# HELPERS
# =========================================================================

read_tree_robust <- function(fn) {
  tr <- tryCatch(read.tree(fn), error = function(e) NULL)
  if (is.null(tr)) tr <- read.nexus(fn)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  tr$edge.length[tr$edge.length <= 0] <- 1e-6
  tr
}

read_phylip_geog <- function(fn) {
  lines    <- readLines(fn, warn = FALSE)
  header   <- strsplit(trimws(lines[1]), "\\s+")[[1]]
  ntax     <- as.integer(header[1])
  narea    <- as.integer(header[2])
  area_str <- sub(".*\\(([^)]*)\\).*", "\\1", lines[1])
  area_names <- strsplit(trimws(area_str), "\\s+")[[1]]

  taxa <- character(ntax)
  mat  <- matrix(0L, nrow = ntax, ncol = narea)
  for (i in 1:ntax) {
    parts    <- strsplit(trimws(lines[1 + i]), "\\s+")[[1]]
    taxa[i]  <- parts[1]
    mat[i, ] <- as.integer(strsplit(parts[2], "")[[1]])
  }
  rownames(mat) <- taxa
  colnames(mat) <- area_names
  mat
}

write_geog_phylip <- function(occ_sub, area_names, fn) {
  lines <- c(
    paste0(nrow(occ_sub), "\t", ncol(occ_sub), "\t(",
           paste(area_names, collapse = " "), ")"),
    paste(rownames(occ_sub),
          apply(occ_sub, 1, paste, collapse = ""),
          sep = "\t")
  )
  writeLines(lines, fn)
}

parse_biogeography_ages <- function(fn, nareas = 14) {
  lines   <- readLines(fn, warn = FALSE)
  age_idx <- grep("^\\s*[0-9]+\\.?[0-9]*\\s*$", lines)
  mats    <- vector("list", length(age_idx))

  for (i in seq_along(age_idx)) {
    s   <- age_idx[i]
    mat <- matrix(NA_real_, nareas, nareas)
    for (j in 1:nareas) {
      parts    <- strsplit(trimws(lines[s + 1 + j]), "\\s+")[[1]]
      mat[j, ] <- as.numeric(parts[-1])
    }
    rownames(mat) <- LETTERS[1:nareas]
    colnames(mat) <- LETTERS[1:nareas]
    mats[[i]] <- mat
  }
  mats
}

# =========================================================================
# 1. LOAD OCCURRENCE MATRIX (PHYLIP)
# =========================================================================

occ_mat        <- read_phylip_geog("croco_occ.txt")
area_names_all <- colnames(occ_mat)

# =========================================================================
# 2. LOAD THE EIGHT DISPERSAL MATRICES
# =========================================================================

mats <- parse_biogeography_ages("biogeography_ages.txt", nareas = 14)
stopifnot(length(mats) == 8)

# =========================================================================
# 3. SINGLE CALIBRATED TREE -> DEC+J run
# =========================================================================

tree_files <- c(
  replace = "calibrated_replace_last_tree_acutus.tre"
)

POLYTOMY_SEED <- 42

for (label in names(tree_files)) {

  tree_file <- tree_files[[label]]

  # ---- All outputs go to the working directory ----
  out_pdf      <- paste0("Eusuchian_Biogeography_DECJ_", label, ".pdf")
  rds_res      <- paste0("results_DECJ_",   label, ".rds")
  rds_tree     <- paste0("tree_",           label, ".rds")
  rds_trng     <- paste0("tipranges_",      label, ".rds")
  geog_matched <- paste0("croco_occ_MATCHED_",         label, ".txt")
  tree_matched <- paste0("calibrated_tree_matched_",   label, ".tre")
  times_file   <- paste0("timeperiods_FIXED_",         label, ".txt")

  # Skip if already done
  if (file.exists(rds_res) && file.exists(rds_tree) && file.exists(rds_trng)) {
    message("=== ", label, ": DEC+J results already exist (", rds_res,
            "), skipping to save time ===")
    next
  }

  message("=== DEC+J: ", label, " (", tree_file, ") -> ", out_pdf, " ===")

  ok <- tryCatch({

    # ---- 3a. Read tree ----
    tr <- read_tree_robust(tree_file)
    tr$edge.length[tr$edge.length <= 0] <- 1e-6

    # ---- 3a-bis. Resolve polytomies ----
    if (!is.binary(tr)) {
      message("  tree has polytomies; resolving with multi2di() (seed = ",
              POLYTOMY_SEED, ")")
      set.seed(POLYTOMY_SEED)
      tr <- multi2di(tr, random = TRUE)
    }
    tr$edge.length[tr$edge.length <= 0] <- 1e-6

    # ---- 3b. Time boundaries ----
    root_age   <- max(node.depth.edgelength(tr))
    boundary_8 <- max(145.001, ceiling(root_age) + 2)
    time_boundaries <- c(5.001, 23.001, 34.001, 56.001,
                         66.001, 100.001, 120.001, boundary_8)
    writeLines(as.character(time_boundaries), times_file)
    message("  root age = ", round(root_age, 2),
            " | boundaries = ", paste(time_boundaries, collapse = " "))

    # ---- 3c. Match taxa ----
    common <- intersect(tr$tip.label, rownames(occ_mat))

    missing_from_tree <- setdiff(rownames(occ_mat), tr$tip.label)
    missing_from_occ  <- setdiff(tr$tip.label,   rownames(occ_mat))
    if (length(missing_from_tree) || length(missing_from_occ)) {
      message("  Name mismatches:")
      if (length(missing_from_tree))
        message("    in occ, not in tree: ", paste(missing_from_tree, collapse = ", "))
      if (length(missing_from_occ))
        message("    in tree, not in occ: ", paste(missing_from_occ, collapse = ", "))
    }

    if (length(common) < length(tr$tip.label)) {
      tr <- drop.tip(tr, setdiff(tr$tip.label, common))
    }
    occ_sub <- occ_mat[common, , drop = FALSE]

    write_geog_phylip(occ_sub, area_names_all, geog_matched)
    write.tree(tr, tree_matched)

    # ---- 3d. Configure BioGeoBEARS (base DEC run object) ----
    run_obj <- define_BioGeoBEARS_run()
    run_obj$trfn    <- tree_matched
    run_obj$geogfn  <- geog_matched
    run_obj$timesfn <- times_file
    run_obj$dispersal_multipliers_fn <- NULL

    run_obj$max_range_size     <- 3
    run_obj$min_branchlength   <- 0.000001
    run_obj$include_null_range <- TRUE

    run_obj$speedup          <- TRUE
    run_obj$use_optimx       <- "GenSA"
    run_obj$num_cores_to_use <- 1
    run_obj$calc_ancprobs    <- TRUE

    # ---- 3d-bis. Turn this into DEC+J ----
    # 'j' is the founder-event speciation parameter.
    run_obj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
    run_obj$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.01
    run_obj$BioGeoBEARS_model_object@params_table["j", "est"]  <- 0.01
    run_obj$BioGeoBEARS_model_object@params_table["j", "min"]  <- 0
    run_obj$BioGeoBEARS_model_object@params_table["j", "max"]  <- 1

    # ---- 3e. Read files & attach tree ----
    run_obj <- readfiles_BioGeoBEARS_run(run_obj)
    run_obj$phy <- tr

    # ---- 3f. Assign the 8 dispersal matrices ----
    tipranges  <- getranges_from_LagrangePHYLIP(lgdata_fn = run_obj$geogfn)
    area_names <- colnames(tipranges@df)

    my_dispersal_mats <- lapply(mats, function(m) {
      rownames(m) <- area_names
      colnames(m) <- area_names
      m
    })
    run_obj$list_of_dispersal_multipliers_mats <- my_dispersal_mats

    # ---- 3g. Section tree ----
    run_obj <- section_the_tree(
      inputs             = run_obj,
      make_master_table  = TRUE,
      plot_pieces        = FALSE,
      cut_fossils        = FALSE,
      fossils_older_than = 0.001
    )

    # ---- 3h. Optimization (DEC+J) ----
    results_DECJ <- bears_optim_run(run_obj)

    # ---- 3h-bis. Save immediately ----
    saveRDS(results_DECJ, rds_res)
    saveRDS(tr,           rds_tree)
    saveRDS(tipranges,    rds_trng)
    message("  saved ", rds_res, " | ", rds_tree, " | ", rds_trng)

    # ---- 3i. Export PDF ----
    pdf(file = out_pdf, width = 12, height = 25)
    tryCatch({
      plot_BioGeoBEARS_results(
        results_object     = results_DECJ,
        analysis_titletxt  = paste("Eusuchian Biogeography (Time-Stratified DEC+J) -", label),
        plotwhat           = "pie",
        label.offset       = 0.2,
        tipcex             = 0.5,
        statecex           = 0.7,
        splitcex           = 0.7,
        titlecex           = 1.0,
        plotsplits         = FALSE,
        include_null_range = TRUE,
        tr                 = run_obj$phy,
        tipranges          = tipranges
      )
    }, finally = {
      if (names(dev.cur()) == "pdf") dev.off()
    })

    TRUE

  }, error = function(e) {
    while (names(dev.cur()) == "pdf") dev.off()
    message("FAILED on ", label, ": ", conditionMessage(e))
    FALSE
  })

  if (isTRUE(ok)) {
    message("SUCCESS -> ", out_pdf)
  } else {
    message("No PDF written for ", label)
  }
}

