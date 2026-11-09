# =========================================================================
# read_rds_make_readable.R
# READ-ONLY. Loads the .rds files from the single-tree analysis and writes
# readable outputs to biogeo/readable/. Nothing is overwritten or modified.
# =========================================================================

library(ape)
library(BioGeoBEARS)

# ---- PATHS ----
biogeo_dir <- "C:/Users/angiu/Documents/pal/TNT/tnt-windows/TNT-bin/trees/biogeo/bio+riomannion/biogeo"
out_dir    <- file.path(biogeo_dir, "readable")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- LABEL ----
# Single-tree analysis, so only one label.
labels <- "biogeographyMIXacutus"

# =========================================================================
# 1. TREES  ->  Newick + node table
# =========================================================================

for (label in labels) {
  tree_rds <- file.path(biogeo_dir, paste0("tree_", label, ".rds"))
  if (!file.exists(tree_rds)) {
    cat("!!! Missing:", tree_rds, "\n")
    next
  }
  tr <- readRDS(tree_rds)
  
  # Newick
  write.tree(tr, file.path(out_dir, paste0("tree_", label, "_readable.tre")))
  
  # Node table: node number, type (tip/internal), label, parent, length
  ntip  <- length(tr$tip.label)
  nnode <- tr$Nnode
  edge  <- tr$edge
  elen  <- tr$edge.length
  
  labels_vec <- c(tr$tip.label,
                  paste0("node_", (ntip + 1):(ntip + nnode)))
  
  node_tbl <- data.frame(
    node      = 1:(ntip + nnode),
    type      = c(rep("tip", ntip), rep("internal", nnode)),
    label     = labels_vec,
    parent    = NA_integer_,
    edge_len  = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(edge))) {
    child <- edge[i, 2]
    node_tbl$parent[child]   <- edge[i, 1]
    node_tbl$edge_len[child] <- elen[i]
  }
  
  write.csv(node_tbl,
            file.path(out_dir, paste0("tree_", label, "_nodes.csv")),
            row.names = FALSE)
  
  cat("--- tree_", label, ".rds ---\n", sep = "")
  cat("  tips:", ntip, " internal nodes:", nnode,
      "  root age:", round(max(node.depth.edgelength(tr)), 3), "\n\n")
}

# =========================================================================
# 2. TIPRANGES  ->  CSV matrix (taxa x areas)
# =========================================================================

for (label in labels) {
  trng_rds <- file.path(biogeo_dir, paste0("tipranges_", label, ".rds"))
  if (!file.exists(trng_rds)) {
    cat("!!! Missing:", trng_rds, "\n")
    next
  }
  tipr <- readRDS(trng_rds)
  
  df <- tipr@df
  write.csv(df,
            file.path(out_dir, paste0("tipranges_", label, "_readable.csv")),
            row.names = TRUE)
  
  cat("--- tipranges_", label, ".rds ---\n", sep = "")
  cat("  taxa:", nrow(df), " areas:", ncol(df), "\n")
  cat("  first 5 taxa:", paste(head(rownames(df), 5), collapse = ", "), "\n")
  cat("  areas:", paste(colnames(df), collapse = " "), "\n\n")
}

# =========================================================================
# 3. RESULTS  ->  parameters, log-likelihood, ancestral probabilities
# =========================================================================

for (label in labels) {
  res_rds <- file.path(biogeo_dir, paste0("results_", label, ".rds"))
  if (!file.exists(res_rds)) {
    cat("!!! Missing:", res_rds, "\n")
    next
  }
  res <- readRDS(res_rds)
  
  cat("=====================================================\n")
  cat("results_", label, ".rds\n", sep = "")
  cat("=====================================================\n")
  cat("Top-level names in results object:\n")
  print(names(res))
  cat("\n")
  
  # --- 3a. Parameter estimates & log-likelihood ---
  opt    <- res$optim_result
  params <- res$params
  
  param_tbl <- data.frame(
    parameter = names(params),
    value     = unlist(params),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  lnL <- NA_real_
  if (!is.null(res$total_loglikelihood))  lnL <- res$total_loglikelihood
  if (is.na(lnL) && !is.null(res$loglik)) lnL <- res$loglik
  if (is.na(lnL) && !is.null(opt$value))  lnL <- -opt$value
  
  param_tbl <- rbind(
    param_tbl,
    data.frame(parameter = "loglikelihood", value = lnL,
               stringsAsFactors = FALSE)
  )
  
  write.csv(param_tbl,
            file.path(out_dir, paste0("results_", label, "_parameters.csv")),
            row.names = FALSE)
  
  cat("Parameters:\n")
  print(param_tbl)
  cat("\n")
  
  # --- 3b. Ancestral state probabilities ---
  P <- res$ML_marginal_prob_each_state_at_branch_top_AT_node
  
  if (is.null(P)) {
    cat("  No ML_marginal_prob_each_state_at_branch_top_AT_node slot found.\n\n")
  } else {
    state_cols <- colnames(P)
    if (is.null(state_cols)) state_cols <- paste0("state_", seq_len(ncol(P)))
    
    P_df <- as.data.frame(P)
    colnames(P_df) <- state_cols
    
    tr    <- readRDS(file.path(biogeo_dir, paste0("tree_", label, ".rds")))
    ntip  <- length(tr$tip.label)
    nnode <- tr$Nnode
    ids   <- c(tr$tip.label, paste0("node_", (ntip + 1):(ntip + nnode)))
    
    if (nrow(P_df) == length(ids)) {
      rownames(P_df) <- ids
    } else {
      rownames(P_df) <- paste0("row_", seq_len(nrow(P_df)))
      cat("  NOTE: nrow(P) =", nrow(P_df),
          "does not equal tips+internals =", length(ids), "\n")
    }
    
    write.csv(P_df,
              file.path(out_dir,
                        paste0("results_", label, "_ancestral_probs_FULL.csv")),
              row.names = TRUE)
    
    # Top-3 summary per node
    top3 <- do.call(rbind, lapply(seq_len(nrow(P_df)), function(i) {
      p   <- as.numeric(P_df[i, ])
      ord <- order(p, decreasing = TRUE)[1:3]
      data.frame(
        node        = rownames(P_df)[i],
        ML_state    = colnames(P_df)[ord[1]],
        ML_prob     = round(p[ord[1]], 4),
        alt_2_state = colnames(P_df)[ord[2]],
        alt_2_prob  = round(p[ord[2]], 4),
        alt_3_state = colnames(P_df)[ord[3]],
        alt_3_prob  = round(p[ord[3]], 4),
        stringsAsFactors = FALSE
      )
    }))
    
    write.csv(top3,
              file.path(out_dir,
                        paste0("results_", label, "_ancestral_probs_top3.csv")),
              row.names = FALSE)
    
    cat("  Wrote ancestral prob matrix:",
        nrow(P_df), "nodes x", ncol(P_df), "states\n")
    cat("  Top-3 summary also written.\n\n")
  }
  
  # --- 3c. Plain-text summary ---
  sink(file.path(out_dir, paste0("results_", label, "_summary.txt")))
  cat("=== results_", label, ".rds ===\n\n", sep = "")
  cat("Parameters:\n"); print(param_tbl)
  cat("\nInputs used:\n")
  if (!is.null(res$inputs$tree)) {
    cat("  tree tips:", length(res$inputs$tree$tip.label), "\n")
  }
  if (!is.null(res$inputs$tipranges)) {
    cat("  tipranges taxa:", nrow(res$inputs$tipranges@df), "\n")
  }
  cat("\nAll top-level names:\n"); print(names(res))
  sink()
  
  cat("--- results_", label, ".rds summary ---\n", sep = "")
  cat("  loglikelihood:", lnL, "\n")
  cat("  parameter rows:", nrow(param_tbl), "\n\n")
}

cat("\nDONE. All readable files are in:\n  ", out_dir, "\n")
cat("Files written:\n")
print(list.files(out_dir))
