# =========================================================================
# plot_pie_reconstructions.R  (DEC+J version)
# Reads from "your_directory_here/readable/" and writes the pie PDF there too.
# =========================================================================

library(ape)
library(RColorBrewer)

biogeo_dir   <- "your_directory_here"
readable_dir <- file.path(biogeo_dir, "readable")
out_dir      <- readable_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

find_file <- function(name) {
  for (d in c(readable_dir, biogeo_dir)) {
    p <- file.path(d, name)
    if (file.exists(p)) return(p)
  }
  stop("File not found: ", name)
}

# ---- 14 distinct colours ----
base14 <- c(brewer.pal(12, "Paired"), "#7F3F00", "#B0B0B0")

area_info <- data.frame(
  code = LETTERS[1:14],
  name = c("North America",
           "Central America / Caribbean",
           "NW South America",
           "South America (remaining)",
           "Europe",
           "North Africa",
           "Sub-Saharan Africa",
           "Madagascar",
           "Middle East / Arabia",
           "Asia (Central / East)",
           "Indian subcontinent",
           "Southeast Asia",
           "Australasia",
           "Antarctica"),
  col  = base14,
  stringsAsFactors = FALSE
)
area_cols <- setNames(area_info$col, area_info$code)
other_col <- "#F2F2F2"

# ---- Dictionary ----
dict_fn <- find_file("node_state_dictionary.csv")
dict <- read.csv(dict_fn, stringsAsFactors = FALSE)
state_to_range <- setNames(dict$range, dict$state_name)

# =========================================================================
# Read the FULL matrix and compute the top-3 states per node
# NOTE: DEC+J output uses the "results_DECJ_<label>_..." prefix.
# =========================================================================
read_top3_from_full <- function(label) {
  fn <- find_file(paste0("results_DECJ_", label, "_ancestral_probs_FULL.csv"))
  probs <- read.csv(fn, row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE)
  
  probs[] <- lapply(probs, function(x) {
    if (is.character(x)) as.numeric(gsub("%", "", x, fixed = TRUE)) / 100 else x
  })
  
  do.call(rbind, lapply(seq_len(nrow(probs)), function(i) {
    p   <- as.numeric(probs[i, ])
    ord <- order(p, decreasing = TRUE)[1:3]
    data.frame(
      node        = rownames(probs)[i],
      ML_state    = colnames(probs)[ord[1]],
      ML_prob     = round(p[ord[1]], 4),
      alt_2_state = colnames(probs)[ord[2]],
      alt_2_prob  = round(p[ord[2]], 4),
      alt_3_state = colnames(probs)[ord[3]],
      alt_3_prob  = round(p[ord[3]], 4),
      stringsAsFactors = FALSE
    )
  }))
}

# =========================================================================
build_pie_matrix <- function(label) {
  
  tree <- readRDS(find_file(paste0("tree_", label, ".rds")))
  ntip <- length(tree$tip.label)
  nnode<- tree$Nnode
  
  tipr <- read.csv(find_file(paste0("tipranges_", label, "_readable.csv")),
                   row.names = 1, check.names = FALSE)
  rownames(tipr) <- gsub("^\ufeff", "", rownames(tipr))
  tipr <- tipr[tree$tip.label, , drop = FALSE]
  
  top3 <- read_top3_from_full(label)
  
  cols_all <- c(LETTERS[1:14], "other")
  pie_mat  <- matrix(0, nrow = ntip + nnode, ncol = length(cols_all),
                     dimnames = list(NULL, cols_all))
  
  # ---- tips ----
  for (i in seq_len(ntip)) {
    present <- colnames(tipr)[tipr[i, ] == 1]
    present <- present[present %in% LETTERS[1:14]]
    if (length(present) > 0) pie_mat[i, present] <- 1 / length(present)
  }
  
  # ---- internal nodes ----
  internal_labels <- paste0("node_", (ntip + 1):(ntip + nnode))
  top3_int <- top3[match(internal_labels, top3$node), ]
  
  safe_range <- function(state) {
    if (is.null(state) || length(state) == 0 || is.na(state)) return("")
    r <- state_to_range[state]
    if (length(r) == 0 || is.na(r)) return("")
    r
  }
  
  for (j in seq_len(nnode)) {
    idx <- ntip + j
    
    s1 <- top3_int$ML_state[j];    s2 <- top3_int$alt_2_state[j]
    p1 <- top3_int$ML_prob[j];     p2 <- top3_int$alt_2_prob[j]
    
    r1 <- safe_range(s1)
    r2 <- safe_range(s2)
    
    a1 <- strsplit(r1, "")[[1]]; a1 <- a1[a1 %in% LETTERS[1:14]]
    a2 <- strsplit(r2, "")[[1]]; a2 <- a2[a2 %in% LETTERS[1:14]]
    
    if (length(a1) > 0) pie_mat[idx, a1] <- pie_mat[idx, a1] + p1 / length(a1)
    if (length(a2) > 0) pie_mat[idx, a2] <- pie_mat[idx, a2] + p2 / length(a2)
    
    used <- sum(pie_mat[idx, LETTERS[1:14]])
    if (used < 1) pie_mat[idx, "other"] <- 1 - used
  }
  
  list(tree = tree, pie_mat = pie_mat)
}

# =========================================================================
make_pie_plot <- function(label) {
  
  obj     <- build_pie_matrix(label)
  tree    <- obj$tree
  pie_mat <- obj$pie_mat
  
  ntip  <- length(tree$tip.label)
  nnode <- tree$Nnode
  
  # Prefix the output with DECJ_ so it does not overwrite the DEC version
  out_pdf <- file.path(out_dir, paste0("DECJ_", label, "_pie_reconstruction_FULL.pdf"))
  
  pdf(out_pdf, width = 18, height = 40)
  par(mar = c(1, 1, 3, 1), xpd = NA)
  
  plot(tree,
       cex            = 1.15,
       label.offset   = 1.2,
       show.tip.label = TRUE,
       lwd            = 0.6,
       main           = paste0("Eusuchian biogeography (", label,
                               ") - DEC+J time-stratified"))
  
  tiplabels(pie       = pie_mat[1:ntip, , drop = FALSE],
            piecol    = c(area_cols, other = other_col),
            pieborder = NA,
            cex       = 0.20)
  
  nodelabels(pie       = pie_mat[(ntip + 1):(ntip + nnode), , drop = FALSE],
             piecol    = c(area_cols, other = other_col),
             pieborder = NA,
             cex       = 0.32)
  
  legend(x = "topleft",
         legend = sprintf("%s  -  %s", area_info$code, area_info$name),
         pt.bg  = area_info$col,
         pch    = 22,
         col    = NA,
         pt.cex = 1.6,
         cex    = 0.85,
         bty    = "n",
         title  = "Geographic areas")
  
  legend(x = "bottomleft",
         legend = c("Tip pie: observed range",
                    "Internal pie: top 2 states",
                    "Light grey: remaining probability"),
         pch    = 21,
         pt.bg  = c("#888888", "#888888", other_col),
         col    = NA,
         pt.cex = 1.6,
         cex    = 0.85,
         bty    = "n",
         title  = "How to read the pies")
  
  dev.off()
  message("Wrote ", out_pdf)
}

# =========================================================================
make_pie_plot("mixed")

message("Done. PDF is in: ", out_dir)
