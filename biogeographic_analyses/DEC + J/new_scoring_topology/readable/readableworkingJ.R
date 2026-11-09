# =========================================================================
# all_states_ordered_DECJ.R  (DEC+J only, new tree)
# Reads from d_j_P/J/readable/, writes tables there too.
# =========================================================================

library(dplyr)
library(tidyr)

readable <- "your_directory_here/readable"
setwd(readable)

threshold <- 0.0001   # keep states with probability > 0.01%

# ---- dictionary ----
dict_csv <- "node_state_dictionary.csv"
if (!file.exists(dict_csv))
  dict_csv <- file.path("..", "node_state_dictionary.csv")
if (!file.exists(dict_csv)) stop("Cannot find node_state_dictionary.csv")

dict <- read.csv(dict_csv, stringsAsFactors = FALSE)
state_to_range <- setNames(dict$range, dict$state_name)

# =========================================================================
# Input files for the DEC+J run
# =========================================================================
find_inputs <- function(label) {
  probs_candidates <- c(
    paste0("results_DECJ_", label, "_ancestral_probs_FULL.csv"),
    paste0("results_",      label, "_ancestral_probs_FULL.csv")
  )
  nodes_candidates <- c(
    paste0("tree_",         label, "_nodes.csv")
  )
  
  probs_file <- probs_candidates[file.exists(probs_candidates)][1]
  nodes_file <- nodes_candidates[file.exists(nodes_candidates)][1]
  
  if (length(probs_file) == 0 || is.na(probs_file))
    stop("No ancestral_probs_FULL file found for ", label,
         " (tried: ", paste(probs_candidates, collapse = ", "), ")")
  if (length(nodes_file) == 0 || is.na(nodes_file))
    stop("No tree nodes file found for ", label,
         " (tried: ", paste(nodes_candidates, collapse = ", "), ")")
  
  list(probs = probs_file, nodes = nodes_file)
}

# =========================================================================
# Main function
# =========================================================================
process_label <- function(label) {
  
  message("Processing ", label, " ...")
  
  inputs <- find_inputs(label)
  message("  using probabilities: ", inputs$probs)
  message("  using tree nodes   : ", inputs$nodes)
  
  out_prefix <- paste0("DECJ_", label)
  
  out_long <- paste0(out_prefix, "_all_states_ordered_long.csv")
  out_wide <- paste0(out_prefix, "_all_states_ordered_wide.csv")
  out_txt  <- paste0(out_prefix, "_all_states_ordered.txt")
  
  # ---- 1. read probability matrix ----
  probs <- read.csv(inputs$probs, row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE)
  
  probs[] <- lapply(probs, function(x) {
    if (is.character(x)) {
      as.numeric(gsub("%", "", x, fixed = TRUE)) / 100
    } else {
      x
    }
  })
  
  # ---- 2. read tree node table ----
  nodes <- read.csv(inputs$nodes, stringsAsFactors = FALSE)
  names(nodes) <- tolower(names(nodes))
  nodes$node   <- as.numeric(nodes$node)
  nodes$parent <- as.numeric(nodes$parent)
  nodes$type   <- as.character(nodes$type)
  nodes$label  <- as.character(nodes$label)
  
  children <- split(nodes$node, nodes$parent)
  children <- children[!is.na(names(children))]
  
  get_desc_tips <- function(id) {
    if (nodes$type[nodes$node == id][1] == "tip") {
      return(nodes$label[nodes$node == id][1])
    }
    ch <- children[[as.character(id)]]
    if (is.null(ch)) return(character(0))
    unlist(lapply(ch, get_desc_tips))
  }
  
  # ---- 3. build long table ----
  internal_nodes <- nodes$node[nodes$type == "internal"]
  long_list <- list()
  
  for (nd in internal_nodes) {
    node_label <- paste0("node_", nd)
    if (!node_label %in% rownames(probs)) next
    
    p <- as.numeric(probs[node_label, ])
    names(p) <- colnames(probs)
    p <- sort(p, decreasing = TRUE)
    p <- p[p > threshold]
    
    if (length(p) == 0) next
    
    desc <- sort(get_desc_tips(nd))
    desc_str <- if (length(desc) <= 5) {
      paste(desc, collapse = ", ")
    } else {
      paste0(paste(head(desc, 5), collapse = ", "),
             ", ... (n=", length(desc), ")")
    }
    
    ranges <- state_to_range[names(p)]
    ranges[is.na(ranges)] <- ""
    
    long_list[[length(long_list) + 1]] <- data.frame(
      node             = nd,
      node_label       = node_label,
      n_tips           = length(desc),
      descendant_tips  = desc_str,
      rank             = seq_along(p),
      state            = names(p),
      range            = ranges,
      probability      = round(as.numeric(p), 6),
      stringsAsFactors = FALSE
    )
  }
  
  long_df <- bind_rows(long_list)
  write.csv(long_df, out_long, row.names = FALSE)
  
  # ---- 4. wide format ----
  wide_df <- long_df %>%
    group_by(node, node_label, n_tips, descendant_tips) %>%
    summarise(
      states_ordered = paste0(range, " (", sprintf("%.4f", probability), ")",
                              collapse = "; "),
      .groups = "drop"
    )
  write.csv(wide_df, out_wide, row.names = FALSE)
  
  # ---- 5. plain-text summary ----
  sink(out_txt)
  for (i in seq_len(nrow(wide_df))) {
    cat("=====================================================\n")
    cat(wide_df$node_label[i], "  (n_tips =", wide_df$n_tips[i], ")\n")
    cat("Tips: ", wide_df$descendant_tips[i], "\n", sep = "")
    cat("States ordered by probability:\n")
    cat(wide_df$states_ordered[i], "\n\n", sep = "")
  }
  sink()
  
  message("  wrote: ", out_long, ", ", out_wide, ", ", out_txt)
  invisible(list(long = long_df, wide = wide_df))
}

# =========================================================================
# Run for the DEC+J tree
# =========================================================================
res_decj <- process_label("new_DECJ")

cat("\nAll done.\n")
cat("new_DECJ: long =", nrow(res_decj$long),
    "rows, wide =", nrow(res_decj$wide), "nodes\n")
