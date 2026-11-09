if (!requireNamespace("ape", quietly = TRUE))       install.packages("ape")
if (!requireNamespace("paleotree", quietly = TRUE)) install.packages("paleotree")
library(ape)
library(paleotree)

# ---- PATHS ----
folder_path <- #YOUR_DIRECTORY_HERE#
age_path    <- paste0(folder_path, "Ages.csv")
tree_files <- c("new_acutus_tree.txt")

# =========================================================================
# Load ages
# =========================================================================
cat(">>> Loading Ages.csv...\n")
age_data <- read.csv(age_path, check.names = FALSE)
rownames(age_data) <- age_data[, 1]
time_matrix <- matrix(
  as.numeric(as.matrix(age_data[, 2:3])),
  ncol = 2,
  dimnames = list(rownames(age_data), c("FAD", "LAD"))
)
cat("    Loaded", nrow(time_matrix), "taxa.\n\n")

# =========================================================================
# Loader: works with TNT "tread ..." files AND with NEXUS tree blocks
# =========================================================================
load_one_tree <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  
  # Drop any line that begins with tread / proc / macro / # / ;
  keep <- !grepl("^\\s*(tread|proc|macro|#|;)", lines, ignore.case = TRUE)
  lines <- lines[keep]
  
  # Find the first line that begins with "(" — that's the start of the Newick
  newick_start <- grep("^\\s*\\(", lines)
  if (length(newick_start) > 0) {
    block <- paste(lines[newick_start[1]:length(lines)], collapse = " ")
  } else {
    ti <- grep("^\\s*tree\\s+\\S+\\s*=", lines, perl = TRUE)
    if (length(ti) > 0) {
      block <- paste(lines[ti[1]:length(lines)], collapse = " ")
    } else {
      stop("No Newick string found")
    }
  }
  
  # Remove any leftover single-quoted strings and square-bracket comments
  block <- gsub("'[^']*'", " ", block)
  block <- gsub("\\[[^]]*\\]", " ", block)
  
  # Repair split-tip errors (space between genus and species)
  block <- gsub("Crocodylus\\s+acutus", "Crocodylus_acutus", block)
  
  # Slice to Newick string
  sp <- regexpr("\\(", block, perl = TRUE)[1]
  ep <- regexpr(";",  block, perl = TRUE)[1]
  if (sp == -1) stop("Could not locate '('")
  if (ep == -1) {
    newick <- substr(block, sp, nchar(block))
    newick <- paste0(newick, ";")
  } else {
    newick <- substr(block, sp, ep)
  }
  
  # Whitespace -> comma
  newick <- gsub("\\s+", ",", newick)
  
  # Fix missing commas between structural elements
  newick <- gsub("\\)(?=\\()", "),", newick, perl = TRUE)
  newick <- gsub("(?<=[^,(;])\\(", ",(", newick, perl = TRUE)
  newick <- gsub("\\)(?=[^,;)])", "),", newick, perl = TRUE)
  
  # Clean up syntax artifacts
  newick <- gsub(",+", ",", newick)
  newick <- gsub("[(],+", "(", newick)
  newick <- gsub(",\\)", ")", newick)
  newick <- gsub(",;", ";", newick)
  
  cat("    Newick length:", nchar(newick), "chars\n")
  
  read.tree(text = newick)
}

# =========================================================================
# Calibrate
# =========================================================================
for (tf in tree_files) {
  cat("==================================================\n")
  cat("Processing:", tf, "\n")
  cat("==================================================\n")
  
  base <- sub("\\.txt$", "", tf)
  tr   <- tryCatch(load_one_tree(paste0(folder_path, tf)),
                   error = function(e) { cat("  !!!", conditionMessage(e), "\n"); NULL })
  if (is.null(tr) || length(tr$tip.label) < 2) { cat("  skip.\n\n"); next }
  
  cat("    Tip count:", length(tr$tip.label), "\n")
  cat("    Nnode    :", tr$Nnode, "\n")
  cat("    First 5 tips:",
      paste(head(tr$tip.label, 5), collapse = " | "), "\n\n")
  
  # Match tips to age table
  matched <- intersect(tr$tip.label, rownames(time_matrix))
  if (length(matched) < 2) {
    fixed <- gsub(" ", "_", tr$tip.label)
    if (length(intersect(fixed, rownames(time_matrix))) > length(matched)) {
      tr$tip.label <- fixed
      matched <- intersect(fixed, rownames(time_matrix))
    }
  }
  cat("    Matched:", length(matched), "of", length(tr$tip.label), "\n")
  
  if (length(matched) < 2) next
  
  unmatched <- setdiff(tr$tip.label, rownames(time_matrix))
  if (length(unmatched) > 0) {
    cat("    Unmatched tips (not in Ages.csv):\n")
    print(unmatched)
  }
  
  pruned_tree <- keep.tip(tr, matched)
  pruned_time <- time_matrix[pruned_tree$tip.label, , drop = FALSE]
  
  calibrated <- tryCatch(
    timePaleoPhy(pruned_tree, pruned_time,
                 type = "mbl", vartime = 2, add.term = TRUE),
    error = function(e) { cat("    timePaleoPhy:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(calibrated)) next
  
  write.tree(calibrated,
             file = paste0(folder_path, "calibrated_", base, ".tre"))
  
  n_tips     <- length(calibrated$tip.label)
  pdf_height <- max(12, n_tips * 0.28)
  pdf_width  <- 16
  
  pdf(paste0(folder_path, "Calibrated_", base, ".pdf"),
      width = pdf_width, height = pdf_height)
  
  par(mar = c(5, 2, 4, 14))
  
  tryCatch({
    plot(calibrated,
         cex          = 0.75,
         edge.width   = 1.2,
         label.offset = 0.3,
         no.margin    = FALSE,
         main         = paste("Calibrated", base))
    axisPhylo()
  }, error = function(e) cat("    plot:", conditionMessage(e), "\n"))
  
  dev.off()
  
  cat("    Wrote Calibrated_", base, ".pdf / calibrated_", base, ".tre\n\n",
      sep = "")
}
cat("DONE.\n")
