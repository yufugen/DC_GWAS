# process_ldblocks.R
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- gsub("_ldblocks.breaks", "_blocks", input_file)

# Define your functions here
load.breaks <- function(prefix, chromosome=NULL) {
  file.name <- paste0(prefix, chromosome, ".breaks")
  if (!file.exists(file.name)) {
    if (!is.null(chromosome) && chromosome == 23) return(load.breaks(prefix, "X"))
    else return(NULL)
  } else return(read.table(file.name, header=T))
}

make.blocks <- function(breaks) {
  N <- dim(breaks)[1]
  blocks <- data.frame(
    start=breaks$POSITION[-N],
    stop=breaks$POSITION[-1]-1,
    size.filt=breaks$INDEX_FILT[-1] - breaks$INDEX_FILT[-N],
    size.all=breaks$INDEX_ALL[-1] - breaks$INDEX_ALL[-N],
    stringsAsFactors=FALSE
  )
  return(blocks)
}

# Load the data
prefix <- gsub(".breaks", "", input_file)
breaks_data <- load.breaks(prefix = prefix)

# Check if data is loaded and then process
if (!is.null(breaks_data)) {
  blocks <- make.blocks(breaks_data)
  # Write to CSV
  write.table(blocks, output_file, sep="\t", row.names = FALSE, quote=FALSE)
} else {
  cat("Failed to load breaks data for:", input_file, "\n")
}

