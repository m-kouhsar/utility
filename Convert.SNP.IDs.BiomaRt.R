#!/usr/bin/env Rscript

message("Loading required libraries...")
suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) < 3) {
  cat("\nUsage: Rscript Convert.SNP.IDs.BiomaRt.R <input_file> <id_col_idx> <genome_ver>\n")
  cat("Example: Rscript Convert.SNP.IDs.BiomaRt.R input.txt 1 hg38\n\n")
  quit(status = 1)
}

Input.File <- arguments[1]
ID.col     <- as.numeric(arguments[2])
Genome.Ver <- arguments[3]

# Validate genome version choice
if (!Genome.Ver %in% c("hg19", "hg38")) stop("Genome.Ver must be 'hg19' or 'hg38'")

message("Reading input data...")
snp_dt <- fread(file = Input.File, header = FALSE, stringsAsFactors = FALSE)

if (ID.col > ncol(snp_dt)) stop("Provided column index exceeds the number of columns in the file.")

# Extract unique target IDs to minimize web query payload size
raw_inputs <- unique(snp_dt[[ID.col]])
raw_inputs <- raw_inputs[!is.na(raw_inputs) & raw_inputs != ""]

# Automate detection: Split inputs into rsids vs coordinates dynamically
is_rsid      <- grepl("^rs\\d+", raw_inputs, ignore.case = TRUE)
rs_inputs    <- raw_inputs[is_rsid]
coord_inputs <- raw_inputs[!is_rsid]

message("Connecting to Ensembl BioMart database...")
if (Genome.Ver == "hg19") {
  snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "https://grch37.ensembl.org")
} else {
  snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "https://www.ensembl.org")
}

snp_attributes <- c("refsnp_id", "chr_name", "chrom_start")
mapping_table  <- data.table(original_id = character(), converted_id = character())

# --- Step 1: Query for inputs that are rsIDs -> Convert to chr:pos ---
if (length(rs_inputs) > 0) {
  message("Processing ", length(rs_inputs), " rsIDs...")
  tryCatch({
    query_res <- getBM(attributes = snp_attributes, filters = "snp_filter", values = rs_inputs, mart = snp_mart)
    if (nrow(query_res) > 0) {
      rs_map <- as.data.table(query_res)
      rs_map[, converted_id := paste(chr_name, chrom_start, sep = ":")]
      setnames(rs_map, "refsnp_id", "original_id")
      mapping_table <- rbind(mapping_table, rs_map[, .(original_id, converted_id)])
    }
  }, error = function(e) { warning("biomaRt query for rsIDs failed: ", e$message) })
}

# --- Step 2: Query for inputs that are chr:pos -> Convert to rsID ---
if (length(coord_inputs) > 0) {
  message("Processing ", length(coord_inputs), " coordinate pairs...")
  
  # Normalize and split coordinates cleanly
  clean_coords <- gsub("^chr", "", coord_inputs, ignore.case = TRUE)
  split_parts  <- tstrsplit(clean_coords, ":")
  
  chrs   <- split_parts[[1]]
  starts <- as.numeric(split_parts[[2]])
  
  values_list <- list(chr_name = chrs, start = starts, end = starts)
  
  tryCatch({
    query_res <- getBM(attributes = snp_attributes, filters = c("chr_name", "start", "end"), values = values_list, mart = snp_mart)
    if (nrow(query_res) > 0) {
      coord_map <- as.data.table(query_res)
      coord_map[, original_id := paste(chr_name, chrom_start, sep = ":")]
      coord_map[, converted_id := refsnp_id]
      
      mapping_table <- rbind(mapping_table, coord_map[, .(original_id, converted_id)])
      
      # Mirror mapping keys if the user's file had explicit "chr" string prefixes
      if (any(grepl("^chr", coord_inputs, ignore.case = TRUE))) {
        coord_map_chr <- copy(coord_map)
        coord_map_chr[, original_id := paste0("chr", original_id)]
        mapping_table <- rbind(mapping_table, coord_map_chr[, .(original_id, converted_id)])
      }
    }
  }, error = function(e) { warning("biomaRt query for coordinates failed: ", e$message) })
}

# --- Step 3: Merge mapping back onto the full dataset ---
message("Merging results...")
mapping_lookup <- unique(mapping_table)

orig_col_name <- colnames(snp_dt)[ID.col]
setnames(snp_dt, orig_col_name, "original_id")

final_dt <- merge(snp_dt, mapping_lookup, by = "original_id", all.x = TRUE)
setnames(final_dt, "original_id", orig_col_name)

# Reposition the new column so it sits next to the parsed input index
col_order <- colnames(final_dt)
target_idx <- which(col_order == orig_col_name)
new_order <- c(col_order[1:target_idx], "converted_id", col_order[(target_idx+1):(length(col_order)-1)])
final_dt <- final_dt[, ..new_order]

# --- Step 4: Write output automatically to the same directory ---
output_file <- paste0(tools::file_path_sans_ext(Input.File), ".ConvSNPID.txt")
message("Writing results to: ", output_file)
fwrite(final_dt, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)
message("Done successfully!")
