#!/usr/bin/env Rscript

# ==============================================================================
# Command Line Script: Bidirectional SNP ID Converter (rsID <-> chr:pos)
# ==============================================================================

# Parse arguments safely
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("\nUsage: Rscript convert_snps.R <file_path> <snp_column> <genome_build> [output_path]\n\n")
  cat("Arguments:\n")
  cat("  <file_path>    Path to the input tab-separated (.txt/.tsv) file.\n")
  cat("  <snp_column>   Name of the column or the 1-based index containing the SNPs.\n")
  cat("  <genome_build> Target genome assembly: 'hg38' or 'hg19'.\n")
  cat("  [output_path]  Optional. Path to save the converted file. If omitted, saves as 'converted_<input_file>'.\n\n")
  cat("Example:\n")
  cat("  Rscript convert_snps.R gwas_data.txt MarkerName hg19\n")
  cat("  Rscript convert_snps.R sumstats.tsv 2 hg38 processed_output.tsv\n\n")
  quit(status = 1)
}

# 1. Load Essential Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

file_path   <- args[1]
snp_col     <- args[2]
genome      <- tolower(args[3])

output_path <- paste0(tools::file_path_sans_ext(file_path), ".ConvSNPID.txt")

# Validate genome version choice
if (!genome %in% c("hg38", "hg19")) {
  stop("Invalid genome build. Must be either 'hg38' or 'hg19'.")
}

# 2. Dynamic package verification based on runtime choices
if (genome == "hg38") {
  if (!requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38", quietly = TRUE)) stop("Missing dependency: SNPlocs.Hsapiens.dbSNP155.GRCh38")
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) stop("Missing dependency: BSgenome.Hsapiens.UCSC.hg38")
  snp_db <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
  genome_db <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else {
  if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37", quietly = TRUE)) stop("Missing dependency: SNPlocs.Hsapiens.dbSNP144.GRCh37")
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) stop("Missing dependency: BSgenome.Hsapiens.UCSC.hg19")
  snp_db <- SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
  genome_db <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
}

# 3. Read Input File
if (!file.exists(file_path)) stop(paste("Input file not found:", file_path))
message("Reading file: ", file_path)
dt <- fread(file_path, sep = "\t")

# Resolve whether snp_col is a name or numerical index position
if (suppressWarnings(!is.na(as.numeric(snp_col)))) {
  col_idx <- as.numeric(snp_col)
  col_name <- colnames(dt)[col_idx]
} else {
  col_name <- snp_col
}

if (!col_name %in% colnames(dt)) stop(paste("Specified SNP column '", snp_col, "' not found in the file headers."))

# Deduplicate target targets to limit lookup latency overhead
raw_snps <- unique(dt[[col_name]])
raw_snps <- raw_snps[!is.na(raw_snps) & raw_snps != ""]

# Separate targets by notation
is_rsid <- grepl("^rs\\d+", raw_snps, ignore.case = TRUE)
rs_inputs <- raw_snps[is_rsid]
coord_inputs <- raw_snps[!is_rsid]

mapping_table <- data.table(input_snp = character(), converted_snp = character())

# --- 4A. Processing rsID -> chr:pos ---
if (length(rs_inputs) > 0) {
  message("Processing ", length(rs_inputs), " rsIDs...")
  clean_rs <- gsub("rs", "", rs_inputs, ignore.case = TRUE)
  
  tryCatch({
    snps_gr <- BSgenome::snpid2grange(snp_db, clean_rs)
    rs_map <- as.data.frame(snps_gr) %>%
      mutate(
        chr = gsub("chr", "", seqnames),
        input_snp = paste0("rs", RefSNP_id),
        converted_snp = paste(chr, start, sep = ":")
      ) %>%
      select(input_snp, converted_snp) %>%
      as.data.table()
    
    mapping_table <- rbind(mapping_table, rs_map)
  }, error = function(e) {
    warning("Error converting rsIDs: ", e$message)
  })
}

# --- 4B. Processing chr:pos -> rsID ---
if (length(coord_inputs) > 0) {
  message("Processing ", length(coord_inputs), " coordinate coordinates...")
  coord_clean <- gsub("^chr", "", coord_inputs, ignore.case = TRUE)
  split_coords <- tstrsplit(coord_clean, ":")
  
  coord_dt <- data.table(
    input_snp = coord_inputs,
    chr = split_coords[[1]],
    pos = as.integer(split_coords[[2]])
  ) %>% filter(!is.na(pos)) %>% as.data.table()
  
  if (nrow(coord_dt) > 0) {
    target_chroms <- if (genome == "hg19") paste0("chr", coord_dt$chr) else coord_dt$chr
    
    tryCatch({
      gr_query <- GenomicRanges::GRanges(
        seqnames = target_chroms,
        ranges = IRanges::IRanges(start = coord_dt$pos, end = coord_dt$pos)
      )
      
      matched_snps <- BSgenome::snpsByOverlaps(snp_db, gr_query)
      coord_map <- as.data.frame(matched_snps) %>%
        mutate(
          chr = gsub("chr", "", seqnames),
          input_snp = paste(chr, start, sep = ":"),
          converted_snp = paste0("rs", RefSNP_id)
        ) %>%
        select(input_snp, converted_snp) %>%
        distinct(input_snp, .keep_all = TRUE) %>% 
        as.data.table()
      
      if (any(grepl(":", gsub("^[^:]:^[^:]:", "", coord_dt$input_snp)))) {
        coord_map <- merge(coord_dt[, .(input_snp, match_key = paste(chr, pos, sep = ":"))], 
                           coord_map, by.x = "match_key", by.y = "input_snp", all.x = TRUE)
        coord_map <- coord_map[, .(input_snp, converted_snp)]
      }
      
      mapping_table <- rbind(mapping_table, coord_map)
    }, error = function(e) {
      warning("Error converting coordinates: ", e$message)
    })
  }
}

# 5. Integrate Results and Relocate Column Next to Input
message("Merging results...")
setnames(mapping_table, "input_snp", col_name)
final_dt <- merge(dt, mapping_table, by = col_name, all.x = TRUE)

orig_idx <- which(colnames(final_dt) == col_name)
cols <- colnames(final_dt)
new_cols <- c(cols[1:orig_idx], "converted_snp", cols[(orig_idx+1):(length(cols)-1)])
final_dt <- final_dt[, ..new_cols]

# 6. Write to Disk
message("Writing output to: ", output_path)
fwrite(final_dt, output_path, sep = "\t", quote = FALSE)
message("Done successfully!")
