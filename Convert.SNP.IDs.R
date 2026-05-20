########################################################################################################
#                                                                                                      #
# Using biomaRt and dbSNP to convert SNP IDs                                                           #
#                                                                                                      #
# Input arguments:                                                                                     #
#                 Input.File: The file contains source IDs (rsid or chr:pos:Ref:Alt format)            #
#                 Out.Prefix: Output files prefix (Directories will be created if it contains path)    #
#                 ID.col: the source ID column name or index in the input file                         #
#                 Genom.Ver: Genome version to use for conversion (GRCh37" or "GRCh38")                #
#                 Source.ID: Source ID type ("rsid" or "chr:pos")                                      #
#                                                                                                      #
########################################################################################################

message("Loading reqiuered libraries...")
library(biomaRt)
library(data.table)
arguments <- commandArgs(T)

message("Reading input data...")
Input.File <- trimws(arguments[1])
Out.Prefix <- trimws(arguments[2])
ID.col <- trimws(arguments[3])
Genome.Ver <- trimws(arguments[4])
Sourc.ID <- trimws(arguments[5])

Out.dir <- dirname(Out.Prefix)
if(!dir.exists(Out.dir)){
  dir.create(Out.dir , recursive = T)
}

Genome.Ver <- match.arg(Genome.Ver , choices = c("GRCh37","GRCh38"))
Sourc.ID <- match.arg(Sourc.ID , choices = c("rsid" , "chr:pos"))

snp_ids = fread(file = Input.File , header = F , stringsAsFactors = F,data.table = F)

ID.col.num <- as.numeric(ID.col)

if(is.na(ID.col.num)){
  if(ID.col %in% names(snp_ids)){
    ID.col.num <- match(ID.col,names(snp_ids))
  }else{
    stop("Can't find column name ",ID.col, " in the input file!")
  }
}

message("Connecting to BioMart database...")
if(Genome.Ver=="GRCh37"){
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp",host = "https://grch37.ensembl.org", verbose = T)
}else if(Genome.Ver=="GRCh38"){
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", verbose = T)
}
snp_attributes = c("refsnp_id", "chr_name", "chrom_start","chrom_end","allele","allele_1")

message("Converting IDs...")
if(Sourc.ID == "rsid" ){
  
  snp_filters = c("snp_filter")
  results = getBM(attributes=snp_attributes, filters=snp_filters, values=snp_ids[,ID.col.num], mart=snp_mart)
  
} else if(Sourc.ID == "chr:pos"){
  
  snp_filters = c("chr_name","start","end")
  results = getBM(attributes=snp_attributes, filters=snp_filters, values=snp_ids[,ID.col.num], mart=snp_mart)
  
}

message("Writing the results...")
fwrite(results, file = paste0(Out.Prefix , "SNP.ID.txt"), row.names = F , sep = "\t",quote = F)
