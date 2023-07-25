# encode_sequence_with_alleles.R
# Given a bed file with sequence (chrom start end AAAAAAAAACCCCCCTTTTGGGGG...)
# and a SNP call file (VCF), generates hdf5
# sample = 253
#   locus = chr1:100-1100
#     encoding = [[], [], [], []]
#     methylated = 10
#     coverage = 20

library(argparse)
library(VariantAnnotation)
library(tidyverse)
library(tidyfast)
library(GenomicRanges)
library(data.table)
library(chunkR)
library(rhdf5)


# Constants ---------------------------------------------------------------

# Simple one-hot encoding
sequence_mapping <- data.frame(
  reference = c("A", "C", "G", "T"),
  A = c(1, 0, 0, 0),
  C = c(0, 1, 0, 0),
  G = c(0, 0, 1, 0),
  `T` = c(0, 0, 0, 1)
)

# Temporary scaffolding to cover ref/alt combinations
template <- data.frame(
  reference = rep(c("A", "C", "G", "T"), times = 4),
  alternate = rep(c("A", "C", "G", "T"),each = 4)
) %>% 
  dplyr::filter(reference != alternate) # Don't care about homo ref

# Reference / Alt case, code both as 0.5
case1 <- template %>%
  dplyr::mutate(variant.call = "0/1") %>%
  dplyr::mutate(
    A = ifelse(reference == "A" | alternate == "A", 0.5, 0),
    C = ifelse(reference == "C" | alternate == "C", 0.5, 0),
    G = ifelse(reference == "G" | alternate == "G", 0.5, 0),
    `T` = ifelse(reference == "T" | alternate == "T", 0.5, 0)
  )

# Alt/Alt case, code alt as 1, everything else as 0
case2 <- template %>%
  dplyr::mutate(variant.call = "1/1") %>%
  dplyr::mutate(
    A = ifelse(alternate == "A", 1, 0),
    C = ifelse(alternate == "C", 1, 0),
    G = ifelse(alternate == "G", 1, 0),
    `T` = ifelse(alternate == "T", 1, 0)
  )

# Funky case (code alt as 0.5, other non-ref as 0.25)
case3 <- template %>%
  dplyr::mutate(variant.call = "1/2") %>%
  dplyr::mutate(
    A = ifelse(alternate == "A", 0.5, ifelse(reference != "A", 0.25, 0)),
    C = ifelse(alternate == "C", 0.5, ifelse(reference != "C", 0.25, 0)),
    G = ifelse(alternate == "G", 0.5, ifelse(reference != "G", 0.25, 0)),
    `T` = ifelse(alternate == "T", 0.5, ifelse(reference != "T", 0.25, 0))
  )

# Join together for merging later
allele_mapping <- rbind(case1, case2, case3)

# Load in SNPs function --------------------------------------------------

# SNPs will be used on each iteration of chunkR
munge_variant_calls <- function(dt) {
  # Takes VCF (read in thru VariantAnnotation::readVcf)
  # and renames columns, adds some derived information,
  # casts to (0, 1] coordinate system
  as.data.frame(geno(dt)$GT) %>%
    rownames_to_column("locus") %>%
    dplyr::rename("variant.call" = 2) %>%
    dplyr::filter(!str_detect(locus, "chrUn"),
                  !str_detect(locus, "random")) %>%
    tidyr::separate(col = locus,
             into = c("chrom", "pos", "reference", "alternate")) %>%
    dplyr::filter(reference %in% c("A", "C", "G", "T")) %>%
    dplyr::mutate(pos = as.numeric(pos)) %>%
    left_join(allele_mapping, by = c("reference", "alternate", "variant.call"))
}




munge_reference_data <- function(dt) {
  # Convert to 1-based coordinate system
  colnames(dt) <- paste0("V", 1:ncol(dt))

  # Calculate interval width
  width <- as.numeric(dt[1, "V3"] - dt[1, "V2"])
  dummy_cols <- paste0("seq", 1:width)

  dt %>%
    # Column names
    dplyr::rename(chrom = "V1", start = "V2", end = "V3",
                  strand = "V4", x = "V5", n = "V6") %>% 
    # Keep track of locus, since these will be written out to final BED files
    dplyr::mutate(locus = paste0(chrom, ":", start, "-", end)) %>%
    # Split the string of letters (AACCTTGG) into however many variables we need
    # It's ok that this produces NAs--there are unevenly sized intervals,
    #which results in empty columns
    dt_separate(col = V7, sep = "", into = dummy_cols) %>%
    tidyr::pivot_longer(cols = tidyr::all_of(dummy_cols),
                        names_to = "seqix",
                        values_to = "reference") %>%
    dplyr::mutate(pos = start +
                  as.numeric(stringr::str_remove(seqix, "seq")) - 1 + 1) %>%
    dplyr::left_join(sequence_mapping, by = "reference") %>%
    dplyr::select(chrom, pos, locus, reference, A, C, G, `T`) %>%
    # empty string returned by separate gives a nonsense row for every chunk
    tidyr::drop_na() %>%
    dplyr::mutate(alternate = NA, variant.call = "0/0")
}

swap_in_snps <- function(ref_dt, snps_dt) {

  # Because we will do this in chunks, we only want to
  # perform substitutions
  ref_dt$is_snp <- FALSE
  snps_dt$is_snp <- TRUE

  snps_in_range_dt <- semi_join(snps_dt, ref_dt,
                                by = c("chrom", "pos", "reference"))

  out <- left_join(ref_dt, snps_dt,
                   by = c("chrom", "pos", "reference")) %>%
         mutate(is_snp = replace_na(is_snp.y, FALSE),
         A = NA, C = NA, G = NA, `T` = NA, variant.call = NA,
         alternate = alternate.y)

  # For rows where there is a SNP (called)
  # we will keep the variant call information. Otherwise,
  # we will keep the reference information
  out[out$is_snp, ] %<>% dplyr::mutate(
                    A = A.y, C = C.y, G = G.y, `T` = `T.y`,
                    variant.call = variant.call.y
                  )

    # Non SNPs
    out[!out$is_snp, ] %<>% dplyr::mutate(
                    A = A.x, C = C.x, G = G.x, `T` = `T.x`,
                    variant.call = variant.call.x
                  )

    swapped <- out %>% 
      transmute(
        chrom, pos, locus, reference, alternate, variant.call,
        A, C, G, `T`)

  return(swapped)
}

write_data_as_bed <- function(data, ofile) {
  if (file.exists(ofile)) {
    fwrite(data, ofile, append = TRUE, sep = "\t")
  } else {
    fwrite(data, ofile, append = FALSE, sep = "\t")
  }
}

# Run the pipeline (functionalzied) -----------------------

parser <- ArgumentParser()

parser$add_argument("--variants", default = "../../data/variant-calls/253.snps.vcf")
parser$add_argument("--reference", default = "../../data/bed-intervals/reference.chunked.bed")
parser$add_argument("--ofile", default = "../../data/training/253.chr1.bed")

args <- parser$parse_args()

# Variant data is derived from a VCF file, all of it can be read in
snps_dt <- munge_variant_calls(VariantAnnotation::readVcf(args$variants))
# snps_dt <- munge_variant_calls(snps_vcf)

# Chunker obkect allows us to iterate more slowly
chunker_obj <- chunker(args$reference, sep = "\t",
                       has_colnames = FALSE,
                       has_rownames = FALSE,
                       chunksize = 10000L)

# Initialize the line counter
lines <- 0

# Important to delete old files
if (file.exists(args$ofile)) {
  unlink(args$ofile)
  cat("Overwriting output file\n")
}

while (next_chunk(chunker_obj)) {
  ref_dt <- munge_reference_data(get_table(chunker_obj))
  swapped_dt <- swap_in_snps(ref_dt, snps_dt)
    
  # Make the data frame "wide" so that we can keep using bed tools operations on it
  swapped_for_bed <- swapped_dt %>%
    group_by(locus) %>%
    summarize(A_encoding = paste(A, collapse=","),
              C_encoding = paste(C, collapse=","),
              G_encoding = paste(G, collapse=","),
              T_encoding = paste(`T`, collapse=",")) %>%
    tidyr::separate(locus, into=c("chrom", "start", "end"))

  write_data_as_bed(swapped_for_bed, args$ofile)

  lines <- lines + nrow(ref_dt)
  cat("Processed ", lines, "lines\n")
}
#END