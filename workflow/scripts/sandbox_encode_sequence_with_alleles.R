#!/usr/bin/env Rscript

# encode_sequence_with_alleles.R
# Given a bed file with sequence (chrom start end AAAAAAAAACCCCCCTTTTGGGGG...)
# and a SNP call file (VCF), generates hdf5
# sample = 253
#   locus = chr1:100-1100
#     encoding = [[], [], [], []]
#     methylated = 10
#     coverage = 20

suppressPackageStartupMessages({
  library(argparse)
  library(VariantAnnotation)
  library(tidyverse)
  library(tidyfast)
  library(GenomicRanges)
  library(data.table)
  library(chunkR)
  library(magrittr)
  library(rhdf5)
})



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


write_data_as_tsv <- function(data, ofile) {

  data2 <- data %>%
    dplyr::transmute(locus, A, C, G, `T`)

  if (file.exists(ofile)) {
    fwrite(data2, ofile, append = TRUE, sep = "\t", col.names = FALSE)
  } else {
    fwrite(data2, ofile, append = FALSE, sep = "\t", col.names = FALSE)
  }
}


construct_h5 <- function(ofile, s) {

  d <- dirname(ofile)

  if (!dir.exists(d)) {
    print("Creating directory")
    dir.create(d)
  }

  if (!file.exists(ofile)) {
    print("Constructing file")
    h5createFile(ofile)
  }

  existing_groups <- h5ls(ofile)$group

  if (!s %in% existing_groups) {
    print("Creating sample group")
    h5createGroup(ofile, sample)
  }
}

write_one_locus_to_h5 <- function(dt, s, l, ofile) {
  #  my_locus <- swapped_dt$locus[1]

  # h5createGroup(ofile, s)
  h5createGroup(ofile, paste0(s, "/", l))
  
  embedding <- dt %>%
    dplyr::filter(locus == l) %>%
    dplyr::select(c(A, C, G, `T`))

  h5write(embedding, ofile, paste0(s, "/", l, "/embedding"))
}



# Run the pipeline (functionalzied) -----------------------

parser <- ArgumentParser()

parser$add_argument("--variants", default = "../data/variant-calls/253.snps.vcf")
parser$add_argument("--reference", default = "../data/bed-intervals/253.chr17.sequence.bed")
parser$add_argument("--odir", default = "../data/training/")
# parser$add_argument("--overwrite", action = "store_true")

args <- parser$parse_args()

# Derived stuff
s <- str_extract(args$variants, "[0-9][0-9][0-9]")
chr <- str_extract(args$reference, "chr[0-9][0-9]")

# unlink(args$ofile)
# Prepare h5
# construct_h5(args$ofile, s)


suffix <- paste0(str_pad(0:500, 3, pad = "0"), ".bed")
file_names <- file.path(args$odir, paste0(s, ".", chr, ".", suffix))


# Variant data is derived from a VCF file, all of it can be read in
snps_dt <- munge_variant_calls(VariantAnnotation::readVcf(args$variants))

# Chunker obkect allows us to iterate more slowly
# Chr1 has ~2,000,000 CpGs, so up to ~200 output files
# chromosome
N_rows <- 10000L

chunker_obj <- chunker(args$reference, sep = "\t",
                       has_colnames = FALSE,
                       has_rownames = FALSE,
                       chunksize = N_rows)

# Initialize the line counter
ix <- 1

while (next_chunk(chunker_obj)) {
  ref_dt <- munge_reference_data(get_table(chunker_obj))
  swapped_dt <- swap_in_snps(ref_dt, snps_dt)

  if (nrow(swapped_dt) == N_rows * 1000) {
    write_data_as_tsv(swapped_dt, file_names[ix])
    ix <- ix + 1
  } else {
    print(nrow(swapped_dt))
  }
 

  # Get unique loci
  # ll_range <- unique(swapped_dt$locus)

  # This needs to be defined anew for every new version of `swapped_dt`
  # write_wrapper <- function(ll) {
  #   write_one_locus_to_h5(swapped_dt, s, ll, args$ofile)
  # }

  # lapply(ll_range, write_wrapper)

  # lines <- lines + nrow(ref_dt)
  # cat("Processed ", lines, "lines\n")
}

# h5closeAll()
#END