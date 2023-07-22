# encode_sequence_with_alleles.R
# Given a bed file with sequence (chrom start end AAAAAAAAACCCCCCTTTTGGGGG...)
# and a SNP call file (VCF), generate 
# Outputs four BED files,

library(argparse)
library(VariantAnnotation)
library(tidyverse)
library(tidyfast)
library(GenomicRanges)
library(data.table)
library(chunkR)


parser <- ArgumentParser()

# by default ArgumentParser will add an help option 
parser$add_argument("--reference", default = "../../data/bed-intervals/reference.chunked.bed")
parser$add_argument("--variants", default = "../../data/variant-calls/253.snps.vcf")
parser$add_argument("--output_prefix", default = "../../data/encoded-data/encoded.sample.")

args <- parser$parse_args()


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
allele.mapping <- rbind(case1, case2, case3)

# Load in SNPs ------------------------------------------------------------
tmp <- VariantAnnotation::readVcf(args$variants)

# Cast to BED standard (0-based, half open)
snps.dt <- 
  geno(tmp)$GT %>% 
  as.data.frame() %>% 
  rownames_to_column("locus") %>% 
  rename("variant.call" = 2) %>%
  separate(col = locus, into = c("chrom", "start", "reference", "alternate")) %>% 
  dplyr::filter(chrom %in% unique(ref.dt$chrom),
                reference %in% c("A", "C", "G", "T")) %>% 
  dplyr::mutate(start = as.numeric(start) - 1,
                end = start + 1) %>% 
  left_join(allele.mapping, by = c("reference", "alternate", "variant.call"))

gc()


# Create function to process reference data (will be read in chunk --------

# Keep everything in BED standard, changes happen to VCF

#--> Load data, this is in bed format (coordinate change is corrected implicitly with
# separate ranging from 0). Keep in BED format (0, 1]
# chr1:12579:12580, but by the end of plumbing, it is chr1:12579
ref.raw.dt <- fread(args$reference)

cast_reference_data_to_snp_format <- function(dt){
  
  # Calculate interval width
  width <- as.numeric(dt[1, "V3"] - dt[1, "V2"])
  dummy.cols <- paste0("seq", 1:width)
  
  dt %>% 
    # Column names
    dplyr::rename(chrom = "V1", start = "V2", end = "V3",
                  strand = "V4", x = "V5", n = "V6") %>% 
    # Keep track of locus, since these will be written out to final BED files
    dplyr::mutate(locus = paste0(chrom, ":", start, "-", end)) %>% 
    # Split the string of letters (AACCTTGG) into however many variables we need
    # It's ok that this produces NAs--there are unevenly sized intervals, which results in 
    # empty columns
    dt_separate(col = V7, sep = "", into = dummy.cols) %>% 
    pivot_longer(cols = all_of(dummy.cols), names_to = "seqix", values_to = "reference") %>% 
    dplyr::mutate(start = start + as.numeric(str_remove(seqix, "seq")) - 1,
                  end = start + 1) %>% 
    left_join(sequence_mapping, by = "reference") %>% 
    dplyr::select(chrom, start, end, locus, reference, A, C, G, `T`) %>% 
    drop_na() %>% # empty string returned by separate gives a nonsense row for every chunk
    dplyr::mutate(alternate = reference, variant.call = "0/0")
}

cast_reference_data_to_snp_format(ref.raw.dt[1:1000, ])

# seq.lengths <- unlist(lapply(ref.raw.dt$V7, str_length))
# summary(seq.lengths)

# Perform the pipeline to pieces at a time.
dim(ref.dt)
 

# May not be able to do this with data.tables after all. How many swaps do we need to make?
# Us GRanges list where each list item is a locus and each 
dim(ref.dt)
dim(snps.dt)

locus.dt <- dplyr::select(ref.dt, c(chrom, start, end, locus))
snps.w.locus.dt <- left_join(snps.dt, locus.dt)


# Try to correct missing values of ref ----------------------------------------------------
merged <- anti_join(ref.dt, snps.w.locus.dt, by = c("chrom", "start", "end", "reference")) %>% 
  bind_rows(snps.w.locus.dt) %>% 
  dplyr::arrange(chrom, start) %>% 
  distinct() 


merged <- bind_rows(snps.w.locus.dt, ref.dt) %>% 
  group_by(chrom, start) %>% 
  dplyr::filter(n() == 1 | variant.call != "0/0")

# What happens when a position is in more than one locus? How to impute that?
head(merged)

# To check whether this works, we need to look at non homozygous ref
merged %>% 
  dplyr::filter(variant.call != "0/0")

my.locus <- unique(merged$locus)[2]
merged %>% 
  dplyr::filter(locus == my.locus,
                variant.call != "0/0") %>% 
  distinct()

nrow(merged)
nrow(ref.dt)
nrow(snps.dt)
