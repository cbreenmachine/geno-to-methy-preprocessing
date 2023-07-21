# encode_sequence_with_alleles.R
# Given a bed file with sequence (chrom start end AAAAAAAAACCCCCCTTTTGGGGG...)
# and a SNP call file (VCF), generate 
# Outputs four BED files,

library(argparse)
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(parallel)

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


#--> Load data, this is in bed format
ref.dt <- fread(args$reference, nrows = 100) %>% 
  dplyr::select(-c(V4, V5, V6)) %>% 
  dplyr::rename(chrom = "V1", start = "V2", end = "V3") %>% 
  dplyr::mutate(locus = paste0(chrom, ":", start, "-", end)) %>% 
  separate(col = V7, sep = "", into = paste0("seq", 0:1000)) %>% 
  pivot_longer(cols = starts_with("seq"), names_to = "seqix", values_to = "reference") %>% 
  dplyr::mutate(start = start + as.numeric(str_remove(seqix, "seq"))) %>% 
  left_join(sequence_mapping, by = "reference") %>% 
  dplyr::select(chrom, start, locus, reference, A, C, G, `T`) %>% 
  drop_na() %>% # empty string returned by separate gives a nonsense row for every chunk
  dplyr::mutate(alternate = reference, variant.call = "0/0")


tmp <- VariantAnnotation::readVcf(args$variants)

# Cast to BED standard (0-based, half open)
snps.dt <- 
  geno(tmp)$GT %>% 
  as.data.frame() %>% 
  rownames_to_column("locus") %>% 
  rename("variant.call" = 2) %>%
  separate(col = locus, into = c("chrom", "start", "reference", "alternate")) %>% 
  dplyr::mutate(start = as.numeric(start) - 1) %>% 
  left_join(allele.mapping, by = c("reference", "alternate", "variant.call")) %>% 
  full_join(snps.dt, dplyr::select(ref.dt, chrom, start, locus), by = c("chrom", "start"))

gc()

# May not be able to do this with data.tables after all. How many swaps do we need to make?
# Us GRanges list where each list item is a locus and each 

#--> Try to correct missing values of ref
dim(ref.dt)
dim(snps.dt)


merged <- anti_join(snps.dt, ref.dt, by = c("chrom", "start", "reference")) %>% 
  bind_rows(ref.dt) 

# What happens when a position is in more than one locus? How to impute that?
head(merged)
merged %>% 
  dplyr::filter(variant.call != "0/0")

# merged <- left_join(ref.dt, snps.dt, by = c("chrom", "start", "reference"))
# drop_na(merged) %>% nrow()

#--> Check that centered on CpG
# tmp <- as.character(ref.gr[1, "V7"])
# L <- nchar(tmp)
# if (substr(tmp, L/2, L/2 + 1) == "CG"){print("Middle characters are CG--good")}

# 
# # One-hot encode ----------------------------------------------------------
# 
# one_hot_encode <- function(x){
#   
#   if (x == "A"){
#     c(1, 0, 0, 0)
#   } else if (x == "C"){
#     c(0, 1, 0, 0)
#   } else if (x == "G"){
#     c(0, 0, 1, 0)
#   } else if (x == "T"){
#     c(0, 0, 0, 1)
#   } else {
#     c(0, 0, 0, 0)
#   }
# }
# 
# # one_hot_encode_v <- Vectorize(one_hot_encode)
# 
# 
# encode_one_chunk <- function(str){
#   v <- base::strsplit(tmp, "")[[1]]
#   as.matrix(do.call(rbind, lapply(X = v, FUN = one_hot_encode)))
# }
# 
# encode_string <- Vectorize(encode_one_chunk)
# 
# 
# ref.gr %>% 
#   dplyr::slice(1:100) %>% 
#   dplyr::group_by(V1, V2) %>% 
#   dplyr::summarize(ohe = encode_string(V7))
# 
# out <- lapply(X=data, FUN=encode_one_chunk)
# str(out)
# 
