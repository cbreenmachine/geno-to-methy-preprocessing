import pandas as pd
import numpy as np
from Bio import SeqIO
# import os
import argparse

print("STARTED HERE")


# Parse parse parse
parser = argparse.ArgumentParser()
parser.add_argument('--ifile') #default = "../../data/reference/by-chrom/chr19.fa")
parser.add_argument('--snp_file') #, default = "../../data/variant-calls/253.snps.vcf")
parser.add_argument('--ofile') #, default = "../../data/encodings/253.chr19.bed")
args = parser.parse_args()




def read_and_filter_vcf(ifile, chrom):
    '''
    Reads and extracts relevant information from VCF file
    '''
    # ifile = "../../data/variant-calls/253.snps.vcf"
    df = pd.read_csv(ifile, 
                     sep='\t', 
                     comment='#', 
                     header=None, 
                     usecols=[0, 1, 3, 4, 9],
                     names=['chrom', 'pos', 'reference', 'alternate', 'extra_info']
                     )
    
    # Filter to shrink the number of comparisons that need to be made
    df = df.loc[df['chrom'] == chrom]
    
    # Pull out the SNP call
    df['variant_call'] = df['extra_info'].str[:3]

    # Get rid of Ns, indicate that ref homozygous
    df = df[df['reference'].isin(['A', 'C', 'G', 'T'])]

    # Makes iteration work
    df = df.drop(columns=['extra_info'])
    df = df.reset_index(drop=True)

    return df




def encode_one_letter(x):
    '''Given A, C, G, T return one-hot encoding'''
    out = np.zeros((4,),dtype = 'float32')

    if x == "A":
        out[0] = 1
    elif x == "C":
        out[1] = 1
    elif x == "G":
        out[2] = 1
    elif x == "T":
        out[3] = 1
    
    return(out)


def encode_snp(row):
    '''
    Given chrom	pos	reference	alternate	variant_call
    encode with fractionals
    '''
    ref, alt, vc = row['reference'], row['alternate'], row['variant_call']

    if vc == "0/1":
        output = (encode_one_letter(ref) + encode_one_letter(alt)) / 2
    elif vc == "1/1":
        output = encode_one_letter(alt)
    elif vc == "1/2":
        alt_split = alt.split(alt)
        non_alt = np.setdiff1d(['A', 'C', 'G', 'T'], [ref])
        output = (encode_one_letter(non_alt[0]) + encode_one_letter(non_alt[1]) + encode_one_letter(non_alt[2])) / 3
    return output


def run_one_hot_encoder(sequence, snp_df):

    l = len(sequence)
    x = np.zeros((l, 4),dtype = 'float32')

    # Remember that i starts at zero, 
    # whereas the positions in the VCF start at 1
    for i, nt in enumerate(sequence):
        # For position
        p = i + 1

        if p in snp_df.index:
            x[i, :] = encode_snp(snp_df.loc[p])
        else:
            x[i, :] = encode_one_letter(nt)
    return x



def main(args):

    print(args)

    # Derive the chromosome 
    chrom = args.ifile.split("/")[-1].split(".")[0]
    

    # Read sequence data (from reference fasta)
    # This is a dictionary with one key. May un-dict this later,
    # but it's safe this way
    seq_data = SeqIO.to_dict(SeqIO.parse(open(args.ifile), 'fasta'))
    # print(seq_data)

    # Read (and filter) called SNPs from VCF
    snp_df = read_and_filter_vcf(args.snp_file, chrom)
    # print(snp_df.head())

    # Workhorse of the script
    out = run_one_hot_encoder(seq_data[chrom], snp_df)

    start = [x for x in range(out.shape[0])]
    out = pd.DataFrame(out, columns=['A','C', 'G', 'T'])

    out['chrom'] = chrom
    out['start'] = start
    out['end'] = out['start'] + 1

    out = out[['chrom', 'start', 'end', 'A', 'C', 'G', 'T']]

    out.to_csv(args.ofile, sep="\t", header=False, index=False)


main(args)