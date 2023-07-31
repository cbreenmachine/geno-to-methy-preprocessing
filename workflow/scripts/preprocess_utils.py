import pandas as pd
import numpy as np
import os

# os.chdir("workflow/scripts")

# Function to read and preprocess VCF file (variant calls)
def read_vcf(ifile):
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
    
    # Pull out the SNP call
    df['variant_call'] = df['extra_info'].str[:3]

    # Get rid of Ns, indicate that ref homozygous
    df = df[df['reference'].isin(['A', 'C', 'G', 'T'])]

    # Makes iteration work
    df.drop(columns=['extra_info'])
    df.reset_index(drop=True, inplace=True)

    return df

def encode_one_letter(x):

    out = np.zeros((1,4),dtype = 'int8')

    if x == "A":
        out[0, 0] = 1
    elif x == "C":
        out[0, 1] = 1
    elif x == "G":
        out[0, 2] = 1
    elif x == "T":
        out[0, 3] = 1
    
    return(out)
        
        
seq_to_encoding = {
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1]
}



def encode_snps(data):
    '''
    Given chrom	pos	reference	alternate	variant_call
    encode with fractionals
    '''

    output = np.zeros(shape=[data.shape[0], 4])

    for i, row in data.iterrows():
        ref, alt, vc = row['reference'], row['alternate'], row['variant_call']

        if vc == "0/1":
            output[i, :] = (encode_one_letter(ref) + encode_one_letter(alt)) / 2
        elif vc == "1/1":
            output[i, :] = encode_one_letter(alt)
        elif vc == "1/2":
            minor2 = np.setdiff1d(['A', 'C', 'G', 'T'], [ref, alt])
            output[i, :] = (encode_one_letter(minor2[0]) + encode_one_letter(minor2[1])) / 4 + encode_one_letter(alt) / 2

    return output


# Encode reference and swap in SNPs on the fly?
# May be slow...


def encode_reference(data):
    encoding = np.zeros(shape=[data.shape[0], 4, 1000])

    for i, row in data.iterrows():

        for j, seq in enumerate(row['sequence']):
            encoding[i, :, j] = encode_one_letter(seq)
    
    # encoding_df = pd.DataFrame(encoding, columns = ['A', 'C', 'G', 'T'])
    # pd.concat([data, encoding], axis=1)
    return(encoding)



def read_sequence_bed(ifile, min_reads=5):
    '''Given a BED file (sample / chrom specific) with the following columns
    chrom start end strand methylated total sequence
    chr1 0 1000 + 10 30 AAAAAAAAAAA...    
    
    '''
    df = pd.read_csv(ifile, 
                     sep='\t', 
                     header=None, 
                     usecols=[0, 1, 2, 4, 5, 6], 
                     names=['chrom', 'start', 'end', 'methylated_counts', 'total_counts', 'sequence'])
    
    df.dropna(inplace = True)

    # Filter low number of reads
    df = df.loc[df['total_counts'] >= 5]

    # Create a locus identifier
    df['locus'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)

    # # Turn AAAAA into five columns of A, A, A,  ...
    ex = pd.DataFrame(df['sequence'].str.split('',expand=True))
    
    # str_split adds an empty column at beginning and end, drop these
    ex = ex.iloc[:, 1:-1]

    my_cols = ["x"+str(i) for i in range(ex.shape[1])]
    ex.columns = my_cols

    # Combine 
    df = pd.concat([df, ex], axis = 1)

    df.wide_to_long(my_cols, stubnames = "x")
    df.reset_index(drop=True, inplace=True)

    return df


snp_df = read_vcf("../../data/variant-calls/253.snps.vcf")
seq_df = read_sequence_bed("../../data/bed-intervals/253.chr17.sequence.bed")

# Want to return reference as well. This will be storage-wise inefficient but makes clear the
# parallel nature of data

# Then iterate through list with zip

# How to efficiently make SNP swap?