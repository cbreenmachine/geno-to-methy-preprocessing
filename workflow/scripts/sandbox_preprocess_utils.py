import pandas as pd
import dask
import dask.dataframe as dd
import numpy as np
import os
from fastparquet import write

# os.chdir("workflow/scripts")

seq_encoding = pd.DataFrame(
    {
    'residue': ['A', 'C', 'G', 'T'],
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1]
    }
)


# Function to read and preprocess VCF file (variant calls)
def read_vcf(ifile, chrom):
    '''
    Reads and extracts relevant information from VCF file
    '''
    # ifile = "../../data/variant-calls/253.snps.vcf"
    df = dd.read_csv(ifile, 
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
     

def encode_snps(row):
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
        # output = (encode_one_letter(minor2[0]) + encode_one_letter(minor2[1])) / 4 + encode_one_letter(alt) / 2
        output = (encode_one_letter(non_alt[0]) + encode_one_letter(non_alt[1]) + encode_one_letter(non_alt[2])) / 3
    return output



def read_sequence_bed(ifile, min_reads=5):
    '''Given a BED file (sample / chrom specific) with the following columns
    chrom start end strand methylated total sequence
    chr1 0 1000 + 10 30 AAAAAAAAAAA...    
    
    '''
    df = dd.read_csv(ifile, 
                     sep='\t', 
                     header=None, 
                     usecols=[0, 1, 2, 4, 5, 6], 
                     names=['chrom', 'start', 'end', 'methylated_counts', 'total_counts', 'sequence'])
    
    df = df.dropna()

    # Filter low number of reads
    df = df.loc[df['total_counts'] >= min_reads]

    # Create a locus identifier
    df['locus'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)

    # Turn AAAAA into five columns of A, A, A,  ...
    ex = df['sequence'].str.split('',expand=True, n=1000, )
    df = df.drop(columns=['sequence'])

    # str_split adds an empty column at beginning, drop it
    ex = ex.iloc[:, 1:]

    my_cols = ["residue"+str(i) for i in range(ex.shape[1])]
    ex.columns = my_cols

    # Combine 
    # Safe to ignore warning since data is derived row-by-row above.
    # In other words, it will match exactly
    df = dd.concat([df, ex], axis = 1)
    df = dd.reshape.melt(df, 
                        id_vars = ['locus', 'chrom', 'start', 'end', 'methylated_counts', 'total_counts', ], 
                        value_vars = my_cols, var_name = "residue_ix", value_name = "residue")

    # Conver 'residue1' --> 1 (int)
    df['residue_ix'] = df['residue_ix'].str.strip('residue').astype(int)
    df['pos'] = df['start'] + df['residue_ix'] + 1
    return df


seq_df = read_sequence_bed("../../data/bed-intervals/253.chr17.sequence.bed")
snp_df = read_vcf("../../data/variant-calls/253.snps.vcf")

# result = snp_df.apply(encode_snps, axis=1, meta = {0: float, 1:float, 2:float, 3:float}, result_type='expand')
# result.columns = ['A', 'C', 'G', 'T']

# snp_df = snp_df.merge(result, left_index=True, right_index=True)



# zz = dd.merge(seq_df, snp_df, on = ['chrom', 'pos'], indicator = True)

# write('2023-07-31-data.parq', zz)
