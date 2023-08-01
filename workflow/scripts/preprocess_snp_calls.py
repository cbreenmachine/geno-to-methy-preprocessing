
import preprocess_utils

# Constants ---------------------------------------------------------------

# Simple one-hot encoding
sequence_mapping = pd.DataFrame({
    'reference': ['A', 'C', 'G', 'T'],
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1]
})

# Temporary scaffolding to cover ref/alt combinations
template = pd.DataFrame({
    'reference': np.repeat(['A', 'C', 'G', 'T'], 4),
    'alternate': np.tile(['A', 'C', 'G', 'T'], 4)
})
template = template[template['reference'] != template['alternate']]  # Don't care about homo ref

# Reference / Alt case, code both as 0.5
case1 = template.copy()
case1['variant.call'] = '0/1'
case1['A'] = np.where((case1['reference'] == 'A') | (case1['alternate'] == 'A'), 0.5, 0)
case1['C'] = np.where((case1['reference'] == 'C') | (case1['alternate'] == 'C'), 0.5, 0)
case1['G'] = np.where((case1['reference'] == 'G') | (case1['alternate'] == 'G'), 0.5, 0)
case1['T'] = np.where((case1['reference'] == 'T') | (case1['alternate'] == 'T'), 0.5, 0)

# Alt/Alt case, code alt as 1, everything else as 0
case2 = template.copy()
case2['variant.call'] = '1/1'
case2['A'] = np.where(case2['alternate'] == 'A', 1, 0)
case2['C'] = np.where(case2['alternate'] == 'C', 1, 0)
case2['G'] = np.where(case2['alternate'] == 'G', 1, 0)
case2['T'] = np.where(case2['alternate'] == 'T', 1, 0)

# Funky case (code alt as 0.5, other non-ref as 0.25)
case3 = template.copy()
case3['variant.call'] = '1/2'
case3['A'] = np.where(case3['alternate'] == 'A', 0.5, np.where(case3['reference'] != 'A', 0.25, 0))
case3['C'] = np.where(case3['alternate'] == 'C', 0.5, np.where(case3['reference'] != 'C', 0.25, 0))
case3['G'] = np.where(case3['alternate'] == 'G', 0.5, np.where(case3['reference'] != 'G', 0.25, 0))
case3['T'] = np.where(case3['alternate'] == 'T', 0.5, np.where(case3['reference'] != 'T', 0.25, 0))

# Join together for merging later
allele_mapping = pd.concat([case1, case2, case3], ignore_index=True)

# Functions --------------------------------------------------------------



# Function to preprocess reference data (bed file)
def munge_reference_data(reference_file):
    dt = pd.read_csv(reference_file, sep='\t', header=None, usecols=[0, 1, 2, 6], names=['chrom', 'start', 'end', 'sequence'])
    dt['locus'] = dt['chrom'] + ':' + dt['start'].astype(str) + '-' + dt['end'].astype(str)
    dt = dt.dropna(subset=['sequence'])
    dt = dt.explode('sequence')
    dt = dt.merge(sequence_mapping, left_on='sequence', right_on='reference')
    dt = dt[['chrom', 'start', 'end', 'locus', 'reference', 'A', 'C', 'G', 'T']]
    dt['alternate'] = np.nan
    dt['variant.call'] = '0/0'
    return dt

# Function to swap SNP information into reference data
def swap_in_snps(ref_dt, snps_dt):
    ref_dt['is_snp'] = False
    snps_dt['is_snp'] = True
    snps_in_range_dt = snps_dt.merge(ref_dt, on=['chrom', 'pos', 'reference'])
    out = ref_dt.merge(snps_dt, on=['chrom', 'pos', 'reference'], how='left', suffixes=('_x', '_y'))
    out['is_snp'] = out['is_snp_y'].fillna(False)
    out['A'] = out.apply(lambda row: row['A_y'] if row['is_snp'] else row['A_x'], axis=1)
    out['C'] = out.apply(lambda row: row['C_y'] if row['is_snp'] else row['C_x'], axis=1)
    out['G'] = out.apply(lambda row: row['G_y'] if row['is_snp'] else row['G_x'], axis=1)
    out['T'] = out.apply(lambda row: row['T_y'] if row['is_snp'] else row['T_x'], axis=1)
    out['variant.call'] = out.apply(lambda row: row['variant.call_y'] if row['is_snp'] else row['variant.call_x'], axis=1)
    out = out[['chrom', 'pos', 'end', 'locus', 'reference', 'alternate', 'variant.call', 'A', 'C', 'G', 'T']]
    return out


def read_methy_data(ifile):
    # ifile = "../../data/bed-intervals/253.chr17.bed"
    methy_data = pd.read_table(
        ifile,
        names = ['chrom', 'start', 'end', 'strand', 'methylated', 'total']
    )
    
    methy_data["locus"] = methy_data["chrom"] + ":" + methy_data["start"].astype(str) + "-" + + methy_data["end"].astype(str)
    
    methy_data = methy_data.astype({ 
        'methylated': 'float32',
        'total': 'float32'
    })
    
    return(methy_data)


def read_encoded_data(ifile):
    # ifile = "../../data/encodings/253.chr17.001.bed"
    encoded_data = pd.read_table(
        ifile, 
        names = ['locus', 'A', 'C', 'G', 'T'],
    )
    
    encoded_data = encoded_data.astype({
        'locus': 'object', 
        'A': 'float32',
        'C': 'float32',
        'G': 'float32',
        'T': 'float32'
    })

    
    return(encoded_data)

def merge_and_separate(methy_data, encoded_data):
    full_dataset = encoded_data.merge(methy_data, on = "locus")
    
    
    
    X = tf.reshape(
        tf.convert_to_tensor(full_dataset[['A', 'C', 'G', 'T']]), 
        shape = [-1, 4, 1000]
    )
    
    
    # Keep the matrices in the right order
    every_nth = list(range(full_dataset.shape[0]))[0::1000]
    y = tf.convert_to_tensor(full_dataset[['methylated', 'total']].iloc[every_nth])
    
    return X, y

# # Function to write data to TSV file
# def write_data_as_tsv(data, ofile):
#     data[['locus', 'A', 'C', 'G', 'T']].to_csv(ofile, sep='\t', header=False, index=False, mode

