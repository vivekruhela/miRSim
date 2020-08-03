#!/usr/bin/python

def update_gff(gff_df,out,ground_truth_filename):
    rna_truth = pd.read_csv(os.path.join(out,ground_truth_filename),sep = '\t')
    rna_truth = rna_truth.set_index('RNA_ID')
    for idx in range(rna_truth.shape[0]):
        chr_start = str(rna_truth.iloc[idx,4])
        chr_end = str(rna_truth.iloc[idx,5])
        gff_df = gff_df[~gff_df.isin([chr_start,chr_end])]
        gff_df = gff_df.dropna()
    return gff_df