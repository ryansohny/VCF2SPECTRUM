#!/usr/bin/env python3

import pandas as pd
import os

def reverse_complement(string):
    try:
        complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        complement_string = "".join([complement_dict[s] for s in string])
    except KeyError:
        raise ValueError("Invalid character other than A,T,G and C")
    return complement_string[::-1]

def load_kmer_normalization_factors(kmer_norm_file):
    """Load k-mer normalization factors from a tab-delimited file"""
    print(f"Loading k-mer normalization factors from {kmer_norm_file}")
    
    # Read the tab-delimited file
    df_norm = pd.read_csv(kmer_norm_file, sep='\t', header=None, names=['kmer', 'norm_factor'])
    
    # Create normalization dictionary for both orientations
    kmer_norm_dict = {}
    
    for _, row in df_norm.iterrows():
        kmer = row['kmer']
        norm_factor = row['norm_factor']
        
        # Store the k-mer as-is
        kmer_norm_dict[kmer] = norm_factor
        
        # Also store the reverse complement if the central base is not C or T
        if kmer[1] not in ["C", "T"]:
            rev_comp_kmer = reverse_complement(kmer)
            kmer_norm_dict[rev_comp_kmer] = norm_factor
    
    return kmer_norm_dict

def extract_trinucleotide(mutation_type):
    """Extract trinucleotide context from SigProfiler SBS96 type (e.g., A[C>A]A -> ACA)"""
    return mutation_type[0] + mutation_type[2] + mutation_type[6]

def main():
    sbs96_file = snakemake.input.sbs96
    kmer_norm_file = snakemake.input.kmer_norm
    output_file = snakemake.output.normalized
    output_dir = snakemake.params.output_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading k-mer normalization factors from {kmer_norm_file}")
    norm_factors = load_kmer_normalization_factors(kmer_norm_file)
    
    print(f"Loading SBS96 data from {sbs96_file}")
    df_sbs96 = pd.read_csv(sbs96_file, sep='\t')
    
    # Calculate original total count
    original_total = df_sbs96.iloc[:, 1].sum()
    print(f"Original total mutation count: {original_total}")
    
    # Apply k-mer normalization
    print("Applying k-mer normalization...")
    normalized_counts = []
    
    for _, row in df_sbs96.iterrows():
        mutation_type = row['MutationType']
        original_count = row.iloc[1]  # Second column contains the count
        
        # Extract trinucleotide context
        trinucleotide = extract_trinucleotide(mutation_type)
        
        # Get normalization factor
        if trinucleotide in norm_factors:
            norm_factor = norm_factors[trinucleotide]
            # Apply normalization by multiplication
            normalized_count = original_count * norm_factor
        else:
            print(f"Warning: Trinucleotide {trinucleotide} not found in normalization data, keeping original count")
            normalized_count = original_count
            
        normalized_counts.append(normalized_count)
    
    # Calculate the scaling factor to preserve total count
    pre_scale_total = sum(normalized_counts)
    if pre_scale_total > 0:
        scaling_factor = original_total / pre_scale_total
        print(f"Pre-scaling total: {pre_scale_total:.2f}")
        print(f"Scaling factor to preserve total: {scaling_factor:.6f}")
        
        # Apply scaling and round to integers
        final_counts = [round(count * scaling_factor) for count in normalized_counts]
    else:
        final_counts = normalized_counts
    
    # Calculate final total count
    final_total = sum(final_counts)
    print(f"Final total mutation count: {final_total}")
    print(f"Count preservation ratio: {final_total/original_total:.6f}")
    
    # Create output dataframe
    df_normalized = df_sbs96.copy()
    df_normalized.iloc[:, 1] = final_counts
    
    # Save normalized data
    print(f"Saving k-mer normalized data to {output_file}")
    df_normalized.to_csv(output_file, sep='\t', index=False)
    
    print("K-mer normalization completed successfully!")
    print(f"Total mutations: {original_total} → {final_total} (ratio: {final_total/original_total:.6f})")

if __name__ == "__main__":
    main()