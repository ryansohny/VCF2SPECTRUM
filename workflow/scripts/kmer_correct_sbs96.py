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

def load_kmer_fractions(kmer_file):
    kmer3_canonical_dict = dict()
    
    with open(kmer_file, 'r') as dfh:
        c = 1
        for line in dfh:
            if c % 2 == 1:
                kmer_count = int(line.strip().lstrip(">"))
            elif c % 2 == 0:
                kmer3_canonical_dict[line.strip()] = kmer_count
            c += 1

    kmer3_canonical_sbs_dict = dict()
    for kmer, count in kmer3_canonical_dict.items():
        if kmer[1] not in ["C", "T"]:
            kmer3_canonical_sbs_dict[reverse_complement(kmer)] = count
        else:
            kmer3_canonical_sbs_dict[kmer] = count

    kmer3_canonical_sbs_dict = dict(sorted(kmer3_canonical_sbs_dict.items()))

    df_canonical_sbs96 = pd.DataFrame(list(kmer3_canonical_sbs_dict.items()), columns=['3mer', 'count'])
    df_canonical_sbs96 = df_canonical_sbs96[["3mer", "count"]].set_index("3mer")
    df_canonical_sbs96_fraction = df_canonical_sbs96.div(df_canonical_sbs96.sum(axis=0), axis=1)
    
    # INFO: Normalize fractions: set minimum fraction to 1 and scale others accordingly
    fraction_dict = df_canonical_sbs96_fraction['count'].to_dict()
    min_fraction = min(fraction_dict.values())
    normalized_fractions = {kmer: frac/min_fraction for kmer, frac in fraction_dict.items()}
    
    return normalized_fractions
    # return df_canonical_sbs96_fraction['count'].to_dict() # NOTE: When normalizing by kmer fractions only (deprecated?)

def extract_trinucleotide(mutation_type):
    """Extract trinucleotide context from SigProfiler SBS96 type (e.g., A[C>A]A -> ACA)"""
    return mutation_type[0] + mutation_type[2] + mutation_type[6]

def main():
    sbs96_file = snakemake.input.sbs96
    kmer_file = snakemake.input.kmer
    output_file = snakemake.output.corrected
    output_dir = snakemake.params.output_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading k-mer fractions from {kmer_file}")
    kmer_fractions = load_kmer_fractions(kmer_file)
    
    print(f"Loading SBS96 data from {sbs96_file}")
    df_sbs96 = pd.read_csv(sbs96_file, sep='\t')
    
    # Calculate original total count
    original_total = df_sbs96.iloc[:, 1].sum()
    print(f"Original total mutation count: {original_total}")
    
    # INFO: Apply k-mer correction
    print("Applying k-mer correction...")
    corrected_counts = []
    
    for _, row in df_sbs96.iterrows():
        mutation_type = row['MutationType']
        original_count = row.iloc[1]  # Second column contains the count
        
        # INFO: Extract trinucleotide context
        trinucleotide = extract_trinucleotide(mutation_type)
        
        # INFO: Get k-mer fraction
        if trinucleotide in kmer_fractions:
            kmer_fraction = kmer_fractions[trinucleotide]
            # Normalize by k-mer fraction
            corrected_count = original_count / kmer_fraction
        else:
            print(f"Warning: Trinucleotide {trinucleotide} not found in k-mer data, keeping original count")
            corrected_count = original_count
            
        corrected_counts.append(corrected_count)
    
    # Calculate the scaling factor to preserve total count
    pre_scale_total = sum(corrected_counts)
    if pre_scale_total > 0:
        scaling_factor = original_total / pre_scale_total
        print(f"Pre-scaling total: {pre_scale_total:.2f}")
        print(f"Scaling factor to preserve total: {scaling_factor:.6f}")
        
        # Apply scaling and round to integers
        final_counts = [round(count * scaling_factor) for count in corrected_counts]
    else:
        final_counts = [round(count) for count in corrected_counts]
    
    # Calculate final total count
    final_total = sum(final_counts)
    print(f"Final total mutation count: {final_total}")
    print(f"Count preservation ratio: {final_total/original_total:.6f}")
    
    # Create output dataframe
    df_corrected = df_sbs96.copy()
    df_corrected.iloc[:, 1] = final_counts
    
    # Save corrected data
    print(f"Saving k-mer corrected data to {output_file}")
    df_corrected.to_csv(output_file, sep='\t', index=False)
    
    print("K-mer correction completed successfully!")
    print(f"Total mutations: {original_total} → {final_total} (ratio: {final_total/original_total:.6f})")

if __name__ == "__main__":
    main()
