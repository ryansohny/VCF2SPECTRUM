#!/usr/bin/env python3
# Author: mhsohny (Min-Hwan "Sohny" Sohn)

import pandas as pd
import os
import sys
from SigProfilerMatrixGenerator import install as genInstall
import sigProfilerPlotting as sigPlt
from SigProfilerAssignment import Analyzer as Analyze

def reverse_complement(string):
    try:
        complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        complement_string = ''.join([complement_dict[s] for s in string])
    except KeyError:
        raise ValueError("Invalid character other than A,T,G and C")
    return complement_string[::-1]

def dbs_context_change(string):
    if len(string) != 5 or string[2] != '>':
        raise ValueError("Input string must be in the format 'NN>NN'")
    
    string_pair = string.split('>')
    new_string = reverse_complement(string_pair[0]) + '>' + reverse_complement(string_pair[1])
    
    return new_string

def read_vcf(path):
    import gzip as gz
    import io
    if path.endswith('.gz'): 
        with gz.open(path, 'rb') as f:
            lines = [l.decode('utf-8') for l in f if not l.startswith(b'##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                       sep='\t'
                       ).rename(columns={'#CHROM': 'CHROM'})
    else:
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str},
                       sep='\t'
                       ).rename(columns={'#CHROM': 'CHROM'})

def create_dbs78_signature(df: pd.core.frame.DataFrame, sampleid: str, output_dir: str):
    dbs78_sigprofiler = ('AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG')
    dbs78_sigprofiler = dict.fromkeys(dbs78_sigprofiler, 0)
    
    df_sorted = df.sort_values(by=["CHROM", "POS"]).reset_index(drop=True)
    
    df_sorted['pos_diff'] = df_sorted.groupby("CHROM")['POS'].diff()
    df_sorted['group_id'] = ((df_sorted['pos_diff'] != 1) | df_sorted['pos_diff'].isna()).cumsum()
    
    for (chrom, group_id), group in df_sorted.groupby(['CHROM', 'group_id']):
        group_sorted = group.sort_values('POS').reset_index(drop=True)
        
        if len(group_sorted) % 2 != 0:
            continue

        for i in range(0, len(group_sorted), 2):
            if i + 1 < len(group_sorted):
                snv1 = group_sorted.iloc[i]
                snv2 = group_sorted.iloc[i + 1]
                
                if snv2['POS'] == snv1['POS'] + 1:
                    ref_dinuc = snv1['REF'] + snv2['REF']
                    alt_dinuc = snv1['ALT'] + snv2['ALT']
                    
                    dbs = f"{ref_dinuc}>{alt_dinuc}"
                
                    if dbs in dbs78_sigprofiler:
                        dbs78_sigprofiler[dbs] += 1
                    elif dbs_context_change(dbs) in dbs78_sigprofiler:
                        dbs78_sigprofiler[dbs_context_change(dbs)] += 1
    
    new_df = pd.DataFrame(dbs78_sigprofiler, index=[0]).T
    new_df.columns = [sampleid]
    new_df.index.name = "MutationType"
    
    output_file = os.path.join(output_dir, f"{sampleid}.DBS78.all")
    new_df.to_csv(output_file, sep="\t")
    
    return output_file

def main():
    # When called from Snakemake, use snakemake object
    if 'snakemake' in globals():
        input_file = snakemake.input[0]
        sample_name = snakemake.params.sample
        output_dir = snakemake.params.output_dir
    else:
        # Command line usage
        if len(sys.argv) != 4:
            print("Usage: python process_dbs78.py <input_vcf_file> <sample_name> <output_dir>")
            sys.exit(1)

        input_file = sys.argv[1]
        sample_name = sys.argv[2]
        output_dir = sys.argv[3]

    # Install GRCh38 reference if needed
    try:
        genInstall.install('GRCh38')
    except FileExistsError:
        # Reference already installed
        pass

    os.makedirs(output_dir, exist_ok=True)

    df = read_vcf(input_file)
    
    output_matrix = create_dbs78_signature(df, sample_name, output_dir)

    sigPlt.plotDBS(
        matrix_path=output_matrix,
        output_path=output_dir,
        project=sample_name,
        plot_type="78",
        savefig_format="pdf",
        percentage=False)

    sigPlt.plotDBS(
        matrix_path=output_matrix,
        output_path=output_dir,
        project=f"{sample_name}.percentage",
        plot_type="78",
        savefig_format="pdf",
        percentage=True)

    Analyze.cosmic_fit(
        output_matrix, 
        output_dir, 
        input_type="matrix", 
        context_type="DINUC", 
        collapse_to_SBS96=False, 
        cosmic_version=3.4, 
        exome=False,
        genome_build="GRCh38", 
        signature_database=None,
        exclude_signature_subgroups=None, 
        export_probabilities=True,
        export_probabilities_per_mutation=False, 
        make_plots=True,
        sample_reconstruction_plots="pdf", 
        verbose=False)

if __name__ == "__main__":
    main()
