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

def trinuc_context_change(string):
    if len(string) != 7 or string[3] != '>':
        raise ValueError("Input string must be in the format 'NNN>NNN'")

    if string[1] not in ['C', 'T']:
        string_pair = string.split('>')
        new_string = reverse_complement(string_pair[0]) + '>' + reverse_complement(string_pair[1])
    else:
        new_string = string
    return new_string

def main():
    # When called from Snakemake, use snakemake object
    if 'snakemake' in globals():
        input_file = snakemake.input.sbs96
        sample_name = snakemake.params.sample
        output_dir = snakemake.params.output_dir
    else:
        # Command line usage
        if len(sys.argv) != 4:
            print("Usage: python process_sbs96.py <input_sbs96_file> <sample_name> <output_dir>")
            sys.exit(1)

        input_file = sys.argv[1]
        sample_name = sys.argv[2]
        output_dir = sys.argv[3]

    # Install GRCh38 reference if needed
    try:
        genInstall.install('GRCh38')
    except (FileExistsError, FileNotFoundError, OSError, PermissionError) as e:
        # Reference already installed or concurrent installation issue
        # Try to wait and retry once
        import time
        time.sleep(2)
        try:
            genInstall.install('GRCh38')
        except Exception:
            # If still fails, continue - reference might already be available
            print(f"Warning: Reference installation issue (continuing): {e}")
            pass

    sbs96_sigprofiler="""A[C>A]A    A[C>A]C    A[C>A]G    A[C>A]T    A[C>G]A    A[C>G]C    A[C>G]G    A[C>G]T    A[C>T]A    A[C>T]C    A[C>T]G    A[C>T]T    A[T>A]A    A[T>A]C    A[T>A]G    A[T>A]T    A[T>C]A    A[T>C]C    A[T>C]G    A[T>C]T    A[T>G]A    A[T>G]C    A[T>G]G    A[T>G]T    C[C>A]A    C[C>A]C    C[C>A]G    C[C>A]T    C[C>G]A    C[C>G]C    C[C>G]G    C[C>G]T    C[C>T]A    C[C>T]C    C[C>T]G    C[C>T]T    C[T>A]A    C[T>A]C    C[T>A]G    C[T>A]T    C[T>C]A    C[T>C]C    C[T>C]G    C[T>C]T    C[T>G]A    C[T>G]C    C[T>G]G    C[T>G]T    G[C>A]A    G[C>A]C    G[C>A]G    G[C>A]T    G[C>G]A    G[C>G]C    G[C>G]G    G[C>G]T    G[C>T]A    G[C>T]C    G[C>T]G    G[C>T]T    G[T>A]A    G[T>A]C    G[T>A]G    G[T>A]T    G[T>C]A    G[T>C]C    G[T>C]G    G[T>C]T    G[T>G]A    G[T>G]C    G[T>G]G    G[T>G]T    T[C>A]A    T[C>A]C    T[C>A]G    T[C>A]T    T[C>G]A    T[C>G]C    T[C>G]G    T[C>G]T    T[C>T]A    T[C>T]C    T[C>T]G    T[C>T]T    T[T>A]A    T[T>A]C    T[T>A]G    T[T>A]T    T[T>C]A    T[T>C]C    T[T>C]G    T[T>C]T    T[T>G]A    T[T>G]C    T[T>G]G    T[T>G]T"""
    sbs96_sigprofiler = sbs96_sigprofiler.split()

    sbs96 = dict()
    for sbs in sbs96_sigprofiler:
        sbs96[f'{sbs[0]}{sbs[2]}{sbs[-1]}>{sbs[0]}{sbs[-3]}{sbs[-1]}'] = sbs

    mutyper_trinuc = [
        'AAA>ACA', 'AAA>AGA', 'AAA>ATA', 'AAC>ACC', 'AAC>AGC', 'AAC>ATC',
        'AAG>ACG', 'AAG>AGG', 'AAG>ATG', 'AAT>ACT', 'AAT>AGT', 'AAT>ATT',
        'ACA>AAA', 'ACA>AGA', 'ACA>ATA', 'ACC>AAC', 'ACC>AGC', 'ACC>ATC',
        'ACG>AAG', 'ACG>AGG', 'ACG>ATG', 'ACT>AAT', 'ACT>AGT', 'ACT>ATT',
        'CAA>CCA', 'CAA>CGA', 'CAA>CTA', 'CAC>CCC', 'CAC>CGC', 'CAC>CTC',
        'CAG>CCG', 'CAG>CGG', 'CAG>CTG', 'CAT>CCT', 'CAT>CGT', 'CAT>CTT',
        'CCA>CAA', 'CCA>CGA', 'CCA>CTA', 'CCC>CAC', 'CCC>CGC', 'CCC>CTC',
        'CCG>CAG', 'CCG>CGG', 'CCG>CTG', 'CCT>CAT', 'CCT>CGT', 'CCT>CTT',
        'GAA>GCA', 'GAA>GGA', 'GAA>GTA', 'GAC>GCC', 'GAC>GGC', 'GAC>GTC',
        'GAG>GCG', 'GAG>GGG', 'GAG>GTG', 'GAT>GCT', 'GAT>GGT', 'GAT>GTT',
        'GCA>GAA', 'GCA>GGA', 'GCA>GTA', 'GCC>GAC', 'GCC>GGC', 'GCC>GTC',
        'GCG>GAG', 'GCG>GGG', 'GCG>GTG', 'GCT>GAT', 'GCT>GGT', 'GCT>GTT',
        'TAA>TCA', 'TAA>TGA', 'TAA>TTA', 'TAC>TCC', 'TAC>TGC', 'TAC>TTC',
        'TAG>TCG', 'TAG>TGG', 'TAG>TTG', 'TAT>TCT', 'TAT>TGT', 'TAT>TTT',
        'TCA>TAA', 'TCA>TGA', 'TCA>TTA', 'TCC>TAC', 'TCC>TGC', 'TCC>TTC',
        'TCG>TAG', 'TCG>TGG', 'TCG>TTG', 'TCT>TAT', 'TCT>TGT', 'TCT>TTT'
    ]

    os.makedirs(output_dir, exist_ok=True)

    # Read mutyper output
    df = pd.read_csv(input_file, sep="\t", header=None).T
    df.columns = ['SBS96_pre', 'Count']

    mutyper_template = pd.DataFrame({
        'SBS96_pre': mutyper_trinuc,
        'Count': 0
    })

    df = mutyper_template.merge(df, on="SBS96_pre", how='left', suffixes=('', '_new'))
    df['Count'] = df['Count_new'].fillna(0)
    df = df.drop(columns=['Count_new'])

    df['SBS96'] = df['SBS96_pre'].apply(trinuc_context_change)
    df['SBS96_SigProfiler'] = df['SBS96'].apply(lambda x: sbs96.get(x, None))

    df['SBS96_SigProfiler'] = pd.Categorical(df['SBS96_SigProfiler'], categories=sbs96_sigprofiler, ordered=True)
    df = df.sort_values(by='SBS96_SigProfiler').reset_index(drop=True)

    df.rename(columns={'SBS96_SigProfiler': 'MutationType', 'Count': sample_name}, inplace=True)
    output_matrix = os.path.join(output_dir, f"{sample_name}.SBS96.all")
    df[['MutationType', sample_name]].to_csv(output_matrix, sep='\t', index=False)

    sigPlt.plotSBS(
        matrix_path=output_matrix, 
        output_path=output_dir, 
        project=sample_name, 
        plot_type="96", 
        savefig_format="pdf",
        percentage=False)

    sigPlt.plotSBS(
        matrix_path=output_matrix, 
        output_path=output_dir, 
        project=f"{sample_name}.percentage", 
        plot_type="96", 
        savefig_format="pdf",
        percentage=True)

    Analyze.cosmic_fit(
        output_matrix, 
        output_dir, 
        input_type="matrix", 
        context_type="96", 
        collapse_to_SBS96=True, 
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
