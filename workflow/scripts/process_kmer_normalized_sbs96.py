#!/usr/bin/env python3
# Author: mhsohny (Min-Hwan "Sohny" Sohn)

import pandas as pd
import os
import sys
from SigProfilerMatrixGenerator import install as genInstall
import sigProfilerPlotting as sigPlt
from SigProfilerAssignment import Analyzer as Analyze

def main():
    # When called from Snakemake, use snakemake object
    if 'snakemake' in globals():
        input_file = snakemake.input.sbs96
        sample_name = snakemake.params.sample
        output_dir = snakemake.params.output_dir
    else:
        # Command line usage
        if len(sys.argv) != 4:
            print("Usage: python process_kmer_normalized_sbs96.py <input_sbs96_file> <sample_name> <output_dir>")
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

    os.makedirs(output_dir, exist_ok=True)

    # The input file is already in SigProfiler format, so we can use it directly
    matrix_file = input_file

    # Generate SBS96 plots
    sigPlt.plotSBS(
        matrix_path=matrix_file, 
        output_path=output_dir, 
        project=sample_name, 
        plot_type="96", 
        savefig_format="pdf",
        percentage=False)

    sigPlt.plotSBS(
        matrix_path=matrix_file, 
        output_path=output_dir, 
        project=f"{sample_name}.percentage", 
        plot_type="96", 
        savefig_format="pdf",
        percentage=True)

    # Perform COSMIC signature analysis
    Analyze.cosmic_fit(
        matrix_file, 
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