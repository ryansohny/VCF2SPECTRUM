#!/usr/bin/env python3
# Author: mhsohny (Min-Hwan "Sohny" Sohn)
import pandas as pd
import os
import sys
import csv
import subprocess
from pathlib import Path
from tqdm import tqdm
from pyfaidx import Fasta
import sigProfilerPlotting as sigPlt
from SigProfilerAssignment import Analyzer as Analyze

def unix_wc(file: Path) -> int:
    wc = subprocess.run(['wc', '-l', str(file)], stdout=subprocess.PIPE)
    return int(wc.stdout.split()[0])

def get_fasta_sequence(fasta_file, chrom: str, start: int, end: int) -> str:
    fasta = Fasta(fasta_file, rebuild=False)
    return fasta[chrom][start-1:end].seq

def count_stretches(sequence: str, nucleotide: str) -> int:
    c = 0
    window = len(nucleotide)

    # in PCAWG, max Repeat size is 5
    max_repeat = 5
    for i in range(max_repeat):
        subseq = sequence[window*i: window*(i+1)]
        if subseq == nucleotide:
            c += 1
        else:
            break
    return c

def get_microhomology(sequence: str, right_flank: str, left_flank: str) -> int:
    """
    :param sequence: deleted sequence from reference
    :param right_flank: right flanking sequence from deleted bases 
    :param left_flank: left flanking sequence from deletd bases
    """
    right_mh = 0
    left_mh = 0

    for i in range(len(sequence)-1):
        if sequence[i] == right_flank[i]:
            right_mh += 1
        else:
            break
    for j in range(1, len(sequence)):
        if sequence[-j] == left_flank[-j]:
            left_mh += 1
        else:
            break

    # INFO: Sanity check
    if max(right_mh, left_mh) == 0:
        raise ValueError("Max(right homology, left homology) cannot be 0 in the code setting! Something's wrong")

    return min(max(right_mh, left_mh), 5)

def id83_generator(inputf: Path, 
                   outputf: Path, 
                   fasta_file: Path) -> dict:
    id83 = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5', '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5', '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5', '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5', '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5', '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5', '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5', '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3', '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5']
    id83_dict = dict.fromkeys(id83, 0)

    with open(inputf, 'r') as dfh, open(outputf, 'w', newline='') as rfh:

        wc = unix_wc(inputf)

        reader = csv.reader(dfh, delimiter='\t')
        writer = csv.writer(rfh, delimiter='\t')

        header = next(reader) # \tCHROM\tPOS\tID
        header.append('ID83')
        writer.writerow(header)

        CHROM_index = header.index("CHROM")
        POS_index = header.index("POS")
        REF_index = header.index("REF")
        ALT_index = header.index("ALT")
        INFO_index = header.index("INFO")

        for row in tqdm(reader, desc="Processing rows", total=wc-1):
            # INFO: check if the record is insertion or deletion
            refseq, altseq = row[REF_index], row[ALT_index]
            reflen, altlen = len(refseq), len(altseq)

            if reflen > altlen:
                deletion = True
                idseq = refseq[1:]
            else:
                deletion = False
                idseq = altseq[1:]

            idsize = abs(reflen - altlen)
            to_pyrimidine = {'A':'T', 'G':'C', 'T':'T', 'C':'C'}

            # INFO: Identifying the sequence context
            # NOTE: Deletion
            if deletion:
                if idsize == 1:
                    seq_context = get_fasta_sequence(fasta_file, 
                                                     row[CHROM_index], 
                                                     int(row[POS_index]) + idsize + 1, 
                                                     int(row[POS_index]) + idsize + (idsize*5)
                                                     )
                    stretches = count_stretches(seq_context, idseq)

                    id83_category = f"{idsize}:Del:{to_pyrimidine[idseq]}:{stretches}"
                    id83_dict[id83_category] += 1
                    row.append(id83_category)
                    writer.writerow(row)
                else:
                    right_seq_context = get_fasta_sequence(fasta_file, 
                                                           row[CHROM_index],
                                                           int(row[POS_index]) + idsize + 1,
                                                           int(row[POS_index]) + idsize + (idsize*5)
                                                           )
                    right_stretches = count_stretches(right_seq_context, idseq)

                    #NOTE: 1. Having repeat units 
                    if right_stretches != 0:
                        id83_category = f"{min(idsize, 5)}:Del:R:{right_stretches}"
                        id83_dict[id83_category] += 1
                        row.append(id83_category)
                        writer.writerow(row)

                    else:
                        left_seq_context = get_fasta_sequence(fasta_file, 
                                                              row[CHROM_index],
                                                              int(row[POS_index]) - (idsize - 2),
                                                              int(row[POS_index])
                                                              )

                        #BUG:check
                        if len(left_seq_context) + 1 != idsize:
                            print(idsize)
                            print(row)
                            print(right_seq_context)
                            print(left_seq_context)
                            raise ValueError("Something's wrong")

                        #NOTE: 2. No microhomology
                        if idseq[0] != right_seq_context[0] and idseq[-1] != left_seq_context[-1]:
                            id83_category = f"{min(idsize, 5)}:Del:R:0"
                            id83_dict[id83_category] += 1
                            row.append(id83_category)
                            writer.writerow(row)
                        #NOTE: 3. Microhomology (NOTE2 can be included in the NOTE3, for efficiency's sake)
                        else:
                            mh = get_microhomology(idseq, right_seq_context, left_seq_context)
                            id83_category = f"{min(idsize, 5)}:Del:M:{mh}"
                            id83_dict[id83_category] += 1
                            row.append(id83_category)
                            writer.writerow(row)

            # NOTE: Insertion
            elif not deletion:
                if idsize == 1:
                    seq_context = get_fasta_sequence(fasta_file, 
                                                     row[CHROM_index], 
                                                     int(row[POS_index]) + 1, 
                                                     int(row[POS_index]) + (idsize*5)
                                                     )
                    stretches = count_stretches(seq_context, idseq)

                    id83_category = f"{idsize}:Ins:{to_pyrimidine[idseq]}:{stretches}"
                    id83_dict[id83_category] += 1
                    row.append(id83_category)
                    writer.writerow(row)

                else:
                    seq_context = get_fasta_sequence(fasta_file, 
                                                     row[CHROM_index],
                                                     int(row[POS_index]) + 1,
                                                     int(row[POS_index]) + (idsize*5)
                                                     )
                    stretches = count_stretches(seq_context, idseq)

                    id83_category = f"{min(idsize, 5)}:Ins:R:{stretches}"
                    id83_dict[id83_category] += 1
                    row.append(id83_category)
                    writer.writerow(row)

            rfh.flush()
    return id83_dict

def main():
    # When called from Snakemake, use snakemake object
    if 'snakemake' in globals():
        input_file = snakemake.input.id83
        sample_name = snakemake.params.sample
        output_dir = snakemake.params.output_dir
        reference_genome = snakemake.input.ref
    else:
        # Command line usage
        if len(sys.argv) != 5:
            print("Usage: python process_id83.py <input_id83_file> <sample_name> <output_dir> <reference_genome>")
            sys.exit(1)

        input_file = sys.argv[1]
        sample_name = sys.argv[2]
        output_dir = sys.argv[3]
        reference_genome = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    temp_output = os.path.join(output_dir, f"{sample_name}_ID83.tsv")
    id83_dict = id83_generator(Path(input_file), Path(temp_output), Path(reference_genome))

    output_matrix = os.path.join(output_dir, f"{sample_name}.ID83.all")
    pd.DataFrame(id83_dict.items(), 
                 columns=['MutationType', sample_name]
                 ).to_csv(output_matrix, sep='\t', index=False)


    sigPlt.plotID(
        matrix_path=output_matrix, 
        output_path=output_dir, 
        project=sample_name, 
        plot_type="83", 
        savefig_format="pdf",
        percentage=False)

    sigPlt.plotID(
        matrix_path=output_matrix, 
        output_path=output_dir, 
        project=f"{sample_name}.percentage", 
        plot_type="83", 
        savefig_format="pdf",
        percentage=True)

    Analyze.cosmic_fit(
        output_matrix, 
        output_dir, 
        input_type="matrix", 
        context_type="83", 
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
