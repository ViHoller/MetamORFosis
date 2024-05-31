import pandas as pd
import re
from Biopython import Seq
import argparse

def read_fasta(f): # update description since fasta file will only contain one sequence
    """
    It reads a fasta file and returns a dictionary with accessions as keys 
    and sequences as values. 
    It raises ValueErrors if file path does not exist.
    this function is not needed in this iteration of the script but kept
    in case it is needed again
    
    Args:
        f (string): File path.
        
    Returns:
        fasta_out (dict): Dictionary with accessions as keys and sequences as 
        values.
    """
    fasta_out = ''
    try:
        with open(f, "r") as h_in:
            for seq in h_in:      
                if seq.startswith(">"):
                    continue
                else:
                    fasta_out += seq
            return fasta_out
    except FileNotFoundError:
        raise ValueError(f"The file path {f} doen't exist.")

def find_all_orfs(genome, threshold):
    stop_list = ['TAA', 'TAG', 'TGA']
    matches = re.finditer('ATG', genome)
    all_orfs = pd.DataFrame(columns=['start_pos', 'stop_pos', 'sequence'])
    for start in matches:
        start_position = start.start()
        position = start_position
        while True:
            if position+3 > len(genome):
                break
            if genome[position:position+3] in stop_list:
                seq = genome[start_position:position+3]
                break
            else:
                position += 3
        if len(seq) > threshold:
            all_orfs.loc[len(all_orfs)] = {'start_pos':start_position, 'stop_pos':position, 'sequence':seq}
    return(all_orfs)

def merge_dfs(df1,df2): # marge dataframes and order them by position
    pass

def write_fasta(df, filename):
    pass

def main(args):
    threshold = args.t
    genome = read_fasta(args.input_filename)
    rev_genome = Seq(genome).reverse_compliment()
    genome_orfs = find_all_ord(genome, threshold)
    rev_genome_orfs = find_all_ord(rev_genome, threshold)
    all_orfs = merge_dfs(genome, ref_genome)
    write_fasta(all_orfs, args.outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    prog="",
    description="",
    )
    parser.add_arguemnt(
    'input_filename',
    help=""     
    )
    parser.add_arguemtn(
    'outfile',
    help=""
    )
    parser.add_argument(
    '--t',
    help="threshold of ORF size" # better description,
    default=123 # ADD THE STANDARD THRESHOLD
    )
    main(parser.parse_argument())

# maybe translation and BLAST
# add later - validation - promoters and codon composition
