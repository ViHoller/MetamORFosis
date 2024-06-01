import pandas as pd
import re
from Bio.Seq import Seq
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

def find_all_orfs(genome, threshold, rev=False):
    '''
    This function searches a sequence for all orfs as a stretch from a start codon to an in frame
    stop codon. This list is passed to remove_overlapping() to remove all overlapping orfs.
    Naming the orfs: >orf_[number]_[start](_rev) 

    Arguments:

        genome:     a string containing the sequence to be analysed

        threshold:  an int, defining the minimum length of an orf

        rev:        boolean, whether it is the reverse strand

    '''

    stop_list = ['TAA', 'TAG', 'TGA']
    matches = re.finditer('ATG', genome)
    all_orfs = pd.DataFrame(columns=['orf_id', 'start_pos', 'stop_pos', 'sequence'])

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
            orf_id = f'>orf_{len(all_orfs)}_{start_position}'
            if rev == True:
                orf_id = orf_id + '_rev'
            all_orfs.loc[len(all_orfs)] = {'orf_id':orf_id, 'start_pos':start_position, 'stop_pos':position, 'sequence':seq}

    all_orfs = remove_overlapping(all_orfs, rev)

    return(all_orfs)

def remove_overlapping(orf_df, rev):
    '''
    this function makes a list of all orfs whose start is inside one of the previous orfs,
    and deletes all rows in the list at the end. In case the orfs are not sorted by start position already,
    the function sorts the dataframe.    

    Arguments:

        orf_df: dataframe containing orfs, still overlapping

        rev:    boolean, defines if reverse strand or not

    '''

    current_stop = 0
    drop_list = []
    orf_df.sort_values(by=['start_pos'], ascending=[True])

    for ind, row in orf_df.iterrows():

        start = row['start_pos']
        
        if start <= current_stop and not rev:
            drop_list.append(ind)
        elif start > current_stop and not rev:
            current_stop = row['stop_pos']
        elif start >= current_stop and rev:
            drop_list.append(ind)
        elif start < current_stop and rev:
            current_stop = row['stop_pos']

    orf_df.drop(index=drop_list, inplace=True)
    orf_df.reset_index(drop=True, inplace=True)

    return(orf_df)

def merge_dfs(df1,df2): # marge dataframes and order them by position
    merged_df = pd.concat([df1,df2], axis=0)
    sorted_orfs = merged_df.sort_values(by='start_pos')
    #print(df1.shape, df2.shape, merged_df.shape, sorted_orfs.shape)
    return sorted_orfs

def write_fasta(df, filename):
    pass

def main(args):
    threshold = args.t
    genome = read_fasta(args.input_filename)
    rev_genome = str(Seq(genome).reverse_complement())
    genome_orfs = find_all_orfs(genome, threshold)
    rev_genome_orfs = find_all_orfs(rev_genome, threshold)
    all_orfs = merge_dfs(genome_orfs, rev_genome_orfs)
    write_fasta(all_orfs, args.outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    prog="",
    description="",
    )
    parser.add_argument(
    'input_filename',
    help=""     
    )
    parser.add_argument(
    'outfile',
    help=""
    )
    parser.add_argument(
    '--t',
    help="threshold of ORF size", # better description
    default=123 # ADD THE STANDARD THRESHOLD
    )
    main(parser.parse_args())

# maybe translation and BLAST
# add later - validation - promoters and codon composition
