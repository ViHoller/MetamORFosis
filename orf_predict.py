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
    Naming the orfs: (rev)_[number]|[start]|[length]
    >filename| added later

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
            orf_id = f'{len(all_orfs)}|start_pos:{start_position}|{len(seq)}_bp'
            if rev == True:
                orf_id = 'rev_' + orf_id
            all_orfs.loc[len(all_orfs)] = {'orf_id':orf_id, 'start_pos':start_position, 'stop_pos':position, 'sequence':seq}

    return(all_orfs)
    
def remove_overlapping(orf_df, rev=False):
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

def get_orfs(genome, threshold):
    '''
    This function is the executive function for getting the orfs. It gets orfs for both forward and reverse strand,
    inverts coordinates for the reverse strand, removes overlapping orfs and concatenates the df into one containing
    all orfs.
    '''
    genome_length = len(genome)
    rev_genome = str(Seq(genome).reverse_complement())
    
    orfs = find_all_orfs(genome, threshold)
    rev_orfs = find_all_orfs(rev_genome, threshold, rev=True)
    
    rev_orfs['start_pos'] = genome_length - rev_orfs['start_pos']
    rev_orfs['stopo_pos'] = genome_length - rev_orfs['stop_pos']

    orfs = remove_overlapping(orfs)
    rev_orfs = remove_overlapping(rev_orfs, rev=True)
    
    all_orfs = merge_dfs(orfs, rev_orfs)
    
    return all_orfs
    


def write_fasta(df, in_filename, out_filename):
    with open(out_filename, 'w') as out_f:
        for index, row in df.iterrows():
            out_f.write(f'>{in_filename}|{row['orf_id']}\n')
            out_f.write(f'{row['sequence']}\n\n')

def main(args):
    threshold = args.t
    genome = read_fasta(args.input_filename)
    all_orfs = get_orfs(genome, threshold)
    write_fasta(all_orfs, args.input_filename, args.outfile)

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
