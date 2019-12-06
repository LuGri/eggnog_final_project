#!/usr/bin/env python

import sys
import argparse
import logging
from collections import Counter, defaultdict
from ete3 import NCBITaxa

# Constants
ncbi = NCBITaxa(dbfile="/apps/etetoolkit/taxa.sqlite")  # location of ete3 database

# Parse and check arguments
parser = argparse.ArgumentParser(description='Parse eggNOG file and determine OGs')
parser.add_argument('-min_occurence', metavar="PERCENT", action="store", type=float, default=0, help="Minimum occurence (percent of genomes where gene is present); should be 0-100, default=0")
parser.add_argument('-min_uniqueness', metavar="PERCENT", action="store", type=float, default=0, help="Minimum uniqueness if present (percent of genomes where gene is present as single-copy); should be 0-100, default=0")
parser.add_argument('-min_occurence_as_singlecopy', metavar="PERCENT", action="store", type=float, default=0, help="Minimum combined occurence+uniqueness (e.g. single-copy in 97%% of all genomes); shouled be 0-100, default=0")
parser.add_argument('-seqids_out', metavar="filename", action="store", type=str, help="Output seqids of proteins in matching OGs'")
parser.add_argument('-missing', metavar="num_missing_OGs", action="store", type=int, default=0, help="If a taxid is lacking at least this number of OGs, output it to standard error (might be interesting)")

args = parser.parse_args()
if (args.min_occurence or args.min_uniqueness) and args.min_occurence_as_singlecopy:
    print(f"""It doesn't seem to make sense to set min. occurence/min. uniqueness
    AND min. occurence+uniqueness at the same time. Pick one.""", file=sys.stderr)
    sys.exit(1)
print(args, file=sys.stderr)

#Create and configure logger 
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

# Initialize variables
all_taxids = set() #after parse_member_file(), contains all taxids occuring in parsed file as a set
good_cogs = set()
cog2taxids = {} #after parse_member_file(), contains dict with all OG-Ids as keys, and corresponding organism taxids in a list as value
cog2seqids = {} #after parse_member_file(), contains dict with all OG-Ids as keys, and corresponding sequence Ids in a list as value
taxid2missing_cogs = defaultdict(list)


def uniques(lst):
    '''Return list of elements which occur only once'''
    counter = Counter(lst)
    return [ el for el,count in counter.items() if count == 1 ]


def calc_stats(taxids, all_taxids):
    '''Calculate occurence and uniqueness of OGs.
    
    Input arguments:
        - taxids of sequences comprising an OG
        - all taxids in this taxonomic group
    Returns: tuple with 3 values
        - % occurence of this OG (with regard to all taxids comprising the tax. group)
        - % uniqueness of this OG (with regard only to the taxids where it occurs)
        - % occurence of this OG as single-copy (with regard to all taxids)
    
    Example:
    
    >>>calc_stats([2,3,4,4], set([2,3,4,5]))
    (75.0, 66.66666666666667, 50.0)
    
    The arguments mean:
    (a) the OG consists of 4 sequences in 3 organisms (taxids): [2,3,4,4]
    (b) the taxonomic group consists of 4 organisms (taxids): [2,3,4,5]
    
    The result means:
    - occurence → 75% (OG is found in 3 out of 4 taxids)
    - uniqueness → 66% (OG is single-copy in 2 out of 3 taxids)
    - occurence_as_singlecopy → 50% (OG occurs as single-copy in 2 out of 4 taxids)
    '''
       
    #taxids: the value of a key/value pair from cog2taxids (which is a list)
    #all_taxids: the "global" variable all_taxids (which is a set)
    

    #could look like this:
    #taxids: ['1226322',  '1262915',  '1519439',  '411467',  '411467',  '445972',  '742738']
    #7 sequences in 6 organisms
    
    #amount of sequences = length of variable "taxids"
    am_seqs = len(taxids)
    #amount of organisms from which sequences are derived = length of unique values of variable "taxids"
    am_taxids_from_this_OG = len(set(taxids))
    #amount of taxids in taxonomic group = length of variable "all_taxids" (must be unique values since it's a set)
    am_all_taxids = len(all_taxids)
    #amount of single copies = length of variable "taxids" with values which only occur once
    am_sc = len(uniques(taxids))
    
    
    #occurrence: amount of unique taxids in this OG divided by amount of all taxids
    occ = (am_taxids_from_this_OG/am_all_taxids) * 100 #multiply with 100 to express as percentage
    
    #uniqueness: amount of taxids in this OG which only occur once divided by amount of unique taxids in this OG
    uniq = (am_sc / am_taxids_from_this_OG) * 100
    
    #occurrence as single copy: amount of taxids in this OG which only occur once divided by amount of all taxids
    occ_as_singlecopy = (am_sc / am_all_taxids) * 100
        
    #return as tuple in this exact order
    return (occ, uniq, occ_as_singlecopy)


def parse_members_file(fin):
    '''e.g. /mirror/eggnog/eggnog_5.0/per_tax_level/1/1_members.tsv.gz
    1	28H59	4	3	565033.GACE_1005,572546.Arcpr_0848,69014.TK0075,69014.TK0091	565033,572546,69014
    -> means: level "1" (root), OG identifier "28H59", 4 sequences in 3 organisms, sequence identifiers, organism taxids 
    '''
    for line in fin:
        level, cog_id, n_seqs, n_taxids, seq_ids, taxids = line.strip().split('\t')

        # Taxids where OG occurs
        taxids = taxids.split(',')
        all_taxids.update(taxids)
        # consistency check
        assert len(taxids) == len(set(taxids)) == int(n_taxids), \
        f'{line}\n {cog_id} -> {len(taxids)}---{len(set(taxids))}---{n_taxids}'

        # Sequences comprising the OG
        seq_ids = seq_ids.split(',')
        assert len(seq_ids) == int(n_seqs)  # consistency check
        seq_ids_taxids = [ el.split('.')[0] for el in seq_ids ]
        
        # Remember the important stuff
        assert cog_id not in cog2taxids  # consistency check
        cog2taxids[cog_id] = seq_ids_taxids
        cog2seqids[cog_id] = seq_ids

        
def missing_taxids(cog): #, taxids):
    '''Determine which OGs are missing from which taxids'''    
    
    #cog = the ID of an COG which fulfilled the initial condition
    #eg. OG which occurs as a single copy in 97% of all taxids
    
    #for every taxid in the input file
    #check if it's not a part of the taxids in which the current cog is present
    #if it isn't, add the current cog-ID to the corresponding taxid in the defaultdict
    
    for taxid in all_taxids:
        if taxid not in cog2taxids[cog]:
            taxid2missing_cogs[taxid].append(cog)
    
    #at the end, defaultdict is filled with
    #keys = all taxids in the input file
    #values = lists of all OGs which fulfilled the initial condition but aren't present in the respective taxid
            
        
def output_seqids(filename, cogs):
    '''Output seqids to file'''
    with open(filename, 'w') as fout:
        for cog in cogs:
            for seqid in cog2seqids[cog]:
                taxid = seqid.split('.')[0]
                print(f'{seqid}\t{cog}', file=fout)


def output_cogs():
    '''Output OGs matching the criteria'''
    print(f'#cog\t%_occurence\tthereof_%_singlecopy\t%_occurence_as_singlecopy', file=sys.stdout)

    good_cogs = set()
    for cog, taxids in cog2taxids.items():
        occurence, uniqueness, occurence_as_singlecopy = calc_stats(taxids, all_taxids)
        logger.debug(f'''{cog} → occurence, uniqueness, occurence_as_singlecopy: 
        {occurence:.3f}, {uniqueness:.3f}, {occurence_as_singlecopy:.3f}''')
        
        #if cog == 'COG1841': print(occurence, uniqueness, occurence_as_singlecopy, file=sys.stderr) 
        if ((occurence < args.min_occurence) or
            (uniqueness < args.min_uniqueness) or
            (occurence_as_singlecopy < args.min_occurence_as_singlecopy)):
            continue
        
        #missing = missing_taxids(cog, taxids) #not sure what I'm supposed to do with the 2nd parameter
        missing_taxids(cog)
        
        good_cogs.add(cog)
        logger.debug(f'{cog} accepted → num good cogs: {len(good_cogs)}')
        
        output = '\t'.join([cog, f'{occurence:.1f}', f'{uniqueness:.1f}', f'{occurence_as_singlecopy:.1f}'])
        print(output, file=sys.stdout)

    # Output sequence identifiers
    if args.seqids_out:
        output_seqids(args.seqids_out, good_cogs)

    # Output taxids with missing cogs
    if args.missing:
        # Your code here
        print("")
        print(f"Taxids which are missing at least {args.missing} OGs which fulfilled the initial condition")
        print("")
        for taxid, cogs in taxid2missing_cogs.items():
            if len(cogs) >= args.missing:
                output = '\t'.join([taxid, str(len(cogs))])
                print(output, file=sys.stdout)
                    

if __name__ == '__main__':
    
    logger.info(f'Reading input...')
    parse_members_file(sys.stdin)
    
    logger.info(f'Outputting result...')
    output_cogs()

