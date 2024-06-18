from collections import defaultdict
from Bio import SeqIO
import argparse
import pysam


parser = argparse.ArgumentParser(
                    prog='mt-paths',
                    description='Given an mtDNA phylogeny (currently phylotree), and a set of mutations, calculate a set of paths through the phylogeny that explains the mutations.',
                    epilog='Text at the bottom of help')
parser.add_argument('phylotree', type=argparse.FileType('rt'))
parser.add_argument('-m', '--muts-positions', nargs='+')
# parser.add_argument('-b', '--bam') ## needs to be implemented!
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag
parser.add_argument('-d', '--debug', action='store_true')  # on/off flag
parser.add_argument("-b", "--bamfile", required=True, help="Input BAM file")
parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
parser.add_argument("-c", "--min-coverage", type=int, default=5, help="Minimum coverage to call consensus") # Need to hook this up to actual functions
parser.add_argument("-a", "--min-agreement", type=float, default=0.8, help="Minimum agreement fraction to call consensus") # Need to hook this up to actual functions
parser.add_argument("--genbank-file", help="Path to the genbank file.")

args = parser.parse_args()


########################3
#### we need to:
##    - add the ability to read a bam file - DONE
##    - identify "foundational" mutations / mutations with some amount of evidence for them ("evidence" will ultimately be defined on the command line) - DONE
##    - have some way (given a position) to get a list of bases at that position: maybe pileup[pos] = {'A' : 6, 'a', 4, etc}? -DONE: base_counts[pos]
##    - somehow distinguish btwn bases that could be deamination: get_support(pos, allele1, allele2)? returns (5, 3), 
##    - I'm currently only matching on "having a mutation at a given position", not based on the actual alleles identified in that position, so we need to do that
##    - Once we have a clear idea of what the program looks like, pull this script apart into multiple files so it stays clean.


######## ~~Jonas' beautiful garden of Zen~~  ########


def extract_protein_coding_regions(genkbank_file):
    """Extract and print protein-coding regions from a GenBank record."""
    regions = []
    record = SeqIO.read(genkbank_file, "genbank")
    for feature in record.features:
        if feature.type == "CDS":
            location = feature.location
            regions.append((location.start, location.end))
    return regions


def is_number_overlapped(ranges, num):
    """Check if a number is overlapped by any range in a sorted list of ranges."""
    # Ensure ranges are sorted by their start positions
    ranges = sorted(ranges, key=lambda x: (x[0], x[1]))

    # Implement a manual binary search
    lo, hi = 0, len(ranges)

    while lo < hi:
        mid = (lo + hi) // 2
        start, end = ranges[mid]

        if num < start:
            hi = mid
        elif num > end:
            lo = mid + 1
        else:
            return True  # num is within the current range

    # Check the previous range to ensure no overlooked overlaps
    if lo > 0:
        start, end = ranges[lo - 1]
        if start <= num <= end:
            return True

    return False


def convert_to_uppercase(case_sensitive_bases: dict):
    '''Converts a base_count dictionary from {'c': 3, 'C': 5, 'T': 2} to {'C': 8, 'T': 2}'''
    converted = defaultdict(int)
    for key, value in case_sensitive_bases.items():
        converted[key.upper()] += value 
    return converted


def passes_filters(bases: dict, min_agreement: float = 0.8, min_coverage: int = 5) -> bool:
    '''Does most of the heavy lifting to filter out unfitting positions. Currently looks at:
    (1) Strand-aware deamination
    (2) Non-majority bases making up at least 1-$min_agreement percent of all bases
    (3) Minimum coverage
    '''

    # Removing likely deamination
    majority_base = max(bases, key=bases.get)
    keys = set(bases.keys())

    # Don't remove if any Gs are on sense strand or any Cs on antisense
    if set(['C', 't']).issubset(keys) and majority_base == 'C' and 'c' not in bases:
        bases.pop('t', None)
    if set(['g', 'A']).issubset(keys) and majority_base == 'g' and 'a' not in bases:
        bases.pop('A', None)

    bases = convert_to_uppercase(bases)

    sum_val = sum(bases.values())
    max_val = max(bases.values())
    if max_val / sum_val > min_agreement:
        return False
    
    # Removing low coverage, note: deamination is already removed at this point
    if sum_val < min_coverage:
        return False
    
    return True


def generate_base_counts(bamfile, reference_fasta):
    bam = pysam.AlignmentFile(bamfile, "rb")
    reference = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
        
    # Iterate over each reference sequence
    for ref_name, ref_seq in reference.items():
        ref_length = len(ref_seq)
        # Create an array to store base counts for each position
        base_counts = [{} for _ in range(ref_length)]
        
        # Iterate over each pileup column
        for pileupcolumn in bam.pileup(ref_name):
            pos = pileupcolumn.pos + 1 
            if pos >= ref_length:
                continue
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]

                    if pileupread.alignment.is_reverse:
                        base = base.lower()

                    base_counts[pos][base] = base_counts[pos].get(base, 0) + 1
    return base_counts


def uncalled_positions(base_counts, protein_coding_regions, min_agreement: float = 0.8, min_coverage: int = 5):
    uncalled_positions = []

    # Create the consensus sequence
    # consensus_seq = []
    for pos in range(len(base_counts)):
        if protein_coding_regions:
            if not is_number_overlapped(protein_coding_regions, pos):
                continue

        if base_counts[pos]:

            bases = base_counts[pos]

            if passes_filters(bases, min_agreement, min_coverage):
                uncalled_positions.append(pos)

            #     consensus_base = max(base_counts[pos], key=base_counts[pos].get)
            #     consensus_seq.append(consensus_base)
            # else:
            #     consensus_seq.append('N')        
        # consensus_seq = ''.join(consensus_seq)
    return uncalled_positions


if args.genbank_file:
    protein_coding_regions = extract_protein_coding_regions(args.genbank_file)
else:
    protein_coding_regions = None

base_counts = generate_base_counts(args.bamfile, args.reference)
# uncalled_positions(base_counts, protein_coding_regions, args.min_agreement, args.min_coverage)
# exit()

####### Ben's eldritch horror #######

tree = []
# pos   = sys.argv[2]

for line in args.phylotree:
    line = line.rstrip().split(',')
    tree += [line]
    pass

#print(len(tree))
#print(len(tree[0]))
#print(len(tree[1]))


tree2 = []
## mapping mutations back to the tree (mut:[line1, line2, line3])
muts_list = defaultdict(list)
## mapping mutation positions back to the tree (mut_pos:[(line1, mut), (line2, mut), (line3, mut)])
muts_pos_list = defaultdict(list)

line_num = 0
for line in tree:

    depth = 0
    haplogroup = ''

    for cell in line:
        if cell != '':
            haplogroup = cell
            break
        depth += 1
        pass

    mutations = line[depth+1].strip().split()

    #print(line)
    #print(depth, haplogroup)
    #print(mutations)

    tree2 += [{'depth':depth,
               'line_num':line_num,
               'haplogroup':haplogroup,
               'mutations':mutations}]
    for m in mutations:
        ## deletions - Deletions are indicated by the letter "d"
        ## following the (range of) position number(s) involved, e.g.
        ## "A249d" or "8281-8289d".
        if 'd' in m: continue

        ## unstable mutations - Mutations between brackets () are
        ## recurrent/unstable within the respective clade, or are yet
        ## uncertain based on current data.
        if '(' in m: continue
        

        ## insertions - Insertions are indicated by the position number
        ## preceding the insertion followed by a dot (.), the relative
        ## insert position, and the inserted base(s), e.g. "2156.1A".											
        if '.' in m: continue

        ## this only works if *all* of the numbers in the mutation tag are the position
        ## i.e. doesn't work for "2156.1A", but we remove those, above
        mut_pos = ''.join(filter(str.isdigit, m))

        print('mutation', mut_pos, m)
        muts_list[m] += [line_num]
        muts_pos_list[mut_pos] += [(line_num,m)]
        pass

    line_num += 1
    pass

## given a position or mutation, go through the tree and find all possible places that this mutation has occurred
## and then return all *paths* leading to those mutations
def find_paths(my_pos):

    ## loop over all mutations (only do this to allow both pos and mut - otherwise could just look it up in the dict)
    for mut,mut_line_nums in muts_list.items():

        
        ## compare both the mutation (A123G) and the position (123)
        ## this allows people to give both
        
        mut_pos = ''.join(filter(str.isdigit, mut))
        if mut != my_pos and mut_pos != my_pos:
            continue

        ## we have found a mutation that matches our position!
        ## and we have all line numbers where it occurs
        print(mut, mut_pos, mut_line_nums)

        all_paths = []
        
        ## for each place where it occurs, walk back up the tree to reconstruct the path
        for mln in mut_line_nums:
            depth = 1000
            # path = []
            path_details = {'mut':mut,
                            'mut_line':mln,
                            'path':[],
                            'mutations':[],
                            'path_lines':[]}

            for line_dict in reversed(tree2[:(mln+1)]):
                next_depth = line_dict['depth']
                if next_depth < depth:
                    print('tree', ' ' * next_depth, line_dict['haplogroup'], ' '.join(line_dict['mutations']))
                    path_details['path'] = [line_dict['haplogroup']] + path_details['path']
                    path_details['mutations'] = [line_dict['mutations']] + path_details['mutations']
                    path_details['path_lines'] = [line_dict['line_num']] + path_details['path_lines']
                    depth = next_depth
                    pass
                pass
            print('path', ' - '.join(path_details['path']))
            all_paths += [path_details]
            pass
        
        pass
    return all_paths

all_pos_paths = {mut:find_paths(mut) for mut in args.muts_positions}
for mut, pos_paths in all_pos_paths.items():
    print()
    print(mut)
    for d in pos_paths:
        print(d['mut'], 'path', ' - '.join(d['path']))
        pass
    pass

## somehow try to compare paths across mutations..
nested_paths = []
for mut, pos_paths in all_pos_paths.items():
    print()
    print('MATCHING', mut)

    ## loop over all paths for this mutation
    for d in pos_paths:
        print('%10s' % d['mut'], 'path    ', ' - '.join(d['path']))
        nest = {'nmuts':1,
                'muts':[mut],
                'paths':[d]}
        
        ## loop over all *other* paths, for other mutations!
        for mut2, pos_paths2 in all_pos_paths.items():
            if mut2 == mut: continue
            for d2 in pos_paths2:

                ## if any of those other paths match, record it (how?)
                if set(d2['path']).issubset(set(d['path'])):
                    print('%10s' % d2['mut'], 'sub-path', ' - '.join(d2['path']))
                    nest['nmuts'] += 1
                    nest['muts'] += [mut2]
                    nest['paths'] += [d2]
                    pass
                pass
            pass

        print(nest['muts'])
        nest['muts']  = [y for _, y in sorted(zip([len(d['path']) for d in nest['paths']], nest['muts']),  key=lambda x : x[0])]
        nest['paths'] = [y for _, y in sorted(zip([len(d['path']) for d in nest['paths']], nest['paths']), key=lambda x : x[0])]
        print(nest['muts'])
        nested_paths += [nest]
        pass
    pass

nested_paths = sorted(nested_paths, key = lambda nest : -nest['nmuts'])
print(nested_paths)

### clean up nested paths (remove any set of paths that is a strict subset of another set)

clean_nested_paths = []
seen_line_sets = []
for nest in nested_paths:
    tmp = set()
    for idx, path_dict in enumerate(nest['paths']):
        print('hey', path_dict)
        tmp.add(path_dict['mut_line'])
        pass
    is_subset = [tmp.issubset(x) for x in seen_line_sets]
    print(is_subset)
    if sum(is_subset) > 0:
        continue
    seen_line_sets += [tmp]
    clean_nested_paths += [nest]
    print('kept')
    pass


seen_muts = []
for nest in clean_nested_paths:
    print(nest['nmuts'], 'mutations:', nest['muts'])
    for idx, path_dict in enumerate(nest['paths']):
        print('**' if nest['muts'][idx] in seen_muts else '  ',
              '%10s' % nest['muts'][idx],
              'sub-path', ' - '.join(path_dict['path']))
        pass
    seen_muts += nest['muts']
    pass


exit()

## old processing
for mut,mut_line_nums in muts_list.items():
    if mut != pos:
        continue
    print(mut, mut_line_nums)


    for mln in mut_line_nums:
        depth = 1000
        path = [mut]
        for line_dict in reversed(tree2[:(mln+1)]):
            next_depth = line_dict['depth']
            if next_depth < depth:
                print('tree', ' ' * next_depth, line_dict['haplogroup'], ' '.join(line_dict['mutations']))
                path = [line_dict['haplogroup']] + path
                depth = next_depth
                pass
            pass
        print('path', ' - '.join(path))
        pass
    
    pass
            
