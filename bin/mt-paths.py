import sys
from collections import defaultdict
import pprint
import argparse

parser = argparse.ArgumentParser(
                    prog='mitoWHAT',
                    description='What the program does',
                    epilog='Text at the bottom of help')
parser.add_argument('phylotree', type=argparse.FileType('rt'))           # positional argument
parser.add_argument('-m', '--muts-positions', nargs='+')      # option that takes a value
parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag
args = parser.parse_args()


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
muts_pos = defaultdict(list)

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
        muts_pos[mut_pos] += [(line_num,m)]
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
            path = []
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
            
