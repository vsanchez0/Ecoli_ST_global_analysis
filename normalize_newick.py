import sys
import re

def normalize(branch):
    '''
    Normalizes the branch length from mutational distance to number of mutations per site

    Parameters:
        length (str): portion of newick denoting branch length, including preceding semicolon

    Return the normalized branch length, including preceding semicolon
    '''
    updated_branch = float(branch.group(0).replace(':', '')) / 5109767
    
    return ':' + format(updated_branch, 'f')

# Import newick string
f = open(sys.argv[1], 'r')
original_newick = f.read()
f.close()

replacement_newick = re.sub(r':\d+(?:\.\d+)?', normalize, original_newick)

# Create output newick with updated branch lengths
o = open(sys.argv[1].split('.')[0] + '_norm.nwk', 'w')
o.write(replacement_newick)
o.close()
