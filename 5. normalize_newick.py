"""
In order for treetime to run as expected, gubbins divergence tree output must be normalized to convey branch lengths that represent
substitution rate / site. This is done by dividing each branch length by genome length.

Ex: python3 5.\ normalize_newick.py --input ./gubbins/st131.final_tree.tre --output ./gubbins/st131_norm.final_tree.tre --gen-length 5109767
"""

import re
import argparse

def normalize(branch, gen_length):
    updated_branch = float(branch.group(0).replace(':', '')) / gen_length
    return ':' + format(updated_branch, 'f')

parser = argparse.ArgumentParser(description='Normalize newick branch lengths by genome length')
parser.add_argument('--input',      required=True,  help='Input newick file')
parser.add_argument('--output',     required=True,  help='Output newick file')
parser.add_argument('--gen-length', required=True,  type=int, help='Genome length to normalize by')
args = parser.parse_args()

with open(args.input, 'r') as f:
    original_newick = f.read()

replacement_newick = re.sub(
    r':\d+(?:\.\d+)?',
    lambda m: normalize(m, args.gen_length),
    original_newick
)

with open(args.output, 'w') as o:
    o.write(replacement_newick)
