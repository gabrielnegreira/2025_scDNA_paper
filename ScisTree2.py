#import packages
import os
os.environ["QT_QPA_PLATFORM"] = "offscreen" #allow Qt to use offscreen backend, so it can run in a non-GUI environment
import sys
import scistree2 as s2
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle, faces, AttrFace, TextFace

#set output dir
out_dir = "inputs/trees"
os.makedirs(out_dir, exist_ok=True)

#get inputs
#af_mat = "inputs/nucleotide_variants/sample_6_both_ad_matrix_for_scistree2.tsv"
af_mat = sys.argv[1] #this will get the file for the matrix as the first argument in the command-line call.
df = pd.read_csv(af_mat, sep="\t", index_col=0) #read it as a pandas dataframe with row names specified in the first column.

#get the name of the input file (needed for naming outputs)
infile_base = os.path.splitext(os.path.basename(af_mat))[0]

#Define a parser that turns "ref,alt" → (int, int), and handles "NA"
def parse_pair(s):
    if pd.isna(s) or s == "NA":
        return (0, 0)    # so no read for reference nor for alternative
    ref, alt = s.split(",")
    return (int(ref), int(alt))

#Apply the `parse_pair` function element‐wise 
#(this converts the ad matrix to the expected format by scistree2)
ad_df = df.map(parse_pair)

#subsample the matrix (to make it run faster, remove later)
#ad_df = ad_df.sample(n = 1000)

#check the object
print(ad_df.head())
print(type(ad_df))

#get cell and site names
cell_names = list(ad_df.columns)
site_names = list(ad_df.index)

#convert it to the expected format ("A three-dimensional NumPy 
#array containing the read counts for each cell at each site")
ad_array = np.array(ad_df.to_numpy().tolist(), dtype=int)

#Convert counts to `scistree.probability.GenotypeProbability` object.
prob = s2.probability.from_reads(ad_array, ado=0.2, seqerr=0.01, posterior=True, af=None, cell_names=cell_names, site_names=site_names)

# SPR local search
caller_spr = s2.ScisTree2(threads=64, max_iter=1000)
tree_spr, imputed_genotype_spr, likelihood_spr = caller_spr.infer(prob)

print('Likelihood of the SPR tree: ', likelihood_spr)
print(type(tree_spr))

# 1) Export the probability matrix as TSV
prob_tsv_file = os.path.join(out_dir, f"{infile_base}_genotype_probabilities.tsv")
prob_matrix = np.array(prob.probs, dtype=float)
pd.DataFrame(prob_matrix, index=prob.site_names, columns=prob.cell_names) \
    .to_csv(prob_tsv_file, sep="\t")
print(f"Wrote genotype probabilities to {prob_tsv_file}")

# Export the imputed genotype matrix as TSV
imputed_tsv_file = os.path.join(out_dir, f"{infile_base}_imputed_genotype.tsv")
pd.DataFrame(imputed_genotype_spr, index=site_names, columns=cell_names) \
    .to_csv(imputed_tsv_file, sep="\t")
print(f"Wrote imputed genotype to {imputed_tsv_file}")

# Export the Newick tree file
## define a function to retrieve the number of mutations in each branch
def get_num_mutations(node):
    return len(node.mutations)

nwk_tree = tree_spr.output(branch_length_func=get_num_mutations)

newick_file = os.path.join(out_dir, f"{infile_base}_inferred_tree.nwk")
#tree_obj.write(format=1, outfile=newick_file)
with open(newick_file, "w") as fh:
    fh.write(nwk_tree)

print(f"Wrote Newick to {newick_file}")

# Parse the Newick string; format=1 understands branch lengths
t = Tree(nwk_tree, format=1)

# Optional: custom layout to show leaf names clearly (good for long labels)
def layout(node):
    if node.is_leaf():
        faces.add_face_to_node(TextFace(node.name, fsize=10), node, column=0, position="aligned")
    # show branch length (mut count) as small label on edges
    if not node.is_root() and node.dist is not None:
        faces.add_face_to_node(TextFace(f"{node.dist:.0f}", fsize=8), node, column=1, position="branch-top")

ts = TreeStyle()
ts.mode = "r"                   # "r" = rectangular; try "c" for circular
ts.show_leaf_name = False       # we draw names via layout() for better control
ts.show_branch_length = False   # we add our own labels above
ts.show_scale = False
ts.layout_fn = layout
ts.branch_vertical_margin = 12  # more space between leaves
ts.scale = 480                  # overall scaling of the drawing

# Output files
png_file = os.path.join(out_dir, f"{infile_base}_inferred_tree.png")
svg_file = os.path.join(out_dir, f"{infile_base}_inferred_tree.svg")

t.render(png_file, tree_style=ts, w=3200, h=2400, units="px")
print(f"Wrote tree PNG to {png_file}")

t.render(svg_file, tree_style=ts)  # SVG is resolution-independent
print(f"Wrote tree SVG to {svg_file}")