#import packages
import os
os.environ["QT_QPA_PLATFORM"] = "offscreen" #allow Qt to use offscreen backend, so it can run in a non-GUI environment
import sys
import scistree2 as s2
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle

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
ad_df = ad_df.sample(n = 1000)

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
caller_spr = s2.ScisTree2(threads=64, max_iter=100000000)
imputed_genotype_spr, tree_spr, likelihood_spr = caller_spr.infer(prob)

print('Likelihood of the SPR tree: ', likelihood_spr)

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
newick_file = os.path.join(out_dir, f"{infile_base}_inferred_tree.nwk")
with open(newick_file, "w") as fh:
    fh.write(tree_spr.rstrip().rstrip(";") + ";\n")
print(f"Wrote Newick to {newick_file}")

# Export a PNG image of the tree
ts = TreeStyle()
ts.show_leaf_name = True   # display cell names
ts.show_branch_length = True
ts.scale = 120

tree_obj = Tree(tree_spr)
png_file = os.path.join(out_dir, f"{infile_base}_inferred_tree.png")
tree_obj.render(png_file, tree_style=ts, w=800, h=600)
print(f"Wrote tree figure to {png_file}")