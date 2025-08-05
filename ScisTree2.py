#import packages
import os
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

#check the object
print(ad_df.head())

#get cell and site names
cell_names = list(ad_df.columns)
site_names = list(ad_df.index)

#convert it to the expected format ("A three-dimensional NumPy 
#array containing the read counts for each cell at each site")
ad_array = np.array(ad_df.to_numpy().tolist(), dtype=int)

print(ad_array.shape)

#Convert counts to `scistree.probability.GenotypeProbability` object.
prob = s2.probability.from_reads(ad_array, ado=0.2, seqerr=0.01, posterior=True, af=None, cell_names=cell_names, site_names=site_names)

# SPR local search
caller_spr = s2.ScisTree2(threads=8, max_iter=10000)
imputed_genotype_spr, tree_spr, likelihood_spr = caller_spr.infer(prob)

print('Likelihood of the NJ tree: ', likelihood_spr)

# 1) Export the imputed genotype matrix as TSV
imputed_tsv_file = os.path.join(out_dir, f"{infile_base}_imputed_genotype_spr.tsv")
if os.path.exists(imputed_tsv_file):
    print(f"⚠️ Warning: Overwriting {imputed_tsv_file}")
pd.DataFrame(imputed_genotype_spr, index=df_subset.index, columns=cell_names) \
    .to_csv(imputed_tsv_file, sep="\t")
print(f"Wrote imputed genotype (SPR) to {imputed_tsv_file}")

# 2) Export the Newick tree file
newick_file = os.path.join(out_dir, f"{infile_base}_inferred_tree_spr.nwk")
with open(newick_file, "w") as fh:
    fh.write(tree_spr.rstrip().rstrip(";") + ";\n")
print(f"Wrote Newick (SPR) to {newick_file}")

# 3) Export a PNG image of the tree
ts = TreeStyle()
ts.show_leaf_name = True   # display cell names
ts.show_branch_length = True
ts.scale = 120

tree_obj = Tree(tree_spr)
png_file = os.path.join(out_dir, f"{infile_base}_inferred_tree_spr.png")
tree_obj.render(png_file, tree_style=ts, w=800, h=600)
print(f"Wrote tree figure (SPR) to {png_file}")

