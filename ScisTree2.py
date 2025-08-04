#import packages
import scistree2 as s2 
import numpy as np 
import pandas as pd
import os 
from ete3 import Tree, TreeStyle

#set output dir
out_dir = "inputs/trees"
os.makedirs(out_dir, exist_ok=True)

#get inputs
af_mat = "inputs/nucleotide_variants/sample_2_ad_matrix_formated.tsv"
df = pd.read_csv(af_mat, sep="\t")

#get the name of the input file (needed for naming outputs)
infile_base = os.path.splitext(os.path.basename(af_mat))[0]

# make a new index by pasting chromosome + "_" + position
df.index = df["chromosome"].astype(str) + "_" + df["position"].astype(str)
# drop the now–redundant columns
df = df.drop(columns=["chromosome", "position"])
#drop columns with only NA values (issue with tsv format)
df = df.dropna(axis=1, how='all')

#Define a parser that turns "ref,alt" → (int, int), and handles "NA"
def parse_pair(s):
    if pd.isna(s) or s == "NA":
        return (0, 0)    # or use None/np.nan if you prefer missing
    ref, alt = s.split(",")
    return (int(ref), int(alt))

# 4)Apply it element‐wise (this converts the ad matrix to the expected format by scistree2)
ad_df = df.map(parse_pair)

#subset it for testing
df_subset = ad_df.sample(n = 10000)

#check the object
print(df_subset.head)

ad_array = np.array(
    [[[*t] for t in row] for row in df_subset.values],
    dtype=int
) 

#calculate probabilities
prob = s2.probability.genotype_probability(ad_array, ado=0.2, seqerr=0.01, posterior=True, af=None)

#construct trees with different methods

# SPR local search
caller_spr = s2.ScisTree2(threads=8)
imputed_genotype_spr, tree_spr, likelihood_spr = caller_spr.infer(prob)

# NNI local search
caller_nni = s2.ScisTree2(nni=True, threads=8)
imputed_genotype_nni, tree_nni, likelihood_nni = caller_nni.infer(prob)

# NJ
caller_nj = s2.ScisTree2(nj=True)
tree_nj= caller_nj.infer(prob)
imputed_genotype_nj, likelihood_nj = caller_nj.evaluate(prob, tree_nj)

#bundle the trees
trees = {
    "spr": tree_spr,
    "nni": tree_nni,
    "nj":  tree_nj,
}


#export the tree plots
ts = TreeStyle()
ts.show_leaf_name   = True   # display cell names at the tips
ts.show_branch_length = True # show branch lengths if you want
ts.scale =  120             # tweak to control the scale bar length

out_dir = "inputs/trees"
os.makedirs(out_dir, exist_ok=True)

for method, newick in trees.items():
    #export the tree file
    fn = f"{infile_base}_inferred_tree_{method}.nwk"
    out_file = os.path.join(out_dir, fn)
    with open(out_file, "w") as fh:
        fh.write(newick.rstrip().rstrip(";") + ";\n")
    print(f"Wrote Newick ({method}) to {out_file}")

    # export the tree plot
    from ete3 import Tree, TreeStyle
    ts = TreeStyle()
    ts.show_leaf_name     = True
    ts.show_branch_length = True
    ts.scale              = 120

    tree_obj = Tree(newick)
    img_fn = f"{infile_base}_inferred_tree_{method}.png"
    img_file = os.path.join(out_dir, img_fn)
    tree_obj.render(img_file, tree_style=ts, w=800, h=600)
    print(f"Wrote figure ({method}) to {img_file}")
