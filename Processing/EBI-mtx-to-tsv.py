import pandas as pd
from scipy.io import mmread
import mygene
import numpy as np
import scanpy as sc
import math

print('Loaded libraries!')

path = '/projects/sysbio/users/nathanmaulding/data/raw/mtxs'

mat = mmread('/projects/sysbio/users/nathanmaulding/data/raw/mtxs/E-MTAB-7407.aggregated_filtered_counts.mtx')
print('got it!') # genes x cells
mat = mat.toarray()
print(mat.shape)

cols = pd.read_csv('/projects/sysbio/users/nathanmaulding/data/raw/mtxs/E-MTAB-7407.aggregated_filtered_counts.mtx_cols', header=None)
#print(cols)
rows = pd.read_csv('/projects/sysbio/users/nathanmaulding/data/raw/mtxs/E-MTAB-7407.aggregated_filtered_counts.mtx_rows', sep='\t', header=None)
#print(rows[0].values)
mg = mygene.MyGeneInfo()
geneinfo = mg.querymany(rows[0].values, scopes='ensembl.gene')
gene_row = []
for i in geneinfo:
    if 'symbol' in i.keys():
        gene_row.append(i['symbol'])
    else:
        gene_row.append(i['query'])


print(cols.shape, rows.shape)

tsv = pd.DataFrame(data=mat, index=gene_row, columns=cols[0].values)
# should be samples as columns and genes as rows
print(tsv)
tsv.to_csv("/projects/sysbio/users/nathanmaulding/data/raw/mtxs/E-MTAB-7407.tsv", sep="\t")



adata = sc.AnnData(tsv)
print("point 1: ", adata)
#adata = adata.transpose()
adata.var_names_make_unique()
print("point 2: ", adata)

### figure ###
# should display genes
#sc.pl.highest_expr_genes(adata, n_top=20, )
# this filters out cells that express less than 200 genes
# and genes that are only expressed in 5% of cells or less
sc.pp.filter_cells(adata, min_genes=200)
#print('test', adata.shape[0], adata.shape[1])
sc.pp.filter_genes(adata, min_cells=math.ceil(adata.shape[0] * 0.05))
print("point 3: ", adata)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True) # removed log1p=False argument because it was failing on bop
### figure ###
#sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

print("point 4: ", adata)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# set min_mean to 0.02 for saelens
sc.pp.highly_variable_genes(adata, min_mean=0.02, max_mean=3, min_disp=0.5)
### figure ###
#sc.pl.highly_variable_genes(adata)
print("point 5: ", adata)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]
print("point 6: ", adata)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
print("point 7: ", adata)
sc.pp.scale(adata, max_value=10)
print("point 8: ", adata)


#print(adata.var_names)
#print(adata.obs_names)
df = pd.DataFrame(adata.T.X)
df = df.set_index(adata.var_names)
df.columns = adata.obs_names
df = np.exp(df)  # this is the antilog of the data since we logged it in pp
print(df)
# output form adata_E-GEOD-#####_exp
df.to_csv(path + 'adata_E-MTAB-7407_exp.tsv', sep='\t')
