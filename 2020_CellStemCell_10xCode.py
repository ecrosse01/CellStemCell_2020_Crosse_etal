import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)
scv.logging.print_version()

scv.settings.set_figure_params('scvelo')

results_file="Jan2020VCS16.h5ad"
adata = sc.read(results_file)

sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

##This qc code has since been updated.
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)
             
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

adata = adata[adata.obs['n_genes'] < 3500, :]
adata = adata[adata.obs['percent_mito'] < 0.05, :]

sc.pp.normalize_total(adata, target_sum=10000)

##This should be the point that you start from .h5ad file
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

cell_cycle_genes = [x.strip() for x in open('/regev_lab_cell_cycle_genes.txt')]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
sc.pp.scale(adata)

adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')

adata.write(results_file)

#########
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

##scv.utils.show_proportions(adata) ##For velocity - not included in the paper.

##2_Scanpy clustering 
scv.pp.moments(adata, n_pcs=40, n_neighbors=8)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

scv.pp.neighbors(adata, n_neighbors=8, n_pcs=40)
scv.tl.umap(adata)
sc.tl.leiden(adata)

#########

sc.pl.umap(adata, color='leiden', legend_loc='on data')
sc.tl.leiden(adata, resolution=0.9, restrict_to=('leiden', ['3', '0', '4',]))
adata.obs['leiden_R'].cat.categories = [str(i) for i in range(len(adata.obs['leiden_R'].cat.categories))]
sc.pl.umap(adata, color='leiden_R', legend_loc='on data')
sc.pl.paga(adata, color=['leiden_R'], edge_width_scale=0.2, threshold=0.25)
adata.obs['clusters'] = adata.obs['leiden_R']

adata.obs['clusters'].cat.categories = [
'1/Endo1', '2/Endo2', '3', '17', '8', '9', '12', '15', '14/Peri', '13', '6', '11/Hem', '5',
'4', '18', '7', '19', '10/HSPC', '16', '20/PGC']

adata.obs['clusters'].cat.categories = [
'1', '2', '3', '17', '8', '9', '12', '15', '14', '13', '6', '11', '5',
'4', '18', '7', '19', '10', '16', '20']

sc.pl.umap(adata, color='clusters', legend_loc='on data')
sc.tl.paga(adata, groups='clusters')
sc.pl.paga(adata, threshold=0.25, edge_width_scale=0.2, layout='fr', random_state=0)
sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color='clusters', legend_loc='on data')
sc.pl.umap(adata, color='clusters', legend_loc='on data', palette = ['mediumseagreen', 'royalblue', 'cornflowerblue', 'rosybrown', 'burlywood', 'indianred', 'peru', 'gold', 'goldenrod', 'orange', 'thistle', 'palevioletred', 'lightsteelblue', 'mediumslateblue', 'grey', 'slateblue', 'silver', 'violet', 'darkorchid', 'forestgreen']) 

sc.pl.umap(adata, color=['CD34', 'CDH5', 'PTPRC', 'RUNX1', 'PDGFRA', 'PDGFRB', 'EPCAM', 'GYPA', 'GJA5', 'HEY2', 'APLNR', 'NT5E'], color_map='magma')
sc.pl.umap(adata, color=['NCOA7', 'SPRR2F', 'ACTA2', 'CYP26A1', 'EMX2', 'PERP', 'PTH', 'TH', 'STMN2', 'POU5F1', 'NANOS3', 'GATA4', 'FGF23', 'DLL4', 'JAG1', 'HEY2'], color_map='magma')
sc.pl.umap(adata, color=['EDN1', 'ECE1', 'EDNRA', 'EDNRB', 'AGTR2', 'ATP6AP2', 'CXCR4', 'CXCL12'], color_map='magma')
sc.pl.paga(adata, threshold=0.35, edge_width_scale=0.4, layout='fr', random_state=0)

scv.tl.velocity_embedding(adata, basis='umap')
scv.pl.velocity_embedding_stream(adata, basis='umap')
scv.tl.terminal_states(adata)
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color=['velocity_pseudotime', 'root_cells', 'end_points'], size=60)

####
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='clusters', legend_loc='on data')

sc.tl.rank_genes_groups(adata, 'clusters', method='wilcoxon')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).to_csv(
'Wilcox.csv')
sc.tl.filter_rank_genes_groups(adata)
pd.DataFrame(adata.uns['rank_genes_groups_filtered']['names']).to_csv(
'WilcoxFilt.csv')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=50, use_raw=True, swap_axes=True, show_gene_labels=False, cmap='coolwarmada', standard_scale='var')

axs = sc.pl.rank_genes_groups_dotplot(adata, groupby='clusters', n_genes=3, dendrogram='dendrogram_clusters', key='rank_genes_groups_filtered')

Genes=['PLAT','ADAMTS1', 'CCL2', 'C1QA', 'EDN1', 'EFNA1', 'CXCL2', 'LEFTY2', 'SPP1', 'NMB', 'TINAGL1', 'FRZB', 'MFAP4',
'ANGPTL6', 'FGFBP3', 'HSD11B1L'] 
gs = sc.pl.matrixplot(adata, Genes, groupby='clusters', dendrogram=True,
                      use_raw=True, cmap='bwr',  swap_axes=True, standard_scale='var')
                      
