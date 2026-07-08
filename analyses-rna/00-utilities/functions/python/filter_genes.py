import pandas as pd
import scanpy as sc
import re


def get_filtered_genes_by_label(adata, celltype_col, min_frac=0.1, output_csv=None):
    """
    Filters genes per label based on minimum fraction of cells expressing the gene.

    Parameters:
    - adata: AnnData object
    - celltype_col: Column in adata.obs to group by (e.g., cell type labels)
    - min_frac: Minimum fraction of cells in group that must express a gene
    - output_csv: If provided, saves the result to this path

    Returns:
    - DataFrame of filtered genes with corresponding label
    """
    filtered_gene = pd.DataFrame()
    # Excludes: mitochondrial (MT-), ribosomal (RPS/RPL/MRPS/MRPL), hemoglobin (HB*), immunoglobulin (IGK/IGL/IGH/JCHAIN)
    pattern = re.compile(r'^(?:MT-|RPS|RPL|MRPS|MRPL|HB[^P]|IGK|IGL|IGH|JCHAIN$)')

    for label in adata.obs[celltype_col].unique():
        adata_subset = adata[adata.obs[celltype_col] == label].copy()
        min_cells = round(adata_subset.shape[0] * min_frac)
        sc.pp.filter_genes(adata_subset, min_cells=min_cells)

        gene_list = pd.DataFrame(adata_subset.var.index, columns=['gene'])

        gene_list = gene_list[~gene_list['gene'].str.contains(pattern, na=False)]
        gene_list[celltype_col] = label
        filtered_gene = pd.concat([filtered_gene, gene_list])

    if output_csv:
        filtered_gene.to_csv(output_csv, index=False)

    return filtered_gene

# filtered_gene_df = get_filtered_genes_by_label(
#     adata,
#     celltype_col='manual.label_l3',
#     min_frac=0.1,
#     output_csv='files/bmmc_l3_filtered_gene_list-no_plasma.csv'
# )
