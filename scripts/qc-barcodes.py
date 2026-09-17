#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import scanpy as sc
import scipy
from pathlib import Path


cbh5 = Path(snakemake.input.cbh5)

# Import cellbender-corrected pooled bam
# Split into three based on cell barcode lists

#DATADIR=Path('/data/CARDPB2/sysbio/singlecell/STUDIES/01_AMP-AD_DiverseCohorts/DATA')
#LIBRARY='VR001'
#SAMPLEDIR = DATADIR / 'FREEMUXLET' / LIBRARY
#singlets=[SAMPLEDIR / f'{N}.singlets' for N in [0,1,2]]

adata = sc.read_10x_h5(cbh5.as_posix())
adata.var_names_make_unique()


adata.layers['counts'] = adata.X.copy()
adata.raw = adata

# Add mitochondrial and ribosomal markers
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)


# Normalize data as CPM
sc.pp.normalize_total(adata)
adata.layers['cpm']=adata.X.copy()

# Logarithmize the data
sc.pp.log1p(adata)
adata.layers['log-norm']=adata.X.copy() 


# Threshold below a given mitochondria percent
adata = adata[adata.obs['pct_counts_mt'] < snakemake.params.mito_max].copy()

# Threshold above a given genes per cell threshold
adata = adata[adata.obs['n_genes_by_counts'] > snakemake.params.genes_min].copy()

# Threshold below a given ribosome threshold
adata = adata[adata.obs['pct_counts_rb'] < snakemake.params.ribo_max].copy()

# Recompute doublet score based on filtered data
sc.pp.scrublet(adata, expected_doublet_rate=(adata.n_obs / 1000) * 0.008, threshold=snakemake.params.doublet_max, n_prin_comps=10)

good_barcodes = list(adata.obs.index)

with open(snakemake.output.good_bcs, 'w') as outfile:
    outfile.write('\n'.join(good_barcodes)+'\n')