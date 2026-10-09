#!/usr/bin/env python3
"""Supplementary Figure 3A: bone-marrow myeloid TNFSF13 (APRIL) expression across HD/MGUS/SMM/MM.

Purpose:      Extends the peripheral-blood myeloid TNFSF13 analysis of Figure 3C into the bone
              marrow, the site where plasma cells survive, using the CD138-depleted marrow
              scRNA-seq of Zavidij et al. (Nat Cancer 2020, GSE124310).

Pipeline:     Load the 32 GSE124310 10x v2 matrices and tag each with its disease group; QC
              (200-5,000 genes, <15% mitochondrial); normalize, log1p, HVG, scale, PCA; Leiden
              cluster and assign lineages from canonical markers; take the per-sample mean
              TNFSF13 in monocytes and dendritic cells, matching the Figure 3C aggregation;
              compare each disease group with HD.

Statistics:   Two-sided Wilcoxon rank-sum with Benjamini-Hochberg correction across the three
              contrasts. The public Zavidij release carries no age or sex metadata, so these
              contrasts are unadjusted.

Inputs:       data/external/zavidij_bm/matrices/ -- one directory per GSE124310 sample, each
              holding matrix.mtx, genes.tsv and barcodes.tsv. Download the supplementary file
              of GSE124310 from GEO and extract it there; the directory names carry the sample
              labels (NBM*, MGUS*, SMMl*, SMMh*, MM*) that set the disease groups.

Outputs:      figures/SupFig3A.png (and PDF + SVG); tables/zavidij_bm_per_sample_summary.csv,
              tables/zavidij_bm_lineage_counts.csv, tables/zavidij_bm_TNFSF13_stats.csv.

Dependencies: Python + scanpy, pandas, numpy, scipy, statsmodels, matplotlib; reads config.py.
"""

import sys, os, warnings, re
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse, stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
# REVISION 2026-10-08: italic statistical symbols, a×10^b instead of e-notation, transparent SVG.
import plot_style as ps

SAMPLES = DATA_DIR / 'external' / 'zavidij_bm' / 'matrices'
TABLES_DIR = REPO_DIR / 'tables'
TABLES_DIR.mkdir(exist_ok=True)
OUT = TABLES_DIR

# Disease-group recoder. Zavidij naming: NBM=healthy; MGUS, SMMh (high-risk), SMMl (low-risk), MM.
def disease_of(name):
    if name.startswith('NBM'): return 'HD'
    if name.startswith('MGUS'): return 'MGUS'
    if name.startswith('SMMh') or name.startswith('SMMl'): return 'SMM'
    if name.startswith('MM'): return 'MM'
    return None

# Per-sample subject ID (strip the 138N / 138N45P suffix)
def subject_of(name):
    return re.sub(r'\.138N(45P)?$', '', name)

# Lineage marker dictionary (canonical PBMC + BM)
LINEAGE_MARKERS = {
    'Mono': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN', 'FCN1'],
    'DC':   ['HLA-DRA', 'CLEC10A', 'CD1C', 'CLEC9A', 'IRF8', 'LILRA4'],
    'T':    ['CD3D', 'CD3E', 'CD3G', 'TRAC'],
    'NK':   ['NKG7', 'GNLY', 'KLRF1', 'KLRD1'],
    'B':    ['CD79A', 'MS4A1', 'CD19', 'CD79B'],
    'Eryth':['HBB', 'HBA1', 'HBA2'],
}

def load_sample(d):
    """Return AnnData for a single sample directory containing matrix.mtx, genes.tsv, barcodes.tsv."""
    X = io.mmread(str(d / 'matrix.mtx')).T.tocsr().astype(np.float32)
    genes = pd.read_csv(d / 'genes.tsv', sep='\t', header=None, names=['ensembl', 'symbol'])
    barcodes = pd.read_csv(d / 'barcodes.tsv', sep='\t', header=None, names=['barcode'])
    # Collapse duplicate symbols -> sum (deterministic; matches scanpy convention)
    if genes['symbol'].duplicated().any():
        # find dup symbols and sum their columns
        dup_mask = genes['symbol'].duplicated(keep=False)
        if dup_mask.any():
            # rare; for our purposes, take the first per symbol
            keep = ~genes['symbol'].duplicated(keep='first')
            X = X[:, keep.values]
            genes = genes[keep].reset_index(drop=True)
    a = sc.AnnData(X=X, obs=barcodes, var=genes)
    a.var.index = a.var['symbol'].astype(str)
    a.var_names_make_unique()
    a.obs.index = a.obs['barcode'].astype(str) + ('_' + d.name)
    a.obs['sample'] = d.name
    a.obs['subject'] = subject_of(d.name)
    a.obs['Disease'] = disease_of(d.name)
    return a

DEPOSITED = DATA_DIR / 'external' / 'zavidij_bm'


def main():
    if not SAMPLES.is_dir():
        # REVISION 2026-10-08: the GSE124310 matrices are not redistributed. Without them the
        # panel is redrawn from the per-sample summary and statistics this script writes, which
        # are deposited in data/external/zavidij_bm/; the plotted values are the same.
        summary = DEPOSITED / 'zavidij_bm_per_sample_summary.csv'
        if not summary.exists():
            raise SystemExit(
                f"GSE124310 sample matrices not found at {SAMPLES}.\n"
                "Download the GSE124310 supplementary file from GEO and extract it there, one "
                "directory per sample containing matrix.mtx, genes.tsv and barcodes.tsv.")
        print(f"GSE124310 matrices not found; redrawing from {DEPOSITED}")
        plot_panel(pd.read_csv(summary), pd.read_csv(DEPOSITED / 'zavidij_bm_TNFSF13_stats.csv'))
        return
    sample_dirs = sorted([p for p in SAMPLES.iterdir() if p.is_dir()])
    print(f"Loading {len(sample_dirs)} samples...")
    adatas = []
    for d in sample_dirs:
        a = load_sample(d)
        adatas.append(a)
        print(f"  {d.name}: {a.n_obs} cells x {a.n_vars} genes  ({a.obs['Disease'].iloc[0]})")
    adata = sc.concat(adatas, axis=0, join='inner', merge='same')
    print(f"\nConcatenated: {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs.groupby('Disease')['sample'].nunique())

    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    pre = adata.n_obs
    adata = adata[(adata.obs['n_genes_by_counts'] >= 200) &
                  (adata.obs['n_genes_by_counts'] <= 5000) &
                  (adata.obs['pct_counts_mt'] < 15)].copy()
    print(f"QC: {pre} -> {adata.n_obs} cells kept")

    # Save raw counts; normalize+log1p; STORE log-norm BEFORE scaling so per-sample TNFSF13
    # is computed on the same scale used in Fig 3C (Mann-Whitney/ANCOVA on log-normalized
    # mean per patient). Scaling is for clustering/PCA only.
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()

    # HVG + scale + PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.6)
    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    # Score lineages: per-cluster mean of marker-gene log-normalized expression
    # Use score_genes per lineage; assign each cluster the lineage with the highest median score.
    for lin, markers in LINEAGE_MARKERS.items():
        present = [g for g in markers if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=present, score_name=f'score_{lin}', use_raw=False)
    cluster_scores = (adata.obs.groupby('leiden')[[f'score_{lin}' for lin in LINEAGE_MARKERS]]
                      .median().reset_index())
    cluster_to_lineage = {}
    for _, r in cluster_scores.iterrows():
        best = max(LINEAGE_MARKERS.keys(),
                   key=lambda lin: r[f'score_{lin}'])
        cluster_to_lineage[r['leiden']] = best
    adata.obs['lineage'] = adata.obs['leiden'].map(cluster_to_lineage)
    print(adata.obs.groupby(['Disease', 'lineage']).size().unstack(fill_value=0))

    # Save lineage counts
    counts = (adata.obs.groupby(['sample', 'Disease', 'lineage']).size()
                       .reset_index(name='n_cells'))
    counts.to_csv(OUT / 'zavidij_bm_lineage_counts.csv', index=False)
    print("Saved lineage counts.")

    # Per-sample mean TNFSF13 in Mono+DC (myeloid)
    myeloid_mask = adata.obs['lineage'].isin(['Mono', 'DC'])
    print(f"Myeloid cells: {myeloid_mask.sum()} / {adata.n_obs}")
    if 'TNFSF13' not in adata.var_names:
        sys.exit("!! TNFSF13 not in gene list after concat")
    april_idx = list(adata.var_names).index('TNFSF13')
    # Use log-normalized layer (NOT the scaled adata.X) for per-cell TNFSF13. This matches the
    # Fig 3C aggregation rule (mean per-patient on log(normalized counts + 1)).
    L = adata.layers['lognorm']
    expr_col = L[:, april_idx]
    expr = np.asarray(expr_col).flatten() if not sparse.issparse(L) \
           else np.asarray(expr_col.todense()).flatten()
    adata.obs['TNFSF13'] = expr

    # Aggregate per sample, restricted to myeloid
    myeloid = adata.obs[myeloid_mask].copy()
    per_sample = (myeloid.groupby(['sample', 'Disease'])
                          .agg(mean_TNFSF13=('TNFSF13', 'mean'),
                               n_myeloid=('TNFSF13', 'size'),
                               pct_TNFSF13_pos=('TNFSF13', lambda v: float((v > 0).mean()*100)))
                          .reset_index())
    per_sample.to_csv(OUT / 'zavidij_bm_per_sample_summary.csv', index=False)
    print("\nPer-sample summary (myeloid):")
    print(per_sample.to_string(index=False))

    # Stats
    DISEASES = ['HD', 'MGUS', 'SMM', 'MM']
    by_d = {d: per_sample[per_sample['Disease'] == d]['mean_TNFSF13'].values for d in DISEASES}
    rows = []
    for d in ['MGUS', 'SMM', 'MM']:
        hd_vals = by_d['HD']; comp = by_d[d]
        if len(comp) == 0:
            rows.append({'cmp': f'HD vs {d}', 'n_HD': len(hd_vals), 'n': 0, 'mean_HD': np.nan,
                         'mean_other': np.nan, 'p': np.nan}); continue
        p = float(stats.mannwhitneyu(hd_vals, comp, alternative='two-sided').pvalue)
        rows.append({'cmp': f'HD vs {d}', 'n_HD': len(hd_vals), 'n': len(comp),
                     'mean_HD': float(hd_vals.mean()), 'mean_other': float(comp.mean()),
                     'p': p})
    pvals = [r['p'] for r in rows if not np.isnan(r['p'])]
    qvals = list(multipletests(pvals, method='fdr_bh')[1]) if pvals else []
    j = 0
    for r in rows:
        if not np.isnan(r['p']):
            r['q_BH'] = qvals[j]; j += 1
        else:
            r['q_BH'] = np.nan
    print("\nStats (MWU + BH across HD-vs-each-group):")
    for r in rows:
        print(f"  {r['cmp']}: n_HD={r['n_HD']} n_other={r['n']} mean_HD={r['mean_HD']:.3f} "
              f"mean_other={r['mean_other']:.3f} p={r['p']:.4g} q_BH={r['q_BH']:.4g}")
    pd.DataFrame(rows).to_csv(OUT / 'zavidij_bm_TNFSF13_stats.csv', index=False)

    plot_panel(per_sample, pd.DataFrame(rows))


def plot_panel(per_sample, stat_df):
    """Draw the panel from the per-sample summary and the HD-vs-group statistics."""
    # REVISION 2026-10-08: plotting moved out of main() so that it can also run from the
    # deposited tables when the GEO matrices are absent.
    DISEASES = ['HD', 'MGUS', 'SMM', 'MM']
    by_d = {d: per_sample[per_sample['Disease'] == d]['mean_TNFSF13'].values for d in DISEASES}
    # Figure: tighter layout for a supplementary panel; matches Fig 3C palette
    diag_colors = {'HD': '#3498db', 'MGUS': '#f1c40f', 'SMM': '#e74c3c', 'MM': '#16a085'}
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    data = [by_d[d] for d in DISEASES]
    bp = ax.boxplot(data, positions=range(len(DISEASES)), widths=0.6,
                    patch_artist=True, showfliers=False)
    for patch, d in zip(bp['boxes'], DISEASES):
        patch.set_facecolor(diag_colors[d]); patch.set_alpha(0.7)
    for el in ['whiskers', 'caps', 'medians']:
        plt.setp(bp[el], color='black', linewidth=1.2)
    for i, (d, vals) in enumerate(zip(DISEASES, data)):
        x = np.linspace(i - 0.12, i + 0.12, max(len(vals), 1))
        ax.scatter(x[:len(vals)], vals, color=diag_colors[d], edgecolor='white',
                   s=48, alpha=0.85, zorder=3, linewidth=0.6)

    # Brackets HD vs MGUS / HD vs SMM / HD vs MM, stacked cleanly with adaptive spacing
    ymax = max([v.max() for v in data if len(v) > 0])
    bracket_order = {'MGUS': 1, 'SMM': 2, 'MM': 3}
    step = max(ymax * 0.18, 0.0025)
    for _, r in stat_df.iterrows():
        g2 = r['cmp'].split(' vs ')[1]; pos = bracket_order[g2]
        y = ymax + step + (pos - 1) * step
        x2 = DISEASES.index(g2)
        ax.plot([0, x2], [y, y], 'k-', lw=1.0)
        ax.plot([0, 0], [y - step * 0.18, y], 'k-', lw=1.0)
        ax.plot([x2, x2], [y - step * 0.18, y], 'k-', lw=1.0)
        qv = r['q_BH']
        lbl = ps.stat('q', qv, digits=2 if qv >= 0.01 else 3)   # REVISION 2026-10-08
        ax.text((0 + x2) / 2.0, y + step * 0.10, lbl, ha='center', va='bottom', fontsize=9)

    ax.set_xticks(range(len(DISEASES)))
    ax.set_xticklabels([f'{d}\n($n$={len(by_d[d])})' for d in DISEASES], fontsize=11)
    # REVISION 2026-10-08: italic gene symbol; values are ln(counts per 10,000 + 1).
    ax.set_ylabel('Mean $\\mathbfit{TNFSF13}$ (APRIL) expression\nin BM myeloid cells [ln(CP10k+1)]',
                  fontsize=10, fontweight='bold')
    ax.set_title('BM myeloid APRIL (Zavidij GSE124310)', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    # y-axis: pad above the highest bracket
    ax.set_ylim(-0.001, ymax + 4 * step)
    plt.tight_layout()
    ps.save(fig, 'SupFig3A', FIGURES_DIR)   # REVISION 2026-10-08: transparent SVG
    plt.close()


if __name__ == '__main__':
    main()
