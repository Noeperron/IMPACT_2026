#!/usr/bin/env python3
"""Supplementary Figure 2: Cell lineage annotation UMAPs and marker heatmaps.

Purpose:      Supplementary Figure 2: per-lineage cell-type annotation UMAPs + canonical marker-gene heatmaps for the integrated scRNA-seq object (B-cell, T-cell, monocyte, dendritic-cell and NK lineages).

Inputs:       Per-lineage subcluster h5ad files in SUBCLUSTER_DIR (bcell, tcell, mono, dc, nk)
              for the UMAP coordinates and subtype labels, and H5AD_NORM for the expression
              values behind the heatmaps; the subcluster objects carry no expression matrix.

Outputs:      figures/SupFig2.png (and PDF + SVG). The UMAP scatter is rasterized inside the
              vector formats, since a fully vector embedding of ~1.08M cells is unusable; all
              text, axes and heatmaps stay vector.

Dependencies: Python + scanpy, matplotlib, numpy, pandas; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.patheffects import withStroke
from matplotlib import font_manager as _fm
import glob as _glob

# Every figure in this paper is set in Arial; fall back to the metrically identical Liberation
# Sans and rewrite the SVG font-family. svg.fonttype="none" keeps text editable.
FONT = "Arial"
_av = {f.name for f in _fm.fontManager.ttflist}
if FONT not in _av:
    for _p in _glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        _fm.fontManager.addfont(_p)
    _av = {f.name for f in _fm.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams['font.family'] = PLOT_FONT
plt.rcParams['font.size'] = 8
plt.rcParams['svg.fonttype'] = 'none'

# REVISION 2026-10-08: italic statistical symbols, a×10^b instead of e-notation, transparent SVG.
import plot_style as ps


def save_figure(basename, dpi=300):
    for ext in ("png", "pdf", "svg"):
        out = FIGURES_DIR / f"{basename}.{ext}"
        # REVISION 2026-10-08: transparent SVG (no page-sized white box in Inkscape).
        plt.savefig(out, dpi=dpi, bbox_inches="tight",
                    **({"transparent": True} if ext == "svg" else {"facecolor": "white"}))
        if ext == "svg" and PLOT_FONT != FONT:
            s = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                s = s.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            s = s.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(s, encoding="utf-8")
        print(f"Saved: {out.name}")

LINEAGE_ORDER = ['bcell', 'tcell', 'mono', 'dc', 'nk']

LINEAGE_CONFIG = {
    'bcell': {
        'file': 'bcell_clean_hmy_pt_subtype.h5ad',
        'annotation_col': 'cell_subtype',
        'cell_types': ['TBC', 'NBC', 'aNBC', 'MBC', 'aMBC', 'IFN+ BC', 'ABC'],
        'markers': [
            'CD19', 'MS4A1', 'CD79A', 'CD79B',
            'CD38', 'IGHM', 'IGHD', 'TCL1A', 'IL4R', 'SELL', 'FCER2', 'CCR7',
            'FOS', 'CD83', 'CD69', 'ZNF331', 'JUNB', 'FOSB', 'NR4A2',
            'CD27', 'TNFRSF13B', 'TCF7',
            'IGHG1', 'IGHG2', 'IGHG3', 'IGHA1', 'IGHA2',
            'IFI44L', 'IFI6', 'ISG15', 'MX1',
            'FCRL5', 'ITGAX', 'TBX21', 'ZEB2'
        ],
        'scale_max': 3, 'title': 'B Cells'
    },
    'tcell': {
        'file': 'tcell_hmy_pt_subtype.h5ad',
        'annotation_col': 'cell_subtype',
        'cell_types': ['CD4+ TN', 'aTN', 'CD4+ TCM', 'Treg', 'CD8+ TN', 'GZMK+ CD8+ TEM', 'GZMB+ CD8+ TEM', 'IFN+ T', 'Tpro', 'Tgd'],
        'markers': [
            'CD3E', 'CD3D', 'CD3G', 'CD4', 'CD8A', 'CD8B',
            'TCF7', 'LEF1', 'SELL', 'CCR7',
            'CD69', 'HLA-DRA', 'NFKBIA',
            'IL7R', 'GATA3', 'CCR4', 'CCR6', 'TNFRSF4', 'RORA', 'KLRB1',
            'IL2RA', 'FOXP3', 'CTLA4', 'IKZF2',
            'GZMK', 'GZMB', 'GZMH', 'PRF1', 'GNLY', 'NKG7',
            'IFIT2', 'IFIT3', 'ISG15',
            'MKI67', 'TOP2A',
            'TRDV2', 'TRGV9', 'TRDC'
        ],
        'scale_max': 3, 'title': 'T Cells'
    },
    'mono': {
        'file': 'mono_hmy_pt_subtype.h5ad',
        'annotation_col': 'cell_subtype',
        'cell_types': ['CD14+ Mono', 'HLA-DRhi Mono', 'IFN+ CD14+ Mono', 'CD16+ Mono'],
        'rename_map': {'HLA-DR-high CD14+ Mono': 'HLA-DRhi Mono'},
        'markers': [
            'CD14', 'FCGR3A', 'CSF1R', 'ITGAM',
            'LYZ', 'S100A8', 'S100A9', 'S100A12', 'VCAN', 'FCN1',
            'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'CD74',
            'ISG15', 'IFI44L', 'IFI6', 'IFIT1', 'IFIT2', 'IFIT3',
            'MS4A7', 'CX3CR1'
        ],
        'scale_max': 3, 'title': 'Monocytes'
    },
    'dc': {
        'file': 'dc_hmy_pt_subtype.h5ad',
        'annotation_col': 'cell_subtype',
        'cell_types': ['cDC1', 'cDC2', 'Cytokine+ cDC2', 'acDC2', 'IFN+ cDC2', 'ccDC2', 'MoDC', 'Cytokine+ MoDC', 'pDC', 'AS-DC'],
        'markers': [
            'ITGAX', 'HLA-DRA', 'HLA-DRB1', 'FLT3',
            'CLEC9A', 'XCR1', 'CADM1', 'BATF3',
            'CLEC10A', 'FCER1A', 'CD1C',
            'IL1B', 'CXCL8', 'CCL3', 'TNF',
            'CD83', 'CCR7', 'LAMP3', 'FOSB',
            'ISG15', 'IFI44L', 'MX1',
            'MKI67',
            'CD14', 'S100A8', 'S100A9',
            'IL3RA', 'LILRA4', 'TCF4', 'IRF7',
            'AXL', 'SIGLEC6'
        ],
        'scale_max': 3, 'title': 'Dendritic Cells'
    },
    'nk': {
        'file': 'nk_hmy_pt_subtype.h5ad',
        'annotation_col': 'cell_subtype',
        'cell_types': ['CD56br NK', 'CD56dim NK', 'aNK', 'IFN+ NK', 'cNK'],
        'markers': [
            'NCAM1', 'KLRD1', 'NCR1', 'NKG7', 'GNLY', 'FCGR3A',
            'SELL', 'IL7R', 'XCL1', 'XCL2', 'GZMK',
            'GZMB', 'GZMH', 'PRF1', 'FGFBP2',
            'CD69', 'IFNG', 'TNF',
            'IFIT2', 'IFIT3', 'ISG15',
            'MKI67'
        ],
        'scale_max': 3, 'title': 'NK Cells'
    }
}

FULL_ADATA = None

def load_full_data():
    """Expression for the marker genes only.

    The heatmaps need about 150 of the 42,090 genes, so only those columns are retained and the
    five per-lineage slices below are then cheap. Note the PEAK is still high, around 45 GB
    observed: anndata reads the backed sparse matrix before subsetting it, so this trades a
    persistent 50 GB object for a transient one. Streaming the CSR row-blocks with h5py and
    accumulating only the wanted columns would remove the peak entirely.
    """
    global FULL_ADATA
    if FULL_ADATA is None:
        print("Loading marker-gene expression...")
        wanted = sorted({m for cfg in LINEAGE_CONFIG.values() for m in cfg['markers']})
        backed = sc.read_h5ad(H5AD_NORM, backed='r')
        keep = [g for g in wanted if g in backed.var_names]
        missing = [g for g in wanted if g not in backed.var_names]
        if missing:
            print(f"  markers absent from the object: {missing}")
        FULL_ADATA = backed[:, keep].to_memory()
        backed.file.close()
        print(f"  {FULL_ADATA.n_obs:,} cells x {FULL_ADATA.n_vars:,} marker genes")
    return FULL_ADATA


def filter_cells(adata, annotation_col):
    """Remove doublets and QC-removed cells."""
    annotations = adata.obs[annotation_col].astype(str)
    mask_doublets = ~annotations.str.startswith('db:')
    mask_qc = ~annotations.isin(['rm', 'PC', 'CLL', 'QC_removed'])
    return adata.obs_names[mask_doublets & mask_qc].tolist()


def get_color_dict(cell_types):
    n_types = len(cell_types)
    cmap = colormaps['tab10'] if n_types <= 10 else colormaps['tab20']
    return {ct: cmap(i % (10 if n_types <= 10 else 20)) for i, ct in enumerate(cell_types)}


def adjust_label_positions(centroids, cell_types, min_dist=0.08):
    """Iteratively repulse overlapping labels."""
    positions = {ct: list(centroids[ct]) for ct in cell_types if ct in centroids}
    if not positions:
        return positions

    all_x = [p[0] for p in positions.values()]
    all_y = [p[1] for p in positions.values()]
    x_range = max(all_x) - min(all_x) if max(all_x) != min(all_x) else 1
    y_range = max(all_y) - min(all_y) if max(all_y) != min(all_y) else 1

    for _ in range(50):
        moved = False
        cts = list(positions.keys())
        for i, ct1 in enumerate(cts):
            for ct2 in cts[i+1:]:
                p1, p2 = positions[ct1], positions[ct2]
                dx = (p2[0] - p1[0]) / x_range
                dy = (p2[1] - p1[1]) / y_range
                dist = np.sqrt(dx**2 + dy**2)
                if 0 < dist < min_dist:
                    push = (min_dist - dist) / 2
                    nx, ny = dx/dist, dy/dist
                    positions[ct1][0] -= nx * push * x_range * 0.5
                    positions[ct1][1] -= ny * push * y_range * 0.5
                    positions[ct2][0] += nx * push * x_range * 0.5
                    positions[ct2][1] += ny * push * y_range * 0.5
                    moved = True
        if not moved:
            break
    return positions


def plot_umap_to_ax(ax, umap_coords, annotations, cell_types, title, show_axes=True):
    color_dict = get_color_dict(cell_types)
    centroids = {}
    for label in np.unique(annotations):
        mask = annotations == label
        centroids[label] = np.median(umap_coords[mask], axis=0)

    for ct in cell_types:
        mask = annotations == ct
        if mask.sum() > 0:
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                      c=[color_dict[ct]], s=0.5, alpha=0.4, rasterized=True)

    adjusted = adjust_label_positions(centroids, cell_types, min_dist=0.12)
    for ct in cell_types:
        if ct in adjusted:
            cx, cy = adjusted[ct]
            ax.text(cx, cy, ct, fontsize=5, fontweight='bold',
                   ha='center', va='center', color='black',
                   path_effects=[withStroke(linewidth=2, foreground='white')])

    ax.set_title(f'{title} ({ps.sym("n", bold=True)}={len(annotations):,})', fontsize=7, fontweight='bold', pad=2)
    if show_axes:
        ax.set_xlabel('UMAP1', fontsize=5, labelpad=1)
        ax.set_ylabel('UMAP2', fontsize=5, labelpad=1)
    else:
        ax.set_xlabel('')
        ax.set_ylabel('')
    ax.tick_params(labelsize=4, pad=1)
    ax.set_aspect('equal', adjustable='box')


def plot_heatmap_to_ax(ax, cell_barcodes, annotations, cell_types, markers, scale_max, title):
    full_adata = load_full_data()
    available_markers = [m for m in markers if m in full_adata.var_names]
    if not available_markers:
        ax.text(0.5, 0.5, 'No markers', ha='center', va='center')
        return

    common_cells = [c for c in cell_barcodes if c in full_adata.obs_names]
    cell_idx = [full_adata.obs_names.get_loc(c) for c in common_cells]
    gene_idx = [full_adata.var_names.get_loc(g) for g in available_markers]

    expr_matrix = full_adata.X[cell_idx, :][:, gene_idx]
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()

    cell_annotations = annotations.loc[common_cells]

    mean_expr = pd.DataFrame(index=cell_types, columns=available_markers, dtype=float)
    for ct in cell_types:
        mask = cell_annotations == ct
        if mask.sum() > 0:
            mean_expr.loc[ct] = np.mean(expr_matrix[mask.values], axis=0)
        else:
            mean_expr.loc[ct] = np.nan

    zscore_expr = mean_expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=0)
    if scale_max is not None:
        zscore_expr = zscore_expr.clip(-scale_max, scale_max)

    cmap = mcolors.LinearSegmentedColormap.from_list('custom', ['steelblue', 'white', 'orangered'], N=256)
    vmax = scale_max if scale_max else zscore_expr.abs().max().max()
    im = ax.imshow(zscore_expr.values, cmap=cmap, aspect='auto', vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(len(available_markers)))
    ax.set_xticklabels(available_markers, rotation=90, fontsize=4, ha='center',
                       fontstyle='italic')
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=5)
    ax.set_title(title, fontsize=7, fontweight='bold', pad=2)

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.01, aspect=15)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label('$z$-score', fontsize=6, labelpad=1)


def load_lineage_data(lineage_name, config):
    h5ad_path = SUBCLUSTER_DIR / config['file']
    print(f"  Loading {config['title']}...")
    adata_sub = sc.read_h5ad(str(h5ad_path))
    valid_cells = filter_cells(adata_sub, config['annotation_col'])
    annotations = adata_sub.obs.loc[valid_cells, config['annotation_col']].astype(str)
    if 'rename_map' in config:
        annotations = annotations.replace(config['rename_map'])
    umap_coords = adata_sub[valid_cells].obsm['X_umap']
    print(f"    {len(valid_cells):,} cells")
    return {'valid_cells': valid_cells, 'annotations': annotations, 'umap_coords': umap_coords}


def main():
    # Load all lineage data
    lineage_data = {}
    for name in LINEAGE_ORDER:
        lineage_data[name] = load_lineage_data(name, LINEAGE_CONFIG[name])

    # Combined figure (US Letter)
    fig = plt.figure(figsize=(8.0, 10.5))
    gs = gridspec.GridSpec(5, 2, figure=fig, width_ratios=[1, 2.8],
                          height_ratios=[1]*5, hspace=0.55, wspace=0.25,
                          top=0.97, bottom=0.06, left=0.06, right=0.98)

    for row, name in enumerate(LINEAGE_ORDER):
        config = LINEAGE_CONFIG[name]
        data = lineage_data[name]
        print(f"  Plotting {config['title']}...")

        ax_umap = fig.add_subplot(gs[row, 0])
        plot_umap_to_ax(ax_umap, data['umap_coords'], data['annotations'].values,
                       config['cell_types'], config['title'], show_axes=False)
        # REVISION 2026-10-08: panel letters A-E, which the legend refers to.
        ps.panel_letter(fig, 0.004, ax_umap.get_position().y1 + 0.006, 'ABCDE'[row], fontsize=10)

        ax_heatmap = fig.add_subplot(gs[row, 1])
        plot_heatmap_to_ax(ax_heatmap, data['valid_cells'], data['annotations'],
                          config['cell_types'], config['markers'],
                          config['scale_max'], config['title'])

    save_figure("SupFig2")
    plt.close()


if __name__ == "__main__":
    main()
