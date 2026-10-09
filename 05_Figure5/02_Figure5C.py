#!/usr/bin/env python3
"""Figure 5C: IL-1B response signature analysis (pre vs post-vx, overall + per-celltype).

Purpose:      Figure 5C: IL-1B response gene signature score (BIOCARTA / Reactome-derived) in peripheral immune cells, paired pre vs post vaccination, in HD/MGUS/SMM. Scores re-normalized from the counts layer on a 2,678-HVG control pool; Wilcoxon signed-rank per group + BH correction.

Inputs:       H5AD_IL1B (the de-identified comprehensive deposit; cells that failed quality control, the Platelets cluster, platelet-containing doublets and plasma cells are excluded); data/hvg_2678_genes.txt control pool; data/il1b_response_genes_human.csv gene set.

Outputs:      figures/Figure5C.png + per-patient IL-1B response-score table.

Dependencies: Python + scanpy, pandas, numpy, matplotlib, seaborn, scipy; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from matplotlib import font_manager

# Every figure in this paper is set in Arial; fall back to the metrically identical Liberation Sans
# and rewrite the SVG font-family afterwards. svg.fonttype="none" keeps text editable.
FONT = "Arial"
_av = {f.name for f in font_manager.fontManager.ttflist}
if FONT not in _av:
    for _p in glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        font_manager.fontManager.addfont(_p)
    _av = {f.name for f in font_manager.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"
# REVISION 2026-10-08: math text ($q$, $r$, $n$) in the plot font, italic, as in the R panels.
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = PLOT_FONT
plt.rcParams["mathtext.it"] = f"{PLOT_FONT}:italic"
import seaborn as sns
from scipy import stats
import gc

# Cell scope. Excluded: cells that failed quality control (QC_removed), the "Platelets" cluster,
# any platelet-containing doublet (db:*+Platelets), and plasma cells (PC). The upstream "MK"
# cluster is platelet contamination (alpha-granule signature PPBP/PF4/TUBB1/NRGN, low UMI/n_genes,
# no CD34/MKI67/RUNX1); it was relabeled "Platelets" in the deposit and the Figure 3B legend
# states that it is excluded from downstream analyses. Plasma cells are excluded so that this
# panel scores the same circulating immune compartment as the rest of the paper; they are 1,348
# of 1,163,643 cells here and their removal leaves every reported statistic unchanged (per-patient
# scores shift by at most 4e-04). Remaining doublet categories are retained.
EXCLUDE_CELLTYPES = ['QC_removed', 'Platelets', 'CLL', 'PC']
EXCLUDE_PATTERNS = []                  # no blanket 'db:' exclusion
EXCLUDE_SUBSTRINGS = ['Platelets']     # catches db:*+Platelets

# Plot colors
COLORS_PREPOST = {'Pre-Vx': '#4878A8', 'Post-Vx': '#E05555'}

# IL-1B response gene list
IL1B_GENES_PATH = DATA_DIR / "il1b_response_genes_human.csv"
HVG_GENES_PATH = DATA_DIR / "hvg_2678_genes.txt"  # 2,678-gene analysis universe (score_genes control pool)


def filter_untreated_smm(df):
    """Keep HD, MGUS as-is; for SMM keep only untreated. Exclude MM."""
    def keep_sample(row):
        if row['Disease'] in ['HD', 'MGUS', 'PROMISE_NEG', 'IgM-MGUS']:
            return True
        elif row['Disease'] == 'SMM':
            return row['Treatment_Status'] == 0.0 or pd.isna(row['Treatment_Status'])
        return False
    return df[df.apply(keep_sample, axis=1)].copy()


def should_exclude_celltype(ct):
    if ct in EXCLUDE_CELLTYPES:
        return True
    for pattern in EXCLUDE_PATTERNS:
        if ct.startswith(pattern):
            return True
    for sub in EXCLUDE_SUBSTRINGS:
        if sub in ct:
            return True
    return False


def paired_effect_size(post, pre):
    """|Z|/sqrt(n_pairs), the same statistic rstatix::wilcox_effsize reports for panels A, B, D-F.
    Z is the tie-corrected standardized signed-rank statistic. Verified to 6 decimals against R."""
    d = np.asarray(post) - np.asarray(pre)
    d = d[d != 0]
    n = len(d)
    if n < 2:
        return np.nan
    rk = stats.rankdata(np.abs(d))
    _, cnt = np.unique(rk, return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24 - (cnt ** 3 - cnt).sum() / 48
    z = (rk[d > 0].sum() - n * (n + 1) / 4) / np.sqrt(var)
    return abs(z) / np.sqrt(n)


def compute_statistics(paired_df, disease_list):
    """Compute paired Wilcoxon statistics for each disease."""
    stats_results = []
    for disease in disease_list:
        subset = paired_df[paired_df['Disease'] == disease]
        n = len(subset)
        if n >= 3:
            pre_vals = subset['Pre-Vx']
            post_vals = subset['Post-Vx']
            try:
                stat, p_value = stats.wilcoxon(post_vals, pre_vals)
            except Exception:
                stat, p_value = np.nan, np.nan
            stats_results.append({
                'Disease': disease, 'n': n,
                'Mean_PreVx': pre_vals.mean(), 'Mean_PostVx': post_vals.mean(),
                'Mean_Diff': post_vals.mean() - pre_vals.mean(),
                'Wilcoxon_stat': stat, 'p_value': p_value,
                'effsize': paired_effect_size(post_vals, pre_vals)
            })
        else:
            stats_results.append({
                'Disease': disease, 'n': n,
                'Mean_PreVx': np.nan, 'Mean_PostVx': np.nan,
                'Mean_Diff': np.nan, 'Wilcoxon_stat': np.nan, 'p_value': np.nan,
                'effsize': np.nan
            })
    out = pd.DataFrame(stats_results)
    # Benjamini-Hochberg across the disease groups tested in this panel.
    pv = out['p_value'].to_numpy(dtype=float)
    ok = ~np.isnan(pv)
    qv = np.full(pv.shape, np.nan)
    if ok.any():
        sel = pv[ok]; m = sel.size; order = np.argsort(sel)
        run = 1.0; tmp = np.empty(m)
        for k in range(m - 1, -1, -1):
            run = min(run, sel[order[k]] * m / (k + 1))
            tmp[order[k]] = run
        qv[ok] = tmp
    out['q_value'] = qv
    return out


def generate_boxplot(paired_df, stats_df, output_path, disease_order=None):
    """Paired pre/post panel matching the layout of Figure 5A, 5B and 5D-F: one facet per disease
    group, boxes for each timepoint, lines connecting each participant, and a single bracket
    carrying the BH-corrected q-value (plus the effect size where the comparison is significant).
    No star or ns glyphs: the value is printed so the reader can judge it."""
    if disease_order is None:
        disease_order = ["HD", "MGUS", "SMM"]
    label_map = {"Healthy": "HD", "HD": "HD", "MGUS": "MGUS", "SMM": "SMM"}
    disease_order = [d for d in disease_order if d in paired_df["Disease"].unique()]

    # R 'steelblue' / 'tomato2' as used by scale_fill_manual() in Figure 5A/B and 5D-F.
    PRE, POST = "#4682B4", "#EE5C42"
    fig, axes = plt.subplots(1, len(disease_order), figsize=(6, 4), dpi=300, sharey=True)
    if len(disease_order) == 1:
        axes = [axes]
    rng = np.random.default_rng(2026)
    positions = [1, 2]

    ymin = min(paired_df["Pre-Vx"].min(), paired_df["Post-Vx"].min())
    ymax = max(paired_df["Pre-Vx"].max(), paired_df["Post-Vx"].max())
    span = max(ymax - ymin, 1e-9)

    for ax, disease in zip(axes, disease_order):
        sub = paired_df[paired_df["Disease"] == disease]
        pre, post = sub["Pre-Vx"].values, sub["Post-Vx"].values
        n = len(sub)

        # Geometry matched to the R panels: geom_boxplot at its default width 0.75, grey20 lines,
        # no caps, median fattened; geom_line size 0.3 alpha 0.5; geom_pointrange shape 21 size 0.5.
        # The R jitter is runif(-0.15, 0.15) drawn independently for every ROW, so a participant's
        # pre and post points sit at different x and the connector is slanted.
        BOX_LINE = "#333333"
        bp = ax.boxplot([pre, post], positions=positions, widths=0.75, showfliers=False,
                        showcaps=False, patch_artist=True,
                        medianprops=dict(color=BOX_LINE, linewidth=2.14),
                        whiskerprops=dict(color=BOX_LINE, linewidth=1.07))
        for patch, colour in zip(bp["boxes"], [PRE, POST]):
            patch.set_facecolor(colour); patch.set_alpha(0.5)
            patch.set_edgecolor(BOX_LINE); patch.set_linewidth(1.07); patch.set_zorder(1)

        jit_pre = rng.uniform(-0.15, 0.15, size=n)
        jit_post = rng.uniform(-0.15, 0.15, size=n)
        for i in range(n):
            ax.plot([positions[0] + jit_pre[i], positions[1] + jit_post[i]], [pre[i], post[i]],
                    color="black", linewidth=0.3, alpha=0.5, zorder=2)
        ax.scatter(positions[0] + jit_pre, pre, s=26, facecolor=PRE, edgecolor="black",
                   linewidth=0.5, zorder=3)
        ax.scatter(positions[1] + jit_post, post, s=26, facecolor=POST, edgecolor="black",
                   linewidth=0.5, zorder=3)

        row = stats_df[stats_df["Disease"] == disease]
        if len(row) and pd.notna(row["q_value"].iloc[0]):
            q_val = float(row["q_value"].iloc[0])
            eff = row["effsize"].iloc[0] if "effsize" in row else np.nan
            p_str = f"$q$={q_val:.3f}" if q_val >= 0.001 else f"$q$={q_val:.1e}"
            # effect size only where the comparison is significant, as in the Olink panels
            if q_val < 0.1 and pd.notna(eff):
                p_str += f", $r$={abs(eff):.2f}"
            bar_y = ymax + 0.08 * span
            ax.plot([positions[0], positions[0], positions[1], positions[1]],
                    [bar_y - 0.02 * span, bar_y, bar_y, bar_y - 0.02 * span],
                    color="black", linewidth=0.9)
            ax.text(1.5, bar_y + 0.03 * span, p_str, ha="center", va="bottom", fontsize=13)

        # Font sizes and the full panel border match theme_bw() in the R panels
        # (axis.text 16, axis.title 18, strip.text 18, panel.border black).
        # REVISION 2026-10-08: larger text (was 18 / 16 / 16 / 13), identical to the R panels.
        ax.set_title(f"{label_map.get(disease, disease)}\n($n$={n})", fontsize=20)
        ax.set_xticks(positions)
        ax.set_xticklabels(["Pre", "Post"], fontsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.set_xlim(0.5, 2.5)
        ax.set_ylim(ymin - 0.05 * span, ymax + 0.30 * span)
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("black")

    # REVISION 2026-10-08: one-line y-axis title at 15 pt, the size used in Figure 5A, B, D-F.
    axes[0].set_ylabel("IL-1\u03b2 response score", fontsize=15)
    plt.tight_layout()
    for ext in ("png", "pdf", "svg"):
        out = str(output_path).replace(".png", f".{ext}")
        # REVISION 2026-10-08: transparent SVG (no white background box in Inkscape).
        plt.savefig(out, dpi=300, bbox_inches="tight", transparent=(ext == "svg"))
        if ext == "svg" and PLOT_FONT != FONT:
            txt = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                txt = txt.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            txt = txt.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(txt, encoding="utf-8")
        print(f"Saved: {out}")
    plt.close()


def main():
    # Load IL1B response genes
    il1b_genes = pd.read_csv(IL1B_GENES_PATH)['gene'].tolist()
    print(f"IL1B response genes: {len(il1b_genes)}")

    # Load h5ad (comprehensive deid object: 42,090 genes + int64 counts layer + UMAP + annotations).
    # MEMORY NOTE: with X stored as csr_matrix(float64) and counts as csr_matrix(int64), the on-disk
    # ~10 GB expands to ~40-50 GB in memory. To make this run on a 58 GB VM, we (i) drop adata.X
    # (the script doesn't use it -- normalization is re-done from the counts layer), then (ii) recast
    # the counts layer to float32, halving its footprint, BEFORE the HVG subset .copy().
    print("Loading h5ad...")
    adata = sc.read_h5ad(H5AD_IL1B)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Free the existing X (it's not used; we re-normalize from counts).
    import scipy.sparse as sp
    adata.X = sp.csr_matrix(adata.X.shape, dtype='float32')
    gc.collect()
    # Recast counts to float32 (matches downstream normalize_total expectation; halves memory).
    adata.layers['counts'] = adata.layers['counts'].astype('float32')
    gc.collect()
    print("Dropped X; counts recast to float32")

    # restrict to the 2,678-gene HVG analysis universe BEFORE normalizing/scoring. This sets
    # score_genes' control-gene pool to the HVG set (scoring on all 42,090 genes degenerates). No-op if
    # the object already carries only the 2,678 genes.
    hvg = [l.strip() for l in open(HVG_GENES_PATH) if l.strip()]
    hvg = [g for g in hvg if g in adata.var_names]
    adata = adata[:, hvg].copy()
    print(f"Restricted to {adata.n_vars:,} HVG genes")

    # Re-normalize from int64 counts (Figure5C recipe; reproduces the published IL-1B result)
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    available_genes = [g for g in il1b_genes if g in adata.var_names]
    print(f"IL1B genes in data: {len(available_genes)}/{len(il1b_genes)}")

    if len(available_genes) < 50:
        print("ERROR: Too few genes found")
        return

    # Score IL1B response. use_raw=False is explicit: score on the 2,678-gene processed matrix (not a
    # .raw layer). the original run inadvertently scored on adata.raw (42,090 genes) via the
    # use_raw=None default; this corrected scoring is reproducible from the deposited object.
    sc.tl.score_genes(adata, gene_list=available_genes, score_name='IL1B_response',
                      ctrl_size=min(100, len(available_genes)), use_raw=False)

    # build the cell table from de-identified obs fields directly
    # (Diagnosis, Timepoint, TreatmentStatus, Deidentified_Patient_ID) instead of joining the
    # original-ID metadata table. Treatment_Status encoded 0=untreated (Never_treated) / 1=treated
    # so the existing filter_untreated_smm() works unchanged.
    timepoint_map = {'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx', 'Post-1st': 'Post-1st'}
    cell_df = pd.DataFrame({
        'PatientID': adata.obs['Deidentified_Patient_ID'].astype(str).values,
        'Annotation_Level_2': adata.obs['Annotation_Level_2'].astype(str).values,
        'IL1B_response': adata.obs['IL1B_response'].values,
        'Disease': adata.obs['Diagnosis'].astype(str).values,
        'Timepoint_Clean': pd.Series(adata.obs['Timepoint'].astype(str).values).map(timepoint_map).values,
        'Treatment_Status': np.where(
            adata.obs['TreatmentStatus'].astype(str).values == 'Never_treated', 0.0, 1.0),
    })
    cell_df = cell_df[cell_df['Timepoint_Clean'].notna()].copy()

    # Restrict to healthy immune cells (myeloid + lymphoid) via EXCLUDE_CELLTYPES. The deposited
    # object retains QC-failed cells, doublets, platelets (not immune) and CLL (excluded
    # paper-wide) under explicit labels, so the filter has to be applied here.
    n_before = len(cell_df)
    keep = ~cell_df['Annotation_Level_2'].map(should_exclude_celltype)
    dropped = cell_df.loc[~keep, 'Annotation_Level_2'].value_counts()
    cell_df = cell_df[keep].copy()
    print(f"Cell-type filter: kept {len(cell_df):,} of {n_before:,} cells "
          f"({100*len(cell_df)/n_before:.1f}%); excluded {n_before-len(cell_df):,}")
    for ct, k in dropped.items():
        print(f"    excluded {ct}: {k:,}")

    del adata
    gc.collect()

    # Overall analysis: aggregate per patient/timepoint
    overall_scores = cell_df.groupby(
        ['PatientID', 'Timepoint_Clean', 'Disease', 'Treatment_Status']
    )['IL1B_response'].mean().reset_index()
    overall_scores.columns = ['Patient_ID', 'Timepoint_Clean', 'Disease', 'Treatment_Status', 'IL1B_score']

    overall_filtered = filter_untreated_smm(overall_scores)
    analysis_df = overall_filtered[
        (overall_filtered['Disease'].isin(['HD', 'MGUS', 'SMM'])) &
        (overall_filtered['Timepoint_Clean'].isin(['Pre-Vx', 'Post-Vx']))
    ].copy()
    # SMM label (data already filtered to untreated only)

    # Paired analysis
    paired_df = analysis_df.pivot_table(
        index=['Patient_ID', 'Disease'], columns='Timepoint_Clean', values='IL1B_score'
    ).reset_index()
    paired_df = paired_df.dropna(subset=['Pre-Vx', 'Post-Vx'])
    print(f"Paired patients: {len(paired_df)}")

    stats_df = compute_statistics(paired_df, ['HD', 'MGUS', 'SMM'])

    for _, row in stats_df.iterrows():
        if pd.notna(row['p_value']):
            print(f"  {row['Disease']} (n={int(row['n'])}): p={row['p_value']:.4f} "
                  f"q={row['q_value']:.4f} r={row['effsize']:.3f}")

    # Plot the PAIRED frame so the panel shows exactly the participants that are tested.
    generate_boxplot(paired_df, stats_df, FIGURES_DIR / "Figure5C.png")
    print("Saved: Figure5C.png")


if __name__ == "__main__":
    main()
