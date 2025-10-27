#/bin/python
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import re
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
import treeswift as ts
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr
import plotly.graph_objects as go
import matplotlib.gridspec as gridspec


OKABE = ['grey', 'white',
    '#D55E00',  # reddish one (vermillion)
    '#E69F00',  # orange one
    '#F0E442',  # yellow one
    '#009E73',  # green/turquoise
    '#56B4E9',  # light blue
    '#0072B2',  # dark blue
    '#882255',  # deep magenta,
    '#CC79A7',   # light purple (magenta),
    '#332288',  # dark blue
]

sns.set_palette("deep")
plt.rcParams.update({"font.size": 14})

def build_summary_df(bm_geno_df: pd.DataFrame, lp_map_df: pd.DataFrame,
                     missing_codes=(-1,)) -> pd.DataFrame:
    bm = bm_geno_df.copy()
    lp = lp_map_df.copy()

    # Ensure numeric keys for matching
    bm["cell_name"] = bm["cell_name"].astype("str").str.strip()
    bm["target_site"] = pd.to_numeric(bm["target_site"], errors="coerce").astype("Int64")

    col_map = {c: int(re.search(r"\d+", str(c)).group()) for c in lp_map_df.columns}
    lp_wide = lp_map_df.rename(columns=col_map)
    lp_wide = lp_wide.copy()
    lp_wide.index = lp_wide.index.astype("str").str.strip()
    lp_wide.index.name = "cell_name"
    lp_long = (
        lp_wide
        .stack()
        .rename("LP_geno")
        .reset_index()
        .rename(columns={"level_1": "target_site"})
    )
    merged = bm.merge(lp_long, on=["cell_name", "target_site"], how="left")
    bm_geno = pd.to_numeric(merged["bM_geno"], errors="coerce")
    lp_geno = pd.to_numeric(merged["LP_geno"], errors="coerce")

    bm_missing = bm_geno.isna() | bm_geno.isin(missing_codes)
    lp_missing = lp_geno.isna() | lp_geno.isin(missing_codes)
    both_present = ~(bm_missing | lp_missing)

    agree = (bm_geno.to_numpy() == lp_geno.to_numpy()) & both_present.to_numpy()
    disagree_true = (~agree) & both_present.to_numpy()
    summary_df = pd.DataFrame({
        "cell_name": merged["cell_name"].to_numpy(),
        "target_site": merged["target_site"].to_numpy(),
        "bM_pmax": merged["bM_pmax"].to_numpy(),
        "bM_geno": bm_geno.to_numpy(),
        "LP_geno": lp_geno.to_numpy(),
        "bm_missing": bm_missing.to_numpy(),
        "lp_missing": lp_missing.to_numpy(),
        "both_present": both_present.to_numpy(),
        "agree": agree,
        "disagree_true": disagree_true
    })
    print(
        f"[build_summary_df] merged rows: {len(summary_df)} | "
        f"LP matches: {(~lp_geno.isna()).sum()} | "
        f"both_present: {both_present.sum()} | "
        f"agree: {summary_df['agree'].sum()} ({summary_df['agree'].sum()/len(summary_df['agree'])}) | "
        f"disagree_true: {summary_df['disagree_true'].sum()} ({summary_df['disagree_true'].sum()/len(summary_df['disagree_true'])})"
    )
    return summary_df

def plot_genotype_confidence(
    summary_df,
    title="baseMemoir probabilities\nLAML-Pro vs baseMemoir genotypes",
    outfile="tmp.pdf",
    bins=40,
    auto_range=True,
    pad_frac=0.03,
    clamp01=True,
    stat="count"
):
    # Split values
    agree_vals    = summary_df.loc[summary_df["agree"],         "bM_pmax"].dropna().to_numpy()
    disagree_vals = summary_df.loc[summary_df["disagree_true"], "bM_pmax"].dropna().to_numpy()

    # Means
    mean_agree    = float(np.mean(agree_vals))    if agree_vals.size    else np.nan
    mean_disagree = float(np.mean(disagree_vals)) if disagree_vals.size else np.nan

    # Tidy df for seaborn
    all_vals = np.concatenate([agree_vals, disagree_vals]) if (agree_vals.size or disagree_vals.size) else np.array([])
    df_plot = pd.DataFrame({
        "bM_pmax": all_vals,
        "Status":  (["Agree"] * agree_vals.size) + (["Disagree"] * disagree_vals.size)
    })

    # Dynamic x-range
    if auto_range and all_vals.size:
        xmin, xmax = np.min(all_vals), np.max(all_vals)
        if xmin == xmax:
            xmin, xmax = xmin - 0.01, xmax + 0.01
        span = xmax - xmin
        xmin -= span * pad_frac
        xmax += span * pad_frac
        if clamp01:
            xmin = max(0.0, xmin)
            xmax = min(1.0, xmax)
        rng = (float(xmin), float(xmax))
    else:
        rng = (0.0, 1.05) if clamp01 else None

    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)

    # Palette: Agree=red, Disagree=blue (so dashed handles are red/blue)
    palette = {"Agree": "C3", "Disagree": "C0"}

    # Hist with gray borders; suppress default legend
    sns.histplot(
        data=df_plot, x="bM_pmax", hue="Status", ax=ax,
        bins=bins, binrange=rng, stat=stat, common_bins=True,
        multiple="layer", palette=palette, alpha=0.5,
        edgecolor="gray", linewidth=0.7, legend=False
    )

    # KDE lines scaled to counts
    if stat == "counts":
        if all_vals.size:
            if rng is None:
                xlim = ax.get_xlim()
                binwidth = (xlim[1] - xlim[0]) / bins
                xs = np.linspace(*xlim, 400)
            else:
                binwidth = (rng[1] - rng[0]) / bins
                xs = np.linspace(*rng, 400)

            if agree_vals.size > 1:
                kde_a = gaussian_kde(agree_vals)
                ax.plot(xs, kde_a(xs) * agree_vals.size * binwidth, color=palette["Agree"], linewidth=2)
            if disagree_vals.size > 1:
                kde_d = gaussian_kde(disagree_vals)
                ax.plot(xs, kde_d(xs) * disagree_vals.size * binwidth, color=palette["Disagree"], linewidth=2)

    # Mean markers (dashed)
    if not np.isnan(mean_agree):
        ax.axvline(mean_agree, linestyle="--", linewidth=2, color=palette["Agree"])
    if not np.isnan(mean_disagree):
        ax.axvline(mean_disagree, linestyle="--", linewidth=2, color=palette["Disagree"])

    # ---- Custom dashed-line legend (main request) ----
    handles = []
    labels  = []
    if not np.isnan(mean_agree):
        handles.append(Line2D([0], [0], color=palette["Agree"], lw=2, ls="--"))
        labels.append(f"Mean Agree ({mean_agree:.2f})")
    if not np.isnan(mean_disagree):
        handles.append(Line2D([0], [0], color=palette["Disagree"], lw=2, ls="--"))
        labels.append(f"Mean Disagree ({mean_disagree:.2f})")
    ax.legend(handles, labels, loc="upper left", title=None, frameon=True)

    ax.set_xlabel("baseMemoir probabilities")
    #ax.set_ylabel("Count")
    ax.set_title(title)

    if rng is not None:
        ax.set_xlim(rng[0], rng[1] + 0.05)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
    return outfile



def clustermap_genos(
    bm_chars: pd.DataFrame,
    lp_chars: pd.DataFrame,
    metric: str = "hamming",
    method: str = "average",
    title: str = "Colony", 
    other: str = "baseMemoir",
    vmin: int = -1, vmax: int = 3,
    outfile: str | None = None,
    n_col_clusters: int = 6,          # how many column clusters to group to the left
    group_blocks_left: bool = True,   # turn grouping on/off
    show_cell_names: bool = True
):
    # 1) align rows/cols
    common_rows = bm_chars.index.intersection(lp_chars.index)
    common_cols = bm_chars.columns.intersection(lp_chars.columns)
    bm = bm_chars.loc[common_rows, common_cols].apply(pd.to_numeric, errors="coerce").fillna(-1)
    lp = lp_chars.loc[common_rows, common_cols].apply(pd.to_numeric, errors="coerce").fillna(-1)
    
    # 2) get LP clustering 
    g_lp = sns.clustermap(
        lp, row_cluster=True, col_cluster=True,
        metric=metric, method=method,
        xticklabels=False, yticklabels=False, cbar=False
    )
    row_order = g_lp.dendrogram_row.reordered_ind
    col_order = g_lp.dendrogram_col.reordered_ind
    Zc = g_lp.dendrogram_col.linkage  # column linkage for grouping
    plt.close(g_lp.fig)               # close the temporary figure

    # reorder matrices by LP order
    lp_ord = lp.iloc[row_order, :].iloc[:, col_order]
    bm_ord = bm.iloc[row_order, :].iloc[:, col_order]

    categories = np.arange(vmin, vmax + 1)
    print("Num categories:", categories)
    cmap   = ListedColormap(OKABE[:len(categories)])
    cmap.set_bad(OKABE[0])  # manually set this to grey                   
    bounds = np.arange(vmin - 0.5, vmax + 1.5, 1.0)
    norm   = BoundaryNorm(bounds, cmap.N)
    
    # side-by-side heatmaps 
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5), dpi=150, constrained_layout=True)
    # BaseMemoir (left)
    if show_cell_names:
        sns.heatmap(
            bm_ord, ax=axes[1], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            #bm_ord, mask=bm_ord==-1, ax=axes[0], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            xticklabels=False, yticklabels=True, cbar=False, rasterized=False,
            # linewidths=0.05, linecolor="white"  
        )
    else:
        sns.heatmap(
            bm_ord, ax=axes[1], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            #bm_ord, mask=bm_ord==-1, ax=axes[0], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            xticklabels=False, yticklabels=False, cbar=False, rasterized=False,
            # linewidths=0.05, linecolor="white"  
        )

    axes[1].set_title(f"{title}: {other} genotypes", fontsize=16)
    axes[1].set_xlabel("Target sites", fontsize=13)
    axes[1].set_ylabel("Cell names", fontsize=13)
    axes[1].tick_params(axis="both", labelsize=11)

    # LAML-Pro (right) 
    sns.heatmap(
        lp_ord, ax=axes[0], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,   
        #lp_ord, mask=lp_ord==-1, ax=axes[1], cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,   
        xticklabels=False, yticklabels=False, cbar=False, rasterized=False,
        # linewidths=0.05, linecolor="white"
    )
    axes[0].set_title(f"{title}: LAML-Pro genotypes", fontsize=16)
    axes[0].set_xlabel("Target sites", fontsize=13)
    axes[0].set_ylabel("") 
    axes[0].tick_params(axis="both", labelsize=11)

    # Black border around each heatmap
    for ax in axes:
        for side in ("left", "right", "top", "bottom"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_linewidth(1.5)
            ax.spines[side].set_edgecolor("black")

    # Shared colorbar
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes, orientation="horizontal", fraction=0.08, pad=0.05
    )
    cbar.set_ticks(categories) #["?/-1"] + [str(x) for x in categories[1:]))
    ticklabels = ["?/-1" if x == -1 else str(x) for x in categories]
    cbar.set_ticklabels(ticklabels)

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=150)

    return fig

def distdict_to_df(distdict):
    # assuming leaf label dict (from treeswift)
    nodes = list(set(distdict.keys()))

    df = pd.DataFrame(np.nan, index=nodes, columns=nodes, dtype=float)
    np.fill_diagonal(df.values, 0.0)

    for u, inner in distdict.items():
        for v, d in inner.items():
            df.loc[u, v] = d
            df.loc[v, u] = d
    return df

def leaf_pairs(df):
    df = df.copy()
    # row order should match col order
    m = df.loc[df.index, df.index].astype(float).copy()
    # ignore self by setting diagonal to pos infty
    np.fill_diagonal(m.values, np.inf)
    nearest = m.idxmin(axis=1)
    dmin = m.min(axis=1)
    return [(i, j, float(dmin.loc[i])) for i, j in zip(m.index, nearest)]

def get_geno_dict(df):
    sites = sorted(df['target_site'].astype(int).unique().tolist())
    d = (df.sort_values('bM_pmax', ascending=False)
           .drop_duplicates(['cell_name', 'target_site']))

    # (cell, site) -> [AA, GG, GA, AG]
    ser = (d.set_index(['cell_name','target_site'])
             [['bMp_AA','bMp_GG','bMp_GA','bMp_AG']]
             .apply(list, axis=1))

    # dict: cell -> vector over sites ([] if missing)
    out = {}
    for cell, sub in ser.groupby(level=0):
        m = sub.droplevel(0).to_dict()
        out[cell] = [m.get(s, []) for s in sites]
    return out

def _as_arr(x):
    # if site is non-missing, return numpy array. else None. 
    return None if (x is None or len(x) == 0) else np.asarray(x, float)

def im_ehd(P, Q):
    # expected hamming distance (ignoring missing)
    dots, count = 0.0, 0
    for p_i, q_i in zip(P, Q):
        p = _as_arr(p_i); q = _as_arr(q_i)
        if (p is not None) and (q is not None):
            dots += 1.0 - float(np.dot(p, q))
            count += 1
    return dots / count if count > 0 else np.nan

def empirical_site_dists(genodict): 
    # if empty everywhere, treat as uniform
    N = max(len(v) for v in genodict.values())
    K = next((len(x) for v in genodict.values() for x in v if x), 4)
    recs = []
    for vec in genodict.values():
        for i, x in enumerate(vec):
            if x and len(x) == K:
                for k, p in enumerate(x):
                    recs.append((i, k, float(p)))
    df = pd.DataFrame(recs, columns=["site", "k", "p"])
    if df.empty:
        base = np.ones(K)/K 
        return [base.copy() for _ in range(N)]

    mat = (df.groupby(["site", "k"])["p"].mean()
             .unstack("k"))
    mat = mat.reindex(range(N))       # ensure all sites present
    mat = mat.fillna(0.0)             # treat unobserved k as 0 prob before renorm

    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums.replace(0, np.nan), axis=0)
    mat = mat.fillna(1.0/K)

    return [mat.loc[i].to_numpy(dtype=float) for i in range(N)]


    any_vec = next(iter(genodict.values()))
    N = len(any_vec)
    hat = []
    for i in range(N):
        obs = []
        for v in genodict.values():
            arr = _as_arr(v[i])
            if arr is not None:
                obs.append(arr)
        if obs:
            p = np.mean(obs, axis=0)
            s = p.sum()
            hat.append(p/s if s>0 else p)

def ehd(P, Q, hat):
    # expected hamming distance (with missing filled in with empirical)
    N = len(hat)
    total = 0.0
    for p_i, q_i, h_i in zip(P, Q, hat):
        p = _as_arr(p_i); q = _as_arr(q_i); h = np.asarray(h_i, float)
        if p is None: p = h
        if q is None: q = h
        total += 1.0 - float(np.dot(p, q))
    return total / N

def sm_ehd(P, Q, hat):
    # shared-missing expected hamming distance
    N = len(hat)
    total = 0.0
    for p_i, q_i, h_i in zip(P, Q, hat):
        p = _as_arr(p_i); q = _as_arr(q_i); h = np.asarray(h_i, float)
        if (p is not None) and (q is not None):
            A = float(np.dot(p, q))
        elif (p is None) and (q is None):
            A = 0.5 + 0.5 * float(np.dot(h, h))
        elif p is None:
            A = float(np.dot(h, q))
        else:  # q is None
            A = float(np.dot(p, h))
        total += (1.0 - A)
    return total / N

def pair_metrics(pairs, G, hat):
    rows = []
    for a, b, dist in pairs:
        a, b = int(a), int(b)
        P, Q = G[a], G[b]
        rows.append((a, b, dist, im_ehd(P, Q), sm_ehd(P, Q, hat), ehd(P, Q, hat)))
    return pd.DataFrame(rows, columns=["leaf1","leaf2","distance","im_ehd","sm_ehd","ehd"])

def _min_max_random(G, hat, metric_func, seed=0):
    rng = np.random.default_rng(seed)
    leaves = [str(k) for k in G.keys()]
    # ensure order and indexable via strings
    GG = {str(k): v for k, v in G.items()}

    mins, maxs, rnds = [], [], []
    for a in leaves:
        Pa = GG[a]
        best_min, best_max = np.inf, -np.inf
        # pick one random partner (not self)
        b_rand = rng.choice([x for x in leaves if x != a])
        Pb = GG[b_rand]
        rnds.append(float(metric_func(Pa, Pb, hat)))

        for b in leaves:
            if b == a:
                continue
            d = float(metric_func(Pa, GG[b], hat))
            if d < best_min: best_min = d
            if d > best_max: best_max = d
        mins.append(best_min)
        maxs.append(best_max)
    return np.array(mins), np.array(maxs), np.array(rnds)

def _per_leaf_from_pairs(df, metric):
    # one value per leaf (use the row where it appears as leaf1)
    s = (df.drop_duplicates('leaf1')
           .set_index('leaf1')[metric]
           .astype(float))
    s.index = s.index.astype(str)
    return s


def plot_concordance_distribution(bm_df, lp_df, genodict, hat, figsize=(12,4), pad_frac=0.06, outfile=None):
    metric_funcs = {
            "EHD":       lambda P, Q, H: ehd(P, Q, H),
            "Ignore-missing EHD": lambda P, Q, H: im_ehd(P, Q),
            "Shared-missing EHD": lambda P, Q, H: sm_ehd(P, Q, H),
        }

    panels = {}
    for name, fn in metric_funcs.items():
        colmap = {"EHD": "ehd", "Ignore-missing EHD": "im_ehd", "Shared-missing EHD": "sm_ehd"}
        col = colmap[name]
        pub = _per_leaf_from_pairs(bm_df, col)
        lpl = _per_leaf_from_pairs(lp_df, col)

        # Min / Largest / Random by scanning all leaves with this metric
        mins, maxs, rnds = _min_max_random(genodict, hat, fn)

        panels[name] = {
            "LAML-Pro":   lpl.reindex(sorted(lpl.index)).to_numpy(),
            "Published":  pub.reindex(sorted(pub.index)).to_numpy(),
            "Min Distance": mins,
            "Largest":      maxs,
            "Random":       rnds
        }

    palette = OKABE[2:]
    order = ["LAML-Pro", "Largest", "Min Distance", "Published", "Random"]
    colors  = dict(zip(order, palette[:len(order)]))

    allvals_global = np.concatenate([
        np.asarray(series[label], float).ravel()
        for series in panels.values() for label in order
    ])
    allvals_global = allvals_global[np.isfinite(allvals_global)]
    lo, hi = (allvals_global.min(), allvals_global.max()) if allvals_global.size else (0.0, 1.0)
    pad = 0.08*(hi - lo) if hi > lo else 0.08        # widen as desired
    xlim = (lo - pad, hi + pad)

    fig, axes = plt.subplots(1, 3, figsize=figsize, dpi=150, sharex=True, sharey=True)
    for ax, (title, series) in zip(axes, panels.items()):
        medians = {}
        for label in order:
            vals = np.asarray(series[label], float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            sns.kdeplot(vals, ax=ax, lw=2, label=label, color=colors[label])
            mval = float(np.median(vals))
            medians[label] = mval
            ax.axvline(mval, ls="--", color=colors[label], alpha=0.9)
        # normalized medians: (median - d_min)/(d_max - d_min)
        d_min = medians.get("Min Distance", np.nan)
        d_max = medians.get("Largest", np.nan)
        denom = (d_max - d_min) if (np.isfinite(d_max) and np.isfinite(d_min)) else np.nan

        def _norm(x):
            return (x - d_min) / denom if (np.isfinite(x) and np.isfinite(denom) and denom > 0) else np.nan

        lp_med  = medians.get("LAML-Pro",  np.nan)
        pub_med = medians.get("Published", np.nan)
        lp_norm  = _norm(lp_med)
        pub_norm = _norm(pub_med)

        # annotate
        ax.text(0.02, 0.96,
                f"med(LP) = {lp_med:.3f}  (norm {lp_norm:.3f})\n"
                f"med(Pub)= {pub_med:.3f}  (norm {pub_norm:.3f})",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=10, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

        ax.set_title(title)
        ax.set_xlabel("Distance") 
        ax.grid(False)

    axes[0].set_ylabel("Density")
    # shared legend centered below
    handles, labels = axes[0].get_legend_handles_labels()
    median_handle = Line2D([0], [0], color="k", lw=1.5, ls="--", label="Median")
    fig.legend(handles + [median_handle], labels + ["Median"],
               ncol=len(order)+1, loc="lower center", bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0,0.06,1,1])
    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.show()
    plt.close(fig)

def plot_concordance_scatterplot(bm_concordance, lp_concordance, figsize=(12,4), pad_frac=0.06, outfile=None):
    # join by leaf1
    b = bm_concordance.drop_duplicates('leaf1').set_index(bm_concordance['leaf1'].astype(str))
    l = lp_concordance.drop_duplicates('leaf1').set_index(lp_concordance['leaf1'].astype(str))
    df = b.join(l, lsuffix="_bm", rsuffix="_lp", how="inner")

    metas = [("ehd", "EHD"),
             ("im_ehd", "Ignore-missing EHD"),
             ("sm_ehd", "Shared-missing EHD")]

    # global limits (shared axes)
    x_all = np.concatenate([df[f"{m}_bm"].to_numpy() for m,_ in metas])
    y_all = np.concatenate([df[f"{m}_lp"].to_numpy() for m,_ in metas])
    mask = np.isfinite(x_all) & np.isfinite(y_all)
    lo, hi = (x_all[mask].min(), x_all[mask].max()) if mask.any() else (0.0, 1.0)
    lo = min(lo, y_all[mask].min()) if mask.any() else lo
    hi = max(hi, y_all[mask].max()) if mask.any() else hi
    pad = pad_frac*(hi-lo) if hi > lo else pad_frac
    xlim = (lo - pad, hi + pad)
    ylim = (lo - pad, hi + pad)

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=figsize, dpi=150)

    for ax, (m, title) in zip(axes, metas):
        x = df[f"{m}_bm"].to_numpy()
        y = df[f"{m}_lp"].to_numpy()

        for xi, yi in zip(x, y):
            if not (np.isfinite(xi) and np.isfinite(yi)): 
                continue
            if yi > xi:  # red vertical
                ax.vlines(xi, xi, yi, color="C3", lw=2)
            elif xi > yi:  # blue horizontal
                ax.hlines(yi, yi, xi, color="C0", lw=2)

        ax.scatter(x, y, s=18, color="C0", edgecolors="k", linewidths=0.3, zorder=3)
        ax.plot(xlim, ylim, "k--", lw=1)  # y = x
        ax.set_xlim(xlim); ax.set_ylim(ylim)

        red_total  = float(np.nansum(np.maximum(y - x, 0)))
        blue_total = float(np.nansum(np.maximum(x - y, 0)))
        ax.text(0.03, 0.97, f"Red total: {red_total:.3f}\nBlue total: {blue_total:.3f}",
                transform=ax.transAxes, ha="left", va="top")
        ax.set_title(title)
        ax.set_xlabel("baseMemoir probabilities:\nPublished sibling")

    axes[0].set_ylabel("baseMemoir probabilities:\nLAML-Pro sibling")
    axes[0].legend([Line2D([0],[0], ls="--", color="k")], ["y = x"], loc="lower right")

    plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=150)
    
    plt.show()
    plt.close(fig)

STATE_ORDER  = [-1, 0, 1, 2, 3]
STATE_LABELS = ['-1', 'AA', 'GG', 'GA', 'AG']

def plot_state_counts(geno_matrix, lp_map_geno_df, title="BaseMemoir vs. LAML-Pro", ylabel="baseMemoir state", outfile=None, use_state_labels=True):
    # 1) align rows/cols
    rows = geno_matrix.index.intersection(lp_map_geno_df.index)
    cols = geno_matrix.columns.intersection(lp_map_geno_df.columns)
    bm = geno_matrix.loc[rows, cols].apply(pd.to_numeric, errors="coerce")
    lp = lp_map_geno_df.loc[rows, cols].apply(pd.to_numeric, errors="coerce")

    # 2) flatten pairs
    pairs = pd.DataFrame({
        "bm": bm.to_numpy().ravel(),
        "lp": lp.to_numpy().ravel()
    }).dropna()
    pairs["bm"] = pairs["bm"].astype(int)
    pairs["lp"] = pairs["lp"].astype(int)

    # 3) counts table with fixed order
    bm_cat = pd.Categorical(pairs["bm"], ordered=True) #categories=STATE_ORDER, ordered=True)
    lp_cat = pd.Categorical(pairs["lp"], ordered=True) #categories=STATE_ORDER, ordered=True)
    counts = pd.crosstab(bm_cat, lp_cat)
    #counts = pd.crosstab(bm_cat, lp_cat).reindex(index=STATE_ORDER, columns=STATE_ORDER, fill_value=0)
    if use_state_labels:
        counts.index   = STATE_LABELS
        counts.columns = STATE_LABELS

    # 4) number of "changes" (off-diagonal, both present)
    true_changes = pairs[(pairs["bm"] != pairs["lp"]) & (pairs["bm"] != -1) & (pairs["lp"] != -1)].shape[0]

    # 5) plot
    plt.figure(figsize=(6.5, 6), dpi=150)
    ax = sns.heatmap(counts, annot=True, fmt='g', cmap=sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True), cbar=True,
                     linewidths=.5, linecolor='white', annot_kws={"fontsize": 8})
    ax.set_title(title, pad=12)
    ax.set_xlabel("LAML-Pro state")
    ax.set_ylabel(ylabel)
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.show()

    return counts


def report_genotype_call_stats(counts, states=['-1','AA','GG','GA','AG'], other_method="baseMemoir", colony=None, report_pretty=True):
    C = counts.reindex(index=states, columns=states, fill_value=0).astype(float)

    N_total = C.values.sum()
    if N_total == 0:
        raise ValueError("Counts matrix is empty.")

    # baseMemoir proportions over all sites
    bm_all = (C.sum(axis=1) / N_total).reindex(states)

    # Split by whether baseMemoir state is missing (-1) or observed (!=-1)
    obs_rows = [s for s in states if s != '-1']
    miss_rows = ['-1']

    # LAML-Pro (obs) over all sites  -> only rows where BM != -1, normalize by N_total
    lp_obs_all = (C.loc[obs_rows].sum(axis=0) / N_total).reindex(states, fill_value=0)

    # LAML-Pro (impute) over all sites -> only rows where BM == -1, normalize by N_total
    lp_imp_all = (C.loc[miss_rows].sum(axis=0) / N_total).reindex(states, fill_value=0)

    # LAML-Pro (obs) over observed sites -> normalize by observed total; '-1' is N/A
    N_obs = C.loc[obs_rows].values.sum()
    lp_obs_over_obs = (C.loc[obs_rows].sum(axis=0) / N_obs).reindex(states, fill_value=np.nan)
    lp_obs_over_obs['-1'] = np.nan  # not applicable

    # LAML-Pro (impute) over missing sites -> normalize by missing total
    N_miss = C.loc[miss_rows].values.sum()
    lp_imp_over_miss = (C.loc[miss_rows].sum(axis=0) / N_miss).reindex(states, fill_value=0)

    rows = pd.DataFrame([
        bm_all, lp_obs_all, lp_imp_all, lp_obs_over_obs, lp_imp_over_miss
    ], index=[
        other_method,
        "LAML-Pro (obs) prop. over all sites",
        "LAML-Pro (impute) prop. over all sites",
        "LAML-Pro (obs) prop. over observed sites",
        "LAML-Pro (impute) prop. over missing sites",
    ])

    if report_pretty:

        cols_pretty = {
            '-1': "Missing/\nSilenced",
            'AA': "Unedited\n(AA)",
            'GG': "GG", 'GA': "GA", 'AG': "AG"
        }
        rows = rows.rename(columns=cols_pretty)

    rows["Total"] = rows.sum(axis=1, skipna=True)

    if colony is not None:
        rows.insert(0, "Colony", colony)

    return rows.round(4)

def save_df_to_pdf(df, filename="table.pdf", title=None, floatfmt="{:.4f}", fontsize=10):
    # format numbers (optional)
    df_disp = df.copy()
    df_disp = df_disp.applymap(lambda x: floatfmt.format(x) if isinstance(x, (int, float, np.floating)) else x)

    n_rows, n_cols = df_disp.shape
    # figure size scales with table size (tweak multipliers if needed)
    fig_w = min(18, max(6, 1.2 * n_cols))
    fig_h = min(24, max(2.5, 0.5 * n_rows + (1.0 if title else 0)))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=fontsize+2, pad=12)

    tbl = ax.table(cellText=df_disp.values,
                   colLabels=df_disp.columns.tolist(),
                   rowLabels=df_disp.index.tolist(),
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1, 1.2)  # widen/heighten cells

    fig.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    return filename



def plot_correlation(centroids_df, bm_distmat, lp_distmat,
                                    title="Colony 2 Leaf-pair distances: phylogenetic vs spatial",
                                    outfile=None, figsize=(10,4), titles=["baseMemoir", "LAML-Pro"]):
    # --- spatial distance matrix from cell centroids ---
    coords = centroids_df.iloc[:, :2].astype(float).copy()
    coords.columns = ["x", "y"]
    labels = coords.index.astype(str)
    spatial = pd.DataFrame(squareform(pdist(coords.values, metric="euclidean")),
                           index=labels, columns=labels)

    # --- make sure distance matrices use string labels & align all three ---
    bm = bm_distmat.copy(); lp = lp_distmat.copy()
    for M in (bm, lp):
        M.index = M.index.astype(str); M.columns = M.columns.astype(str)

    common = spatial.index.intersection(bm.index).intersection(lp.index)
    spatial = spatial.loc[common, common]
    bm = bm.loc[common, common]
    lp = lp.loc[common, common]

    # use only unique unordered pairs (upper triangle without diagonal)
    n = len(common)
    iu = np.triu_indices(n, k=1)
    y = spatial.to_numpy()[iu]
    x_bm = bm.to_numpy()[iu]
    x_lp = lp.to_numpy()[iu]

    # shared limits
    xs = np.concatenate([x_bm, x_lp]); ys = y
    lo_x, hi_x = xs.min(), xs.max()
    lo_y, hi_y = ys.min(), ys.max()
    pad_x = 0.05*(hi_x-lo_x) if hi_x > lo_x else 0.05
    pad_y = 0.07*(hi_y-lo_y) if hi_y > lo_y else 0.07

    # --- plotting ---
    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=figsize, dpi=150)
    for ax, (xvals, subtitle) in zip(axes, [(x_bm, titles[0]), (x_lp, titles[1])]):
        sns.regplot(x=xvals, y=y, ax=ax,
                    scatter_kws=dict(s=22, alpha=0.6),
                    line_kws=dict(lw=2, alpha=0.9))
        r, p = pearsonr(xvals, y)
        rho, ps = spearmanr(xvals, y)
        ax.text(0.03, 0.97,
                f"Pearson r = {r:.3f} (p={p:.1e})\n"
                f"Spearman ρ = {rho:.3f} (p={ps:.1e})",
                transform=ax.transAxes, va="top", ha="left",
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
        ax.set_title(subtitle)
        ax.set_xlabel("Phylogenetic distance")

    axes[0].set_ylabel("Spatial distance")
    axes[0].set_xlim(lo_x - pad_x, hi_x + pad_x)
    axes[0].set_ylim(lo_y - pad_y, hi_y + pad_y)

    if title:
        fig.suptitle(title, y=1.02, fontsize=16)
        plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=150)
    plt.show()
    return fig, axes

def plot_bl_variance(bm_tree_df, lp_tree_df, title, outfile):
    fig = plt.figure(figsize=(3.6, 6), dpi=150)
    gs  = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=0.05)
    ax_hist = fig.add_subplot(gs[0])
    ax      = fig.add_subplot(gs[1], sharex=ax_hist)

    # scatter (x=branch length, y=cumulative length)
    ax.scatter(bm_tree_df["branch_len"], bm_tree_df["cum_len"], s=22, color="C0", alpha=0.85,
               label=f"baseMemoir (var={np.var(bm_tree_df['cum_len']):.4g})")
    ax.scatter(lp_tree_df["branch_len"], lp_tree_df["cum_len"], s=22, color="C1", alpha=0.85,
               label=f"LAML-Pro (var={np.var(lp_tree_df['cum_len']):.4g})")

    ax.set_xlabel("Branch length")
    ax.set_ylabel("Cumulative branch length")
    ax.legend(loc="lower right", fontsize=8, frameon=False)

    bins = 50
    ax_hist.hist(bm_tree_df["branch_len"], bins=bins, color="C0", alpha=0.5, edgecolor="none")
    ax_hist.hist(lp_tree_df["branch_len"], bins=bins, color="C1", alpha=0.5, edgecolor="none")
    ax_hist.tick_params(axis="x", labelbottom=False)
    ax_hist.set_ylabel("Count")
    ax_hist.spines["right"].set_visible(False)
    ax_hist.spines["top"].set_visible(False)

    # shared x-limits (pad a bit)
    xmin = min(bm_tree_df["branch_len"].min(), lp_tree_df["branch_len"].min())
    xmax = max(bm_tree_df["branch_len"].max(), lp_tree_df["branch_len"].max())
    pad  = 0.05 * (xmax - xmin) if xmax > xmin else 0.5
    ax.set_xlim(xmin - pad, xmax + pad)
    if title:
        fig.suptitle(title, y=1.02, fontsize=16)
        plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=150)
    return fig

def branch_table(tree):
    rows = []
    for u in tree.traverse_preorder():
        bl = u.get_edge_length()
        height = tree.distance_between(tree.root, u)
        rows.append({"node": u.label, "branch_len": bl, "cum_len": height})
    return pd.DataFrame(rows)


def add_internal_labels(tree):
    i = 0
    for node in tree.traverse_postorder():
        if node.label is None:
            node.label = f"internal_{i}"
            i += 1
    return tree



def plot_tree_3d(ts_tree, merged_df, x_col="x", y_col="y", title="Tree in 3D (x, y, cumulative branch length)",
                 outfile=None):
    """
    ts_tree: a treeswift.Tree object
    merged_df: pandas DataFrame with index=leaf names, columns [x_col, y_col]
    """

    # ---- 1) Basic sanity checks
    if not isinstance(merged_df.index, pd.Index):
        raise ValueError("merged_df must have an index of leaf names.")
    for c in (x_col, y_col):
        if c not in merged_df.columns:
            raise KeyError(f"Column '{c}' not found in merged_df")

    # ---- 2) Utilities for TreeSwift that are robust to attribute name variants
    def is_leaf(node):
        return getattr(node, "is_leaf", lambda: len(get_children(node)) == 0)()

    def get_children(node):
        return getattr(node, "children", getattr(node, "get_children", lambda: [])()) or []

    def get_parent(node):
        return getattr(node, "parent", getattr(node, "get_parent", lambda: None)())

    def get_label(node):
        return getattr(node, "label", getattr(node, "get_label", lambda: None)())

    def get_edge_length(node):
        # edge length is from node to its parent; root may have None
        el = getattr(node, "edge_length", getattr(node, "get_edge_length", lambda: None)())
        return 0.0 if el is None else float(el)

    root = getattr(ts_tree, "root", getattr(ts_tree, "get_root", lambda: None)())
    if root is None:
        raise ValueError("Could not determine root of the tree.")

    # ---- 3) Gather leaf coords from merged_df
    # Allow numeric strings vs ints: coerce index to string for matching
    leaf_xy = {}
    df_index_as_str = merged_df.copy()
    df_index_as_str.index = df_index_as_str.index.map(str)
    for node in ts_tree.traverse_preorder():
        if is_leaf(node):
            lab = get_label(node)
            if lab is None:
                continue
            lab = str(lab)
            if lab in df_index_as_str.index:
                x = float(df_index_as_str.loc[lab, x_col])
                y = float(df_index_as_str.loc[lab, y_col])
                leaf_xy[node] = (x, y)

    if not leaf_xy:
        raise ValueError("No leaf labels from the tree matched merged_df.index. "
                         "Make sure the leaf names equal the DataFrame index.")

    # ---- 4) Compute descendant-leaf means for every node (postorder)
    # For each node, keep a running sum and count of descendant leaf coords
    sum_counts = {}  # node -> (sum_x, sum_y, count)
    for node in ts_tree.traverse_postorder():
        if is_leaf(node) and node in leaf_xy:
            x, y = leaf_xy[node]
            sum_counts[node] = [x, y, 1]
        else:
            sx, sy, cnt = 0.0, 0.0, 0
            for ch in get_children(node):
                if ch in sum_counts:
                    sx += sum_counts[ch][0]
                    sy += sum_counts[ch][1]
                    cnt += sum_counts[ch][2]
            # If a clade has no leaves in df, skip (won't plot that node)
            if cnt > 0:
                sum_counts[node] = [sx, sy, cnt]

    xy = {}  # node -> (x, y)
    for node, (sx, sy, cnt) in sum_counts.items():
        xy[node] = (sx / cnt, sy / cnt)

    # ---- 5) Compute cumulative branch lengths (z) from the root (preorder)
    z = {root: 0.0}
    for node in ts_tree.traverse_preorder():
        for ch in get_children(node):
            z[ch] = z.get(node, 0.0) + get_edge_length(ch)

    # ---- 6) Build edge segments (lines) and node markers
    # - one lines trace with NaN separators
    # - one markers trace for leaves (blue)
    # - one markers trace for internals (red)

    line_x, line_y, line_z = [], [], []
    leaf_x, leaf_y, leaf_z = [], [], []
    internal_x, internal_y, internal_z = [], [], []

    leaf_labels = []
    for node in ts_tree.traverse_preorder():
        # skip nodes that don't have xy (e.g., clades with no matching leaves)
        if node not in xy or node not in z:
            continue

        node_x, node_y, node_z = xy[node][0], xy[node][1], z[node]

        if is_leaf(node):
            leaf_x.append(node_x); leaf_y.append(node_y); leaf_z.append(node_z)
            leaf_labels.append(node.label)
        else:
            internal_x.append(node_x); internal_y.append(node_y); internal_z.append(node_z)

        par = get_parent(node)
        if par is not None and par in xy and par in z:
            px, py = xy[par]
            pz = z[par]
            # segment parent -> child
            line_x += [px, node_x, None]
            line_y += [py, node_y, None]
            line_z += [pz, node_z, None]

    # ---- 7) Create Plotly figure
    lines = go.Scatter3d(
        x=line_x, y=line_y, z=line_z,
        mode="lines",
        line=dict(width=3),
        name="Edges",
        hoverinfo="skip"
    )

    leaves = go.Scatter3d(
        x=leaf_x, y=leaf_y, z=leaf_z,
        text=leaf_labels,
        textposition="top center",
        textfont=dict(size=12),
        mode="markers+text",
        marker=dict(size=5, color="green"),
        name="Leaves",
        hoverinfo="skip"
    )

    internals = go.Scatter3d(
        x=internal_x, y=internal_y, z=internal_z,
        mode="markers",
        marker=dict(size=5, color="blue"),
        name="Internal nodes"
    )

    fig = go.Figure(data=[lines, internals, leaves])
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title=x_col,
            yaxis_title=y_col,
            zaxis_title="Cumulative branch length"
        ),
        legend=dict(itemsizing="constant")
    )

    # ---- 8) Open in browser
    if outfile:
        from pathlib import Path
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(outfile, include_plotlyjs=True, full_html=True)
        print(f"Saved interactive HTML to: {outfile}")
    try:
        from plotly.offline import plot as _offline_plot
        _offline_plot(fig, auto_open=True)
    except Exception:
        # Fallback: display in notebook environments (if applicable)
        fig.show()

    return fig

def edge_ratio_table(ts_tree, merged_df, x_col="x", y_col="y"):
    """
    Build a table with, for every node:
      - (x,y) position (leaf-based descendant mean for internals)
      - parent (x,y)
      - edge length (to parent)
      - Euclidean distance in (x,y) to parent
      - ratio = edge_length / xy_distance

    Parameters
    ----------
    ts_tree : treeswift.Tree
    merged_df : pd.DataFrame
        index = leaf names (string/int ok), columns include x_col, y_col
    """
    # ---- helpers for TreeSwift (robust to attribute name variants)
    def is_leaf(n):       return getattr(n, "is_leaf", lambda: len(get_children(n)) == 0)()
    def get_children(n):  return getattr(n, "children", getattr(n, "get_children", lambda: [])()) or []
    def get_parent(n):    return getattr(n, "parent", getattr(n, "get_parent", lambda: None)())
    def get_label(n):     return getattr(n, "label", getattr(n, "get_label", lambda: None)())
    def edge_len(n):
        el = getattr(n, "edge_length", getattr(n, "get_edge_length", lambda: None)())
        return 0.0 if el is None else float(el)
    root = getattr(ts_tree, "root", getattr(ts_tree, "get_root", lambda: None)())
    if root is None:
        raise ValueError("Cannot determine tree root.")

    # ---- gather leaf coordinates from merged_df (match on stringified labels)
    df_idx = merged_df.copy()
    df_idx.index = df_idx.index.map(str)
    leaf_xy = {}
    for n in ts_tree.traverse_preorder():
        if is_leaf(n):
            lab = get_label(n)
            if lab is None:
                continue
            lab = str(lab)
            if lab in df_idx.index:
                leaf_xy[n] = (float(df_idx.loc[lab, x_col]), float(df_idx.loc[lab, y_col]))

    if not leaf_xy:
        raise ValueError("No leaf labels matched merged_df.index.")

    # ---- descendant-mean (x,y) for every node (postorder)
    sum_counts = {}
    for n in ts_tree.traverse_postorder():
        if is_leaf(n) and n in leaf_xy:
            x, y = leaf_xy[n]
            sum_counts[n] = [x, y, 1]
        else:
            sx = sy = 0.0
            cnt = 0
            for ch in get_children(n):
                if ch in sum_counts:
                    sx += sum_counts[ch][0]
                    sy += sum_counts[ch][1]
                    cnt += sum_counts[ch][2]
            if cnt > 0:
                sum_counts[n] = [sx, sy, cnt]

    xy = {n: (sx / cnt, sy / cnt) for n, (sx, sy, cnt) in sum_counts.items()}

    # ---- cumulative z (not required for the ratio, but handy for debugging/exports)
    z = {root: 0.0}
    for n in ts_tree.traverse_preorder():
        for ch in get_children(n):
            z[ch] = z.get(n, 0.0) + edge_len(ch)

    # ---- assemble rows
    rows = []
    for n in ts_tree.traverse_preorder():
        lab = get_label(n)
        lab_str = str(lab) if lab is not None else ""
        par = get_parent(n)

        if n not in xy:
            # (subtree had no matched leaves) — skip row entirely
            continue

        x, y = xy[n]
        zz = z.get(n, np.nan)

        # parent values
        if par is not None and par in xy:
            px, py = xy[par]
            pz = z.get(par, np.nan)
            e = edge_len(n)
            d = math.hypot(x - px, y - py)
            ratio = (e / d) if d > 0 else (np.nan if e == 0 else np.inf)
            par_label = get_label(par)
            par_label_str = str(par_label) if par_label is not None else ""
        else:
            px = py = pz = np.nan
            e = d = ratio = np.nan
            par_label_str = ""

        rows.append({
            "node_label": lab_str,
            "is_leaf": is_leaf(n),
            "x": x, "y": y, "z": zz,
            "parent_label": par_label_str,
            "parent_x": px, "parent_y": py, "parent_z": pz,
            "edge_length": e,
            "xy_distance_to_parent": d,
            "bl_over_xy": ratio
        })

    df = pd.DataFrame(rows)
    # (Optional) bring leaves first, then internals
    df = df.sort_values(by=["is_leaf"], ascending=False).reset_index(drop=True)
    return df

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

def _palette_for(n, base=None):
    """
    Build a list of at least n distinct colors, starting with Okabe–Ito,
    then topping up from matplotlib's tab20 if needed.
    """
    base = list(base) if base is not None else _OKABE_ITO
    if len(base) >= n:
        return base[:n]
    # top up with tab20
    import matplotlib.cm as cm
    tab20 = [cm.get_cmap("tab20")(i) for i in range(20)]
    # convert RGBA to hex
    def rgba2hex(rgba):
        r,g,b,a = rgba
        return "#{:02X}{:02X}{:02X}".format(int(r*255), int(g*255), int(b*255))
    extras = [rgba2hex(c) for c in tab20]
    out = base + extras
    if len(out) < n:  # if still short, repeat extras
        reps = (n - len(out) + len(extras) - 1) // len(extras)
        out += (extras * reps)
    return out[:n]

def plot_genotype_panels(
    mats,
    labels=None,
    title=None,
    outfile=None,
    row_order=None,
    col_order=None,
    vmin=-1,
    vmax=8,
    missing_code=-1,
    okabe_colors=None,
    figsize_per_panel=(4.0, 4.5),
    dpi=150,
    show_ticks=False,
    rasterized=False,
):
    """
    Plot multiple genotype matrices as side-by-side heatmaps with a shared discrete colorbar.

    Parameters
    ----------
    matrices : list[pd.DataFrame]
        Each matrix is cells (rows) x target sites (cols). Values are integers in [vmin, vmax]
        with `missing_code` (default -1) indicating missing and will be masked.
        All matrices should be aligned to the same index/columns or be alignable via reindex.
    labels : list[str] | None
        Panel titles (one per matrix). If None, panels are named 'Matrix 1', 'Matrix 2', ...
    title : str | None
        Figure-level title prefix. If provided, each panel title is f"{title}: {labels[i]}".
    outfile : str | None
        If provided, the figure is saved to this path.
    row_order, col_order : array-like | None
        Optional explicit ordering of rows/cols to apply to all matrices.
        If provided, they are used to reindex (dropping absent items).
    vmin, vmax : int
        Minimum and maximum *coded* state values (inclusive). Typically -1..8.
    missing_code : int
        Value to treat as missing and mask out of the heatmap (default: -1).
    okabe_colors : list[str] | None
        Custom color list to use (falls back to Okabe–Ito + extras).
        Colors are mapped to integer categories from vmin..vmax (inclusive).
    figsize_per_panel : tuple(float, float)
        Size of each panel; overall width scales with number of matrices.
    dpi : int
        Figure DPI.
    show_ticks : bool
        If True, show axis tick labels; otherwise hide them (default False).
    rasterized : bool
        If True, rasterize the heatmap artists (useful for PDF size).

    Returns
    -------
    (fig, axes)
    """
    n_panels = len(mats)

    reordered_mats = [M.iloc[row_order, :].iloc[:, col_order] for M in mats]
    mats = reordered_mats

    # Build discrete categories & colormap
    categories = np.arange(vmin, vmax + 1, dtype=int)  # includes -1 if vmin=-1
    colors = _palette_for(len(categories), base=OKABE)
    cmap = ListedColormap(colors)
    # Make missing cells transparent via mask; but keep a "bad" color for safety
    cmap.set_bad(colors[0])  # first color (often used for -1) in case any NaNs slip through

    bounds = np.arange(vmin - 0.5, vmax + 1.5, 1.0)  # bin edges centered on ints
    norm = BoundaryNorm(bounds, cmap.N)

    # Figure & axes
    fig_w = figsize_per_panel[0] * n_panels
    fig_h = figsize_per_panel[1]
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, fig_h), dpi=dpi, constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    # Plot each panel
    for i, (ax, M) in enumerate(zip(axes, mats)):
        # mask = (M.values == missing_code) | np.isnan(M.values)
        sns.heatmap(
            M,  ax=ax, # mask=mask,
            cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            xticklabels=show_ticks, yticklabels=show_ticks, cbar=False,
            rasterized=rasterized
        )
        # Titles/labels
        panel_title = labels[i]
        if title:
            ax.set_title(f"{title}: {panel_title}", fontsize=16)
        else:
            ax.set_title(panel_title, fontsize=16)

        ax.set_xlabel("Target sites", fontsize=13)
        ax.set_ylabel("Cell names" if i == 0 else "", fontsize=13)
        ax.tick_params(axis="both", labelsize=11)

        # Draw a clean box around each panel
        for side in ("left", "right", "top", "bottom"):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_linewidth(1.5)
            ax.spines[side].set_edgecolor("black")

    # Shared colorbar
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes, orientation="horizontal", fraction=0.08, pad=0.05
    )
    cbar.set_ticks(categories)
    ticklabels = [f"?/{missing_code}" if x == missing_code else str(x) for x in categories]
    cbar.set_ticklabels(ticklabels)
    cbar.ax.tick_params(labelsize=11)

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", dpi=dpi)
    return fig, axes

