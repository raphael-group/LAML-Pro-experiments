from scipy import stats
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import pandas as pd
import re

plt.rcParams.update({"font.size": 9})

def _qq_stats_text(x):
    x = x[np.isfinite(x)]
    if x.size == 0: 
        return "n=0"
    # Shapiro reliable up to ~5000; subsample if larger
    xs = x if x.size <= 5000 else np.random.default_rng(0).choice(x, 5000, replace=False)
    W, p_sw = stats.shapiro(xs)
    D, p_ks = stats.kstest(x, 'norm')              # vs N(0,1) since already z-scored
    A = stats.anderson(x, dist='norm')
    # 5% critical value + reject flag
    idx5 = np.where(A.significance_level == 5)[0]
    crit5 = float(A.critical_values[idx5[0]]) if idx5.size else np.nan
    #rej5  = A.statistic > crit5 if np.isfinite(crit5) else False
    return f"SW W={W:.3f}, p={p_sw:.1e}\nKS D={D:.3f}, p={p_ks:.1e}\nAD A²={A.statistic:.3f}, 5%>{crit5:.3f}" #[{ 'rej' if rej5 else 'ok' }]"

def scree(X, out=None, target=0.90, title="basememoir over all states and colonies"):
    # PCA
    vr = PCA().fit(X).explained_variance_ratio_
    cv = np.cumsum(vr); k = np.arange(1, cv.size + 1)

    # Plot
    plt.figure(figsize=(8, 4.5), dpi=150)
    plt.plot(k, cv, marker='o', linewidth=2.5)
    plt.axhline(target, ls='--', color='red', linewidth=1.5, label=f'{int(target*100)}% variance')
    plt.title(title, fontsize=16)
    plt.xlabel('Number of Principal Components', fontsize=14)
    plt.ylabel('Cumulative Explained Variance', fontsize=14)
    plt.xticks(k, fontsize=12); plt.yticks(fontsize=12)
    plt.ylim(0, 1.05); plt.xlim(0.55, k[-1]+0.55)
    #plt.grid(True, alpha=0.35)
    plt.legend(loc='lower right', fontsize=12) #, frameon=True)
    plt.tight_layout()
    if out: plt.savefig(out, bbox_inches='tight')
    plt.show()
    return vr, cv

def learn_kde(data, var_type=None):
    X = np.asarray(data, float)
    mask = np.isfinite(X).all(axis=1)
    X = X[mask]
    if var_type is None:
        var_type = 'c' * X.shape[1]
    return KDEMultivariate(data=X, var_type=var_type, bw='normal_reference')

def score_with_kdes(X_pca, kde_by_class, eps=1e-300):
    X_pca = np.asarray(X_pca, float)
    classes = np.array(list(kde_by_class.keys()))
    logpdf = np.column_stack([
        np.log(np.maximum(kde_by_class[c].pdf(X_pca), eps))
        for c in classes
    ])
    y_pred = classes[np.argmax(logpdf, axis=1)]

    argmax_idx = logpdf.argmax(axis=1)
    preds = pd.DataFrame({
        "pred_label": classes[argmax_idx],
        "argmax_idx": argmax_idx,
        "argmax_logpdf": logpdf[np.arange(len(argmax_idx)), argmax_idx],
    })
    scores = pd.DataFrame(logpdf, columns=[f'state{x}_prob' for x in classes])

    return preds, scores 
    

def reshape_to_cmat(tbl):
    # tbl should have columns: target_site, cell_name, pred_label
    state_wide = (
        tbl.pivot_table(
            index='target_site',
            columns='cell_name',
            values='pred_label',
            aggfunc='first'   
        )
        .sort_index()                
        .sort_index(axis=1)          
        .fillna(-1)
        .astype(int)
    )
    return state_wide

def _site_slice(dfw, site, id_cols, codebook_df):
    """Extract one site's row-block with mapped features & labels using a codebook."""
    out = dfw[id_cols].copy()
    out["target_site"] = site

    # Predictions / labels tied to this site (if present)
    if site in dfw.columns:
        out["pet_state"] = dfw[site]
    if f"{site}_actual" in dfw.columns:
        out["seq_state"] = dfw[f"{site}_actual"]
    if f"{site}_brightest" in dfw.columns:
        out["brightest_state"] = dfw[f"{site}_brightest"]
    if f"{site}_prob" in dfw.columns:
        out["pet_prob"] = dfw[f"{site}_prob"]

    # --- Use the codebook to determine the 9 r-bit columns for this site ---
    cb_site = codebook_df.query("site == @site").copy()

    # Put edits first, then the site's 'unedited' bit last (if present)
    cb_site["is_unedited"] = (cb_site["edit"].astype(str).str.lower() == "unedited")
    cb_site = cb_site.sort_values(["is_unedited", "bit"], ascending=[True, True])

    cols_for_site = cb_site["bit"].tolist()

    # Expect exactly 9 bits per site (8 edits + 1 unedited)
    missing = [c for c in cols_for_site if c not in dfw.columns]
    if missing:
        raise ValueError(
            f"Missing r-bit columns in dfw for site={site}: {missing}. "
            f"Expected from codebook: {cols_for_site}"
        )
    if len(cols_for_site) != 9:
        raise ValueError(
            f"Expected 9 bits for site={site} from codebook, got {len(cols_for_site)}: {cols_for_site}"
        )

    # --- Map to feature_0..feature_8 (numeric) ---
    block_vals = dfw[cols_for_site].apply(pd.to_numeric, errors="coerce").to_numpy()

    for j in range(9):
        out[f"feature_{j}"] = block_vals[:, j]

    # (Optional but handy) carry the bit names and edit labels for reference
    # These are constant per site → repeat so they travel with the rows.
    #for j, (bit, edit) in enumerate(zip(cb_site["bit"].tolist(), cb_site["edit"].tolist())):
    #    out[f"feature_{j}_bit"] = bit
    #    out[f"feature_{j}_edit"] = edit

    return out
