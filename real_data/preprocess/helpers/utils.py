from scipy import stats
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import pandas as pd


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

def scree(X, out=None, target=0.90):
    # PCA
    vr = PCA().fit(X).explained_variance_ratio_
    cv = np.cumsum(vr); k = np.arange(1, cv.size + 1)

    # Plot
    plt.figure(figsize=(8, 4.5), dpi=150)
    plt.plot(k, cv, marker='o', linewidth=2.5)
    plt.axhline(target, ls='--', color='red', linewidth=1.5, label=f'{int(target*100)}% variance')
    plt.title('baseMemoir over all states and colonies', fontsize=16)
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
