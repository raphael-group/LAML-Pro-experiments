import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import json
import os, re

palette = list(plt.colormaps['tab20'].colors)

def _common_suffix(strings):
    rev = [s[::-1] for s in strings]
    suf_rev = os.path.commonprefix(rev)
    return suf_rev[::-1]

def _next_boundary(s, L):
    """
    Return the smallest index >= L that lands on a token boundary.
    Tokens are: [A-Za-z0-9_]+ (keep words and numbers together) or a single non-word char.
    """
    if L >= len(s):
        return len(s)
    # Precompute ends of tokens once per string
    ends = getattr(_next_boundary, "_cache", None)
    if ends is None or ends.get("key") != s:
        tokens = list(re.finditer(r'[A-Za-z0-9_]+|[^A-Za-z0-9_]', s))
        ends_set = {m.end() for m in tokens}
        _next_boundary._cache = {"key": s, "ends": sorted(ends_set)}
    for e in _next_boundary._cache["ends"]:
        if e >= L:
            return e
    return len(s)

def minimal_unique_labels(labels):
    if not labels:
        return labels
    pref = os.path.commonprefix(labels)
    suf  = _common_suffix(labels)
    # Trim shared prefix/suffix
    trimmed = [s[len(pref): len(s)-len(suf) if suf else None] or s for s in labels]

    # Shortest unique prefix per trimmed label, cut only at token boundaries
    uniq = []
    for i, s in enumerate(trimmed):
        got = None
        for L in range(1, len(s)+1):
            Lb = _next_boundary(s, L)
            cand = s[:Lb]
            if sum(t.startswith(cand) for t in trimmed) == 1:
                got = cand
                break
        uniq.append(got or f"{s}_{i}")  # fallback if needed

    # Final de-dupe just in case
    for i, s in enumerate(uniq):
        if uniq.count(s) > 1:
            uniq[i] = f"{s}_{i}"
    return uniq

# Use custom style
plt.style.use(os.path.expanduser("~/plotting/paper.mplstyle"))

# Parse arguments
parser = argparse.ArgumentParser(description='Plot log likelihoods from multiple results.json files')
parser.add_argument('results', nargs='+', type=str, help='Paths to JSON files with log_likelihoods list')
parser.add_argument('--labels', nargs='*', help='Optional labels for each curve')
args = parser.parse_args()

# Set up plot
plt.figure(figsize=(6, 4))

filepaths = args.results
raw_labels = [
    (args.labels[i] if args.labels and i < len(args.labels)
     else os.path.basename(filepath).split('.')[1])
    for i, filepath in enumerate(filepaths)  # filepaths = list of your input files
]

run_labels = minimal_unique_labels(raw_labels)


# Iterate through each results file
for i, filepath in enumerate(args.results):
    color = palette[i % len(palette)]
    with open(filepath, 'r') as f:
        data = json.load(f)

    log_likelihoods = data['log_likelihoods']
    iterations = list(range(len(log_likelihoods)))
    print(f"{filepath}: {len(log_likelihoods)} iterations")

    label = run_labels[i]
    df = pd.DataFrame({
        'Iteration': iterations,
        'Log Likelihood': log_likelihoods,
        'Run': label,
    })

    plt.plot(df['Iteration'], df['Log Likelihood'], label=df['Run'].iloc[0], linewidth=1.0, color=color)

    plt.scatter(
        df['Iteration'].iloc[-1],
        df['Log Likelihood'].iloc[-1],
        marker='x',
        color='black',
        zorder=5
    )



# Finalize plot
plt.title('Log-density vs Iteration')
plt.xlabel('Simulated Annealing Iterations')
plt.ylabel('Log-density')
plt.legend(frameon=True, facecolor='white', edgecolor='black', loc='lower right')

plt.grid(True)
plt.tight_layout()
plt.show()

