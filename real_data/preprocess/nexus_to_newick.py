#!/usr/bin/env python3
import re, sys, pathlib

# Usage: python nexus_to_newick.py baseMemoir.pos5.nexus > out.newick

p = pathlib.Path(sys.argv[1]).read_text()

# 1) TRANSLATE map (if present)
tmatch = re.search(r'(?is)translate\s+(.*?);', p)
trans = {}
if tmatch:
    for part in re.split(r',\s*', tmatch.group(1).strip()):
        if not part:
            continue
        m = re.match(r'^\s*([^,\s=]+)\s*=?\s*(.+?)\s*$', part)
        if not m:
            continue
        k, v = m.group(1), m.group(2).rstrip(',')
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
            v = v[1:-1]
        trans[k] = v

# 2) First TREE’s Newick
nmatch = re.search(r'(?im)^\s*tree\s+[^=]*=\s*(.*?);', p)
if not nmatch:
    sys.exit("No TREE found.")
newick = nmatch.group(1).strip()

# 3) Retranslate labels (avoid touching branch lengths after ':')
if trans:
    keys = sorted(trans.keys(), key=len, reverse=True)
    pat = re.compile(r'(?<!:)\b(?:' + '|'.join(map(re.escape, keys)) + r')\b')
    DELIMS = set("(),:;[]")
    def needs_quotes(s): return any(c.isspace() or c in DELIMS for c in s)
    def repl(m):
        lab = trans[m.group(0)]
        return "'" + lab.replace("'", "''") + "'" if needs_quotes(lab) else lab
    full_tree = pat.sub(repl, newick)
else:
    full_tree = newick

# helper: strip square-bracket comments while respecting quotes
def strip_bracket_comments(s: str) -> str:
    out = []
    i, n = 0, len(s)
    in_sq = False  # inside single-quoted label
    while i < n:
        c = s[i]
        if in_sq:
            out.append(c)
            # handle escaped '' inside quotes
            if c == "'":
                if i + 1 < n and s[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                in_sq = False
            i += 1
            continue
        # not in quotes: start of comment?
        if c == '[':
            # skip until the next ']' (non-nesting per Newick/Nexus comments)
            j = s.find(']', i + 1)
            if j == -1:
                break  # malformed; drop remainder
            i = j + 1
            continue
        if c == "'":
            in_sq = True
            out.append(c)
            i += 1
            continue
        out.append(c)
        i += 1
    return ''.join(out)

stripped_tree = strip_bracket_comments(full_tree).strip()

# 4) Print both, each ending with ';'
def ensure_semicolon(x: str) -> str:
    x = x.strip()
    return x if x.endswith(';') else x + ';'

print(ensure_semicolon(full_tree))
print(ensure_semicolon(stripped_tree))
