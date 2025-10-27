#!/usr/bin/env python3
# newick_to_digraph.py
from ete3 import Tree
import networkx as nx
from pathlib import Path
import argparse
from typing import Tuple, Optional

def newick_to_digraph(
    newick_string: str,
    *,
    length_attr: str = "length",
    time_attr: str = "time",
    preserve_internal_labels: bool = True,
    internal_prefix: str = "I",
    ete_format: int = 1
) -> Tuple[nx.DiGraph, str]:
    """
    Convert a Newick string into a rooted NetworkX DiGraph (arborescence).

    - Edges are directed parent->child.
    - Edge attribute `{length_attr}` stores branch length (float; 0.0 if missing).
    - Node attribute `{time_attr}` stores cumulative distance from root (depth).
    - Node identifiers:
        * Leaves: their Newick labels (must be unique in input).
        * Internal nodes:
            - If the Newick has internal labels and preserve_internal_labels=True, use them.
            - Otherwise, assign stable names: I1, I2, ...

    Returns
    -------
    (G, root_name)
        G: nx.DiGraph arborescence
        root_name: identifier of the root node in G
    """
    # Parse with ETE. format=1 supports internal node names if present.
    ete_tree = Tree(newick_string.strip(), format=ete_format)

    G = nx.DiGraph()
    internal_counter = 0
    name_cache = {}

    def get_node_name(node) -> str:
        nonlocal internal_counter
        if node.is_leaf():
            if not node.name:
                # Leaves really should be named in most phylo data;
                # fall back to a stable generated name if missing.
                internal_counter += 1
                return f"{internal_prefix}Leaf{internal_counter}"
            return node.name
        # Internal node
        if preserve_internal_labels and node.name:
            return node.name
        if node in name_cache:
            return name_cache[node]
        internal_counter += 1
        gen = f"{internal_prefix}{internal_counter}"
        name_cache[node] = gen
        return gen

    # First pass: create nodes and compute times via preorder traversal
    root_name = get_node_name(ete_tree)
    G.add_node(root_name, **{time_attr: 0.0})

    # Stack holds (ete_node, parent_name)
    stack = [(ete_tree, None)]
    while stack:
        cur, parent_name = stack.pop()
        cur_name = get_node_name(cur)

        # Ensure node exists (might have been added for root)
        if cur_name not in G:
            G.add_node(cur_name)

        # Compute children times and add edges
        for child in cur.children:
            child_name = get_node_name(child)
            # ETE uses node.dist as the length of the edge from node.parent to node
            # For a child, the edge length is child.dist
            bl = float(child.dist) if child.dist is not None else 0.0

            # Parent time
            parent_time = G.nodes[cur_name].get(time_attr, 0.0)
            child_time = parent_time + bl

            G.add_node(child_name, **{time_attr: child_time})
            G.add_edge(cur_name, child_name, **{length_attr: bl})

            stack.append((child, cur_name))

    # Sanity check: arborescence?
    # (A well-formed Newick should yield an arborescence.)
    if not nx.is_arborescence(G):
        # This typically means duplicated names or multifurcations with name collisions.
        # NetworkX arborescence allows multifurcations; name collisions are the usual culprit.
        raise ValueError("Resulting graph is not an arborescence; check for duplicate node labels.")

    return G, root_name


# --- CLI wrapper (optional) ---
def main():
    ap = argparse.ArgumentParser(description="Convert Newick to NetworkX DiGraph with time (depth).")
    ap.add_argument("newick_file", help="Path to .nwk/.newick file")
    ap.add_argument("-o", "--output", help="Output GraphML/GEXF/GPickle/JSON path (by extension). "
                                           "Default: <input>.graphml")
    ap.add_argument("--length-attr", default="length", help="Edge attribute name for branch length (default: length)")
    ap.add_argument("--time-attr", default="time", help="Node attribute name for cumulative depth (default: time)")
    ap.add_argument("--no-preserve-internal-labels", action="store_true",
                    help="Do not use internal labels from Newick even if present")
    ap.add_argument("--internal-prefix", default="I", help="Prefix for generated internal node names (default: I)")
    ap.add_argument("--ete-format", type=int, default=1, help="ETE Newick parse format (default: 1)")
    args = ap.parse_args()

    in_path = Path(args.newick_file)
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")

    newick_text = in_path.read_text()
    G, root = newick_to_digraph(
        newick_text,
        length_attr=args.length_attr,
        time_attr=args.time_attr,
        preserve_internal_labels=not args.no_preserve_internal_labels,
        internal_prefix=args.internal_prefix,
        ete_format=args.ete_format
    )

    # Basic stats
    leaves = [n for n in G.nodes if G.out_degree(n) == 0]
    internals = [n for n in G.nodes if G.out_degree(n) > 0]
    print(f"Root: {root}")
    print(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")
    print(f"Leaves: {len(leaves)} | Internal nodes: {len(internals)}")

    # Decide output path/format by extension
    out_path = Path(args.output) if args.output else in_path.with_suffix(".graphml")
    ext = out_path.suffix.lower()
    if ext == ".graphml":
        nx.write_graphml(G, out_path)
    elif ext == ".gexf":
        nx.write_gexf(G, out_path)
    elif ext == ".gpickle":
        nx.write_gpickle(G, out_path)
    elif ext == ".json":
        from networkx.readwrite import json_graph
        out_path.write_text(
            __import__("json").dumps(json_graph.node_link_data(G))
        )
    else:
        raise SystemExit(f"Unsupported output format by extension: {ext}\n"
                         f"Use one of: .graphml, .gexf, .gpickle, .json")

    print(f"Saved: {out_path}")

if __name__ == "__main__":
    main()

