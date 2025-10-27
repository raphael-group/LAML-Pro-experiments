import numpy as np
from ete3 import Tree
import argparse
import pandas as pd
import math
class Square_Parsimony():
    def __init__(self, tree, leaf_states, use_branch_length=True, epsilon=1e-6):
        """
        Parameters:
          tree (str or ete3.Tree): If a string, assumed to be Newick format.
          leaf_states (dict): Dictionary mapping leaf node names to state vectors 
                              (e.g. np.array of shape (dim,)).
          use_branch_length (bool, optional): Whether to use branch lengths. If False, every branch is 1.
          epsilon (float, optional): Small value to avoid division by zero.
        Raises:
          ValueError: If the set of leaf names in the tree does not match the keys in leaf_states.
        """
        if isinstance(tree, str):
            self.tree = Tree(tree, format=1)
        else:
            self.tree = tree
            
        # Ensure the root has a name.
        if self.tree.get_tree_root().name == "":
            self.tree.get_tree_root().name = "root"
            
        self.use_branch_length = use_branch_length
        # Extract the leaf names from the tree.
        leaf_names = [node.name for node in self.tree.get_leaves()]
        if set(leaf_names) != set(leaf_states.keys()):
            print("tree leaf names", set(leaf_names))
            print("centroids leaf names", set(leaf_states.keys()))
            raise ValueError("leaf_states keys do not match tree leaf names")
            
        self.leaf_states = leaf_states
        # Infer the dimension (assume all state vectors have the same length).
        self.dim = len(next(iter(leaf_states.values())))
        print("Dimension:", self.dim)
        self.epsilon = epsilon
        if not self.use_branch_length:
            for node in self.tree.traverse():
                node.dist = 1
            self.epsilon = 0  # if branch lengths are not used, set epsilon to 0

        # To store per-dimension minimum squared sums
        self.min_square_sum = np.zeros(self.dim)
        self.min_square_sum_total = 0

        # Set up additional attributes for dynamic programming:
        self.set_parameters()
        self.compute_min_square_sum()
        self.infer_ancestral()
        # for node in self.tree.traverse():
        #     print(node.name, node.p_array, node.q_array, node.inferred_state)
        # Finally, update the "state" attribute for internal nodes from the computed inferred state.
        for node in self.tree.traverse():
            if not node.is_leaf():
                # Set the final inferred vector (make a copy to avoid accidental aliasing)
                node.state = node.inferred_state.copy()
            else:
                # For leaves, use the provided state.
                node.state = self.leaf_states[node.name]

    def set_parameters(self):
        """
        For each node in the tree, add three features:
          - p_array: a (dim x 3) NumPy array holding the coefficients for the postorder pass.
          - q_array: a (dim x 3) NumPy array for the preorder pass.
          - inferred_state: a length-d vector (NumPy array) for the inferred state in each dimension.
        """
        for node in self.tree.traverse():
            node.add_feature("p_array", np.zeros((self.dim, 3)))
            node.add_feature("q_array", np.zeros((self.dim, 3)))
            node.add_feature("inferred_state", np.zeros(self.dim))

    def compute_min_square_sum(self):
        """
        Postorder pass: For each dimension, process leaves then internal nodes in postorder.
        For leaves, initialize the quadratic cost for the observed state. For internal nodes,
        sum over the contributions of each child (using a recurrence that, for internal children,  
        uses the formula:
        
           new_term = [ p1/(vp*p1 + 1), p2/(vp*p1 + 1), - p2**2/(4*(p1 + 1/vp)) + p3 ]
        
        where vp = child.dist + epsilon.
        At the root, compute the inferred state as the minimizer of the quadratic, and record the  
        minimum squared sum.
        """
        # Process nodes in postorder (children before parent)
        nodes = list(self.tree.traverse(strategy="postorder"))
        var_penalty = 0

        for i in range(self.dim):
            # Initialize leaf nodes.
            num_edges = 0
            for node in self.tree.get_leaves():
                vp = node.dist + self.epsilon
                p_val = self.leaf_states[node.name][i]
                node.inferred_state[i] = p_val
                # Set quadratic coefficients: [a, b, c]
                node.p_array[i] = np.array([1/vp, -2*p_val/vp, (p_val**2)/vp])
            # Process each internal node.
            for node in nodes:
                num_edges += 1
                vp = node.dist + self.epsilon
                var_penalty += np.log(1/math.sqrt(2 * math.pi * vp))
                if node.is_leaf():
                    continue
                else:
                    for child in node.children:
                        vp = child.dist + self.epsilon
                        if child.is_leaf():
                            node.p_array[i] += child.p_array[i]
                        else:
                            p1, p2, p3 = child.p_array[i]
                            # Compute the updated coefficients according to the recurrence.
                            new_a = p1 / (vp * p1 + 1)
                            new_b = p2 / (vp * p1 + 1)
                            new_c = - (p2 ** 2) / (4 * (p1 + 1/vp)) + p3
                            node.p_array[i] += np.array([new_a, new_b, new_c])
                    if node.is_root():
                        node.inferred_state[i] = - node.p_array[i][1] / (2 * node.p_array[i][0])
                        # Initialize q_array at the root as a copy of p_array.
                        node.q_array[i] = node.p_array[i].copy()
                        min_sq = (node.p_array[i][0] * node.inferred_state[i]**2 +
                                  node.p_array[i][1] * node.inferred_state[i] +
                                  node.p_array[i][2])
                        self.min_square_sum[i] = min_sq
                        self.min_square_sum_total += min_sq
        print("Total min square sum:", self.min_square_sum_total)
        print("Total maximum likelihood under Brownian motion:", var_penalty - 0.5 * self.min_square_sum_total)
        print("Avg. min square sum (over branches):", self.min_square_sum_total /num_edges)
        print("Avg. maximum likelihood under Brownian motion (over branches):", (var_penalty - 0.5 * self.min_square_sum_total)/num_edges)
        

    def infer_ancestral(self):
        """
        Preorder pass: For each dimension, update internal nodes (except root and leaves) by  
        aggregating information from the parent and children. The update step uses the parent's  
        q_array and adds terms from the node itself and its children. Then the node's inferred state  
        is updated as the minimizer of the updated quadratic.
        """
        nodes = list(self.tree.traverse(strategy="preorder"))
        for i in range(self.dim):
            for node in nodes:
                if node.is_leaf() or node.is_root():
                    continue
                else:
                    q1,q2,q3 = node.up.q_array[i]
                    vp = node.dist + self.epsilon
                    node.q_array[i] += np.array([1/vp, -2 * node.up.inferred_state[i] / vp, 0])
                    for child in node.children:
                        vp_child = child.dist + self.epsilon
                        if child.is_leaf():
                            node.q_array[i] += child.p_array[i]
                        else:
                            c_p1, c_p2, _ = child.p_array[i]
                            node.q_array[i] += np.array([c_p1 / (vp_child * c_p1 + 1),
                                                         c_p2 / (vp_child * c_p1 + 1),
                                                         0])
                    # Update the inferred state in dimension i.
                    node.inferred_state[i] = - node.q_array[i][1] / (2 * node.q_array[i][0])

    def mean_squared_error(self, ground_truth):
        """
        ground_truth: dict[node_name]: ground_truth_data
        Calculate the mean squared error between the predicted states and the actual 
        internal node states.
        """
        mse = 0.0
        count = 0
        for node in self.tree.traverse():
            if node.is_leaf():
                continue
            if node.name in ground_truth:
                mse += np.sum((node.state - ground_truth[node.name]) ** 2)
                count += 1
        return mse / count if count > 0 else float('inf')

    def get_inferred_states(self, output_file=None):
        """
        Return a list of internal node states as a list of lists: [node.name, state_vector...].
        Optionally, if output_file is provided, write a CSV file with the inferred states.
        """
        inferred_states = []
        for node in self.tree.traverse():
            if not node.is_leaf():
                inferred_states.append([node.name] + list(node.state))
        node_names = [state[0] for state in inferred_states]
        state_values = [state[1:] for state in inferred_states]
        df = pd.DataFrame(state_values, index=node_names, columns=[f"Dim_{i}" for i in range(self.dim)])
        if output_file is not None:
            df.to_csv(output_file, index=True)
        return df

# ------------------------------
# Example of using the class:
# ------------------------------

def simple_test_euclidean():
    newick_str = "((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);"
    tree = Tree(newick_str, format=1)

    # Create a dummy leaf_states dictionary.
    # Let’s say our states are 2-dimensional numeric vectors.
    leaf_states = {
        "A": np.array([5.0, 2.0]),
        "B": np.array([7.0, 3.0]),
        "C": np.array([3.0, 4.0]),
        "D": np.array([9.0, 1.0])
    }

    # Instantiate the Square_Parsimony solver.
    # (use_branch_length=False forces branch lengths to 1)
    sp = Square_Parsimony(tree, leaf_states, use_branch_length=False, epsilon=1e-6)

    # After running the inference, the tree nodes will have their "state" attribute assigned.
    # Print out the inferred state at each internal node.
    for node in sp.tree.traverse():
        if not node.is_leaf():
            print("Node", node.name, "inferred state:", node.state)

    # Optionally, get a list of inferred states.
    inferred = sp.get_inferred_states()
    print("Inferred internal node states:")
    for entry in inferred:
        print(entry)

def test_simple_1d():
    # Path to the Newick tree file
    test_tree_path = "/n/fs/ragr-data/users/viola/mouse_dev/scripts/c_elegans/test_continous_tree.nwk"
    
    # Load the tree from the Newick file
    tree = Tree(test_tree_path, format=1)
    
    # Create a DataFrame with leaf states
    test_df = pd.DataFrame([[0.0], [0.0], [0.0], [0.98], [1.2], [1.3], [0.0]], 
                            index=["a", "b", "c", "d", "e", "f", "g"])
    
    # Convert the DataFrame to a dictionary for leaf_states
    leaf_states = {index: row.values for index, row in test_df.iterrows()}
    
    # Instantiate the Square_Parsimony solver
    model = Square_Parsimony(tree, leaf_states, use_branch_length=False, epsilon=1e-6)
    
    # Print inferred states for internal nodes
    for node in model.tree.traverse():
        if not node.is_leaf():
            print("Node", node.name, "inferred state:", node.state)
    
    # Optionally, get a list of inferred states  
    inferred = model.get_inferred_states()
    print("Inferred internal node states:")
    for entry in inferred:
        print(entry)


def test_simple_4_leaves_1D():
    newick = "((A:1,B:1)N1:1,(C:1,D:1)N2:1)Root;"
    tree = Tree(newick, format=1)
    Omega = np.array([0.001 * i for i in range(1000)])
    leaf_states = {"A": [0.9], "B": [0.99], "C": [0.2], "D": [0.01]}
    # for leaf, state in leaf_states.items():
    #     Omega = np.vstack([Omega, state])
    model = Square_Parsimony(tree, leaf_states, use_branch_length=False, epsilon=1e-6)
    # for node in model.tree.traverse():
    #     if not node.is_leaf():
    #         print("Node", node.name, "inferred state:", node.state)

    # Optionally, get a list of inferred states.
    inferred = model.get_inferred_states()
    print("Inferred internal node states:")
    for entry in inferred:
        print(entry)

def main():
    parser = argparse.ArgumentParser(description="Infer internal node states using Square Parsimony.")
    parser.add_argument("--tree", "-t",type=str, required=True, help="Path to Newick tree file.")
    parser.add_argument("--leaf_df", "-l",type=str, required=True, help="CSV file with leaf states (index: leaf names, columns: vector dimensions).")
    parser.add_argument("--output", "-o", type=str, required=True, help="Output CSV file for inferred internal node states.")
    parser.add_argument("--use_branch_length", action="store_true", help="Use branch lengths from the tree.")
    parser.add_argument("--epsilon", type=float, default=1e-6, help="Small value to avoid division by zero.")

    args = parser.parse_args()

    # Load tree
    tree = Tree(args.tree, format=1)

    # Load leaf states DataFrame
    df = pd.read_csv(args.leaf_df, index_col=0, header=None, sep="\t")
    leaf_states = {str(idx): row.values for idx, row in df.iterrows()}

    # Run inference
    model = Square_Parsimony(tree, leaf_states, use_branch_length=args.use_branch_length, epsilon=args.epsilon)

    # Output inferred internal node states
    inferred = model.get_inferred_states(output_file=args.output)
    print(f"Wrote inferred internal node states to {args.output}")

if __name__ == "__main__":
    main()
