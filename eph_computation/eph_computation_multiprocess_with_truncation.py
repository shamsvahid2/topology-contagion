import numpy as np
import pandas as pd
import scipy.stats
import gudhi as gd
import networkx as nx
import random
from tqdm import tqdm
import os
from concurrent.futures import ProcessPoolExecutor  # ADDED (multiprocessing): for parallel execution

# --- Persistent Homology Computation ---
def compute_persistent_homology(graph, filtration):
    """
    Compute the persistent homology of a graph with a given filtration.

    :param graph: A networkx graph.
    :param filtration: A list of node filtration values, corresponding to each node in the graph.
    :return: A list of persistent homology generators of dimension one.
    """
    # Initialize the simplex tree for Gudhi
    st = gd.SimplexTree()

    # Insert nodes with filtration values
    for i, node in enumerate(graph.nodes()):
        st.insert([node], filtration=filtration[i])

    # Insert edges with filtration equal to the max of filtration of its vertices
    for edge in graph.edges():
        st.insert(list(edge), filtration=max(
            filtration[list(graph.nodes()).index(edge[0])],
            filtration[list(graph.nodes()).index(edge[1])]
        ))

    # Expand to clique complex
    st.extend_filtration()

    # Remove null cycles
    persistence_blocks = st.extended_persistence(min_persistence=1e-5)
    tmp = []
    tmp.extend(persistence_blocks[0])
    tmp.extend(persistence_blocks[1])
    tmp.extend(persistence_blocks[2])
    tmp.extend(persistence_blocks[3])
    persistence = tmp

    # Ensure birth <= death
    tmp = []
    for gen in persistence:
        tmp.append((gen[0], (min(gen[1][0], gen[1][1]), max(gen[1][0], gen[1][1]))))

    persistence = tmp

    # Filter for dimension one homology
    dim1_res = []
    for gen in persistence:
        if gen[0] == 1:
            if gen[1][0] != gen[1][1]:
                dim1_res.append(gen)
    return dim1_res

# --- Worker for one row (used in multiprocessing) ---
def _compute_eph_corr_row(args):  # ADDED (multiprocessing)
    """
    Worker for a single row: compute EPH and corr for one simulation.
    This is used in multiprocessing; all arguments must be picklable.

    Parameters
    ----------
    args : tuple
        (order_prl_row, graph, node_list, deg_list)

    Returns
    -------
    eph_value : float
    corr_value : float
    """
    order_prl_row, graph, node_list, deg_list = args

    # Indices of nodes that actually got infected (not NaN in order_prl_row)
    not_nan_ind = np.where(~np.isnan(order_prl_row))[0]
    sgraph_nodes = node_list[not_nan_ind]
    sgraph = graph.subgraph(sgraph_nodes)
    filt = order_prl_row[not_nan_ind]
    sgraph = nx.Graph(sgraph)

    persistent_homology = compute_persistent_homology(
        graph=sgraph,
        filtration=-filt  # keep sign convention
    )

    # Lifetime of generators
    life_set = []
    for gen in persistent_homology:
        birth, death = gen[1]
        life_set.append(death - birth)
    eph_value = np.nanmean(life_set)

    # Compute PRL (Spearman correlation between degree and infection order)
    order_curr = order_prl_row  # includes NaNs for non-infected
    try:
        rho, _ = scipy.stats.spearmanr(deg_list, order_curr, nan_policy='omit')
    except Exception:
        # Preserve original behavior: on failure, set rho = 0
        # (printing is omitted in worker to avoid cluttering output)
        rho = 0.0

    return eph_value, rho

# --- Add EPH and Correlation ---
def add_eph_and_corr_to_csv(csv_path, graph, up_to_step_infection=None):  # ADDED: up_to_step_infection
    """
    Add EPH and correlation (corr) columns to the given CSV file based on the graph.

    :param csv_path: Path to the CSV file.
    :param graph: The corresponding graph (networkx.Graph).
    :param up_to_step_infection: int or None.
        If not None, nodes whose infection step is strictly greater than this
        value are treated as never infected (set to NaN) before computing
        EPH and corr.
    """
    # Choose output name; encode the cutoff (if any) in the filename  # ADDED
    if up_to_step_infection is None:  # ADDED
        output_csv = csv_path.replace('.csv', '_EPH-v3.csv')
    else:  # ADDED
        output_csv = csv_path.replace('.csv', f'_EPH-v3_step-{up_to_step_infection}.csv')

    if os.path.exists(output_csv):
        print(output_csv, "already exists, continue")
        return -1

    sims = pd.read_csv(csv_path)

    # Identify columns that correspond to nodes (names are integers, possibly negative)
    order_columns = [col for col in sims.columns
                     if str(col).isdigit() or (str(col).lstrip('-').isdigit())]
    order = sims[order_columns].copy(deep=True)
    order_prl = order.copy(deep=True)

    # --- Infection-step cutoff, as requested ---  # ADDED
    if up_to_step_infection is not None:  # ADDED
        def _mask_late_infections(row):  # ADDED
            # Values > cutoff are set to NaN and will be treated as "never infected"  # ADDED
            return row.where(row <= up_to_step_infection, np.nan)  # ADDED
        order = order.apply(_mask_late_infections, axis=1)      # ADDED
        order_prl = order_prl.apply(_mask_late_infections, axis=1)  # ADDED

    # Original behavior: push NaNs to (row.max() + 1) so "never infected" are ranked latest
    # All NaNs (including the ones we just introduced for > up_to_step_infection) are handled here.
    order = order.apply(lambda row: row.fillna(row.max() + 1), axis=1)

    # Precompute arrays for multiprocessing  # ADDED (multiprocessing)
    order_prl_np = order_prl.to_numpy()
    node_list = np.array(list(graph.nodes()))
    deg_list = np.array(list(dict(graph.degree()).values()))

    EPH = []
    corr = []

    # Build argument list for each row  # ADDED (multiprocessing)
    args_list = []
    for row_index in range(order_prl_np.shape[0]):
        args_list.append((
            order_prl_np[row_index, :],  # row of infection times (with NaNs)
            graph,
            node_list,
            deg_list
        ))

    # Use multiprocessing to compute EPH and corr per row  # ADDED (multiprocessing)
    with ProcessPoolExecutor() as executor:
        results_iter = executor.map(_compute_eph_corr_row, args_list)
        for eph_value, rho in tqdm(results_iter, total=len(args_list)):
            EPH.append(eph_value)
            corr.append(rho)

    sims['EPH'] = EPH
    sims['corr'] = corr

    sims.to_csv(output_csv, index=False)
    print(f"Updated CSV saved to {output_csv}")

# --- Main driver ---
def main():
    # Network and CSV details
    network_list = ['email_eu_modified']


    # Simple contagion parameters
    simple_betas = [0.02, 0.03, 0.04, 0.05]

    # # Complex contagion parameters
    # complex_thetas = [2, 3, 4, 5, 6]
    # complex_probs = [-1, 0.02, 0.03, 0.04]
    #

    complex_thetas = [4,5,6]
    complex_probs = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
    # Infection-step cutoffs to consider (treat later infections as "never infected")  # ADDED
    up_to_step_infections = [ 2, 3, 4, 5, 6, 7, 8, 9, 10]  # ADDED

    network_dir = '../networks/unweighted'
    results_dir = '../result_email_26feb/unweighted'


    for network_name in network_list:
        # Load and preprocess graph
        graph_path = os.path.join(network_dir, f'G_unweighted_{network_name}.graphml')
        G = nx.read_graphml(graph_path)

        # GraphML nodes are strings; convert to integer labels for consistency
        mapping = {node: int(float(node)) for node in G.nodes()}
        G = nx.relabel_nodes(G, mapping)

        # # Process simple contagion CSVs
        # for beta in simple_betas:
        #     csv_path = os.path.join(
        #         results_dir,
        #         f'unweighted_simple_{network_name}_beta_{beta}.csv'
        #     )
        #     if os.path.exists(csv_path):
        #         for cutoff in up_to_step_infections:  # ADDED
        #             add_eph_and_corr_to_csv(csv_path, G, up_to_step_infection=cutoff)  # ADDED

        # Process complex contagion CSVs
        for theta in complex_thetas:
            for prob_complex in complex_probs:
                csv_path = os.path.join(
                    results_dir,
                    f'unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv'
                )
                if os.path.exists(csv_path):
                    for cutoff in up_to_step_infections:  # ADDED
                        add_eph_and_corr_to_csv(csv_path, G, up_to_step_infection=cutoff)  # ADDED

if __name__ == "__main__":
    main()
