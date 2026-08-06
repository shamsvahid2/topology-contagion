import numpy as np
import pandas as pd
import scipy.stats
import gudhi as gd
import networkx as nx
import random
from tqdm import tqdm
import os

# --- Persistent Homology Computation ---
def compute_persistent_homology(graph, filtration):
    """
    Compute the persistent homology of a graph with a given filtration.

    :param graph: A networkx graph.
    :param filtration: A list of node filtration values, corresponding to each node in the graph.
    :return: Persistent homology of dimension one.
    """
    # Create simplex tree
    st = gd.SimplexTree()

    # Add vertices with filtration values
    for i, node in enumerate(graph.nodes):
        st.insert([node], filtration=filtration[i])

    # Add edges
    for edge in graph.edges:
        st.insert(list(edge), filtration=max(
            filtration[list(graph.nodes()).index(edge[0])],
            filtration[list(graph.nodes()).index(edge[1])]
        ))
    # Expand to clique complex
    st.extend_filtration()

    # Remove null cycles
    persistence = st.extended_persistence(min_persistence=1e-5)
    tmp = []
    tmp.extend(persistence[0])
    tmp.extend(persistence[1])
    tmp.extend(persistence[2])
    tmp.extend(persistence[3])
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

# --- Add EPH and Correlation ---
def add_eph_and_corr_to_csv(csv_path, graph, partial_ratio):
    """
    Add EPH and correlation (corr) columns to the given CSV file based on the graph.

    :param csv_path: Path to the CSV file.
    :param graph: The corresponding graph (networkx.Graph).
    """
    output_csv = csv_path.replace('.csv', f'limited-node_{partial_ratio}_EPH-v2.csv')
    if os.path.exists(output_csv):
        print(output_csv, " already exist continue")
        return -1

    sims = pd.read_csv(csv_path)

    # Remove extra columns based on expected CSV format
    order_columns = [col for col in sims.columns if str(col).isdigit() or (str(col).lstrip('-').isdigit())]

    #select observable nodes
    partial_n_nodes = int(partial_ratio*len(order_columns))
    seed_curr = np.random.randint(1e8)
    np.random.seed(seed_curr)
    order_columns_obs = np.random.choice(order_columns,partial_n_nodes,replace=False)

    order = sims[order_columns_obs].copy(deep=True)
    order_prl = order.copy(deep=True)
    order = order.apply(lambda row: row.fillna(row.max() + 1), axis=1)




    EPH = []
    corr = []

    for row_index in tqdm(range(len(order))):

#         print(1)
        not_nan_ind = np.where(~order_prl.iloc[row_index].isna())[0]
#         print(2)
        g_node_list = np.array(list(graph.nodes()))
#         print(3)
        sgraph = graph.subgraph(g_node_list[not_nan_ind])
#         print(4)
        filt = order_prl.iloc[row_index][not_nan_ind]
#         print(filt)
#         print(len(sgraph.nodes()))
#         print(5)
        sgraph = nx.Graph(sgraph)
        persistent_homology = compute_persistent_homology(graph=sgraph, filtration=-filt)  # FIXED: sign was missing (should be -filt, matching the other EPH scripts) — see README "Corrections applied"
#         print(6)


#         # Define filtration
#         filt = list(-order.iloc[row_index].values)

#         # Compute EPH
#         persistent_homology = compute_persistent_homology(graph=graph, filtration=filt)

        # Lifetime of generators
        life_set = []
        for gen in persistent_homology:
            life_set.append(gen[1][1] - gen[1][0])
        EPH.append(np.nanmean(life_set))

        # Compute PRL (Spearman correlation)
        deg_list = list(dict(graph.degree).values())
        order_curr = list(order_prl.iloc[row_index].values)

        try:
            rho, _ = scipy.stats.spearmanr(deg_list, order_curr, nan_policy='omit')
        except Exception as e:
            print(f'PRL method could not find rho for row {row_index}: {e}')
            rho = 0

        corr.append(rho)

    sims['EPH'] = EPH
    sims['corr'] = corr
    sims['seed_obs'] = seed_curr

    output_csv = csv_path.replace('.csv', f'limited-node_{partial_ratio}_EPH-v2.csv')
    sims.to_csv(output_csv, index=False)
    print(f"Updated CSV saved to {output_csv}")

# --- Main Function ---
def main():
    # Network and CSV details
    network_list = ['conf', 'school']
    network_list = ['email_eu_modified']

    simple_betas = [0.02, 0.03, 0.04, 0.05]
#     simple_betas = [0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.021, 0.03, 0.04, 0.05]
#     complex_thetas = [2, 3, 4, 5, 6, 7, 8, 9]
    complex_thetas = [2, 3, 4]
    complex_probs = [-1]
    observable_ratio = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13]
    observable_ratio = [0.11,0.12,0.13]
    results_dir = '../result_85stop_10%init_nodes_limited-v4/unweighted/'
    network_dir = '../networks/unweighted/'

    for network_name in network_list:
        # Load and preprocess graph
        graph_path = os.path.join(network_dir, f'G_unweighted_{network_name}.graphml')
        G = nx.read_graphml(graph_path)

        mapping = {node: int(float(node)) for node in G.nodes()}
        G = nx.relabel_nodes(G, mapping)

        # Process simple contagion CSVs
        for obs_ratio in observable_ratio:
            for beta in simple_betas:
                csv_path = os.path.join(results_dir, f'unweighted_simple_{network_name}_beta_{beta}.csv')
                if os.path.exists(csv_path):
                    add_eph_and_corr_to_csv(csv_path, G, obs_ratio)

        # Process complex contagion CSVs
        for theta in complex_thetas:
            for prob_complex in complex_probs:
                for obs_ratio in observable_ratio:
                    csv_path = os.path.join(results_dir, f'unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv')
                    if os.path.exists(csv_path):
                        add_eph_and_corr_to_csv(csv_path, G, obs_ratio)


if __name__ == "__main__":
    main()
