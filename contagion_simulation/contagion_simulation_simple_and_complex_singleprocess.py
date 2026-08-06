import networkx as nx
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import copy
import os

# --- Configuration: Define base directories ---
NETWORK_DIR = '../networks/unweighted'  # relative to this script's location; see networks/README.md for where to point this
RESULT_DIR = '../result_email_26feb/unweighted'

# --- Graph Loading and Preprocessing ---
def load_graph(file_path):
    """
    Load a graph from a GraphML file and preprocess node labels.
    """
    G = nx.read_graphml(file_path)
    # GraphML node labels are strings; here we convert them to ints for consistency.
    mapping = {node: int(float(node)) for node in G.nodes()}
    G = nx.relabel_nodes(G, mapping)
    return G

# --- Contagion Transfer Function (Simple & Complex) ---
def perform_contagion_transfers_unweighted_version(
    graph,
    infected_list,
    mode='simple',
    threshold_complex=2,
    beta_simple=0.5,
    prob_complex=-1
):
    """
    Perform one step of contagion process (simple or complex) on the given graph.

    Parameters
    ----------
    graph : networkx.Graph
        Underlying unweighted network.
    infected_list : list
        Current list of infected nodes.
    mode : {'simple', 'complex'}
        Determines which contagion rule to use.
    threshold_complex : int
        Threshold for complex contagion (number of infected neighbors).
    beta_simple : float
        Per-edge infection probability for simple contagion.
    prob_complex : float
        Sub-threshold infection probability for complex contagion.
        If -1, the probabilistic (below-threshold) channel is disabled.

    Returns
    -------
    new_infection_marker : bool
        True if at least one new node was infected in this step.
    infected_list : list
        Updated list of infected nodes (in-place extended).
    n_simple_infection : int
        Number of *below-threshold* (probabilistic) infections in this step.
    n_complex_infection : int
        Number of *above-threshold* (deterministic) infections in this step.
    infection_type_this_step : dict
        # ADDED: maps node -> code for this step (only for complex mode):
        #         1 = above-threshold adoption
        #         2 = below-threshold adoption
        #         (empty dict for simple mode)
    """
    new_infection_marker = False
    new_infected_nodes = []
    all_nodes = list(graph.nodes())
    n_complex_infection = 0
    n_simple_infection = 0

    infection_type_this_step = {}  # ADDED: per-step record of adoption type for newly infected nodes

    for node in all_nodes:
        if node in infected_list:
            # Skip already infected nodes; contagion is SI (no recovery).
            continue

        # Count how many neighbors are currently infected.
        infected_neighbors_count = sum(
            neighbor in infected_list for neighbor in graph.neighbors(node)
        )

        if mode == 'simple':
            # Simple contagion: independent transmission on each infected edge.
            infection_probability = 1 - (1 - beta_simple) ** infected_neighbors_count
            if random.uniform(0, 1) <= infection_probability:
                new_infection_marker = True
                new_infected_nodes.append(node)

        elif mode == 'complex':
            # Complex contagion:
            # Above-threshold adoption: deterministic once threshold is reached.
            if infected_neighbors_count >= threshold_complex:
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_complex_infection += 1
                infection_type_this_step[node] = 1  # ADDED: 1 = above-threshold adoption

            # Below-threshold adoption: probabilistic "simple-like" spillover.
            elif (
                0 < infected_neighbors_count
                and random.uniform(0, 1) <= prob_complex
            ):
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_simple_infection += 1
                infection_type_this_step[node] = 2  # ADDED: 2 = below-threshold adoption

    # Update the infected list in-place with new infections from this step.
    infected_list.extend(new_infected_nodes)

    return (
        new_infection_marker,
        infected_list,
        n_simple_infection,
        n_complex_infection,
        infection_type_this_step,  # ADDED: return per-step infection type info
    )

# --- Simple Contagion Process ---
def perform_simple_contagion_unweighted(beta_simple, graph, initially_infected, stop_ratio):
    """
    Simulate the simple contagion process until stop_ratio of nodes are infected.

    Returns
    -------
    step_of_infection_per_node : np.ndarray
        For each node (in the order of graph.nodes()), the time step at which
        it became infected (0 for initial seeds), or NaN if never infected.
    """
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    target_infected = int(stop_ratio * total_nodes)
    time_current = 0

    # Stores time of infection for each node (aligned with list(graph.nodes())).
    step_of_infection_per_node = [None] * len(graph.nodes)

    # Initial seeds are considered infected at time step 0.
    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current

    # Evolve until the desired fraction of nodes becomes infected.
    while len(infected_list) < target_infected:
        time_current += 1
        # For simple contagion, we ignore the extra return (infection type).
        new_infections, infected_list, _, _, _ = perform_contagion_transfers_unweighted_version(  # ADDED: unpack extra return
            graph,
            infected_list,
            mode='simple',
            beta_simple=beta_simple
        )

        # Update infection times for nodes that newly became infected in this step.
        for node_index in infected_list:
            node_pos = list(graph.nodes).index(node_index)
            if step_of_infection_per_node[node_pos] is None:
                step_of_infection_per_node[node_pos] = time_current

    return np.array(step_of_infection_per_node)

# --- Complex Contagion Process ---
def perform_complex_contagion_unweighted(threshold, graph, initially_infected, prob_complex, stop_ratio=None):
    """
    Simulate complex contagion with thresholding and optional probabilistic component.

    Parameters
    ----------
    threshold : int
        Complex contagion threshold (number of infected neighbors required
        for deterministic adoption).
    graph : networkx.Graph
    initially_infected : list
        Initial seed nodes.
    prob_complex : float
        Probability of below-threshold adoption.
        If -1, below-threshold channel is disabled.
    stop_ratio : float or None
        If not None and prob_complex != -1, stop once this fraction of nodes
        is infected. If prob_complex == -1, the process runs until no new
        infections occur.

    Returns
    -------
    step_of_infection_per_node : np.ndarray
        Time of infection for each node.
    simple_infection_rate : float
        Ratio (# below-threshold infections) / (# above-threshold infections),
        or 0 if no above-threshold infections occurred.
    infection_type_per_node : np.ndarray
        # ADDED: For each node (aligned with graph.nodes()):
        #         0 = initial seed or never infected via neighbors
        #         1 = infected via above-threshold (complex) adoption
        #         2 = infected via below-threshold (probabilistic) adoption
    """
    complex_contagion_count = 0
    simple_contagion_count = 0
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    target_infected = int(stop_ratio * total_nodes) if stop_ratio and prob_complex != -1 else None
    time_current = 0

    # Time of infection per node.
    step_of_infection_per_node = [None] * len(graph.nodes)
    # ADDED: Store infection type per node (0 = seed/never infected, 1 = above, 2 = below).
    infection_type_per_node = [0] * len(graph.nodes)  # ADDED

    # Initial seeds: infected at time 0, with infection type = 0 (external seeding).
    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current
        # infection_type_per_node[node_pos] remains 0 for seeds (no adoption event).

    if prob_complex == -1:
        # Pure threshold model (no below-threshold noise).
        while True:
            time_current += 1
            (
                new_infections,
                infected_list,
                n_simp,
                n_comp,
                infection_type_this_step,  # ADDED
            ) = perform_contagion_transfers_unweighted_version(  # ADDED: unpack infection_type_this_step
                graph,
                infected_list,
                mode='complex',
                threshold_complex=threshold,
                prob_complex=prob_complex
            )
            complex_contagion_count += n_comp
            simple_contagion_count += n_simp  # will be 0 in this regime

            # ADDED: record infection type for nodes newly infected in this step
            for node_index, infection_code in infection_type_this_step.items():  # ADDED
                node_pos = list(graph.nodes).index(node_index)  # ADDED
                if infection_type_per_node[node_pos] == 0:  # ADDED: only set if not set before
                    infection_type_per_node[node_pos] = infection_code  # ADDED

            if not new_infections:
                # No more new infections; process has reached a fixed point.
                break

            # Update infection times for newly infected nodes.
            for node_index in infected_list:
                node_pos = list(graph.nodes).index(node_index)
                if step_of_infection_per_node[node_pos] is None:
                    step_of_infection_per_node[node_pos] = time_current
    else:
        # Threshold + probabilistic below-threshold adoption, with a stopping ratio.
        while len(infected_list) < target_infected:
            time_current += 1
            (
                new_infections,
                infected_list,
                n_simp,
                n_comp,
                infection_type_this_step,  # ADDED
            ) = perform_contagion_transfers_unweighted_version(  # ADDED: unpack infection_type_this_step
                graph,
                infected_list,
                mode='complex',
                threshold_complex=threshold,
                prob_complex=prob_complex
            )
            complex_contagion_count += n_comp
            simple_contagion_count += n_simp

            # ADDED: record infection type for nodes newly infected in this step
            for node_index, infection_code in infection_type_this_step.items():  # ADDED
                node_pos = list(graph.nodes).index(node_index)  # ADDED
                if infection_type_per_node[node_pos] == 0:  # ADDED
                    infection_type_per_node[node_pos] = infection_code  # ADDED

            # Update infection times for newly infected nodes.
            for node_index in infected_list:
                node_pos = list(graph.nodes).index(node_index)
                if step_of_infection_per_node[node_pos] is None:
                    step_of_infection_per_node[node_pos] = time_current

    # Compute ratio of below-threshold to above-threshold infections.
    simple_infection_rate = 0
    if complex_contagion_count > 0:
        simple_infection_rate = simple_contagion_count / complex_contagion_count

    return (
        np.array(step_of_infection_per_node),
        simple_infection_rate,
        np.array(infection_type_per_node),  # ADDED: per-node infection type
    )

# --- Simulation Runner ---
def run_simulation(network_name, contagion_type, params, iterations, output_file):
    """
    Run multiple contagion simulations for a given network and save results.

    For 'simple' contagion, only the infection time matrix is saved.
    For 'complex' contagion, an additional CSV is saved that encodes, for
    each node and simulation, whether its adoption was above or below
    threshold.

    The additional CSV has the same rows as the main result:
      * columns 0..N-1: infection type code per node
      * 'seed' : simulation index
      * '<contagion_type>_param' : beta or theta
      * 'simple_infection_ratio' : as in the main CSV
    """
    graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
    G = load_graph(graph_path)

    order_collection_simulations = []        # Infection timing per simulation
    seed_sim = []                            # Seed index for each simulation
    contagion_param_sim = []                 # Beta (simple) or theta (complex)
    simple_ratio_infection_list = []         # Ratio of below-threshold to above-threshold infections

    infection_type_collection_simulations = []  # ADDED: per-simulation infection type matrices for complex contagion

    for sample_it in tqdm(range(iterations)):
        # Initialize seed nodes for this simulation (variance reduction via fixed RNG seed).
        list_inf = init_nodes_var_reduction(
            sample_it,
            list(G.nodes()),
            params.get("initial_infected", 10)
        )

        if contagion_type == 'simple':
            beta = params['beta']
            stop_ratio = params['stop_ratio']
            order = perform_simple_contagion_unweighted(
                beta_simple=beta,
                graph=copy.deepcopy(G),
                initially_infected=copy.deepcopy(list_inf),
                stop_ratio=stop_ratio
            )
            contagion_param_sim.append(beta)
            simple_ratio_infection_list.append(1)  # For simple contagion, treat ratio as 1 by convention.

            # For simple contagion, we do NOT record infection_type_collection_simulations
            # (logic preserved; extra CSV only for complex).

        elif contagion_type == 'complex':
            threshold = params['threshold']
            prob_complex = params['prob_complex']
            stop_ratio = params.get('stop_ratio', None)
            order, simple_rate, infection_type = perform_complex_contagion_unweighted(  # ADDED: infection_type returned
                threshold=threshold,
                graph=copy.deepcopy(G),
                initially_infected=copy.deepcopy(list_inf),
                prob_complex=prob_complex,
                stop_ratio=stop_ratio
            )
            contagion_param_sim.append(threshold)
            simple_ratio_infection_list.append(simple_rate)
            infection_type_collection_simulations.append(infection_type)  # ADDED: store per-node infection types

        order_collection_simulations.append(order)
        seed_sim.append(sample_it)

    # Main result: infection timing per node and simulation.
    result = pd.DataFrame(order_collection_simulations)
    result['seed'] = seed_sim
    result[contagion_type + '_param'] = contagion_param_sim
    result['simple_infection_ratio'] = simple_ratio_infection_list

    result.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

    # ADDED: For complex contagion, also save an infection-type CSV.
    if contagion_type == 'complex' and len(infection_type_collection_simulations) > 0:  # ADDED
        infection_type_df = pd.DataFrame(infection_type_collection_simulations)  # ADDED
        infection_type_df['seed'] = seed_sim  # ADDED: match meta columns with main CSV
        infection_type_df[contagion_type + '_param'] = contagion_param_sim  # ADDED
        infection_type_df['simple_infection_ratio'] = simple_ratio_infection_list  # ADDED

        # Name: original file name with '_infection_type' suffix.
        infection_type_output_file = output_file.replace('.csv', '_infection_type.csv')  # ADDED
        infection_type_df.to_csv(infection_type_output_file, index=False)  # ADDED
        print(f"Infection type details saved to {infection_type_output_file}")  # ADDED

# --- Helper for Node Initialization ---
def init_nodes_var_reduction(iter_cnt, arr, k):
    """
    Initialize k random seed nodes using a reproducible RNG based on iter_cnt.

    Parameters
    ----------
    iter_cnt : int
        Simulation index; used as RNG seed for variance reduction.
    arr : list
        List of candidate node IDs.
    k : int
        Number of initial infected nodes.

    Returns
    -------
    list
        List of chosen seed nodes.
    """
    rng = np.random.RandomState(iter_cnt)
    return list(rng.choice(arr, k, replace=False))

# --- Main Simulation Loop ---
def main():
    # Choose your networks here
    network_list = [
        'email_eu_modified'
    ]

    iterations = 100
    stop_ratio = 0.85

    # Simple contagion parameters
    simple_betas = [0.02, 0.03, 0.04, 0.05]

    # Complex contagion parameters
    complex_thetas = [2, 3, 4, 5, 6]
    complex_probs = [0.2,0.3,0.4,0.5,0.6]

    # # --- Run Simple Contagion Simulations ---
    # for network_name in network_list:
    #     graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
    #     G = load_graph(graph_path)
    #
    #     for beta in simple_betas:
    #         output_file = os.path.join(
    #             RESULT_DIR,
    #             f'unweighted_simple_{network_name}_beta_{beta}.csv'
    #         )
    #         if os.path.exists(output_file):
    #             # Preserve existing behavior: if the timing CSV exists, skip re-running.
    #             print(output_file, "already exists, skipping.")
    #             continue
    #         simple_params = {
    #             'beta': beta,
    #             'stop_ratio': stop_ratio,
    #             'initial_infected': int(np.round(len(G.nodes()) * 0.01))
    #         }
    #         run_simulation(network_name, 'simple', simple_params, iterations, output_file)

    # --- Run Complex Contagion Simulations ---
    for network_name in network_list:
        graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
        G = load_graph(graph_path)

        for theta in complex_thetas:
            for prob_complex in complex_probs:
                output_file = os.path.join(
                    RESULT_DIR,
                    f'unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv'
                )
                if os.path.exists(output_file):
                    # NOTE: This check is deliberately unchanged to preserve logic:
                    # if the main timing CSV exists, we skip the simulation, and
                    # therefore *also* skip generating the extra infection-type CSV.
                    print(output_file, "already exists, skipping.")
                    continue

                complex_params = {
                    'threshold': theta,
                    'prob_complex': prob_complex,
                    'stop_ratio': None if prob_complex == -1 else stop_ratio,
                    'initial_infected': int(np.round(len(G.nodes()) * 0.01))
                }
                run_simulation(network_name, 'complex', complex_params, iterations, output_file)

if __name__ == "__main__":
    main()
