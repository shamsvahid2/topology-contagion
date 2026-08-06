import networkx as nx
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import copy
import os
import math  # relative threshold support (phi -> ceil(phi * deg))

# --- Configuration: Define base directories ---
NETWORK_DIR = '../networks/unweighted'  # relative to this script's location; see networks/README.md for where to point this
RESULT_DIR = '../result_small_world/unweighted'

# --- Graph Loading and Preprocessing ---
def load_graph(file_path):
    """
    Load a graph from a GraphML file and preprocess node labels.
    """
    G = nx.read_graphml(file_path)
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

    For complex thresholds:
      - int >= 1  : absolute threshold (needs that many infected neighbors)
      - float in (0,1]: relative threshold φ (needs ceil(φ * degree(node)))
    Returns:
      new_infection_marker (bool),
      updated_infected_list (list),
      n_simple_infection (int),
      n_complex_infection (int),
      causes_dict (dict: node -> list of infected neighbors at adoption),
      adoption_type_dict (dict: node -> 1 if at/above threshold; 0 if below-threshold)
        * For mode='simple', adoption_type_dict is empty.
    """
    new_infection_marker = False
    new_infected_nodes = []
    causes_dict = {}          # node -> list[int] (its infected neighbors at adoption step)
    adoption_type_dict = {}   # node -> 1 (complex) or 0 (below threshold)
    all_nodes = list(graph.nodes())
    n_complex_infection = 0
    n_simple_infection = 0

    infected_set = set(infected_list)

    for node in all_nodes:
        if node in infected_set:
            continue

        # Count infected neighbors *at this step*, and keep the set for "causes"
        inf_neigh = [nbr for nbr in graph.neighbors(node) if nbr in infected_set]
        infected_neighbors_count = len(inf_neigh)

        if mode == 'simple':
            # Standard independent cascade per-step probability
            infection_probability = 1 - (1 - beta_simple) ** infected_neighbors_count
            if random.uniform(0, 1) <= infection_probability:
                new_infection_marker = True
                new_infected_nodes.append(node)
                causes_dict[node] = inf_neigh[:]  # record who was infected when this node adopted

        elif mode == 'complex':
            deg = graph.degree[node]

            # Interpret threshold as relative if (0,1], else absolute
            if isinstance(threshold_complex, float) and 0 < threshold_complex <= 1:
                if deg == 0:
                    threshold_count = deg + 1  # impossible
                else:
                    threshold_count = math.ceil(threshold_complex * deg)
            else:
                threshold_count = int(threshold_complex)

            if infected_neighbors_count >= max(threshold_count, 1):
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_complex_infection += 1
                causes_dict[node] = inf_neigh[:]
                adoption_type_dict[node] = 1  # adopted at/above threshold
            elif infected_neighbors_count > 0 and random.uniform(0, 1) <= prob_complex:
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_simple_infection += 1
                causes_dict[node] = inf_neigh[:]
                adoption_type_dict[node] = 0  # adopted below threshold (probabilistic)

    infected_list.extend(new_infected_nodes)
    return (
        new_infection_marker,
        infected_list,
        n_simple_infection,
        n_complex_infection,
        causes_dict,
        adoption_type_dict,
    )

# --- Simple Contagion Process ---
def perform_simple_contagion_unweighted(beta_simple, graph, initially_infected, stop_ratio):
    """
    Simulate the simple contagion process until stop_ratio of nodes are infected.
    Returns:
      step_of_infection_per_node: np.array[int or None]
      causes_per_node: list[str or None] serialized as "u|v|w" at adoption; "" for seeds
    """
    nodes_order = list(graph.nodes())
    node_index = {n: i for i, n in enumerate(nodes_order)}

    infected_list = initially_infected[:]
    infected_set = set(infected_list)
    total_nodes = len(nodes_order)
    target_infected = int(stop_ratio * total_nodes)
    time_current = 0
    step_of_infection_per_node = [None] * total_nodes
    causes_per_node = [None] * total_nodes

    # Seeds
    for n in initially_infected:
        pos = node_index[n]
        step_of_infection_per_node[pos] = time_current
        causes_per_node[pos] = ""  # seeds have no prior causes

    while len(infected_list) < target_infected:
        time_current += 1
        (
            new_infections,
            infected_list,
            _,
            _,
            causes_dict,
            _,
        ) = perform_contagion_transfers_unweighted_version(
            graph, infected_list, mode='simple', beta_simple=beta_simple
        )

        # Mark adoption time/causes only for genuinely-new nodes in this step
        for n in causes_dict.keys():
            pos = node_index[n]
            if step_of_infection_per_node[pos] is None:
                step_of_infection_per_node[pos] = time_current
                # store as pipe-joined sorted IDs for CSV readability
                cause_str = "|".join(str(x) for x in sorted(causes_dict[n]))
                causes_per_node[pos] = cause_str

        if not new_infections:
            # no new adoptions; break to avoid infinite loop (shouldn't happen with target ratio)
            break

    return np.array(step_of_infection_per_node), causes_per_node

# --- Complex Contagion Process ---
def perform_complex_contagion_unweighted(threshold, graph, initially_infected, prob_complex, stop_ratio=None):
    """
    Simulate complex contagion with thresholding (absolute or relative) and optional probabilistic component.

    threshold:
      - int >= 1  -> absolute threshold (needs that many infected neighbors)
      - float in (0,1] -> relative threshold φ (needs ceil(φ * degree(node)) infected neighbors)

    Returns:
      step_of_infection_per_node: np.array[int or None]
      simple_infection_rate: float (below-threshold/at-or-above-threshold)
      adoption_type_per_node: list[0/1/None] (1=at/above threshold OR seed; 0=below-threshold)
      causes_per_node: list[str or None], pipe-joined neighbors at adoption; "" for seeds
    """
    nodes_order = list(graph.nodes())
    node_index = {n: i for i, n in enumerate(nodes_order)}

    complex_contagion_count = 0
    simple_contagion_count = 0

    infected_list = initially_infected[:]
    total_nodes = len(nodes_order)
    target_infected = int(stop_ratio * total_nodes) if (stop_ratio and prob_complex != -1) else None
    time_current = 0
    step_of_infection_per_node = [None] * total_nodes
    adoption_type_per_node = [None] * total_nodes  # 1 complex or seed; 0 below-threshold; None never infected
    causes_per_node = [None] * total_nodes

    # Seeds (mark as complex-type=1 "otherwise" per your instruction)
    for n in initially_infected:
        pos = node_index[n]
        step_of_infection_per_node[pos] = time_current
        adoption_type_per_node[pos] = 1
        causes_per_node[pos] = ""

    def do_one_wave():
        nonlocal complex_contagion_count, simple_contagion_count, time_current, infected_list
        time_current += 1
        (
            new_infections,
            infected_list,
            n_simp,
            n_comp,
            causes_dict,
            adoption_type_dict,
        ) = perform_contagion_transfers_unweighted_version(
            graph, infected_list, mode='complex', threshold_complex=threshold, prob_complex=prob_complex
        )
        complex_contagion_count += n_comp
        simple_contagion_count += n_simp

        # Update first-adoption timestamps, types, and causes
        for n, tflag in adoption_type_dict.items():
            pos = node_index[n]
            if step_of_infection_per_node[pos] is None:
                step_of_infection_per_node[pos] = time_current
                adoption_type_per_node[pos] = int(tflag)
                cause_str = "|".join(str(x) for x in sorted(causes_dict.get(n, [])))
                causes_per_node[pos] = cause_str

        return new_infections

    if prob_complex == -1:
        # keep going until no new adoptions
        while True:
            new_infections = do_one_wave()
            if not new_infections:
                break
    else:
        # stop when reaching target ratio
        while len(infected_list) < target_infected:
            do_one_wave()

    simple_infection_rate = 0.0
    if complex_contagion_count > 0:
        simple_infection_rate = simple_contagion_count / complex_contagion_count

    return (
        np.array(step_of_infection_per_node),
        simple_infection_rate,
        adoption_type_per_node,
        causes_per_node,
    )

# --- Helper for Node Initialization ---
def init_nodes_var_reduction(iter_cnt, arr, k):
    """
    Initialize k random seed nodes using reproducible RNG.
    """
    rng = np.random.RandomState(iter_cnt)
    return list(rng.choice(arr, k, replace=False))

# --- Simulation Runner ---
def run_simulation(network_name, contagion_type, params, iterations, output_file):
    """
    Run multiple contagion simulations for a given network and save results.

    Also writes:
      - For complex: *_adoption_type.csv (0 below-threshold, 1 otherwise/seed)
      - For simple & complex: *_causes.csv (pipe-joined list of infected neighbors at adoption)
    """
    graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
    G = load_graph(graph_path)

    order_collection_simulations = []
    seed_sim = []
    contagion_param_sim = []
    simple_ratio_infection_list = []

    # new collectors
    adoption_type_rows = []  # complex only
    causes_rows = []         # both modes

    for sample_it in tqdm(range(iterations)):
        list_inf = init_nodes_var_reduction(sample_it, list(G.nodes()), params.get("initial_infected", 10))

        if contagion_type == 'simple':
            beta = params['beta']
            stop_ratio = params['stop_ratio']
            order, causes = perform_simple_contagion_unweighted(
                beta_simple=beta,
                graph=copy.deepcopy(G),
                initially_infected=copy.deepcopy(list_inf),
                stop_ratio=stop_ratio
            )
            contagion_param_sim.append(beta)
            simple_ratio_infection_list.append(1)
            causes_rows.append(causes)

        elif contagion_type == 'complex':
            threshold = params['threshold']  # int (absolute) or float in (0,1] (relative)
            prob_complex = params['prob_complex']
            stop_ratio = params.get('stop_ratio', None)
            order, simple_rate, adoption_type, causes = perform_complex_contagion_unweighted(
                threshold=threshold,
                graph=copy.deepcopy(G),
                initially_infected=copy.deepcopy(list_inf),
                prob_complex=prob_complex,
                stop_ratio=stop_ratio
            )
            contagion_param_sim.append(threshold)
            simple_ratio_infection_list.append(simple_rate)
            adoption_type_rows.append(adoption_type)
            causes_rows.append(causes)

        order_collection_simulations.append(order)
        seed_sim.append(sample_it)

    # --- Main result (unchanged) ---
    result = pd.DataFrame(order_collection_simulations)
    result['seed'] = seed_sim
    result[contagion_type + '_param'] = contagion_param_sim
    result['simple_infection_ratio'] = simple_ratio_infection_list
    result.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

    # --- Causes CSV (both modes) ---
    causes_df = pd.DataFrame(causes_rows)
    causes_df['seed'] = seed_sim
    causes_df[contagion_type + '_param'] = contagion_param_sim
    causes_df['simple_infection_ratio'] = simple_ratio_infection_list
    causes_path = output_file.replace('.csv', '_causes.csv')
    causes_df.to_csv(causes_path, index=False)
    print(f"Causes saved to {causes_path}")

    # --- Adoption-type CSV (complex only) ---
    if contagion_type == 'complex':
        adoption_df = pd.DataFrame(adoption_type_rows)
        adoption_df['seed'] = seed_sim
        adoption_df['complex_param'] = contagion_param_sim
        adoption_df['simple_infection_ratio'] = simple_ratio_infection_list
        adoption_path = output_file.replace('.csv', '_adoption_type.csv')
        adoption_df.to_csv(adoption_path, index=False)
        print(f"Adoption types saved to {adoption_path}")

# --- Main Simulation Loop ---
def main():
    # Choose your networks here
    network_list = []
    for seed in range(20):
        network_list.append(f'smallworld_p0.1_seed{seed}')
        network_list.append(f'smallworld_p0.2_seed{seed}')
        network_list.append(f'smallworld_p0.3_seed{seed}')
        network_list.append(f'smallworld_p0.4_seed{seed}')
        network_list.append(f'smallworld_p0.5_seed{seed}')

    iterations = 200
    stop_ratio = 0.85

    # Simple contagion parameters
    simple_betas = [0.02, 0.03, 0.04, 0.05]

    # Complex contagion parameters
    # Use ints for absolute or floats in (0,1] for relative thresholds
    complex_thetas = [0.1, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5]
    complex_probs = [0.02]

    # --- Run Simple Contagion Simulations ---
    for network_name in network_list:
        graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
        G = load_graph(graph_path)

        for beta in simple_betas:
            output_file = os.path.join(RESULT_DIR, f'unweighted_simple_{network_name}_beta_{beta}.csv')
            if os.path.exists(output_file):
                print(output_file, "already exists, skipping.")
                continue
            simple_params = {
                'beta': beta,
                'stop_ratio': stop_ratio,
                'initial_infected': int(np.round(len(G.nodes()) * 0.01))
            }
            run_simulation(network_name, 'simple', simple_params, iterations, output_file)

    # --- Run Complex Contagion Simulations ---
    for network_name in network_list:
        graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
        G = load_graph(graph_path)

        for theta in complex_thetas:
            for prob_complex in complex_probs:
                output_file = os.path.join(RESULT_DIR, f'unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv')
                if os.path.exists(output_file):
                    print(output_file, "already exists, skipping.")
                    continue

                complex_params = {
                    'threshold': theta,  # int (absolute) or float in (0,1] (relative)
                    'prob_complex': prob_complex,
                    'stop_ratio': None if prob_complex == -1 else stop_ratio,
                    'initial_infected': int(np.round(len(G.nodes()) * 0.01))
                }
                run_simulation(network_name, 'complex', complex_params, iterations, output_file)

if __name__ == "__main__":
    main()
