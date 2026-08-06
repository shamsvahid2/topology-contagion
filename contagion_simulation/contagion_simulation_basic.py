import networkx as nx
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import copy
import os

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
def perform_contagion_transfers_unweighted_version(graph, infected_list, mode='simple', threshold_complex=2, beta_simple=0.5, prob_complex=-1):
    """
    Perform one step of contagion process (simple or complex) on the given graph.
    """
    new_infection_marker = False
    new_infected_nodes = []
    all_nodes = list(graph.nodes())
    n_complex_infection = 0
    n_simple_infection = 0

    for node in all_nodes:
        if node in infected_list:
            continue

        infected_neighbors_count = sum(neighbor in infected_list for neighbor in graph.neighbors(node))

        if mode == 'simple':
            infection_probability = 1 - (1 - beta_simple)**infected_neighbors_count
            if random.uniform(0, 1) <= infection_probability:
                new_infection_marker = True
                new_infected_nodes.append(node)

        elif mode == 'complex':
            if infected_neighbors_count >= threshold_complex:
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_complex_infection += 1
            elif 0 < infected_neighbors_count and random.uniform(0, 1) <= prob_complex:
                new_infection_marker = True
                new_infected_nodes.append(node)
                n_simple_infection += 1

    infected_list.extend(new_infected_nodes)
    return new_infection_marker, infected_list, n_simple_infection, n_complex_infection

# --- Simple Contagion Process ---
def perform_simple_contagion_unweighted(beta_simple, graph, initially_infected, stop_ratio):
    """
    Simulate the simple contagion process until stop_ratio of nodes are infected.
    """
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    target_infected = int(stop_ratio * total_nodes)
    time_current = 0
    step_of_infection_per_node = [None] * len(graph.nodes)

    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current

    while len(infected_list) < target_infected:
        time_current += 1
        new_infections, infected_list, _, _ = perform_contagion_transfers_unweighted_version(
            graph, infected_list, mode='simple', beta_simple=beta_simple)

        for node_index in infected_list:
            node_pos = list(graph.nodes).index(node_index)
            if step_of_infection_per_node[node_pos] is None:
                step_of_infection_per_node[node_pos] = time_current

    return np.array(step_of_infection_per_node)

# --- Complex Contagion Process ---
def perform_complex_contagion_unweighted(threshold, graph, initially_infected, prob_complex, stop_ratio=None):
    """
    Simulate complex contagion with thresholding and optional probabilistic component.
    """
    complex_contagion_count = 0
    simple_contagion_count = 0
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    target_infected = int(stop_ratio * total_nodes) if stop_ratio and prob_complex != -1 else None
    time_current = 0
    step_of_infection_per_node = [None] * len(graph.nodes)

    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current

    if prob_complex == -1:
        while True:
            time_current += 1
            new_infections, infected_list, n_simp, n_comp = perform_contagion_transfers_unweighted_version(
                graph, infected_list, mode='complex', threshold_complex=threshold, prob_complex=prob_complex)
            complex_contagion_count += n_comp
            simple_contagion_count += n_simp
            if not new_infections:
                break
            for node_index in infected_list:
                node_pos = list(graph.nodes).index(node_index)
                if step_of_infection_per_node[node_pos] is None:
                    step_of_infection_per_node[node_pos] = time_current
    else:
        while len(infected_list) < target_infected:
            time_current += 1
            new_infections, infected_list, n_simp, n_comp = perform_contagion_transfers_unweighted_version(
                graph, infected_list, mode='complex', threshold_complex=threshold, prob_complex=prob_complex)
            complex_contagion_count += n_comp
            simple_contagion_count += n_simp
            for node_index in infected_list:
                node_pos = list(graph.nodes).index(node_index)
                if step_of_infection_per_node[node_pos] is None:
                    step_of_infection_per_node[node_pos] = time_current

    simple_infection_rate = 0
    if complex_contagion_count > 0:
        simple_infection_rate = simple_contagion_count / complex_contagion_count

    return np.array(step_of_infection_per_node), simple_infection_rate

# --- Simulation Runner ---
def run_simulation(network_name, contagion_type, params, iterations, output_file):
    """
    Run multiple contagion simulations for a given network and save results.
    """
    graph_path = os.path.join(NETWORK_DIR, f'G_unweighted_{network_name}.graphml')
    G = load_graph(graph_path)

    order_collection_simulations = []
    seed_sim = []
    contagion_param_sim = []
    simple_ratio_infection_list = []

    for sample_it in tqdm(range(iterations)):
        list_inf = init_nodes_var_reduction(sample_it, list(G.nodes()), params.get("initial_infected", 10))

        if contagion_type == 'simple':
            beta = params['beta']
            stop_ratio = params['stop_ratio']
            order = perform_simple_contagion_unweighted(beta_simple=beta, graph=copy.deepcopy(G),
                                                        initially_infected=copy.deepcopy(list_inf), stop_ratio=stop_ratio)
            contagion_param_sim.append(beta)
            simple_ratio_infection_list.append(1)

        elif contagion_type == 'complex':
            threshold = params['threshold']
            prob_complex = params['prob_complex']
            stop_ratio = params.get('stop_ratio', None)
            order, simple_rate = perform_complex_contagion_unweighted(threshold=threshold, graph=copy.deepcopy(G),
                                                                      initially_infected=copy.deepcopy(list_inf),
                                                                      prob_complex=prob_complex, stop_ratio=stop_ratio)
            contagion_param_sim.append(threshold)
            simple_ratio_infection_list.append(simple_rate)

        order_collection_simulations.append(order)
        seed_sim.append(sample_it)

    result = pd.DataFrame(order_collection_simulations)
    result['seed'] = seed_sim
    result[contagion_type + '_param'] = contagion_param_sim
    result['simple_infection_ratio'] = simple_ratio_infection_list

    result.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# --- Helper for Node Initialization ---
def init_nodes_var_reduction(iter_cnt, arr, k):
    """
    Initialize k random seed nodes using reproducible RNG.
    """
    rng = np.random.RandomState(iter_cnt)
    return list(rng.choice(arr, k, replace=False))

# --- Main Simulation Loop ---
def main():
    # Choose your networks here
    network_list = [
 'smallworld_p0.1_seed0']

    iterations = 200
    stop_ratio = 0.85

    # Simple contagion parameters
    simple_betas = [0.02, 0.03, 0.04, 0.05]

    # Complex contagion parameters
    complex_thetas = [2, 3, 4, 5, 6, 7, 8, 9]
    complex_probs = [-1, 0.02, 0.03, 0.04]

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
                    'threshold': theta,
                    'prob_complex': prob_complex,
                    'stop_ratio': None if prob_complex == -1 else stop_ratio,
                    'initial_infected': int(np.round(len(G.nodes()) * 0.01))
                }
                run_simulation(network_name, 'complex', complex_params, iterations, output_file)

if __name__ == "__main__":
    main()
