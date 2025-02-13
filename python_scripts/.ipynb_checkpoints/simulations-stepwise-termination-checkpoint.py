import networkx as nx
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import copy
import os


# --- Graph Loading and Preprocessing ---
def load_graph(file_path):
    """
    Load a graph from a GraphML file and preprocess it.

    Parameters:
    - file_path (str): Path to the GraphML file.

    Returns:
    - graph (networkx.Graph): Loaded and preprocessed graph.
    """
    G = nx.read_graphml(file_path)
    mapping = {node: int(float(node)) for node in G.nodes()}
    G = nx.relabel_nodes(G, mapping)
    return G

# --- Contagion Transfer ---
def perform_contagion_transfers_unweighted_version(graph, infected_list, mode='simple', threshold_complex=2, beta_simple=0.5, prob_complex=-1):
    """
    Simulates contagion transfers within a network for both simple and complex contagions.

    Parameters:
    - graph (networkx.Graph): The network graph over which the contagion spreads.
    - infected_list (list): List of nodes that are currently infected.
    - mode (str, optional): Specifies the mode of contagion ('simple' or 'complex'). Defaults to 'simple'.
    - threshold_complex (int, optional): Threshold for complex contagion. Defaults to 2.
    - beta_simple (float, optional): Transmission probability for simple contagion. Defaults to 0.5.
    - prob_complex (float, optional): Probability of infection in complex mode if neighbors are below the threshold. Defaults to -1.

    Returns:
    - (bool, list): A tuple containing a boolean for new infections and the updated list of infected nodes.
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
    return new_infection_marker, infected_list,n_simple_infection,n_complex_infection

# --- Simple Contagion ---
def perform_simple_contagion_unweighted(beta_simple, graph, initially_infected, stop_ratio):
    """
    Perform simple contagion process on a graph.

    Parameters:
    - beta_simple (float): Transmission probability for simple contagion.
    - graph (networkx.Graph): The network graph.
    - initially_infected (list): List of nodes initially infected.
    - stop_ratio (float): Portion of nodes required to be infected for termination.

    Returns:
    - np.array: Steps at which each node was infected.
    """
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    time_current = 0
    new_infections = True
    step_of_infection_per_node = [None] * len(graph.nodes)

    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current

    while time_current < stop_ratio:
        time_current += 1
        new_infections, infected_list, _,_ = perform_contagion_transfers_unweighted_version(
            graph, infected_list, mode='simple', beta_simple=beta_simple)

        for node_index in infected_list:
            node_pos = list(graph.nodes).index(node_index)
            if step_of_infection_per_node[node_pos] is None:
                step_of_infection_per_node[node_pos] = time_current

    return np.array(step_of_infection_per_node)

# --- Complex Contagion ---
def perform_complex_contagion_unweighted(threshold, graph, initially_infected, prob_complex, stop_ratio=None):
    """
    Perform complex contagion process on a graph.

    Parameters:
    - threshold (int): Threshold for complex contagion.
    - graph (networkx.Graph): The network graph.
    - initially_infected (list): List of nodes initially infected.
    - prob_complex (float): Probability of infection if neighbors are below the threshold.
    - stop_ratio (float, optional): Portion of nodes required to be infected for termination (applies if prob_complex != -1).

    Returns:
    - np.array: Steps at which each node was infected.
    """
    #for counting simple and complex infection
    complex_contagion_count = 0
    simple_contagion_count = 0 
    
    infected_list = initially_infected[:]
    total_nodes = len(graph.nodes())
    time_current = 0
    new_infections = True
    step_of_infection_per_node = [None] * len(graph.nodes)

    for node_index in initially_infected:
        node_pos = list(graph.nodes).index(node_index)
        step_of_infection_per_node[node_pos] = time_current
        

    if prob_complex == -1:
        
        while new_infections:
            time_current += 1
            new_infections, infected_list, n_simp, n_comp = perform_contagion_transfers_unweighted_version(
                graph, infected_list, mode='complex', threshold_complex=threshold, prob_complex=prob_complex)

            complex_contagion_count += n_comp
            simple_contagion_count += n_simp
            
            for node_index in infected_list:
                node_pos = list(graph.nodes).index(node_index)
                if step_of_infection_per_node[node_pos] is None:
                    step_of_infection_per_node[node_pos] = time_current
    elif prob_complex > 0:
        
        while time_current < stop_ratio:
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
        simple_infection_rate = simple_contagion_count/complex_contagion_count 

    return np.array(step_of_infection_per_node),simple_infection_rate

# --- Simulation Runner ---
def run_simulation(network_name, contagion_type, params, iterations, output_file):
    """
    Run contagion simulation and save the results.

    Parameters:
    - network_name (str): Name of the network.
    - contagion_type (str): 'simple' or 'complex'.
    - params (dict): Parameters for the contagion simulation.
    - iterations (int): Number of simulation iterations.
    - output_file (str): File to save the simulation results.

    Returns:
    - None
    """
    # Load graph
    graph_path = f'../networks/unweighted/G_unweighted_{network_name}.graphml'
    G = load_graph(graph_path)

    # Store results
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
            stop_ratio = params.get('stop_ratio', None)  # Include stop_ratio for complex contagion
            order,simple_rate = perform_complex_contagion_unweighted(threshold=threshold, graph=copy.deepcopy(G),
                                                         initially_infected=copy.deepcopy(list_inf), prob_complex=prob_complex, stop_ratio=stop_ratio)
            simple_ratio_infection_list.append(simple_rate)
            contagion_param_sim.append(threshold)

        order_collection_simulations.append(order)
        seed_sim.append(sample_it)

    result = pd.DataFrame(order_collection_simulations)
    result['seed'] = seed_sim
    result[contagion_type + '_param'] = contagion_param_sim
    result['simple_infection_ratio'] = simple_ratio_infection_list

    # Save results
    result.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


# --- Utilities ---
def init_nodes_var_reduction(iter_cnt, arr, k):
    """
    Initialize a random selection of nodes.

    Parameters:
    - iter_cnt (int): Iteration count for seed generation.
    - arr (list): List of nodes.
    - k (int): Number of nodes to select.

    Returns:
    - list: List of selected nodes.
    """
    rng = np.random.RandomState(iter_cnt)
    return list(rng.choice(arr, k, replace=False))

# --- Main Function ---
def main():
    # Define simulation parameters
    network_list = ['email_eu_modified', 'conf', 'school']
    network_list = ['email_eu_modified', 'conf', 'school']
#     network_list = ['barabasi_albert']
    iterations = 200  # Number of simulations per combination

    # Parameters for simple contagion
    simple_betas = [0.02, 0.03, 0.04, 0.05]
    simple_betas = [0.02, 0.03, 0.04, 0.05]
    
    stop_ratio = 3  # after 3 loops of infection

    # Parameters for complex contagion
    complex_thetas = [2, 3, 4, 5, 6, 7, 8, 9]
    complex_probs = [-1, 0.02, 0.03, 0.04]  # Include -1 for threshold-based contagion
    
    # Run simple contagion simulations
    for network_name in network_list:
        # Load graph
        graph_path = f'../networks/unweighted/G_unweighted_{network_name}.graphml'
        G = load_graph(graph_path)
        for beta in simple_betas:
            output_file = f'../result_3loopsstop_10%init_nodes-v3/unweighted/unweighted_simple_{network_name}_beta_{beta}.csv'
            if os.path.exists(output_file):
                print(output_file, "   already exist continue")
                continue 
            simple_params = {
                'beta': beta,
                'stop_ratio': stop_ratio,
                'initial_infected': int(np.round(len(G.nodes())*1/100))
            }
            run_simulation(
                network_name,
                contagion_type='simple',
                params=simple_params,
                iterations=iterations,
                output_file=output_file
            )

    # Run complex contagion simulations
    for network_name in network_list:
        graph_path = f'../networks/unweighted/G_unweighted_{network_name}.graphml'
        G = load_graph(graph_path)
        for theta in complex_thetas:
            for prob_complex in complex_probs:
                if prob_complex > 0:
                    output_file = f'../result_3loopsstop_10%init_nodes-v3/unweighted/unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv'
                    if os.path.exists(output_file):
                        print(output_file, "   already exist continue")
                        continue 
                    complex_params = {
                        'threshold': theta,
                        'prob_complex': prob_complex,
                        'stop_ratio': stop_ratio,  # Only applies when prob_complex != -1
                        'initial_infected': int(np.round(len(G.nodes())*1/100))
                    }
                    run_simulation(
                        network_name,
                        contagion_type='complex',
                        params=complex_params,
                        iterations=iterations,
                        output_file=output_file
                    )
                if prob_complex == -1:
                    output_file = f'../result_3loopsstop_10%init_nodes-v3/unweighted/unweighted_complex_{network_name}_theta_{theta}_prob_{prob_complex}.csv'
                    if os.path.exists(output_file):
                        print(output_file, "   already exist continue")
                        continue 
                    complex_params = {
                        'threshold': theta,
                        'prob_complex': prob_complex,
                        'stop_ratio': None,
                        'initial_infected': int(np.round(len(G.nodes())*1/100))
                    }
                    run_simulation(
                        network_name,
                        contagion_type='complex',
                        params=complex_params,
                        iterations=iterations,
                        output_file=output_file
                    )

if __name__ == "__main__":
    main()
