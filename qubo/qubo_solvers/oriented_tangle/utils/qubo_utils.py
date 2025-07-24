import numpy as np
import networkx as nx
from itertools import product
from math import floor
from qubo_solvers.logging import get_logger

logger = get_logger(__name__)


def edge2node_qubo_matrix_from_graph(graph: nx.DiGraph, alpha: float | None=None, penalties: list | None=None) -> tuple[np.ndarray, float, int, int]:
    """Constructs the QUBO matrix corresponding to a graph. Also returns the offset of the model, the max time and the number of nodes.

    Args:
        graph (nx.DiGraph): the node-weighted graph describing the problem.
        alpha (float, optional): the proportion of extra time allowed to paths over the maximum weight.

    Returns:
        tuple[np.ndarray, float, int, int]: qubo_matrix, offset, T_max, V
    """
    nodes = list(graph.nodes)
    V = int(len(nodes) / 2)
    total_weight = int(sum(graph.nodes[node]["weight"] for node in nodes))
    
    # T_max = total weight + "a bit"
    if alpha is None:
        alpha = 1.1
    T_max = floor(total_weight * alpha)
    logger.info(f'V: {V}, T: {T_max}')
    
    if penalties is None:
        # Penalty Values
        lambda_t = 100
        lambda_g = 50
        lambda_w = 1    
    else:
        lambda_t = penalties[0]
        lambda_g = penalties[1]
        lambda_w = penalties[2]
    logger.info(f'Penalties. t: {lambda_t}, g: {lambda_g}, w: {lambda_w}')
    
    # Note: we add an end node with parity 0 and 1, we only want 1 of them. We will delete the other at the end.
    qubo_matrix = np.zeros((T_max, V + 1, 2, T_max, V + 1, 2), dtype=np.int16)
    
    # Path constraint
    for t in range(T_max):
        for i, b in product(range(V), range(2)):
            qubo_matrix[t, i, b, t, i, b] -= lambda_t
            qubo_matrix[t, V, 0, t, i, b] += 2 * lambda_t
                
        qubo_matrix[t, V, 0, t, V, 0] -= lambda_t
        
        for i, j, bi, bj in product(range(V), range(V), range(2), range(2)):
            if not (i == j and bi == bj):
                qubo_matrix[t, i, bi, t, j, bj] += lambda_t
    
    # Set start/end nodes
    start_nodes=set()
    end_nodes= set()
    for node, val in dict(graph.nodes.data('start')).items():
        if val == 'start':
            print(f'Found start node:{node}')
            node_index_in_qubo = floor(nodes.index(node)/ 2)
            start_nodes.add(node_index_in_qubo)
        if val == 'end':
            print(f'Found end node:{node}')
            node_index_in_qubo = floor(nodes.index(node)/ 2)
            end_nodes.add(node_index_in_qubo)
    
    start_nodes = list(start_nodes)        
    end_nodes = list(end_nodes)
    print(f'Start nodes: {start_nodes}, End nodes: {end_nodes}')
    exist_start_nodes = len(start_nodes) > 0
    exist_end_nodes = len(end_nodes) > 0
    
    # Graph step constraints
    for t in range(T_max - 1):
        for i, j, bi, bj in product(range(V), range(V), range(2), range(2)):
            if (nodes[2 * i + bi], nodes[2 * j + bj]) not in graph.edges:
                qubo_matrix[t, i, bi, t+1, j, bj] += lambda_g
        for i, bi in product(range(V), range(2)):
            qubo_matrix[t, V, 0, t+1, i, bi] += 5 * lambda_g
            if exist_end_nodes:
                if i not in end_nodes:
                    qubo_matrix[t, i, bi, t+1, V, 0] += lambda_g
    if exist_start_nodes:
        start_node = start_nodes[0]
        for b in range(2):
            qubo_matrix[0, start_node, b, 0, start_node, b] -= lambda_g
            qubo_matrix[0, start_node, b, 0, start_nodes, 1 - b] += lambda_g
        
                
    # Weights constraints
    # ( sum_t ( x_{t,i,b} ) - w[i, b] ) ^2 = sum_t ( (1 - 2 * w[i, b]) x_{t,i,b} ) + sum_{t1 /= t2} x_{t1,i,b} x_{t2,i,b} + w[i,b]^2
    for i in range(V):
        for b in range(2):  
            for t in range(T_max):
                qubo_matrix[t, i, b, t, i, b] -= (2 * graph.nodes[nodes[2 * i + b]]["weight"] - 1) * lambda_w
        
            for t1, t2 in product(range(T_max), range(T_max)):
                if not (t1 == t2):
                    qubo_matrix[t1, i, b, t2, i, b] += lambda_w

    qubo_matrix = qubo_matrix.reshape((T_max * (V+1) * 2), (T_max * (V+1) * 2))
    qubo_matrix = 0.5 * (qubo_matrix + qubo_matrix.T)

    # Delete rows and columns corresponding to the extra end node we do not need
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 0)
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 1)
    
    offset = lambda_t * T_max  + lambda_w * int(sum(graph.nodes[nodes[2 * i + b]]["weight"] ** 2 for i in range(V) for b in range(2))) + (1 if exist_start_nodes else 0)  * lambda_g
    
    # normalisation = np.max(np.abs(qubo_matrix))
    # qubo_matrix = qubo_matrix / normalisation
    # offset = offset / normalisation
    
    return qubo_matrix, offset, T_max, V


def qubo_matrix_from_graph(graph: nx.DiGraph, alpha: float | None=None, penalties: list | None=None) -> tuple[np.ndarray, float, int, int]:
    """Constructs the QUBO matrix corresponding to a graph. Also returns the offset of the model, the max time and the number of nodes.

    Args:
        graph (nx.DiGraph): the node-weighted graph describing the problem.
        alpha (float, optional): the proportion of extra time allowed to paths over the maximum weight.

    Returns:
        tuple[np.ndarray, float, int, int]: qubo_matrix, offset, T_max, V
    """
    nodes = list(graph.nodes)
    V = int(len(nodes) / 2)
    total_weight = int(sum(graph.nodes[node]["weight"] for node in nodes) / 2)
    
    # T_max = total weight + "a bit"
    if alpha is None:
        alpha = 1.1
    T_max = floor(total_weight * alpha)
    logger.info(f'V: {V}, T: {T_max}')
    
    if penalties is None:
        # Penalty Values
        lambda_t = 100
        lambda_g = 50
        lambda_w = 1    
    else:
        lambda_t = penalties[0]
        lambda_g = penalties[1]
        lambda_w = penalties[2]
    logger.info(f'Penalties. t: {lambda_t}, g: {lambda_g}, w: {lambda_w}')
    
    # Note: we add an end node with parity 0 and 1, we only want 1 of them. We will delete the other at the end.
    qubo_matrix = np.zeros((T_max, V + 1, 2, T_max, V + 1, 2), dtype=np.int16)
    
    # Path constraint
    for t in range(T_max):
        for i, b in product(range(V), range(2)):
            qubo_matrix[t, i, b, t, i, b] -= lambda_t
            qubo_matrix[t, V, 0, t, i, b] += 2 * lambda_t
                
        qubo_matrix[t, V, 0, t, V, 0] -= lambda_t
        
        for i, j, bi, bj in product(range(V), range(V), range(2), range(2)):
            if not (i == j and bi == bj):
                qubo_matrix[t, i, bi, t, j, bj] += lambda_t
    
    # Set start/end nodes
    start_nodes=set()
    end_nodes= set()
    for node, val in dict(graph.nodes.data('start')).items():
        if val == 'start':
            print(f'Found start node:{node}')
            node_index_in_qubo = floor(nodes.index(node)/ 2)
            start_nodes.add(node_index_in_qubo)
        if val == 'end':
            print(f'Found end node:{node}')
            node_index_in_qubo = floor(nodes.index(node)/ 2)
            end_nodes.add(node_index_in_qubo)
    
    start_nodes = list(start_nodes)        
    end_nodes = list(end_nodes)
    print(f'Start nodes: {start_nodes}, End nodes: {end_nodes}')
    exist_start_nodes = len(start_nodes) > 0
    exist_end_nodes = len(end_nodes) > 0
    
    # Graph step constraints
    for t in range(T_max - 1):
        for i, j, bi, bj in product(range(V), range(V), range(2), range(2)):
            if (nodes[2 * i + bi], nodes[2 * j + bj]) not in graph.edges:
                qubo_matrix[t, i, bi, t+1, j, bj] += lambda_g
        for i, bi in product(range(V), range(2)):
            qubo_matrix[t, V, 0, t+1, i, bi] += 5 * lambda_g
            if exist_end_nodes:
                if i not in end_nodes:
                    qubo_matrix[t, i, bi, t+1, V, 0] += lambda_g
    if exist_start_nodes:
        start_node = start_nodes[0]
        for b in range(2):
            qubo_matrix[0, start_node, b, 0, start_node, b] -= lambda_g
            qubo_matrix[0, start_node, b, 0, start_nodes, 1 - b] += lambda_g
        
                
    # Weights constraints
    for i in range(V):
        for t in range(T_max):
            for b in range(2):
                qubo_matrix[t, i, b, t, i, b] -= (2 * graph.nodes[nodes[2 * i]]["weight"] - 1) * lambda_w
        
        for t1, t2 in product(range(T_max), range(T_max)):
            for b1, b2 in product(range(2), range(2)):
                if not (t1 == t2 and b1 == b2):
                    qubo_matrix[t1, i, b1, t2, i, b2] += lambda_w

    qubo_matrix = qubo_matrix.reshape((T_max * (V+1) * 2), (T_max * (V+1) * 2))
    qubo_matrix = 0.5 * (qubo_matrix + qubo_matrix.T)

    # Delete rows and columns corresponding to the extra end node we do not need
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 0)
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 1)
    
    offset = lambda_t * T_max  + lambda_w * int(sum(graph.nodes[nodes[2 * i]]["weight"] ** 2 for i in range(V))) + (1 if exist_start_nodes else 0)  * lambda_g
    
    # normalisation = np.max(np.abs(qubo_matrix))
    # qubo_matrix = qubo_matrix / normalisation
    # offset = offset / normalisation
    
    return qubo_matrix, offset, T_max, V