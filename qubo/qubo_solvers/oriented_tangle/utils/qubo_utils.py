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
    total_weight = sum(graph.nodes[node]["weight"] for node in nodes)
    
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
    qubo_matrix = np.zeros((T_max, V + 1, 2, T_max, V + 1, 2), dtype=float)
    
    # Path constraint: for each t
    # ( sum_{i,b} ( x_{t,i,b} ) - 1 ) ^2 
    # = sum_{(i1,b1), (i2,b2)} ( x_{t,i1,b1} x_{t,i2,b2} ) + sum_{i,b} ( - 2 * x_{t,i,b} )  +  1        
    T_qubo_matrix = lambda_t * (
        np.ones((V+1, 2, V+1, 2), dtype=float) # + sum_{(i1,b1),(i2,b2)} x_{t,i1,b1} x_{t,i2,b2}
        - 2 * np.array([1]+ ([0]*2*(V+1)+[1]) * (2*(V+1)-1), dtype=float).reshape((V+1,2,V+1,2)) # -2 * sum_{i,b} ( - x_{t,i,b} )
    )
    for t in range(T_max):
        qubo_matrix[t, :, :, t, :, :] += T_qubo_matrix
    
    
    # Graph step constraints: for each (t, t+1)
    # 1 - sum_{i1, b1} ( x_{t,i1,b1} * ( sum_{i2,b2 neighbour i1, b1} x_{t+1,i2,b2} ) )
    G_qubo_matrix = np.zeros((V+1, 2, V+1, 2), dtype=float)
    G_qubo_matrix[:, :, V, 0] = 0
    for u, v in graph.edges:
        u_idx = nodes.index(u)
        i, bi = u_idx // 2, u_idx % 2
        v_idx = nodes.index(v)
        j, bj = v_idx // 2, v_idx % 2
        G_qubo_matrix[i, bi, j, bj] = 0
    for t in range(T_max - 1):
        qubo_matrix[t, :, :, t+1, :, :] = G_qubo_matrix
            
                
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
    
    offset = lambda_t * T_max  + lambda_w * sum(graph.nodes[nodes[2 * i + b]]["weight"] ** 2 for i in range(V) for b in range(2))
    
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
    total_weight = sum(graph.nodes[node]["weight"] for node in nodes) / 2
    
    # T_max = total weight + "a bit"
    if alpha is None:
        alpha = 1.1
    T_max = max(floor(total_weight * alpha), 1)
    logger.info(f'V: {V}, T: {T_max}')
    
    if penalties is None:
        # Penalty Values
        lambda_t = 200
        lambda_g = 50
        lambda_w = 1    
    else:
        lambda_t = penalties[0]
        lambda_g = penalties[1]
        lambda_w = penalties[2]
    logger.info(f'Penalties. t: {lambda_t}, g: {lambda_g}, w: {lambda_w}')
    
        # Note: we add an end node with parity 0 and 1, we only want 1 of them. We will delete the other at the end.
    qubo_matrix = np.zeros((T_max, V + 1, 2, T_max, V + 1, 2), dtype=float)
    
    # Path constraint: for each t
    # ( sum_{i,b} ( x_{t,i,b} ) - 1 ) ^2 
    # = sum_{(i1,b1), (i2,b2)} ( x_{t,i1,b1} x_{t,i2,b2} ) + sum_{i,b} ( - 2 * x_{t,i,b} )  +  1        
    T_qubo_matrix = lambda_t * (
        np.ones((V+1, 2, V+1, 2), dtype=float) # + sum_{(i1,b1),(i2,b2)} x_{t,i1,b1} x_{t,i2,b2}
        - 2 * np.array([1]+ ([0]*2*(V+1)+[1]) * (2*(V+1)-1), dtype=float).reshape((V+1,2,V+1,2)) # -2 * sum_{i,b} ( - x_{t,i,b} )
    )
    for t in range(T_max):
        qubo_matrix[t, :, :, t, :, :] += T_qubo_matrix
    
    
    # Graph step constraints: for each (t, t+1)
    # 1 - sum_{i1, b1} ( x_{t,i1,b1} * ( sum_{i2,b2 neighbour i1, b1} x_{t+1,i2,b2} ) ) - 0.8 *  x_{t, end 0} * ( sum_{i, b} x_{t, i , b} )
    # The middle term means the "1" is cancelled if x_{t,i1,b1} == 1 and x_{t, i2, b2} == 1 for a valid graph step 
    # The last term allows new contigs to be started at low cost
    G_qubo_matrix = np.zeros((V+1, 2, V+1, 2), dtype=float)
    G_qubo_matrix[:, :, V, 0] = -1 * lambda_g # End node is always a neighbour
    for u, v in graph.edges:
        u_idx = nodes.index(u)
        i, bi = u_idx // 2, u_idx % 2
        v_idx = nodes.index(v)
        j, bj = v_idx // 2, v_idx % 2
        G_qubo_matrix[i, bi, j, bj] = -1 * lambda_g # - sum_{i1, b1} ( x_{t,i1,b1} * ( sum_{i2,b2 neighbour i1, b1} x_{t+1,i2,b2} ) )
    # for i, b in product(range(V), range(2)):
    #     G_qubo_matrix[V, 0, i, b] = -0.8 * lambda_g # - 0.8 *  x_{t, end 0} * ( sum_{i, b} x_{t, i , b} )
    
    for t in range(T_max - 1):
        qubo_matrix[t, :, :, t+1, :, :] = G_qubo_matrix 
    
     
    # Weight constraint
    # for each i
    # ( sum_{t,b} x_{t,i,b} - w(i) ) * ( sum_{t,b} x_{t,i,b} - w(i) + 0.5) 
    # = w(i)(w(i) - 0.5) - (2w(i) - 0.5) * sum_{t,b} ( x_{t,i,b} ) + sum_{t1,b1,t2,b2} ( x_{t1,i,b1} x_{t2,i,b2})
    
    # UPDATE: multiply by log(len(i)) to reflect the fact that we have more confidence in the weights of longer nodes
    def W_qubo_matrix(weight, length):
        return lambda_w * (
                np.ones((T_max, 2, T_max, 2), dtype=float) # sum_{t1,b1,t2,b2} ( x_{t1,i,b1} x_{t2,i,b2})
                - (2 * weight - 0.5) * np.array([1]+ ([0]*2*(T_max)+[1]) * (2*(T_max)-1), dtype=float).reshape((T_max,2,T_max,2)) # - (2w(i) - 0.5) * sum_{t,b} ( x_{t,i,b} )
            )
        # return np.log10(length) * lambda_w * (
        #         np.ones((T_max, 2, T_max, 2), dtype=float) # sum_{t1,b1,t2,b2} ( x_{t1,i,b1} x_{t2,i,b2})
        #         - (2 * weight - 0.5) * np.array([1]+ ([0]*2*(T_max)+[1]) * (2*(T_max)-1), dtype=float).reshape((T_max,2,T_max,2)) # - (2w(i) - 0.5) * sum_{t,b} ( x_{t,i,b} )
        #     )
        
        
    for i in range(V):
        node = graph.nodes[nodes[2*i]]
        qubo_matrix[:, i, :, :, i, :] += W_qubo_matrix(node["weight"], node["length"])
        

    qubo_matrix = qubo_matrix.reshape((T_max * (V+1) * 2), (T_max * (V+1) * 2))
    qubo_matrix = 0.5 * (qubo_matrix + qubo_matrix.T)

    # Delete rows and columns corresponding to the extra end node we do not need
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 0)
    qubo_matrix = np.delete(qubo_matrix, [np.ravel_multi_index((t, V, 1), dims=(T_max, V+1, 2)) for t in range(T_max)], 1)
    
    offset = (
        lambda_t * T_max  
        + lambda_g * (T_max - 1)
        + lambda_w * sum([
            (graph.nodes[nodes[2 * i]]["weight"] ** 2 - 0.5 * graph.nodes[nodes[2 * i]]["weight"]) # * np.log10(graph.nodes[nodes[2 * i]]["length"])
        for i in range(V)])
    )
    
    # normalisation = np.max(np.abs(qubo_matrix))
    # qubo_matrix = qubo_matrix / normalisation
    # offset = offset / normalisation
    
    return qubo_matrix, offset, T_max, V