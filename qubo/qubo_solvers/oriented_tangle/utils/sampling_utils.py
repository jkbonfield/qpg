import numpy as np
import networkx as nx
import re
import gurobipy as gp
import subprocess
from qubo_solvers.definitions import QuboDescription
from qubo_solvers.logging import get_logger
from gurobipy import GRB
from math import floor

logger = get_logger(__name__)

rng = np.random.default_rng()

def mqlib_sample_qubo(qubo_description: QuboDescription):
    input_filepath = f"{qubo_description.data_dir}/mqlib_input_{qubo_description.filename}.txt"

    paths = {}
    for time_limit in qubo_description.time_limits:
        paths[time_limit] = []
        
        for _ in range(qubo_description.jobs):
            if qubo_description.T < 5:
                logger.info(f'Small problem, T = {qubo_description.T}. Setting time limit to <=5')
                time_limit = min(time_limit, 5)
            elif qubo_description.T < 10:
                logger.info(f'Small problem, T = {qubo_description.T}. Setting time limit to <=10')
                time_limit = min(time_limit, 10)
            # Run the MQLib solver and capture output
            process = subprocess.run(["MQLib", "-fQ", input_filepath, "-h", "PALUBECKIS2004bMST2", "-r", str(time_limit), "-s", str(rng.integers(0, 65535)), "-ps"], capture_output=True)

            out = process.stdout.decode("utf-8")

            try:
                # First line of output includes run data. 3rd line contains the solution.
                out_data = [x for x in out.split('\n') if len(x) > 0]
                solution = out_data[2].split()
                solution = [int(x) for x in solution]
                logger.info(out_data[0].split(','))
                logger.info(out_data[0].split(',')[-1])
                solution_energy = float(out_data[0].split(',')[3])
            except (ValueError, IndexError):
                logger.error('Could not parse mqlib data')
                logger.error(out)
                paths[time_limit].append(([], np.inf, []))
                continue
            energy = qubo_description.offset - solution_energy
            path = sample_list_to_path(solution, qubo_description.graph, qubo_description.T, qubo_description.V)
            paths[time_limit].append((solution, energy, path))
            
    return paths


def dwave_sample_qubo(qubo_description: QuboDescription) -> dict[int, tuple]:
    """Perform a batch of annealing on a given Binary Quadratic Model.

    Args:
        qubo_description (QuboDescription): a description of the problem
        
    Returns:
        (dict): Returns the best sample, energy and path for each job run.
    """
    
    from dimod import BQM
    from dwave.system import LeapHybridSampler
    bqm = BQM(qubo_description.Q, 'BINARY')
    bqm.offset = qubo_description.offset
    sampler = LeapHybridSampler()
    
    paths = {}
    for time_limit in qubo_description.time_limits:
        paths[time_limit] = []
        for _ in range(qubo_description.jobs):
            sampleset = sampler.sample(bqm, time_limit, label=f'oriented_{qubo_description.filename}')
            try:
                logger.info(f"D-Wave access time: {round(sampleset.info['run_time'] / 10 ** 6)}")
            except KeyError:
                pass
            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            path = sample_list_to_path(np.array(list(best_sample.values())), qubo_description.graph, qubo_description.T, qubo_description.V)
            paths[time_limit].append((best_sample, best_energy, path))
            
    return paths


def gurobi_sample_qubo(qubo_description: QuboDescription):
    total_weight = int(sum(qubo_description.graph.nodes[node]["weight"] for node in list(qubo_description.graph.nodes)) / 2)
    
    logger.info(f'Offset: {qubo_description.offset}')
    logger.info(f'Total weight: {total_weight}')
    logger.info(f'T_max: {qubo_description.T}')
    
    paths = {}
    Q = np.array(qubo_description.Q)
    with gp.Env() as env, gp.Model(env=env) as model:
        model_vars = model.addMVar(shape=Q.shape[0], vtype=GRB.BINARY, name="x")
        model.setObjective(model_vars @ Q @ model_vars, GRB.MINIMIZE)
        model.Params.BestObjStop = - qubo_description.offset
        
        for time_limit in qubo_description.time_limits:
            paths[time_limit] = []
            model.Params.TimeLimit = time_limit
            for _ in range(qubo_description.jobs):
                model.Params.Seed = rng.integers(0, 100000)
                model.optimize()
                energy = model.ObjVal + qubo_description.offset
                path = sample_list_to_path(model_vars.X, qubo_description.graph, qubo_description.T, qubo_description.V)
                paths[time_limit].append((model_vars.X, energy, path))
                model.reset()
    
    return paths


def sample_array_to_path(sample_array: np.ndarray, nodes: list, V: int):
    nz = np.nonzero(sample_array == 1)
    return [
        (
            int(nz[0][i]), 
            nodes[nz[1][i] * 2 + nz[2][i]] if nz[1][i] in range(V) else 'end'
        ) for i in range(nz[0].shape[0])
    ]


def sample_list_to_path(sample_list: np.ndarray, graph: nx.Graph, T_max: int, V: int):
    for idx in [t * (V + 1) * 2 + V * 2 + 1 for t in range(T_max)]:
        sample_list = np.insert(sample_list, idx, 0)
    sample_array = sample_list.reshape((T_max, V + 1, 2))
    return sample_array_to_path(sample_array, list(graph.nodes), V)
    

def print_path(path: list):
    """Pretty print a path"""
    num_per_line = 6
    if len(path) < num_per_line:
        print(path)
        return
    
    for i in range(floor(len(path) / num_per_line)):
        print(path[i * num_per_line: (i + 1) * num_per_line])
    if not (i + 1) * num_per_line == len(path):
        print(path[(i + 1)*num_per_line:])
        
        
def get_original_vertex_name(vertex_name):
    pattern = r'(.+)_([\+\-])+'
    match = re.search(pattern, vertex_name)
    if match is None:
        raise Exception('Could not retrieve vertex name')
    else:
        return match.group(1)
        

def validate_path(path: list, graph: nx.Graph):
    """Checks the constraints for a path on a graph.
    
    In particular:
     - does the path go along graph edges at each time step
     - is each node visited the correct number of times
     - is exactly one node visited per time step

    Args:
        path (list): _description_
        graph (nx.Graph): _description_
    """
    logger.info("Best path:")
    print_path(path)
    if not len(path):
        return
    
    end_nodes = set()
    start_nodes = set()
    for node, val in dict(graph.nodes.data('start')).items():
        if val == 'end':
            end_nodes.add(node)
        elif val == 'start':
            start_nodes.add(node)
    if len(end_nodes) > 0:
        end_nodes.add('end')
    
    if len(start_nodes) > 0 and path[0][1] not in start_nodes:
        logger.info('Did not start at start')
    
    time_offset = 0
    i = 0
    while i < len(path):
        if i + time_offset == path[i][0]:
            i += 1
            continue
        if path[i][0] < i + time_offset:
            logger.info(f'Extra node at time {path[i][0]}')
            time_offset -= 1
            i += 1
            continue
        if path[i][0] > i + time_offset:
            logger.info(f'Skipped time {path[i][0] - 1}')
            time_offset += 1
            i += 1
            continue
    
    node_dict = {node: 0 for node in graph.nodes}
    node_dict['end'] = 0
    
    for x in range(len(path) - 1):
        v1 = path[x][1]
        node_dict[v1] += 1
        v2 = path[x + 1][1]            
        if v1 == 'end' and not v2 == 'end':
            logger.info(f'Left the end node at path entry {x}')
        elif (not v1 == 'end') and (not v2 == 'end') and ((v1, v2) not in graph.edges):
            logger.info(f'Broke graph edge at path entry {x}')
        elif len(end_nodes) > 0 and (v2 == 'end') and (v1 not in end_nodes):
            logger.info(f'Went to end node illegally at path entry {x}')
    if len(path) > 1:
        node_dict[v2] += 1
    if len(path) == 1:
        node_dict[path[0][1]] += 1
    
    nodes = list(graph.nodes)
    for i in range(int(len(nodes) / 2)):
        visits = node_dict[nodes[2 * i]] + node_dict[nodes[2 * i + 1]]
        missing_visits = graph.nodes[nodes[2 * i]]["weight"] - visits
        if  missing_visits != 0:
            logger.info(f'Did not meet node weight for node: {get_original_vertex_name(nodes[2 * i])}. Missing visits: {missing_visits}')
            
            
            
def validate_edge2node_path(path: list, graph: nx.Graph):
    """Checks the constraints for a path on a graph.
    
    In particular:
     - does the path go along graph edges at each time step
     - is each node visited the correct number of times
     - is exactly one node visited per time step

    Args:
        path (list): _description_
        graph (nx.Graph): _description_
    """
    logger.info("Best path:")
    print_path(path)
    if not len(path):
        return
    
    end_nodes = set()
    start_nodes = set()
    for node, val in dict(graph.nodes.data('start')).items():
        if val == 'end':
            end_nodes.add(node)
        elif val == 'start':
            start_nodes.add(node)
    if len(end_nodes) > 0:
        end_nodes.add('end')
    
    if len(start_nodes) > 0 and path[0][1] not in start_nodes:
        logger.info('Did not start at start')
    
    time_offset = 0
    i = 0
    while i < len(path):
        if i + time_offset == path[i][0]:
            i += 1
            continue
        if path[i][0] < i + time_offset:
            logger.info(f'Extra node at time {path[i][0]}')
            time_offset -= 1
            i += 1
            continue
        if path[i][0] > i + time_offset:
            logger.info(f'Skipped time {path[i][0] - 1}')
            time_offset += 1
            i += 1
            continue
    
    node_dict = {node: 0 for node in graph.nodes}
    node_dict['end'] = 0
    
    for x in range(len(path) - 1):
        v1 = path[x][1]
        node_dict[v1] += 1
        v2 = path[x + 1][1]            
        if v1 == 'end' and not v2 == 'end':
            logger.info(f'Left the end node at path entry {x}')
        elif (not v1 == 'end') and (not v2 == 'end') and ((v1, v2) not in graph.edges):
            logger.info(f'Broke graph edge at path entry {x}')
        elif len(end_nodes) > 0 and (v2 == 'end') and (v1 not in end_nodes):
            logger.info(f'Went to end node illegally at path entry {x}')
    if len(path) > 1:
        node_dict[v2] += 1
    
    nodes = list(graph.nodes)
    for i in range(int(len(nodes))):
        visits = node_dict[nodes[i]]
        missing_visits = graph.nodes[nodes[i]]["weight"] - visits
        if  missing_visits != 0:
            logger.info(f'Did not meet node weight for node: {nodes[i]}. Missing visits: {missing_visits}')

