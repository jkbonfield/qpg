import gfapy
import networkx as nx
from qubo_solvers.logging import get_logger

logger = get_logger(__name__)


def edge2node_oriented_graph(filename: str, copy_numbers: list[str]):
    """Reads a .gfa file into an oriented graph, where each node has a positive and negative version.

    Args:
        filename (str): filepath to read.
        copy_numbers (list[str]): list of copy numbers for positive and negative nodes

    Returns:
        nx.Graph: corresponding oriented graph.
    """
    gfa = gfapy.Gfa.from_file(filename, vlevel=0)
    graph = nx.DiGraph()
    for index, segment_line in enumerate(gfa.segments):
        graph.add_node(f'{segment_line.name}_+', weight=copy_numbers[2*index], start=segment_line.st)
        graph.add_node(f'{segment_line.name}_-', weight=copy_numbers[2*index+1], start=segment_line.st)
    for edge_line in gfa.edges:
        v1 = edge_line.sid1
        v2 = edge_line.sid2
        graph.add_edges_from([
            (f'{v1.name}_{v1.orient}', f'{v2.name}_{v2.orient}'),
        ])
        v1.invert()
        v2.invert()
        graph.add_edges_from([
            (f'{v2.name}_{v2.orient}', f'{v1.name}_{v1.orient}'),
        ])
    return graph


def oriented_graph_with_copy_numbers(filename, copy_numbers: list[float] | None ):
    """Reads a .gfa file into an oriented graph, where each node has a positive and negative version.

    Args:
        filename (str): filepath to read.

    Returns:
        nx.Graph: corresponding oriented graph.
    """
    gfa = gfapy.Gfa.from_file(filename, vlevel=0)
    
    if copy_numbers is None:
        copy_numbers = [segment_line.SC for segment_line in gfa.segments]
    
    graph = nx.DiGraph()
    for index, segment_line in enumerate(gfa.segments):
        graph.add_node(f'{segment_line.name}_+', weight=copy_numbers[index], start=segment_line.st)
        graph.add_node(f'{segment_line.name}_-', weight=copy_numbers[index], start=segment_line.st)
    for edge_line in gfa.edges:
        v1 = edge_line.sid1
        v2 = edge_line.sid2
        graph.add_edges_from([
            (f'{v1.name}_{v1.orient}', f'{v2.name}_{v2.orient}'),
        ])
        v1.invert()
        v2.invert()
        graph.add_edges_from([
            (f'{v2.name}_{v2.orient}', f'{v1.name}_{v1.orient}'),
        ])
    return graph
