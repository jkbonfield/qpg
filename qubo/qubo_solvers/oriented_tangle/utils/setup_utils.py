import pickle
import os
import argparse
from qubo_solvers.definitions import DATA_DIR, OUT_DIR, Solver, COVERAGE_SUFFIX, QuboDescription
from qubo_solvers.oriented_tangle.utils.graph_utils import oriented_graph_with_copy_numbers

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath', default=f'{DATA_DIR}/test.gfa')
parser.add_argument('-t', '--times', help='delimited list input', 
    type=lambda s: [int(item) for item in s.split(',')])
parser.add_argument('-j', '--jobs', type=int)
parser.add_argument('-s', '--solver', required=True)
parser.add_argument('-d', '--data-dir', default=f'{OUT_DIR}/oriented')


def setup() -> QuboDescription:
    args = parser.parse_args()

    if args.solver in set(item.value for item in Solver):
        solver = Solver(args.solver)
    else:
        raise Exception(f'Solver {args.solver} not implemented yet.')

    
    filename = os.path.basename(args.filepath)
        
    with open(f'{args.data_dir}/{filename}.{COVERAGE_SUFFIX}', 'r') as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise Exception(f'Could not read copy numbers from {args.data_dir}/{filename}.{COVERAGE_SUFFIX}')
    copy_numbers = [int(x) for x in lines[2].split()]
        
    graph = oriented_graph_with_copy_numbers(args.filepath, copy_numbers)
    
    with open(f'{args.data_dir}/qubo_data_{filename}.pkl', 'rb') as f:
        data = pickle.load(f)
    
    
    return QuboDescription(filename=filename, data_dir=args.data_dir, graph=graph, time_limits=args.times, jobs=args.jobs,
                        Q=data['qubo_matrix'], offset=data['offset'], T=data['T_max'], V=data['V'], solver=solver)
