import pickle
from collections import Counter
from datetime import datetime
import os
import numpy as np
import argparse
import re
from qubo_solvers.oriented_tangle.utils.sampling_utils import dwave_sample_qubo, mqlib_sample_qubo, gurobi_sample_qubo, validate_path, validate_edge2node_path
from qubo_solvers.definitions import Solver, QuboDescription
from qubo_solvers.logging import get_logger


logger = get_logger(__name__)


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath')
parser.add_argument('-t', '--times', help='delimited list input', 
    type=lambda s: [int(item) for item in s.split(',')])
parser.add_argument('-j', '--jobs', type=int)
parser.add_argument('-s', '--solver', required=True)
parser.add_argument('-d', '--data-dir')
parser.add_argument('-o', '--output')
parser.add_argument('--edge2node', action='store_true')


def setup(args) -> QuboDescription:
    if args.solver in set(item.value for item in Solver):
        solver = Solver(args.solver)
    else:
        raise Exception(f'Solver {args.solver} not implemented yet.')

    filename = os.path.basename(args.filepath)
            
    try:
        with open(f'{args.data_dir}/qubo_data_{filename}.pkl', 'rb') as f:
            data = pickle.load(f)
    except FileNotFoundError:
        raise Exception('Run build_oriented_qubo_matrix first!')
    
    return QuboDescription(filename=filename, data_dir=args.data_dir, graph=data['graph'], time_limits=args.times, jobs=args.jobs,
                        Q=data['Q'], offset=data['offset'], T=data['T_max'], V=data['V'], solver=solver)


def main():
    args = parser.parse_args()
    qubo_description = setup(args)
    logger.info(f'Start sampling using: {qubo_description.solver.name}')
    if qubo_description.solver == Solver.DWAVE:
        paths = dwave_sample_qubo(qubo_description)
    elif qubo_description.solver == Solver.MQLIB:
        paths = mqlib_sample_qubo(qubo_description)
    elif qubo_description.solver == Solver.GUROBI:
        paths = gurobi_sample_qubo(qubo_description)
    else:
        raise Exception(f'Unrecognised solver: {qubo_description.solver}')

    for time_limit in qubo_description.time_limits:
        for i in range(qubo_description.jobs):
            if args.edge2node:
                logger.info('Validating edge2node')
                validate_edge2node_path(paths[time_limit][i][2], qubo_description.graph)
            else:
                validate_path(paths[time_limit][i][2], qubo_description.graph)
            logger.info(f'Energy of path: {paths[time_limit][i][1]}')

    now = datetime.now().strftime('%d%m%Y_%H%M')
    save_file = qubo_description.data_dir + f'/{qubo_description.solver.value}_{qubo_description.filename}_{now}'   
        
    with open(save_file, 'wb') as f:
        pickle.dump(paths, f)
        
    compile_path = qubo_description.data_dir + f'/{qubo_description.solver.value}.{qubo_description.filename}.compiled.txt'
    counts = {
        time_limit: Counter([float(paths[time_limit][i][1]) for i in range(len(paths[time_limit]))]) for time_limit in qubo_description.time_limits
    }
    with open(compile_path, 'a') as f:
        for time_limit in qubo_description.time_limits:
            f.write(f'{time_limit}: {counts[time_limit]},')

    if args.output is not None:
        for time_limit in qubo_description.time_limits:
            for i in range(qubo_description.jobs):         
                path_fragments = []
                fragment = [] 
                for path_step in paths[time_limit][i][2]:
                    node = path_step[1]
                    if node == 'end':
                        path_fragments.append(fragment)
                        fragment = []
                    else:
                        match: re.Match | None = re.search(r'([!-)+-<>-~][!-~]*)\_([+-])', node)
                        if match is None:
                            raise Exception(f'Could not parse node: {node}')
                        fragment.append(f'{">" if match.group(2) == "+" else "<"}{match.group(1)}\n')
                path_fragments.append(fragment)

                # TODO: handle multiple contigs?
                # TODO: need to check that a fragment doesn't have broken graph steps?
                longest_frag_index = np.argmax([len(fragment) for fragment in path_fragments])        
                with open(f'{args.output}.{time_limit}.{i}', 'w') as f:
                    for node in path_fragments[longest_frag_index]:
                        f.write(node)

    return 0


if __name__ == "__main__":
    exit(main())