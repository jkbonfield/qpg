import numpy as np
import pickle
import os
import argparse
from pathlib import Path
from qubo_solvers.oriented_tangle.utils.graph_utils import oriented_graph_with_copy_numbers
from qubo_solvers.oriented_tangle.utils.qubo_utils import qubo_matrix_from_graph
from qubo_solvers.logging import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filepath')
    parser.add_argument('-c', '--copy-numbers', help='delimited list input', 
        type=lambda s: [float(item) for item in s.split(',') if len(item)])
    parser.add_argument('-p', '--penalties', help='delimited list input', 
        type=lambda s: [int(item) for item in s.split(',') if len(item)])
    parser.add_argument('-d', '--data-dir')

    args = parser.parse_args()

    logger.info(f'Building qubo matrix for {args.filepath}')
    filename = os.path.basename(args.filepath)
    Path(args.data_dir).mkdir(exist_ok=True, parents=True)

    if args.copy_numbers is None:
        logger.info('No copy numbers provided')
    else:
        logger.info(f'Copy numbers provided: {args.copy_numbers}')
    
    logger.info(f'Getting graph from {args.filepath}')
    graph = oriented_graph_with_copy_numbers(args.filepath, args.copy_numbers)
    Q, offset, T_max, V = qubo_matrix_from_graph(graph, penalties=args.penalties)
    
    logger.info('Saving data')
    to_save = {
        'Q': Q.tolist(),
        'offset': offset,
        'T_max': T_max,
        'V': V,
        'graph': graph
    }
    with open(f'{args.data_dir}/qubo_data_{filename}.pkl', 'wb') as file:
        pickle.dump(to_save, file)


    logger.info('Writing to MQLib format')
    mqlib_filepath = f'{args.data_dir}/mqlib_input_{filename}.txt'

    ut_qubo_matrix = np.triu(Q)
    non_zero = np.nonzero(ut_qubo_matrix)
    non_zero_count = int(non_zero[0].shape[0])

    with open(mqlib_filepath, 'w') as f:
        f.write(f'{Q.shape[0]} {non_zero_count}\n')
        to_write = ''
        for i in range(len(non_zero[0])):
            to_write += f'{non_zero[0][i] + 1} {non_zero[1][i] + 1} {-Q[non_zero[0][i], non_zero[1][i]]} \n'
            if i % 500 == 0:
                f.write(to_write)
                to_write = ''

        f.write(to_write)
    logger.info('Finished building oriented qubo matrix')


if __name__ == "__main__":
    exit(main())