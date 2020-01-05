#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

from sklearn import preprocessing

from redun import Ligand

# available scoring methods
scores = [
    'rfscore_v1', 'rfscore_v2', 'rfscore_v3', 'nnscore', 'vina_affinity'
]

# application ascii art, typical
app_info = '''
    ███████╗ ██████╗ ██████╗ ██████╗ ██╗   ██╗███╗   ██╗
    ██╔════╝██╔════╝██╔═══██╗██╔══██╗██║   ██║████╗  ██║
    ███████╗██║     ██║   ██║██████╔╝██║   ██║██╔██╗ ██║
    ╚════██║██║     ██║   ██║██╔══██╗██║   ██║██║╚██╗██║
    ███████║╚██████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║
    ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝

Analyse similarities between complexes ligand-protein by docking them against each other's receptors
'''


def main():
    parser = argparse.ArgumentParser(
        description=app_info,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', nargs='?', default='-', type=str,
                        help='smiles database filename')
    parser.add_argument('ligands', nargs='+', type=int, help='ligands indexes')
    parser.add_argument('-o', '--output', type=str, help='output filename')
    args = parser.parse_args()

    # read from stdin if no filename provided
    if args.database == '-':
        # write to .temp file because Chem needs filename
        with open('.temp', 'w') as temp:
            temp.writelines([line for line in sys.stdin])
        args.database = '.temp'

    # parse all ligands in smiles database
    with open(args.database) as db:
        ligands = [Ligand(line) for line in db]

    # holds ligands order, helpful when displaying results
    results_order = []
    # every ligand scores
    results = {}
    for ligand_index in args.ligands:
        for ligand_sec_index in args.ligands:
            # ligand key displayed in results
            result_key = f'{ligand_index}->{ligand_sec_index}'
            # preserve order
            results_order.append(result_key)
            # result file should be named like this when docked with script 'dock'
            result_file = f'docking_ligand{ligand_index}_protein{ligand_sec_index}_rescored.mol2'
            ret = 0
            # if file doesn't exists, run 'dock' script
            if not os.path.exists(result_file):
                ret = subprocess.call(['./dock', args.database, str(ligand_index), str(ligand_sec_index)])
            # if returned something other than 0, error occured
            if ret:
                print(f'''[-] Docking {ligands[ligand_index].short_info} to protein from {ligands[
                    ligand_sec_index].short_info} ended with error''')
                sys.exit(1)
            # open result file
            with open(result_file, 'r') as result_file:
                current_score = {}
                for line in result_file:
                    # break loop after one molecule, get only the best one
                    if '@<TRIPOS>MOLECULE' in line:
                        break
                    # add scores
                    for score in scores:
                        if score in line:
                            line_elements = line.split()
                            # check if score is number
                            if len(line_elements) < 3 or not line_elements[2].replace('.', '', 1).lstrip(
                                    '-+').isdigit():
                                print(f'[-] Wrong file format, error in line: {line}')
                                sys.exit(1)
                            # save score under proper key
                            current_score[score] = float(line_elements[2])
                # add all scores under ligand key
                results[result_key] = current_score
    # convert dict to list
    results_list = [[ligand_score[score] for score in scores] for ligand_score in results.values()]
    # normalize scores results
    results_list = preprocessing.normalize(results_list)
    # absolute value because of vina affinity
    results_list = [list(map(lambda x: abs(x), ligand_scores)) for ligand_scores in results_list]
    # sum all scores
    results_list = [sum(ligand_scores) for ligand_scores in results_list]
    # use results_order to get ligands ids back
    results_sum = {results_order[index]: ligand_score for index, ligand_score in enumerate(results_list)}

    # remove output file if exists
    if args.output and os.path.exists(args.output):
        print(f'[*] Overwriting output file = {args.output}')
        os.remove(args.output)

    # convert dict to string, score with 8 places after dot
    results_to_print = '\n'.join([f'{key}\t{score:.8f}' for key, score in results_sum.items()])
    print('Results format:\n\tligand->receptor\tscore\n')
    print(results_to_print)
    # save to file if specified
    if args.output:
        with open(args.output, 'w') as out:
            out.write(results_to_print)


if __name__ == '__main__':
    main()
