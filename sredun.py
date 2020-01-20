#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

index = 1
count = 1


class TermColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class Receptor:
    # global index, used to connecting data between plots and stdout
    index = 0

    def __init__(self, pdb: str, directory: str = None):
        split_pdb = pdb.splitlines()
        if len(split_pdb) < 1 or len(split_pdb[0].split()) < 2:
            # wrong smiles format
            self.name = None
            self.index = None
            self.info = None
            self.short_info = None
            self.pdb = None
            self.directory = None
            self.filename = None
            self.cf_filename = None
            return

        # receptor name/id
        self.name = split_pdb[0].split()[1]

        # receptor pdb file
        self.pdb = pdb

        # name to lowercase
        self.name = self.name.lower()

        # get rid of _POCKET suffix
        if self.name.endswith('_pocket'):
            self.name = self.name[:-len('_pocket')]

        # receptor short info as string
        self.short_info = f'[{self.index}]: {self.name}'

        # receptor full info as string
        self.info = f'[{self.index}]:\t{self.name}'

        # receptor global index
        self.index = f'{Receptor.index}'

        # directory
        self.directory = directory

        # receptor filename
        self.filename = f'{self.name}_pocket.pdb' if not directory else f'{directory}/{self.name}_pocket.pdb'

        # receptor chemical features filename
        self.cf_filename = f'{self.name}_pocket-cf.pdb' if not directory else f'{directory}/{self.name}_pocket-cf.pdb'

        Receptor.index += 1


def prepare_files(receptor: Receptor) -> bool:
    # save files for computation
    print(f'[*] [{receptor.name}] Checking files...', end='')

    # pdb file
    pdb_file_path = receptor.filename
    if not Path(pdb_file_path).exists():
        with open(pdb_file_path, 'w') as rec_file:
            rec_file.write(receptor.pdb)

    if receptor.pdb and Path(pdb_file_path).exists():
        print(f'\t{TermColors.OKGREEN}{pdb_file_path}{TermColors.ENDC}', end='')
    else:
        print(f'\t{TermColors.FAIL}{pdb_file_path}{TermColors.ENDC}', end='')
        return False

    # cf pdb file
    ret = 0
    cf_pdb_file_path = receptor.cf_filename
    if not Path(cf_pdb_file_path).exists():
        # here we provide pdb filename because it's input
        ret = subprocess.call(['java', 'AssignChemicalFeatures', pdb_file_path])

    if ret or not Path(cf_pdb_file_path).exists():
        print(f'\t{TermColors.FAIL}{cf_pdb_file_path}{TermColors.ENDC}')
        return False
    else:
        print(f'\t{TermColors.OKGREEN}{cf_pdb_file_path}{TermColors.ENDC}')
        return True


def get_ga_score(glosa_output: str) -> float:
    lines = glosa_output.splitlines()
    for line in lines:
        if line.startswith('GA-score'):
            score = line.split(':')[1]
            score = float(score.strip())
            return score
    return -1.0


def receptors_similarity(receptor1: Receptor, receptor2: Receptor) -> float:
    global index, count
    output_filename = f'{receptor1.name}_{receptor2.name}.out' if not receptor1.directory else f'{receptor1.directory}/{receptor1.name}_{receptor2.name}.out'
    output_filename_opt = f'{receptor2.name}_{receptor1.name}.out' if not receptor1.directory else f'{receptor1.directory}/{receptor2.name}_{receptor1.name}.out'
    cmd = ['glosa', '-s1', receptor1.filename, '-s1cf', receptor1.cf_filename, '-s2', receptor2.filename, '-s2cf',
           receptor2.cf_filename]

    print(f'[*] [{receptor1.name}] [{index % count:3} / {count:3}] Comparing to: {receptor2.name}\tscore:\t', end='')
    index += 1

    if not Path(output_filename).exists() and not Path(output_filename_opt).exists():
        glosa_output = subprocess.check_output(cmd).decode(sys.stdout.encoding).strip()
        if not glosa_output:
            print(f'{TermColors.FAIL}ERROR{TermColors.ENDC}\tCMD = {" ".join(cmd)}')
            return 0.0
        else:
            with open(output_filename, 'w') as out_file:
                out_file.write(glosa_output)
            with open(output_filename_opt, 'w') as out_file:
                out_file.write(glosa_output)
    else:
        if Path(output_filename).exists():
            with open(output_filename, 'r') as out_file:
                glosa_output = out_file.read()
        else:
            with open(output_filename_opt, 'r') as out_file:
                glosa_output = out_file.read()

    score = get_ga_score(glosa_output)
    if score >= 0.0:
        if score > 0.8:
            print(f'{TermColors.WARNING}{score}{TermColors.ENDC}')
        else:
            print(f'{TermColors.OKGREEN}{score}{TermColors.ENDC}')
    else:
        print(f'{TermColors.FAIL}ERROR{TermColors.ENDC}')
        Path(output_filename).unlink()

    return score if score >= 0.0 else 0.0


def receptor_compare(receptor: Receptor, receptors: list) -> list:
    return [receptors_similarity(receptor, sec_receptor) for sec_receptor in receptors]


# application modes, read help
app_modes = [
    'all', 'sim', 'dist', 'map', 'pro'
]

# application ascii art, typical
app_info = '''
    ███████╗██████╗ ███████╗██████╗ ██╗   ██╗███╗   ██╗
    ██╔════╝██╔══██╗██╔════╝██╔══██╗██║   ██║████╗  ██║
    ███████╗██████╔╝█████╗  ██║  ██║██║   ██║██╔██╗ ██║
    ╚════██║██╔══██╗██╔══╝  ██║  ██║██║   ██║██║╚██╗██║
    ███████║██║  ██║███████╗██████╔╝╚██████╔╝██║ ╚████║
    ╚══════╝╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚═╝  ╚═══╝

Helping finding similarities in protein database
'''


def main():
    parser = argparse.ArgumentParser(description=app_info, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', nargs='?', default='-', type=str,
                        help='receptors database filename (database as merged all .pdb files)')
    parser.add_argument('-m', '--mode', type=str, default='all', choices=app_modes,
                        help=
                        '''
                            application mode:
                                sim = similarities graph,
                                dist = similarity distance graph,
                                map = show index: protein mapping,
                                pro = show protein similar receptors
                        ''')
    parser.add_argument('-p', '--protein', type=int, default=-1,
                        help='select protein to show its similarities')
    parser.add_argument('-c', '--concurrency', action='store_true',
                        help='parallel computation')
    parser.add_argument('-t', '--threshold', type=float, default=0.0,
                        help='similarity threshold inclusive (0.0 - 1.0)')
    parser.add_argument('-d', '--dir', type=str, default='analysis',
                        help='working directory, all database computation will be done there')
    parser.add_argument('-o', '--output', type=str,
                        help='output filename, graphs will be saved with suffix _sim.png/_dist.png')
    args = parser.parse_args()

    # read from stdin if no filename provided
    if args.database == '-':
        with open('.temp', 'w') as temp:
            temp.writelines([line for line in sys.stdin])
        args.database = '.temp'

    # parse all receptors in merged pdb database
    with open(args.database) as db:
        pdbs = db.read().split('END')
        pdbs = [f'{pdb.strip()}\nTER' for pdb in pdbs if pdb.strip()]
        receptors = [Receptor(pdb, args.dir) for pdb in pdbs]
        if not all([prepare_files(receptor) for receptor in receptors]):
            print(f'[-] Preparing files failed')
            sys.exit(1)

    # create directory
    Path(args.dir).mkdir(parents=True, exist_ok=True)

    global index, count
    index = 1
    count = len(receptors)

    if args.concurrency:
        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count())

        # using scoring function, compare all vs all
        raw_similarities = pool.starmap(receptor_compare, [(rec, receptors) for rec in receptors])

        pool.close()
    else:
        # using scoring function, compare all vs all
        raw_similarities = [receptor_compare(rec, receptors) for rec in receptors]

    # filter scores with proper threshold
    similarities = [list(map(lambda x: round(x, 4) if x >= args.threshold else 0.0, sim)) for sim in raw_similarities]

    # remove output file if exists
    if args.output and os.path.exists(args.output):
        print(f'[*] Overwriting output file = {args.output}')
        os.remove(args.output)

    # receptors mapping
    if args.mode in ['all', 'map', 'sim', 'dist']:
        section_header = '========= RECEPTORS MAPPING ========='
        # get mapping as string and print to stdout
        mapping = '\n'.join([receptor.info for receptor in receptors])
        print(section_header)
        print(mapping)
        # if output filename was provided, write to file
        if args.output:
            with open(args.output, 'a') as out:
                out.write(f'{section_header}\n')
                out.write(mapping)

    # selected receptor similar proteins
    if args.mode in ['all', 'pro', 'sim', 'dist']:
        section_header = '========= RECEPTOR SIMILARITIES ========='
        if 0 <= args.protein < len(receptors):
            # get raw similarities for receptors
            receptor_similarities = raw_similarities[args.protein]
            # convert to dict where index: similarity
            receptor_sim_dict = {k: v for k, v in enumerate(receptor_similarities)}
            # sort dict by similarities where highest scores are at the bottom
            receptor_sim_dict = {k: v for k, v in
                                 sorted(receptor_sim_dict.items(), key=lambda item: item[1], reverse=True)}
            # round to 4 places
            receptor_sim_dict = {k: round(v, 4) for k, v in receptor_sim_dict.items()}
            # get similarities as string
            receptor_sims = '\n'.join([f'{sim:.4f}\t{receptors[key].info}' for key, sim in receptor_sim_dict.items()])
            header = f'''Selected receptor:\n\t{receptors[args.protein].info}\n\nscore\tindex\tid\n'''
            print(section_header)
            print(header)
            print(receptor_sims)
            if args.output:
                with open(args.output, 'a') as out:
                    out.write(f'{section_header}\n')
                    out.write(header)
                    out.write(receptor_sims)
        elif args.protein > 0:
            print(f'[-] Wrong receptor selected, available indexes: (0 - {len(receptors) - 1})')

    # similarities graph
    if args.mode in ['all', 'sim']:
        plt.figure()
        plt.grid()
        plt.title(f'Receptors similarity map')
        # show similarities between receptos as hot map
        plt.imshow(similarities, cmap='hot')
        plt.colorbar()
        # if output filename was provided, save .png to file
        if args.output:
            plt.savefig(f'{args.output}_sim.png')

    # dist graph
    if args.mode in ['all', 'dist']:
        plt.figure()
        plt.title(f'Receptors similarity graph where distance is similarity')
        # create distance array from list of lists
        dist_mat = np.array(similarities)
        # create matrix
        G = nx.from_numpy_matrix(dist_mat)
        # calculate what is the step to color nodes based on their order
        step = 1.0 / len(receptors)
        # calculate distances between nodes based on distance matrix based on similarities
        pos = nx.spring_layout(G)
        # draw labels and nodes without edges
        nx.draw_networkx_nodes(G, pos, node_size=100, cmap='hot', alpha=0.8,
                               node_color=[step * i for i in range(len(receptors))])
        nx.draw_networkx_labels(G, pos, labels={receptors.index(receptor): receptor.index for receptor in receptors})
        # if output filename was provided, save .png to file
        if args.output:
            plt.savefig(f'{args.output}_dist.png')

    # if mode required plots, show them now
    if args.mode in ['all', 'sim', 'dist']:
        plt.show()


if __name__ == '__main__':
    main()
