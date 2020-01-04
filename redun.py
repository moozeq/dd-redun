#!/usr/bin/env python3
import argparse
import os
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

from DD_Scaffolds import scaffolds


# structure used to store all info about ligand
class Ligand:
    # global index, used to connecting data between plots and stdout
    index = 0

    def __init__(self, smi: str):
        split_smi = smi.split()
        if len(split_smi) < 2:
            # wrong smiles format
            self.smile = None
            self.name = None
            self.index = None
            self.info = None
            return

        # ligand smile
        self.smile = split_smi[0]

        # ligand name/id
        self.name = split_smi[1]

        # get rid of _ligand suffix
        if self.name.endswith('_ligand'):
            self.name = self.name[:-len('_ligand')]

        # ligand info as string
        self.info = f'[{self.index}]: {self.name} {self.smile}'

        # ligand global index
        self.index = f'{Ligand.index}'
        Ligand.index += 1


# compering methods used in scaffolding
comp_modes = [
    'RINGS_WITH_LINKERS_1', 'RINGS_WITH_LINKERS_2', 'MURCKO_1', 'MURCKO_2', 'OPREA_1', 'OPREA_2', 'OPREA_3',
    'SCHUFFENHAUER_1', 'SCHUFFENHAUER_2', 'SCHUFFENHAUER_3', 'SCHUFFENHAUER_4', 'SCHUFFENHAUER_5'
]

# shortcuts for compering methods
comp_short_modes = [
    'rwl1', 'rwl2', 'mur1', 'mur2', 'opr1', 'opr2', 'opr3', 'sch1', 'sch2', 'sch3', 'sch4', 'sch5'
]

# application modes, read help
app_modes = [
    'all', 'sim', 'scaffs', 'dist', 'map'
]

# application ascii art, typical
app_art = '''
    ██████╗ ███████╗██████╗ ██╗   ██╗███╗   ██╗
    ██╔══██╗██╔════╝██╔══██╗██║   ██║████╗  ██║
    ██████╔╝█████╗  ██║  ██║██║   ██║██╔██╗ ██║
    ██╔══██╗██╔══╝  ██║  ██║██║   ██║██║╚██╗██║
    ██║  ██║███████╗██████╔╝╚██████╔╝██║ ╚████║
    ╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚═╝  ╚═══╝

'''


def main():
    parser = argparse.ArgumentParser(description=app_art + 'Finding similarities in smiles database',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', nargs='?', default='-', type=str,
                        help='smiles database filename')
    parser.add_argument('-m', '--mode', type=str, default='all', choices=app_modes,
                        help=
                        '''
                            application mode:
                                sim = similarities graph,
                                scaffs = scaffolds,
                                dist = similarity distance graph,
                                map = show index: ligand mapping
                        ''')
    parser.add_argument('-s', '--sim', type=str, default='dice', choices=['dice', 'tanimoto'],
                        help='fingerprint similarity scoring function')
    parser.add_argument('-c', '--comp', type=str, default='rwl1', choices=comp_short_modes,
                        help=
                        '''
                            inhibitors compering method:
                                rwl = rings with linkers,
                                mur = murcko,
                                opr = oprea,
                                sch = schuffenhauer
                        ''')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='similarity threshold inclusive (0.0 - 1.0)')
    parser.add_argument('-o', '--output', type=str,
                        help='output filename, graphs will be saved with suffix _sim.png/_dist.png')
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

    # load ligands to Chem module
    smiles = Chem.SmilesMolSupplier(args.database, titleLine=False)

    # use proper similarity scoring function and fingerprinting
    if args.sim == 'dice':
        fps = [AllChem.GetMorganFingerprint(smile, 2) for smile in smiles]
        sim_scoring = DataStructs.DiceSimilarity
    elif args.sim == 'tanimoto':
        fps = [FingerprintMols.FingerprintMol(smile) for smile in smiles]
        sim_scoring = DataStructs.FingerprintSimilarity
    else:
        print(f'[-] Wrong similiarity scoring function')
        sys.exit(1)

    # using scoring function, compare all vs all
    similarities = [[round(sim_scoring(fp, cfp), 4) for cfp in fps] for fp in fps]
    # filter scores with proper threshold
    similarities = [list(map(lambda x: x if x >= args.threshold else 0.0, sim)) for idx, sim in enumerate(similarities)]

    # remove output file if exists
    if args.output and os.path.exists(args.output):
        print(f'[*] Overwriting output file = {args.output}')
        os.remove(args.output)

    # scaffolding
    if args.mode in ['all', 'scaffs']:
        # strip smiles to scaffolds
        scaffolds.strip(args.database)
        # merge ligands with same scaffolds
        scaffs = scaffolds.merge(args.database, args.comp)
        # show results
        results = scaffolds.show_results(scaffs)
        # if output filename was provided, write to file
        if args.output:
            with open(args.output, 'a') as out:
                out.write('========= SCAFFOLDS =========')
                out.write(results)

    # ligands mapping
    if args.mode in ['all', 'map', 'sim', 'dist']:
        # get mapping as string and print to stdout
        mapping = '\n'.join([ligand.info for ligand in ligands])
        print(mapping)
        # if output filename was provided, write to file
        if args.output:
            with open(args.output, 'a') as out:
                out.write('========= LIGANDS MAPPING =========')
                out.write(mapping)

    # similarities graph
    if args.mode in ['all', 'sim']:
        plt.figure()
        plt.title(f'Ligands {args.sim} similarity map')
        # show similarities between ligands as hot map
        plt.imshow(similarities, cmap='hot')
        plt.colorbar()
        # if output filename was provided, save .png to file
        if args.output:
            plt.savefig(f'{args.output}_sim.png')

    # dist graph
    if args.mode in ['all', 'dist']:
        plt.figure()
        plt.title(f'Ligands {args.sim} similarity graph where distance is similarity')
        # create distance array from list of lists
        dist_mat = np.array(similarities)
        # create matrix
        G = nx.from_numpy_matrix(dist_mat)
        # calculate what is the step to color nodes based on their order
        step = 1.0 / len(fps)
        # calculate distances between nodes based on distance matrix based on similarities
        pos = nx.spring_layout(G)
        # draw labels and nodes without edges
        nx.draw_networkx_nodes(G, pos, node_size=100, cmap='hot', alpha=0.8,
                               node_color=[step * i for i in range(len(fps))])
        nx.draw_networkx_labels(G, pos, labels={ligands.index(ligand): ligand.index for ligand in ligands})
        # if output filename was provided, save .png to file
        if args.output:
            plt.savefig(f'{args.output}_dist.png')

    # if mode required plots, show them now
    if args.mode in ['all', 'sim', 'dist']:
        plt.show()


if __name__ == "__main__":
    main()
