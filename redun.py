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
        self.info = f'[{self.index}]:\t{self.name}\t{self.smile}'

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
    'all', 'sim', 'scaffs', 'dist', 'map', 'lig'
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
                                map = show index: ligand mapping,
                                lig = show ligand similar molecules
                        ''')
    parser.add_argument('-s', '--sim', type=str, default='dice', choices=['dice', 'tanimoto'],
                        help='fingerprint similarity scoring function')
    parser.add_argument('-l', '--ligand', type=int, default=-1,
                        help='select ligand to show its similarities')
    parser.add_argument('-c', '--comp', type=str, default='rwl1', choices=comp_short_modes,
                        help=
                        '''
                            inhibitors compering method:
                                rwl = rings with linkers,
                                mur = murcko,
                                opr = oprea,
                                sch = schuffenhauer
                        ''')
    parser.add_argument('-t', '--threshold', type=float, default=0.0,
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
    raw_similarities = [[sim_scoring(fp, cfp) for cfp in fps] for fp in fps]
    # filter scores with proper threshold
    similarities = [list(map(lambda x: round(x, 4) if x >= args.threshold else 0.0, sim)) for idx, sim in
                    enumerate(raw_similarities)]

    # remove output file if exists
    if args.output and os.path.exists(args.output):
        print(f'[*] Overwriting output file = {args.output}')
        os.remove(args.output)

    # scaffolding
    if args.mode in ['all', 'scaffs']:
        section_header = '========= SCAFFOLDS ========='
        # strip smiles to scaffolds
        scaffolds.strip(args.database)
        # merge ligands with same scaffolds
        scaffs = scaffolds.merge(args.database, args.comp)
        # show results
        print(section_header)
        results = scaffolds.show_results(scaffs)
        # if output filename was provided, write to file
        if args.output:
            with open(args.output, 'a') as out:
                out.write(f'{section_header}\n')
                out.write(results)

    # ligands mapping
    if args.mode in ['all', 'map', 'sim', 'dist']:
        section_header = '========= LIGANDS MAPPING ========='
        # get mapping as string and print to stdout
        mapping = '\n'.join([ligand.info for ligand in ligands])
        print(section_header)
        print(mapping)
        # if output filename was provided, write to file
        if args.output:
            with open(args.output, 'a') as out:
                out.write(f'{section_header}\n')
                out.write(mapping)

    # selected ligand similar molecules
    if args.mode in ['all', 'lig', 'sim', 'dist']:
        section_header = '========= LIGAND SIMILARITIES ========='
        if 0 <= args.ligand < len(ligands):
            # get raw similarities for ligand
            ligand_similarities = raw_similarities[args.ligand]
            # convert to dict where index: similarity
            ligand_sim_dict = {k: v for k, v in enumerate(ligand_similarities)}
            # sort dict by similarities where highest scores are at the bottom
            ligand_sim_dict = {k: v for k, v in sorted(ligand_sim_dict.items(), key=lambda item: item[1], reverse=True)}
            # round to 4 places
            ligand_sim_dict = {k: round(v, 4) for k, v in ligand_sim_dict.items()}
            # get similarities as string
            ligand_sims = '\n'.join([f'{sim:.4f}\t{ligands[key].info}' for key, sim in ligand_sim_dict.items()])
            header = f'''Selected ligand:\n\t{ligands[
                args.ligand].info}\nSimilarity scoring mode:\n\t{args.sim}\n\nscore\tindex\tid\tsmiles\n'''
            print(section_header)
            print(header)
            print(ligand_sims)
            if args.output:
                with open(args.output, 'a') as out:
                    out.write(f'{section_header}\n')
                    out.write(header)
                    out.write(ligand_sims)
        elif args.ligand > 0:
            print(f'[-] Wrong ligand selected, available indexes: (0 - {len(ligands) - 1})')

    # similarities graph
    if args.mode in ['all', 'sim']:
        plt.figure()
        plt.grid()
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
