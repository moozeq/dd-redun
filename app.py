#!/usr/bin/env python3
import argparse
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

from DD_Scaffolds import scaffolds


class Ligand:
    index = 0

    def __init__(self, smi: str):
        split_smi = smi.split()
        if len(split_smi) < 2:
            self.smile = None
            self.name = None
            return
        self.smile = split_smi[0]
        self.name = split_smi[1]
        if self.name.endswith('_ligand'):
            self.name = self.name[:-len('_ligand')]
        self.index = f'{Ligand.index}'
        Ligand.index += 1


comp_modes = [
    'RINGS_WITH_LINKERS_1', 'RINGS_WITH_LINKERS_2', 'MURCKO_1', 'MURCKO_2', 'OPREA_1', 'OPREA_2', 'OPREA_3',
    'SCHUFFENHAUER_1', 'SCHUFFENHAUER_2', 'SCHUFFENHAUER_3', 'SCHUFFENHAUER_4', 'SCHUFFENHAUER_5']

app_modes = [
    'all', 'sim', 'scaff', 'dist', 'map'
]


def main():
    parser = argparse.ArgumentParser(description='Finding similarities in smiles database')
    parser.add_argument('database', nargs='?', default='-', type=str, help='smiles database filename')
    parser.add_argument('-m', '--mode', type=str, default='all', choices=app_modes,
                        help=
                        '''
                            application mode:
                                sim = similarities graph,
                                scaff = scaffolds,
                                dist = similarity distance graph,
                                map = show index: ligand mapping
                        ''')
    parser.add_argument('-s', '--sim', type=str, default='dice', choices=['dice', 'tanimoto'],
                        help='fingerprint similarity scoring function')
    parser.add_argument('-c', '--comp', type=str, help='inhibitors compering method', default='RINGS_WITH_LINKERS_1',
                        choices=comp_modes)
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help='similarity threshold')
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()
    if args.database == '-':
        with open('.temp', 'w') as temp:
            temp.writelines([line for line in sys.stdin])
        args.database = '.temp'

    with open(args.database) as db:
        ligands = [Ligand(line) for line in db]
    ligands_map = {ligand.index: ligand.name for ligand in ligands}

    smiles = Chem.SmilesMolSupplier(args.database, titleLine=False)
    fps = [AllChem.GetMorganFingerprint(smile, 2) for smile in smiles]

    sim_scoring = DataStructs.FingerprintSimilarity if args.sim == 'tanimoto' else DataStructs.DiceSimilarity
    similarities = [[round(sim_scoring(fp, cfp), 4) for cfp in fps] for fp in fps]
    similarities = [list(map(lambda x: x if x > args.threshold else 0.0, sim)) for idx, sim in enumerate(similarities)]

    # ligands mapping
    if args.mode in ['all', 'map']:
        for index in ligands_map:
            print(f'{index}: {ligands_map[index]} {ligands[int(index)].smile}')

    # similarities graph
    if args.mode in ['all', 'sim']:
        plt.title(f'Ligands {args.sim} similarity map')
        plt.imshow(similarities, cmap='hot')
        plt.colorbar()
        plt.show()

    # scaffolding
    if args.mode in ['all', 'scaff']:
        scaffolds.strip(args.database)
        scaffs = scaffolds.merge(args.database, args.comp)
        scaffolds.show_results(scaffs)

    # dist graph
    if args.mode in ['all', 'dist']:
        dist_mat = np.array(similarities)
        G = nx.from_numpy_matrix(dist_mat)
        step = 1.0 / len(fps)
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos, node_size=100, cmap='hot', alpha=0.8,
                               node_color=[step * i for i in range(len(fps))])
        nx.draw_networkx_labels(G, pos, labels={ligands.index(ligand): ligand.index for ligand in ligands})
        plt.title(f'Ligands {args.sim} similarity graph where distance is similarity')
        plt.show()


if __name__ == "__main__":
    main()
