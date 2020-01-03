#!/usr/bin/env python3
import argparse
import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def main():
    parser = argparse.ArgumentParser(description='Finding similarities in smiles database')
    parser.add_argument('database', nargs='?', default='-', type=str, help='smiles database filename')
    parser.add_argument('-s', '--sim', type=str, default='dice', choices=['dice', 'tanimoto'],
                        help='fingerprint similarity scoring function')
    parser.add_argument('-p', '--precision', type=int, default=2,
                        help='similarity scoring precision')
    args = parser.parse_args()
    if args.database == '-':
        with open('.temp', 'w') as temp:
            temp.writelines([line for line in sys.stdin])
        args.database = '.temp'

    smiles = Chem.SmilesMolSupplier(args.database, titleLine=False)
    fps = [AllChem.GetMorganFingerprint(smile, args.precision) for smile in smiles]

    sim_scoring = DataStructs.FingerprintSimilarity if args.sim == 'tanimoto' else DataStructs.DiceSimilarity
    similarities = [[round(sim_scoring(fp, cfp), args.precision) for cfp in fps] for fp in fps]
    distances = [[round(1.0 - sim, args.precision) for sim in sim_list] for sim_list in similarities]

    dist_mat = np.array(distances)
    g = nx.from_numpy_matrix(dist_mat)
    nx.draw(g)
    plt.show()


if __name__ == "__main__":
    main()
