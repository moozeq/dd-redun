# Description
Finding similarities and redundancy in chemical data sets.

# Installation

```bash
git clone https://github.com/moozeq/DD_Redun.git

cd DD_Redun
git submodule update --init --recursive
pip3 install -r requirements.txt
```

## Prepare database

### Ligands

1. Download PDBBind database (e.g. [CASF-2016](http://www.pdbbind.org.cn/casf.asp))
2. You may need to remove **4mme** complex, because ligand from this complex is causing error when creating fingerprints
3. In *CASF-2016/coreset* run script from below (simply getting smiles and id for each ligand):
    ```bash
    for f in *; do obabel -imol2 ${f}/${f}_ligand.mol2 -osmi | awk '{print $1" "$2}' >> db.smi; done
    ```
4. Database should be at *CASF-2016/coreset/db.smi*, copy it into DD_Redun folder
5. For docking functionality you need to also copy *coreset* folder into DD_Redun folder

### Receptors

6. In *CASF-2016/coreset* run script from below (simply merging all x_pocket.pdb files into one file database):
    ```bash
    for f in *; do cat ${f}/${f}_pocket.pdb >> prots.pdb; done
    ```
7. Database should be at *CASF-2016/coreset/prots.pdb*, copy it into DD_Redun folder
8. G-LoSA should be build used clang or g++:
    ```bash
    g++ glosa.cpp -o glosa
    ```

## Requirements
### Main functionality
- [openbabel](http://openbabel.org/wiki/Main_Page)

### Scaffolds
- [strip-it](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/strip-it/1.0.2/strip-it.html)

### Docking
- [AutoDock Vina](http://vina.scripps.edu/)
- [ODDT](https://pythonhosted.org/oddt/)
- PDBBind coreset

### Receptors
- [G-LoSA](https://compbio.lehigh.edu/GLoSA/index.html)

# Usage

## Help
```bash
./redun.py -h
./scorun.py -h
./sredun.py -h
```

## Database file
```bash
./redun.py db.smi
./scorun.py db.smi [ints]
./sredun.py prots.pdb
```

## In pipeline
```bash
cat db.smi | ./redun.py
cat db.smi | ./scorun.py [ints]
cat prots.pdb | ./sredun.py
```