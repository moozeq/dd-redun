# Description
Finding similarities and redundancy in chemical data sets.

# Installation

```bash
git clone https://github.com/moozeq/DD_Redun.git

cd DD_Redun
git submodule update --init --recursive
pip3 install -r requirements.txt
```

## Requirements

- [openbabel](http://openbabel.org/wiki/Main_Page)
- [strip-it](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/strip-it/1.0.2/strip-it.html)

# Usage

## Help
```bash
./redun.py -h
```

## In pipeline
```bash
cat db.smi | ./redun.py
```

## Database file
```bash
./redun.py db.smi
```

## Prepare database
1. Download PDBBind database (e.g. CASF-2016)
2. In **CASF-2016/coreset** run script from below (simply getting smiles and id for each ligand):
    ```bash
    for f in *; do obabel -imol2 ${f}/${f}_ligand.mol2 -osmi | awk '{print $1" "$2}' >> results.smi; done
    ```
3. Database should be at **CASF-2016/coreset/results.smi**, copy it to DD_Redun folder