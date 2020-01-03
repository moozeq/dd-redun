# Description
Finding similarities and redundancy in chemical data sets.

# Installation

```bash
git clone https://github.com/moozeq/DD_Redun.git

cd DD_Redun
git submodule update --init --recursive
pip3 install -r requirements.txt
```

# Usage

## Help
```bash
./app.py -h
```

## In pipeline
```bash
cat db.smi | ./app.py
```

## Database file
```bash
./app.py db.smi
```

## Prepare database
Required:
- open babel

1. Download PDBBind database (e.g. CASF-2016)
2. In __CASF-2016/coreset__ run script from below (simply getting smiles and id for each ligand):
    ```bash
    for f in *; do obabel -imol2 ${f}/${f}_ligand.mol2 -osmi | awk '{print $1" "$2}' >> results.smi; done
    ```
3. Database should be at __CASF-2016/coreset/results.smi__