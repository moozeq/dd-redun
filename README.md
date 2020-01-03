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