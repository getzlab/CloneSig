# CloneSig

CloneSig for identification of genes with convergent mutations across multiple cancer clones.

Implemented in Convergent evolution within and across patients reveals mechanisms of cancer drug resistance.

## Usage

```bash
python Run_CloneSig.py input.maf genes_list.txt [--n_perm N_PERMUTATIONS] [--step {permutations,genes,all}]
```

## Installation Requirements
1. Clone this repository:
```bash
git clone https://github.com/getzlab/CloneSig.git
cd CloneSig
```

2. Install CloneSig and its dependencies:
```bash
pip install -e .
```

3. Install CurveBall:
```bash
# Option 1: Install as an extra
pip install -e .[curveball]

# Option 2: Manual installation
git clone https://github.com/getzlab/CurveBall.git
cd CurveBall
pip install -e .
```


# Concepts

# Inputs

# Outputs
