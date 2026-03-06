#!/bin/bash
#SBATCH --job-name=fix_dgl_sp2d
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --partition=g24
#SBATCH --output=/home/yangl_pacagen_com/fix_dgl_sp2d_%j.out
#SBATCH --error=/home/yangl_pacagen_com/fix_dgl_sp2d_%j.err

set -e

echo "=== Fix DGL + SP2D at $(date) ==="

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate drug_lamp_test

echo "=== Upgrading DGL for torch 2.1 + cu121 ==="
pip install "dgl>=2.1.0" -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html

echo "=== Patching torch_geometric x_map for SP2D ==="
SMILES_PY=$(python -c "from torch_geometric.utils import smiles; print(smiles.__file__)")
echo "Patching: $SMILES_PY"

python -c "
from torch_geometric.utils.smiles import x_map
print('Before:', x_map['hybridization'])
"

# Insert SP2D before OTHER in the x_map
python -c "
import ast, re

smiles_file = '${SMILES_PY}'
with open(smiles_file, 'r') as f:
    content = f.read()

# Replace the hybridization list to include SP2D
old = \"'SP3D2', 'OTHER'\"
new = \"'SP3D2', 'SP2D', 'OTHER'\"
if old in content and 'SP2D' not in content:
    content = content.replace(old, new)
    with open(smiles_file, 'w') as f:
        f.write(content)
    print('Patched successfully')
else:
    print('Already patched or pattern not found')
"

echo "=== Verifying ==="
python -c "
import torch
print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
print('cuda available:', torch.cuda.is_available())

import torch_geometric
print('torch_geometric:', torch_geometric.__version__)

from torch_geometric.utils.smiles import x_map
print('hybridization types:', x_map['hybridization'])
print('SP2D supported:', 'SP2D' in x_map['hybridization'])

import dgl
print('dgl:', dgl.__version__)

import dgllife
print('dgllife:', dgllife.__version__)

print('All checks passed!')
"

echo "=== Done at $(date) ==="
