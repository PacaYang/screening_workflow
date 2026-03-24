#!/bin/bash
#SBATCH --job-name=reinstall_deps
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --partition=g24
#SBATCH --output=/home/yangl_pacagen_com/reinstall_deps_%j.out
#SBATCH --error=/home/yangl_pacagen_com/reinstall_deps_%j.err

set -e

echo "=== Starting dependency reinstall at $(date) ==="
echo "Host: $(hostname)"

source /home/yangl_pacagen_com/miniconda3/etc/profile.d/conda.sh
conda activate drug_lamp_test

echo "=== Before upgrade ==="
python -c "import torch; print('torch:', torch.__version__, 'cuda:', torch.version.cuda)"
pip show torch-geometric dgl dgllife 2>/dev/null | grep -E "^(Name|Version)"

echo "=== Upgrading PyTorch to 2.1.2+cu121 ==="
pip install torch==2.1.2 torchvision==0.16.2 torchaudio==2.1.2 --index-url https://download.pytorch.org/whl/cu121

echo "=== Upgrading torch_geometric ==="
pip install torch_geometric==2.5.3

echo "=== Upgrading DGL ==="
pip install dgl -f https://data.dgl.ai/wheels/torch-2.1/cu121/repo.html

echo "=== Verifying installation ==="
python -c "
import torch
print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
print('cuda available:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('gpu:', torch.cuda.get_device_name(0))

import torch_geometric
print('torch_geometric:', torch_geometric.__version__)

from torch_geometric.utils.smiles import x_map
print('hybridization types:', x_map['hybridization'])
print('SP2D supported:', 'SP2D' in x_map['hybridization'])

import dgl
print('dgl:', dgl.__version__)

import dgllife
print('dgllife:', dgllife.__version__)
"

echo "=== Done at $(date) ==="
