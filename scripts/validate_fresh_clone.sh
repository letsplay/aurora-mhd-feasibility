#!/bin/bash
set -e
echo "=== AURORA-MHD Fresh Clone Validation ==="

echo "1. Creating temp directory..."
TMPDIR=$(mktemp -d)
cd "$TMPDIR"

echo "2. Cloning repo..."
git clone https://github.com/letsplay/aurora-mhd-feasibility.git
cd aurora-mhd-feasibility

echo "3. Setting up venv..."
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt --quiet

echo "4. Running tests..."
pytest tests/ -v

echo "5. Executing notebooks..."
jupyter nbconvert --execute --to notebook notebooks/01_magnet_trade.ipynb
jupyter nbconvert --execute --to notebook notebooks/02_energy_closure.ipynb
jupyter nbconvert --execute --to notebook notebooks/03_operating_envelope.ipynb

echo "6. Checking outputs..."
test -f results/figures/d1_trade.png && echo "  ✓ d1_trade.png"
test -f results/figures/d2_balance.png && echo "  ✓ d2_balance.png"
test -f results/figures/d3_envelope.png && echo "  ✓ d3_envelope.png"

echo ""
echo "=== ALL PASSED ==="
rm -rf "$TMPDIR"