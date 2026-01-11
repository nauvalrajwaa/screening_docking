#!/bin/bash

echo "======================================"
echo "   MetaScreener Setup Script"
echo "======================================"

# 1. Check Python
if ! command -v python3 &> /dev/null; then
    echo "[ERROR] Python 3 is not installed."
    exit 1
fi

echo "[*] Installing Python dependencies..."
pip install -r requirements.txt

# 2. Check Homebrew (Mac)
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! command -v brew &> /dev/null; then
        echo "[WARNING] Homebrew not found. Skipping system package checks."
    else
        echo "[*] Checking system dependencies via Homebrew..."
        
        # Check OpenBabel
        if ! command -v obabel &> /dev/null; then
             echo "[+] Installing OpenBabel..."
             brew install open-babel
        else
             echo "[OK] OpenBabel found."
        fi

        # Check Vina
        if ! command -v vina &> /dev/null; then
             echo "[+] Installing AutoDock Vina..."
             brew install vina
        else
             echo "[OK] Vina found."
        fi
    fi
else
    echo "[!] Non-Mac OS detected. Please install 'openbabel' and 'vina' manually."
fi

echo "======================================"
echo "   Setup Complete!"
echo "======================================"
echo "Run example:"
echo "python main.py --help"
