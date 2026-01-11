#!/bin/bash

echo "======================================"
echo "   MetaScreener Setup Script (Linux)"
echo "======================================"

# 1. Check Python
if ! command -v python3 &> /dev/null; then
    echo "[ERROR] Python 3 is not installed."
    exit 1
fi

echo "[*] Installing Python dependencies..."
pip install -r requirements.txt

# 2. Check System Dependencies (Linux/Debian/Ubuntu)
if command -v apt-get &> /dev/null; then
    echo "[*] checking system dependencies..."
    
    # Check OpenBabel
    if ! command -v obabel &> /dev/null; then
            echo "[+] Installing OpenBabel..."
            sudo apt-get update
            sudo apt-get install -y openbabel
    else
            echo "[OK] OpenBabel found."
    fi

    # Check Vina
    if ! command -v vina &> /dev/null; then
            echo "[+] Installing AutoDock Vina..."
            sudo apt-get install -y autodock-vina
    else
            echo "[OK] Vina found."
    fi
else
    echo "[!] 'apt-get' not found. This script supports Debian/Ubuntu-based Linux systems."
    echo "    Please install 'openbabel' and 'autodock-vina' manually using your package manager."
fi

echo "======================================"
echo "   Setup Complete!"
echo "======================================"
echo "Run example:"
echo "python main.py --help"
