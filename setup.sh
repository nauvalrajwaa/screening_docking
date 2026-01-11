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

    # Check AutoDock-GPU
    if ! command -v autodock-gpu &> /dev/null; then
        echo "[+] Installing AutoDock-GPU (CUDA)..."
        # Install build dependencies
        sudo apt-get install -y git build-essential

        # Clone and compile
        echo "[*] Cloning AutoDock-GPU..."
        if [ -d "AutoDock-GPU" ]; then
            rm -rf AutoDock-GPU
        fi
        git clone https://github.com/ccsb-scripps/AutoDock-GPU.git
        
        cd AutoDock-GPU
        echo "[*] Compiling AutoDock-GPU (DEVICE=CUDA)..."
        # Assuming CUDA paths are set or standard. The notebook set env vars, but we might rely on system defaults or user env.
        # Exporting notebook-like paths just in case, though rigorous setup might need more checks.
        export GPU_INCLUDE_PATH=/usr/local/cuda/include
        export GPU_LIBRARY_PATH=/usr/local/cuda/lib64
        
        make DEVICE=CUDA NUMWI=64
        
        # Install to local bin or move binary
        if [ -f "bin/autodock-gpu" ]; then
            echo "[+] Installation successful. Moving binary to /usr/local/bin..."
            sudo mv bin/autodock-gpu /usr/local/bin/autodock-gpu
        else
            echo "[!] Compilation failed or binary not found."
        fi
        
        cd ..
    else
        echo "[OK] AutoDock-GPU found."
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
