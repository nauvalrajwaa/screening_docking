import subprocess
import os

            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(input=smiles)
        
        if process.returncode != 0:
            print(f"Obabel conversion failed: {stderr}")
            return False
            
        return os.path.exists(output_path)
    except Exception as e:
        print(f"Conversion error: {e}")
        return False
