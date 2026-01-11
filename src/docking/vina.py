import subprocess
import os
import re
from .base import Docker

class VinaDocker(Docker):
    def dock(self, ligand_path, receptor_path, center, size, output_path, **kwargs):
        """
        Run Vina docking.
        """
        log_path = output_path.replace('.pdbqt', '.log')
        
        cmd = [
            self.binary_path,
            '--receptor', receptor_path,
            '--ligand', ligand_path,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(size[0]),
            '--size_y', str(size[1]),
            '--size_z', str(size[2]),
            '--out', output_path,
            '--cpu', str(kwargs.get('cpu', 4))
        ]
        
        print(f"Running Vina: {' '.join(cmd)}")
        try:
            # Vina writes log to stdout if --log is not specified or compatible
            with open(log_path, 'w') as log_file:
                subprocess.run(cmd, check=True, stdout=log_file, stderr=subprocess.STDOUT, text=True)
            
            return self._parse_score(log_path)
        except subprocess.CalledProcessError as e:
            print(f"Vina failed: {e}")
            return None
        except Exception as e:
            print(f"Vina execution error: {e}")
            return None

    def _parse_score(self, log_path):
        """Parse the best affinity score from Vina log."""
        if not os.path.exists(log_path):
            return None
            
        best_score = None
        with open(log_path, 'r') as f:
            for line in f:
                # Mode |   affinity | dist from best mode
                #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
                # -----+------------+----------+----------
                #    1         -9.3      0.000      0.000
                if line.strip().startswith('1'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            best_score = float(parts[1])
                            return best_score
                        except ValueError:
                            pass
        return best_score
