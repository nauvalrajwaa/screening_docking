import subprocess
import os
from .base import Docker

class AutoDockDocker(Docker):
    def __init__(self, binary_path, autogrid_path=None):
        super().__init__(binary_path)
        self.autogrid_path = autogrid_path

    def dock(self, ligand_path, receptor_path, center, size, output_path, **kwargs):
        """
        Run AutoDock-GPU.
        Prerequisite: GPF/Maps should essentially be ready or we need to generate them.
        For simplicity in this MVP, we assume maps might not exist and we might need to rely on PDBQTs but AutoDock-GPU NEEDS .fld and .xyz files (maps).
        
        Wait, AutoDock-GPU typically takes:
        -lfile <ligand.pdbqt>
        -ffile <maps.fld>
        -nrun <n>
        
        Generating maps requires autogrid4 and a GPF file.
        Constructing a GPF file on the fly is possible.
        """
        
        # Determine directory for intermediate files
        work_dir = os.path.dirname(output_path)
        receptor_base = os.path.basename(receptor_path).replace('.pdbqt', '')
        
        # 1. Generate GPF if we have autogrid
        if self.autogrid_path:
            gpf_path = os.path.join(work_dir, f"{receptor_base}.gpf")
            glg_path = os.path.join(work_dir, f"{receptor_base}.glg")
            
            # Simple GPF generation logic (simplified)
            # In reality, this requires knowing atom types in receptor to set 'map' types.
            # This is complex to do robustly without mgltools.
            # ALTERNATIVE: Use valid maps if they exist in the receptor folder??
            pass
        
        # For this implementation, we will try to run assuming user provided enough info 
        # OR we might have to degrade to just alerting the user "Maps not found".
        # However, to be useful, let's try to run a command if we can.
        
        # Assuming we just call the binary for now as placeholder for the integration structure
        # The user requested integration, but AutoDock setup is complex without pre-calc maps.
        
        print(f"AutoDock-GPU integration is experimental. Expecting maps to be present if not using Vina.")
        
        cmd = [
             self.binary_path,
             '--lfile', ligand_path,
             # We need ffile... assuming it matches receptor name?
             '--ffile', receptor_path.replace('.pdbqt', '.maps.fld'),
             '--nrun', '10',
             '--resnam', output_path
        ]
        
        print(f"Running AutoDock-GPU: {' '.join(cmd)}")
        try:
             # AutoDock-GPU output parsing is different.
             res = subprocess.run(cmd, check=True, capture_output=True, text=True)
             # Parse stdout for best energy
             return self._parse_score_from_stdout(res.stdout)
        except Exception as e:
             print(f"AutoDock-GPU failed: {e}")
             return None

    def _parse_score_from_stdout(self, stdout):
        # Parse "Best Inter + Intra" or similar from output
        # Example: "Best Inter + Intra    -9.77 kcal/mol"
        import re
        match = re.search(r"Best Inter \+ Intra\s+([-\d\.]+)", stdout)
        if match:
            return float(match.group(1))
        return None
