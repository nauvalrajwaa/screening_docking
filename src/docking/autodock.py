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
        
        # 1. Check if maps exist (AutoDock-GPU needs .maps.fld and .maps.xyz or .map files)
        # We generally check for .maps.fld
        fld_path = receptor_path.replace('.pdbqt', '.maps.fld')
        
        if not os.path.exists(fld_path) and self.autogrid_path:
             print(f"[AutoDock] Maps not found. Generating using {self.autogrid_path}...")
             if not self._generate_maps(receptor_path, center, size, work_dir):
                 print("[AutoDock] Map generation failed.")
                 return None
        
        if not os.path.exists(fld_path):
             print(f"[AutoDock] Error: Map file {fld_path} missing and could not be generated.")
             return None

        cmd = [
             self.binary_path,
             '--lfile', ligand_path,
             '--ffile', fld_path,
             '--nrun', '10',
             '--resnam', output_path
        ]
        
        if kwargs.get('gpu_id'):
            cmd.extend(['-D', str(kwargs.get('gpu_id'))])
        
        
        print(f"Running AutoDock-GPU: {' '.join(cmd)}")
        try:
             # AutoDock-GPU output parsing is different.
             res = subprocess.run(cmd, check=True, capture_output=True, text=True)
             # Parse stdout for best energy
             return self._parse_score_from_stdout(res.stdout)
        except subprocess.CalledProcessError as e:
             print(f"AutoDock-GPU failed with exit code {e.returncode}")
             print(f"[AutoDock-GPU] --- Error Output ---")
             print(e.stderr if e.stderr else e.stdout)
             print(f"[AutoDock-GPU] --- End Error ---")
             return None
        except Exception as e:
             print(f"AutoDock-GPU failed: {e}")
             return None

    def _generate_maps(self, receptor_path, center, size, work_dir):
        """
        Generate GPF and run AutoGrid4.
        """
        try:
            receptor_base = os.path.basename(receptor_path).replace('.pdbqt', '')
            gpf_path = os.path.join(os.path.dirname(receptor_path), f"{receptor_base}.gpf")
            glg_path = os.path.join(os.path.dirname(receptor_path), f"{receptor_base}.glg")
            
            # 1. Identify Atom Types in Receptor
            atom_types = set()
            with open(receptor_path, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        parts = line.split()
                        if len(parts) >= 12:  # PDBQT format has atom type at the end
                            atype = parts[-1]
                            if atype and len(atype) <= 2:
                                atom_types.add(atype)
            
            # Standard AutoDock atom types for ligands
            ligand_types = ['C', 'A', 'N', 'NA', 'OA', 'S', 'SA', 'HD', 'H']
            
            # 2. Write GPF (AutoDock GPF format)
            spacing = 0.375
            npts = [int(s / spacing) for s in size]
            
            receptor_abs = os.path.abspath(receptor_path)
            receptor_dir = os.path.dirname(receptor_abs)
            
            with open(gpf_path, 'w') as gpf:
                gpf.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
                gpf.write(f"gridfld {receptor_base}.maps.fld\n")
                gpf.write(f"spacing {spacing}\n")
                gpf.write(f"receptor_types {' '.join(sorted(atom_types))}\n")
                gpf.write(f"ligand_types {' '.join(ligand_types)}\n")
                gpf.write(f"receptor {receptor_base}.pdbqt\n")
                gpf.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
                gpf.write("smooth 0.5\n")
                
                # Map files (relative paths)
                for at in ligand_types:
                    gpf.write(f"map {receptor_base}.{at}.map\n")
                
                gpf.write(f"elecmap {receptor_base}.e.map\n")
                gpf.write(f"dsolvmap {receptor_base}.d.map\n")
                gpf.write("dielectric -0.1465\n")
            
            # 3. Run AutoGrid4 (from receptor directory)
            print(f"[AutoDock] Running AutoGrid4 on {gpf_path}...")
            cmd = [self.autogrid_path, '-p', os.path.basename(gpf_path), '-l', os.path.basename(glg_path)]
            
            result = subprocess.run(
                cmd, 
                cwd=os.path.dirname(receptor_abs),
                capture_output=True, 
                text=True
            )
            
            if result.returncode != 0:
                print(f"[AutoDock] AutoGrid4 failed with exit code {result.returncode}")
                print(f"[AutoDock] Log file: {glg_path}")
                if os.path.exists(glg_path):
                    with open(glg_path, 'r') as log:
                        print(f"[AutoDock] --- AutoGrid4 Log ---")
                        print(log.read())
                        print(f"[AutoDock] --- End Log ---")
                return False
            
            return True
            
        except Exception as e:
            print(f"[AutoDock] GPF/AutoGrid generation failed: {e}")
            if 'glg_path' in locals() and os.path.exists(glg_path):
                with open(glg_path, 'r') as log:
                    print(f"[AutoDock] --- AutoGrid4 Log ---")
                    print(log.read())
                    print(f"[AutoDock] --- End Log ---")
            return False

    def _parse_score_from_stdout(self, stdout):
        # Parse "Best Inter + Intra" or similar from output
        # Example output line: "Best Inter + Intra    -9.77 kcal/mol"
        # The score is typically on a line that looks like a table row
        import re
        
        # Look for pattern with actual numeric value, not dashes
        match = re.search(r'Best Inter \+ Intra\s+([-\d\.]+)\s+kcal/mol', stdout)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                pass
        
        # Fallback: search for any line mentioning "best" and a score
        for line in stdout.split('\n'):
            if 'best inter + intra' in line.lower() and 'kcal/mol' in line.lower():
                # Extract the first float-like number
                numbers = re.findall(r'[-+]?\d+\.?\d*', line)
                for num in numbers:
                    try:
                        val = float(num)
                        if -50 < val < 50:  # Reasonable binding affinity range
                            return val
                    except ValueError:
                        continue
        
        return None
