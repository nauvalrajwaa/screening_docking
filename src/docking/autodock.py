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
        except Exception as e:
             print(f"AutoDock-GPU failed: {e}")
             return None

    def _generate_maps(self, receptor_path, center, size, work_dir):
        """
        Generate GPF and run AutoGrid4.
        """
        try:
            receptor_base = os.path.basename(receptor_path).replace('.pdbqt', '')
            gpf_path = os.path.join(os.path.dirname(receptor_path), f"{receptor_base}.gpf") # Save GPF next to receptor
            glg_path = os.path.join(os.path.dirname(receptor_path), f"{receptor_base}.glg")
            
            # 1. Identify Atom Types in Receptor
            atom_types = set()
            with open(receptor_path, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # PDBQT atom type is usually last column (or near last)
                        # Example: ATOM      1  N   MET A  41      26.960 -10.590  -9.088  1.00  0.00     0.063 N
                        parts = line.split()
                        atype = parts[-1]
                        # Valid AutoDock types usually: C, A, N, O, S, H, HD, HS, NA, OA, SA, Cl, Br, F, I, Mg, Mn, Zn, Ca, Fe
                        # We just take what is there assuming PDBQT is valid
                        if atype and len(atype) <= 2:
                            atom_types.add(atype)
            
            # Common ligand types to support in map (Union of Receptor + Ligand possibilities)
            # Simplification: standard set
            common_types = ['A', 'C', 'HD', 'N', 'NA', 'OA', 'S', 'SA', 'Cl', 'F', 'Br', 'I'] 
            # Merge with receptor types (only valid ones) to ensure we cover binding site
            # Actually, AutoGrid needs 'map' lines for LIGAND atom types.
            # The receptor types are used for the 'receptor' parameter but maps are generated for probes.
            # We should include all standard types to be safe.
            
            # 2. Write GPF
            # npts calculation (spacing 0.375 default)
            spacing = 0.375
            npts = [int(s / spacing) for s in size]
            
            with open(gpf_path, 'w') as gpf:
                gpf.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
                gpf.write(f"gridfld {receptor_path.replace('.pdbqt', '.maps.fld')}\n")
                gpf.write(f"spacing {spacing}\n")
                gpf.write(f"receptor_types {' '.join(sorted(atom_types))}\n") # Types present in receptor
                gpf.write(f"ligand_types {' '.join(common_types)}\n") # Types in ligand
                gpf.write(f"receptor {receptor_path}\n")
                gpf.write(f"gridcenter {center[0]} {center[1]} {center[2]}\n")
                gpf.write("smooth 0.5\n")
                
                # Map files
                for at in common_types:
                    gpf.write(f"map {receptor_path.replace('.pdbqt', f'.{at}.map')}\n")
                
                gpf.write(f"elecmap {receptor_path.replace('.pdbqt', '.e.map')}\n")
                gpf.write(f"dsolvmap {receptor_path.replace('.pdbqt', '.d.map')}\n")
                gpf.write("dielectric -0.1465\n")
            
            # 3. Run AutoGrid4
            print(f"[AutoDock] Running AutoGrid4 on {gpf_path}...")
            cmd = [self.autogrid_path, '-p', gpf_path, '-l', glg_path]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            
            return True
            
        except Exception as e:
            print(f"[AutoDock] GPF/AutoGrid generation failed: {e}")
            return False

    def _parse_score_from_stdout(self, stdout):
        # Parse "Best Inter + Intra" or similar from output
        # Example: "Best Inter + Intra    -9.77 kcal/mol"
        import re
        match = re.search(r"Best Inter \+ Intra\s+([-\d\.]+)", stdout)
        if match:
            return float(match.group(1))
        return None
