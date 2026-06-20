***Carazolol Redocking Validation in the β2-Adrenergic Receptor (PDB 2RH1)***
**Overview**

This project validates a molecular-docking workflow by redocking the crystallographic ligand carazolol (ligand code: CAU) into the β2-adrenergic receptor structure PDB 2RH1.

The objective was to determine whether AutoDock Vina could reproduce the experimentally observed ligand-binding pose. Docked poses were compared with the crystallographic ligand using heavy-atom root-mean-square deviation (RMSD).

**Tools**
AutoDock Vina 1.2.7
RDKit
UCSF ChimeraX
Python
NumPy
NetworkX

**Docking System**
Receptor: β2-adrenergic receptor
PDB structure: 2RH1
Ligand: Carazolol
Ligand code: CAU
Docking Parameters
Grid center: (-30.262, 9.439, 6.943)
Grid size: 21.390 × 18.353 × 14.183 Å
Exhaustiveness: 32

The docking box was centered on the crystallographic carazolol-binding site.

**Workflow**
Prepared the 2RH1 receptor structure.
Extracted and prepared the crystallographic carazolol ligand.
Generated receptor and ligand PDBQT files.
Defined a docking box around the native binding site.
Performed molecular redocking with AutoDock Vina.
Converted docked poses for structural analysis.
Calculated heavy-atom RMSD values between the crystallographic and docked poses.
Visualized the native and redocked ligand conformations in ChimeraX.
