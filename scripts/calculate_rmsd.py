#!/usr/bin/env python3
"""
Calculate heavy-atom RMSD between a crystallographic ligand and each docked pose.

This script compares coordinates in the original receptor reference frame.
It does not translate or rotate the molecules before calculating RMSD.
Equivalent atom mappings are identified from molecular connectivity, and the
mapping with the lowest RMSD is reported.

Example:
    python calculate_rmsd.py \
        --native ligand/CAU_native.pdb \
        --poses results/CAU_redocked.sdf \
        --csv results/rmsd_results.csv

Requirements:
    pip install rdkit numpy networkx
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import networkx as nx
import numpy as np
from rdkit import Chem


def load_native_ligand(path: Path) -> Chem.Mol:
    """Load the crystallographic ligand from a PDB file."""
    mol = Chem.MolFromPDBFile(
        str(path),
        removeHs=False,
        sanitize=False,
    )
    if mol is None:
        raise ValueError(f"Could not read native ligand: {path}")

    # RMSD is calculated over heavy atoms only.
    return Chem.RemoveHs(mol, sanitize=False)


def load_docked_poses(path: Path) -> list[Chem.Mol]:
    """Load all docked poses from a multi-molecule SDF file."""
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=False,
    )

    poses: list[Chem.Mol] = []
    for index, mol in enumerate(supplier, start=1):
        if mol is None:
            raise ValueError(f"Could not read pose {index} from {path}")
        poses.append(Chem.RemoveHs(mol, sanitize=False))

    if not poses:
        raise ValueError(f"No docked poses were found in: {path}")

    return poses


def molecule_to_graph(mol: Chem.Mol) -> nx.Graph:
    """
    Convert an RDKit molecule to a connectivity graph.

    Node matching uses atomic number. Bond order is intentionally ignored
    because PDB files may not preserve ligand bond orders reliably.
    """
    graph = nx.Graph()

    for atom in mol.GetAtoms():
        graph.add_node(
            atom.GetIdx(),
            atomic_number=atom.GetAtomicNum(),
        )

    for bond in mol.GetBonds():
        graph.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
        )

    return graph


def calculate_best_rmsd(native: Chem.Mol, pose: Chem.Mol) -> float:
    """
    Calculate the minimum heavy-atom RMSD across valid graph isomorphisms.

    Coordinates remain in their original receptor coordinate frame; no
    least-squares alignment is performed.
    """
    if native.GetNumAtoms() != pose.GetNumAtoms():
        raise ValueError(
            "Native ligand and docked pose have different heavy-atom counts: "
            f"{native.GetNumAtoms()} vs {pose.GetNumAtoms()}"
        )

    native_graph = molecule_to_graph(native)
    pose_graph = molecule_to_graph(pose)

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        native_graph,
        pose_graph,
        node_match=lambda left, right: (
            left["atomic_number"] == right["atomic_number"]
        ),
    )

    native_xyz = np.asarray(native.GetConformer().GetPositions(), dtype=float)
    pose_xyz = np.asarray(pose.GetConformer().GetPositions(), dtype=float)

    best_rmsd = math.inf
    mapping_found = False

    # Mapping direction: native atom index -> pose atom index
    for mapping in matcher.isomorphisms_iter():
        mapping_found = True
        native_indices = sorted(mapping)
        pose_indices = [mapping[index] for index in native_indices]

        squared_distances = np.sum(
            (native_xyz[native_indices] - pose_xyz[pose_indices]) ** 2,
            axis=1,
        )
        rmsd = float(np.sqrt(np.mean(squared_distances)))
        best_rmsd = min(best_rmsd, rmsd)

    if not mapping_found:
        raise ValueError(
            "No atom mapping was found between the native ligand and pose."
        )

    return best_rmsd


def extract_vina_score(mol: Chem.Mol) -> float | None:
    """Extract the Vina free-energy score from Meeko metadata when available."""
    if not mol.HasProp("meeko"):
        return None

    try:
        metadata = json.loads(mol.GetProp("meeko"))
        score = metadata.get("free_energy")
        return float(score) if score is not None else None
    except (json.JSONDecodeError, TypeError, ValueError):
        return None


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    """Write pose number, Vina score, and RMSD to a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["pose", "vina_score_kcal_mol", "rmsd_angstrom"],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate heavy-atom RMSD values for docked ligand poses "
            "relative to a crystallographic ligand."
        )
    )
    parser.add_argument(
        "--native",
        type=Path,
        default=Path("ligand/CAU_native.pdb"),
        help="Path to the crystallographic ligand PDB file.",
    )
    parser.add_argument(
        "--poses",
        type=Path,
        default=Path("results/CAU_redocked.sdf"),
        help="Path to the multi-pose docked SDF file.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("results/rmsd_results.csv"),
        help="Output CSV path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        native = load_native_ligand(args.native)
        poses = load_docked_poses(args.poses)

        rows: list[dict[str, object]] = []

        print("Pose  Vina score (kcal/mol)  RMSD (Ã…)")
        print("----  ---------------------  --------")

        for pose_number, pose in enumerate(poses, start=1):
            rmsd = calculate_best_rmsd(native, pose)
            score = extract_vina_score(pose)

            score_text = f"{score:.3f}" if score is not None else "N/A"
            print(f"{pose_number:>4}  {score_text:>21}  {rmsd:>8.3f}")

            rows.append(
                {
                    "pose": pose_number,
                    "vina_score_kcal_mol": (
                        f"{score:.3f}" if score is not None else ""
                    ),
                    "rmsd_angstrom": f"{rmsd:.3f}",
                }
            )

        write_csv(rows, args.csv)
        print(f"\nSaved CSV results to: {args.csv}")
        return 0

    except (FileNotFoundError, ValueError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
