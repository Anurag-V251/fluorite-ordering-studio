#!/usr/bin/env python3
"""
Interactive harmonic free-energy calculations with:
    NEP model -> Calorine CPUNEP -> Phonopy
For every selected structure, this script writes:
    free_energy_<structure_name>.csv
    free_energy_<structure_name>_metadata.txt
Energy normalization in the CSV:
    E_static, F_vib, F_total : eV/atom
Definitions:
    F_total(T) = E_static + F_vib(T)
This is a fixed-cell harmonic calculation. It is NOT quasi-harmonic:
there is no thermal expansion or temperature-dependent equilibrium volume.
"""
# FIX 1: restored mangled dunder name (**future** -> __future__)
from __future__ import annotations

from pathlib import Path
import glob
import re
import sys
import numpy as np
from ase.io import read
from ase.units import kJ, mol
from calorine.calculators import GPUNEP,CPUNEP
from calorine.tools import get_force_constants, relax_structure


# ----------------------- User-interface helpers ----------------------- #

def ask_float(prompt: str, default: float) -> float:
    """Ask for a float, accepting a default value on empty input."""
    answer = input(f"{prompt} [{default}]: ").strip()
    return default if answer == "" else float(answer)


def ask_int(prompt: str, default: int) -> int:
    """Ask for an integer, accepting a default value on empty input."""
    answer = input(f"{prompt} [{default}]: ").strip()
    return default if answer == "" else int(answer)


def ask_yes_no(prompt: str, default: bool = False) -> bool:
    """Ask a yes/no question."""
    default_text = "Y/n" if default else "y/N"
    answer = input(f"{prompt} [{default_text}]: ").strip().lower()
    if answer == "":
        return default
    if answer in ("y", "yes"):
        return True
    if answer in ("n", "no"):
        return False
    raise ValueError("Please answer y/yes or n/no.")


def safe_stem(filename: str) -> str:
    """Produce a safe filename stem."""
    stem = Path(filename).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)


def choose_files(pattern: str, title: str, allow_multiple: bool) -> list[str]:
    """Interactively select file(s) found with glob."""
    files = sorted(glob.glob(pattern))
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)
    if files:
        for i, filename in enumerate(files, start=1):
            print(f"  {i:2d}: {filename}")
    else:
        print(f"  No files found with pattern: {pattern}")

    if allow_multiple:
        prompt = (
            "\nEnter comma-separated indices, filenames, or glob patterns "
            "(for example: 1,2 or POSCAR_*): "
        )
    else:
        prompt = "\nEnter one index, filename, or glob pattern: "

    selection = input(prompt).strip()
    if not selection:
        raise ValueError("No selection was provided.")

    selected = []
    for item in selection.split(","):
        item = item.strip()
        if not item:
            continue
        if item.isdigit() and 1 <= int(item) <= len(files):
            selected.append(files[int(item) - 1])
        else:
            matches = sorted(glob.glob(item))
            if matches:
                selected.extend(matches)
            elif Path(item).is_file():
                selected.append(item)
            else:
                raise FileNotFoundError(f"Cannot find: {item}")

    # Preserve order and remove duplicates.
    selected = list(dict.fromkeys(selected))

    if not allow_multiple and len(selected) != 1:
        raise ValueError("Please select exactly one NEP model file.")

    # FIX 5: guard against an empty selection list
    if not selected:
        raise ValueError("No valid files were selected.")

    return selected


# ----------------------- Calculation ----------------------- #

def calculate_free_energy(
    structure_file: str,
    nep_file: str,
    supercell: list[int],
    mesh: int,
    t_min: float,
    t_max: float,
    t_step: float,
    relax: bool,
    fmax: float,
) -> tuple[dict, str]:
    """
    Calculate E_static, F_vib(T), and F_total(T) for one structure.

    Returns
    -------
    data : dict
        Arrays and metadata for writing to CSV.
    metadata : str
        Text report for a companion metadata file.
    """
    atoms = read(structure_file)

    # A fresh calculator is used for each structure.
    #calc = GPUNEP(nep_file,command="/opt/GPUMD/src/gpumd")
    calc=CPUNEP(nep_file)
    atoms.calc = calc

    if relax:
        print(f"  Relaxing {structure_file} with {nep_file} with gpu ...")
        relax_structure(atoms, fmax=fmax,minimizer='lbfgs')
        print("done relaxing")
    #print("changing to CPUNEP")
    #cal=CPUNEP(nep_file)
    #atoms.calc=cal
    # Static NEP energy, normalized to eV/atom.
    e_static_eV_atom = atoms.get_potential_energy() / len(atoms)

    # Force constants from finite displacements evaluated with the NEP.
    print(f"  Building force constants: supercell = {supercell}")
    phonon = get_force_constants(atoms, calc, supercell)

    # Phonopy mesh sampling — capture the Mesh result object (Phonopy v4 API).
    # FIX 4: run_mesh() returns a Mesh result object; store it to access
    #         frequencies without relying on the transitional phonon.mesh
    #         property, which is planned for deprecation in Phonopy v5.
    print(f"  Running Phonopy mesh: {mesh}")
    mesh_result = phonon.run_mesh(mesh)

    # FIX 3 & 6: use the current Phonopy v4 result-object API.
    #   run_thermal_properties() returns a ThermalProperties object whose
    #   attributes are the stable public surface:
    #     .temperatures  [K]
    #     .free_energy   [kJ/mol per primitive cell]
    #     .entropy       [J/(mol·K) per primitive cell]
    #     .heat_capacity [J/(mol·K) per primitive cell]
    #   The legacy get_thermal_properties_dict() is deprecated.
    print("  Computing thermal properties ...")
    tp = phonon.run_thermal_properties(
        t_min=t_min,
        t_max=t_max,
        t_step=t_step,
    )

    temperature_K                = np.asarray(tp.temperatures,   dtype=float)
    f_vib_kJ_mol_primitive       = np.asarray(tp.free_energy,    dtype=float)
    entropy_J_molK_primitive     = np.asarray(tp.entropy,        dtype=float)
    cv_J_molK_primitive          = np.asarray(tp.heat_capacity,  dtype=float)

    # Important: normalize relative to the primitive cell selected by Phonopy,
    # not necessarily the input POSCAR cell.
    n_primitive_atoms = len(phonon.primitive)

    # kJ/mol per primitive cell -> eV/primitive cell -> eV/atom.
    # ASE units: kJ ~ 0.01036 eV, mol ~ 1/6.022e23, so kJ/mol ~ 1.036e-2 eV.
    # Confirmed against ase.units documentation: 100 * kJ/mol == 1.0364 eV.
    f_vib_eV_atom = (
        f_vib_kJ_mol_primitive * kJ / mol / n_primitive_atoms
    )

    # Entropy and Cv per atom, useful for downstream analysis.
    entropy_J_molK_atom = entropy_J_molK_primitive / n_primitive_atoms
    cv_J_molK_atom      = cv_J_molK_primitive      / n_primitive_atoms

    # Total harmonic Helmholtz free energy.
    f_total_eV_atom = e_static_eV_atom + f_vib_eV_atom

    # Diagnostic: a stable harmonic reference should not have appreciable
    # negative (imaginary) phonon frequencies.
    # FIX 4 (continued): read frequencies from the captured Mesh result object
    #   instead of the deprecated phonon.mesh.frequencies accessor.
    frequencies_THz       = np.asarray(mesh_result.frequencies)
    minimum_frequency_THz = float(np.min(frequencies_THz))
    n_negative            = int(np.count_nonzero(frequencies_THz < -0.05))

    data = {
        "temperature_K":       temperature_K,
        "E_static_eV_atom":    np.full_like(temperature_K, e_static_eV_atom, dtype=float),
        "F_vib_eV_atom":       f_vib_eV_atom,
        "F_total_eV_atom":     f_total_eV_atom,
        "entropy_J_molK_atom": entropy_J_molK_atom,
        "Cv_J_molK_atom":      cv_J_molK_atom,
    }

    metadata = f"""Harmonic free-energy calculation
================================
Structure input file: {structure_file}
NEP model file:       {nep_file}
Relaxed with NEP:     {relax}
Relaxation fmax:      {fmax} eV/Angstrom
Input-cell atoms:     {len(atoms)}
Phonopy primitive-cell atoms: {n_primitive_atoms}
Phonon supercell:     {supercell}
Phonopy mesh density: {mesh}
Temperature grid:     {t_min} to {t_max} K; step = {t_step} K
Static NEP energy:    {e_static_eV_atom:.12f} eV/atom
F_vib at first T:     {f_vib_eV_atom[0]:.12f} eV/atom
F_total at first T:   {f_total_eV_atom[0]:.12f} eV/atom
Minimum mesh frequency: {minimum_frequency_THz:.8f} THz
Modes below -0.05 THz:  {n_negative}
Interpretation:
  F_total(T) = E_static + F_vib(T)
WARNING:
  If significant imaginary modes are present, the harmonic free energy
  is not a reliable thermodynamic free energy for that structure.
"""
    return data, metadata


def write_results(output_stem: str, data: dict, metadata: str) -> None:
    """Write CSV data and metadata report."""
    csv_file  = Path(f"free_energy_{output_stem}.csv")
    meta_file = Path(f"free_energy_{output_stem}_metadata.txt")

    header = (
        "temperature_K,"
        "E_static_eV_atom,"
        "F_vib_eV_atom,"
        "F_total_eV_atom,"
        "entropy_J_molK_atom,"
        "Cv_J_molK_atom"
    )
    table = np.column_stack([
        data["temperature_K"],
        data["E_static_eV_atom"],
        data["F_vib_eV_atom"],
        data["F_total_eV_atom"],
        data["entropy_J_molK_atom"],
        data["Cv_J_molK_atom"],
    ])
    np.savetxt(
        csv_file,
        table,
        delimiter=",",
        header=header,
        comments="",
        fmt="%.12e",
    )
    meta_file.write_text(metadata, encoding="utf-8")
    print(f"  Wrote: {csv_file}")
    print(f"  Wrote: {meta_file}")


# ----------------------- Main program ----------------------- #

def main() -> None:
    print("\nNEP + Calorine + Phonopy: harmonic free-energy calculation")
    print("Run this script in the directory containing structure and NEP files.")

    # Select structures. The default catches common POSCAR-style filenames.
    structure_pattern = input(
        "\nStructure-file pattern [POSCAR*]: "
    ).strip() or "POSCAR*"
    structure_files = choose_files(
        structure_pattern,
        "Available structure files",
        allow_multiple=True,
    )

    # FIX 5 (continued): explicit post-selection guard.
    if not structure_files:
        raise ValueError("No structure files were selected. Aborting.")

    # Select exactly one NEP file.
    nep_pattern = input("\nNEP-file pattern [*.txt]: ").strip() or "*.txt"
    nep_file = choose_files(
        nep_pattern,
        "Available NEP model files",
        allow_multiple=False,
    )[0]

    print("\nCalculation settings")
    t_min  = ask_float("Minimum temperature, K",  0.0)
    t_max  = ask_float("Maximum temperature, K",  1000.0)
    t_step = ask_float("Temperature step, K",     10.0)
    if t_min < 0 or t_max <= t_min or t_step <= 0:
        raise ValueError("Require 0 <= Tmin < Tmax and Tstep > 0.")

    sc_text  = input("Phonon supercell, e.g. 5 5 5 [5 5 5]: ").strip() or "5 5 5"
    supercell = [int(x) for x in sc_text.replace(",", " ").split()]
    if len(supercell) != 3 or any(x < 1 for x in supercell):
        raise ValueError("The supercell must contain three positive integers.")

    mesh  = ask_int("Phonopy q-point mesh density", 100)
    relax = ask_yes_no(
        "Relax atomic positions/cell with the NEP before phonons?",
        default=False,
    )
    fmax = ask_float("Relaxation force threshold, eV/Angstrom", 1e-4) if relax else 1e-4

    print("\nSelected structures:")
    for filename in structure_files:
        print(f"  - {filename}")
    print(f"Selected NEP model: {nep_file}\n")

    for structure_file in structure_files:
        print("=" * 72)
        print(f"Processing: {structure_file}")
        print("=" * 72)
        try:
            data, metadata = calculate_free_energy(
                structure_file=structure_file,
                nep_file=nep_file,
                supercell=supercell,
                mesh=mesh,
                t_min=t_min,
                t_max=t_max,
                t_step=t_step,
                relax=relax,
                fmax=fmax,
            )
            write_results(safe_stem(structure_file), data, metadata)
            print(
                "  E_static = "
                f"{data['E_static_eV_atom'][0]:.8f} eV/atom"
            )
            print(
                "  F_total at Tmin = "
                f"{data['F_total_eV_atom'][0]:.8f} eV/atom"
            )
        except Exception as exc:
            print(f"\nERROR while processing {structure_file}: {exc}")
            print("Continuing with the next selected structure.\n")

    print("\nFinished.")


# FIX 1: restored mangled dunder name (**name** -> __name__)
if __name__ == "__main__":
    main()
