import numpy as np
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGS, BFGS


def relax_structure(structure, fmax=1e-3, steps=500, minimizer="bfgs", constant_cell=False, **kw):
    ucf = structure if constant_cell else FrechetCellFilter(structure)
    kw["logfile"] = kw.get("logfile", None)
    (LBFGS if minimizer == "lbfgs" else BFGS)(ucf, **kw).run(fmax=fmax, steps=steps)


def get_force_constants(structure, calculator, supercell_matrix):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    ph_at = PhonopyAtoms(symbols=list(structure.symbols), cell=np.array(structure.cell),
                         scaled_positions=structure.get_scaled_positions(), masses=structure.get_masses())
    ph = Phonopy(ph_at, np.diag(supercell_matrix) if np.ndim(supercell_matrix) == 1 else supercell_matrix)
    ph.generate_displacements()
    forces = []
    for sc in ph.supercells_with_displacements:
        a = Atoms(symbols=sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        calc = type(calculator)(__file__)            # fresh wells around the (undisplaced) supercell
        ref = Atoms(symbols=ph.supercell.symbols, cell=ph.supercell.cell,
                    scaled_positions=ph.supercell.scaled_positions, pbc=True)
        calc.ref = ref.get_scaled_positions(wrap=False)
        a.calc = calc
        forces.append(a.get_forces().copy())
    ph.forces = forces
    ph.produce_force_constants()
    return ph
