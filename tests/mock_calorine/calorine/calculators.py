import os
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

K = {"Zr": 10.0, "Ce": 8.0, "O": 5.0}


class CPUNEP(Calculator):
    implemented_properties = ["energy", "forces", "stress", "free_energy"]

    def __init__(self, model_filename, **kw):
        super().__init__(**kw)
        if not os.path.isfile(model_filename):
            raise FileNotFoundError(model_filename)
        self.ref = None

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        sym = np.array(self.atoms.get_chemical_symbols())
        fail = os.environ.get("MOCK_NEP_FAIL")                 # "6" or "6,16": structures with that many Ce fail
        if fail and str(int((sym == "Ce").sum())) in fail.split(","):
            raise RuntimeError("mock NEP failure (MOCK_NEP_FAIL)")
        f = self.atoms.get_scaled_positions(wrap=False)
        if self.ref is None or len(self.ref) != len(f):        # wells fixed when first used
            self.ref = f.copy()
            self.ref[sym == "O", 0] -= 0.01
        k = np.array([K[s] for s in sym])
        if os.environ.get("MOCK_NEP_SOFT"):
            k[sym == "O"] = -1.0
        df = f - self.ref
        df -= np.round(df)
        d = df @ np.array(self.atoms.cell)
        e = 0.5 * float((k[:, None] * d ** 2).sum())
        v = self.atoms.get_volume()
        s = np.einsum("i,ia,ib->ab", k, d, d) / v          # d(E)/d(strain) / V
        self.results = {"energy": e, "free_energy": e, "forces": -k[:, None] * d,
                        "stress": np.array([s[0, 0], s[1, 1], s[2, 2], s[1, 2], s[0, 2], s[0, 1]])}
