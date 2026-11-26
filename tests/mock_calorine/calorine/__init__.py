"""TEST-ONLY stand-in for calorine (no NEP model needed).  NOT physics.

CPUNEP(model) is an Einstein crystal: every atom sits in a harmonic well of stiffness k (eV/Å^2, per species)
around its starting fractional position, displaced by -0.01 (fractional, along x) for O so that a relaxation
has something to do.  The harmonic free energy is then known analytically, which the campaign test checks.
Set MOCK_NEP_SOFT=1 to give O a negative stiffness (-> imaginary modes) for the flagging test.
"""
