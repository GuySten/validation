"""One-group analytical transport benchmarks from

    B. D. Ganapol, "Analytical Benchmarks for Nuclear Engineering
    Applications: Case Studies in Neutron Transport Theory,"
    NEA/DB/DOC(2008)1, OECD Nuclear Energy Agency (2008).

All problems use a total cross section of 1 cm^-1 so that lengths are in mean
free paths, and isotropic scattering with c secondaries per collision. Only
the converged (most accurate) column of each qualification table is kept.

"""

from dataclasses import dataclass, field


# Benchmark 3.1, Table 3.1.1: scalar flux from a normally incident plane beam
# source at x = 0 in an infinite medium with c = 0.9. Values are (x, flux)
# for a source of unit strength.
INFINITE_MEDIUM_BEAM = {
    'c': 0.9,
    'flux': [
        (-5.0, 1.02506e-01),
        (-4.0, 1.73596e-01),
        (-3.0, 2.94474e-01),
        (-2.0, 5.01677e-01),
        (-1.0, 8.66764e-01),
        (-0.01, 1.66724e+00),
        (0.01, 2.71018e+00),
        (1.0, 2.17754e+00),
        (2.0, 1.42338e+00),
        (3.0, 8.83953e-01),
        (4.0, 5.36683e-01),
        (5.0, 3.22109e-01),
    ],
}

# Benchmark 3.2, Table 3.2.1: scalar flux in a half-space x > 0 with c = 0.99
# for a beam of unit strength normally incident on the free surface.
HALF_SPACE_BEAM = {
    'c': 0.99,
    'flux': [
        (1.0e-5, 2.47287e+00),
        (0.5, 3.23075e+00),
        (1.0, 3.31086e+00),
        (1.5, 3.21262e+00),
        (2.0, 3.04188e+00),
        (2.5, 2.84333e+00),
        (3.0, 2.63836e+00),
        (3.5, 2.43757e+00),
        (4.0, 2.24612e+00),
        (4.5, 2.06633e+00),
        (5.0, 1.89899e+00),
    ],
}

# Benchmark 3.3, Table 3.3.1: exiting angular fluxes (integrated over the
# azimuth) for a slab of thickness 1 with c = 0.9 and a normally incident beam
# of unit strength on the left surface. Values are (|mu|, flux).
SLAB_BEAM = {
    'c': 0.9,
    'thickness': 1.0,
    'reflected': [
        (1.0, 4.20002e-01),
        (0.9, 4.47770e-01),
        (0.8, 4.78944e-01),
        (0.7, 5.13947e-01),
        (0.6, 5.53090e-01),
        (0.5, 5.96303e-01),
        (0.4, 6.42534e-01),
        (0.3, 6.88481e-01),
        (0.2, 7.26613e-01),
        (0.1, 7.45339e-01),
    ],
    'transmitted': [
        (1.0, 3.80529e-01),
        (0.9, 4.01409e-01),
        (0.8, 4.23766e-01),
        (0.7, 4.47267e-01),
        (0.6, 4.71072e-01),
        (0.5, 4.93358e-01),
        (0.4, 5.10404e-01),
        (0.3, 5.15163e-01),
        (0.2, 4.96655e-01),
        (0.1, 4.48066e-01),
    ],
}

# Benchmark 3.4, Table 3.4.3(b): critical radii R of a bare infinite cylinder
# as a function of c. (The table repeats the c = 1.36 entry in place of
# c = 1.26, so c = 1.26 is not included.)
CRITICAL_CYLINDER_RADII = {
    1.01: 1.3125516e+01, 1.02: 9.0432549e+00, 1.03: 7.2314062e+00,
    1.04: 6.1499466e+00, 1.05: 5.4112883e+00, 1.06: 4.8657524e+00,
    1.07: 4.4416675e+00, 1.08: 4.0998327e+00, 1.09: 3.8167572e+00,
    1.10: 3.5773913e+00, 1.12: 3.1922162e+00, 1.14: 2.8934088e+00,
    1.16: 2.6531096e+00, 1.18: 2.4546128e+00, 1.20: 2.2872093e+00,
    1.22: 2.1436738e+00, 1.24: 2.0189307e+00, 1.28: 1.8120184e+00,
    1.30: 1.7250029e+00, 1.32: 1.6466161e+00, 1.34: 1.5755658e+00,
    1.36: 1.5108142e+00, 1.38: 1.4515168e+00, 1.40: 1.3969786e+00,
    1.42: 1.3466215e+00, 1.44: 1.2999606e+00, 1.46: 1.2565852e+00,
    1.48: 1.2161455e+00, 1.50: 1.1783408e+00, 1.52: 1.1429117e+00,
    1.54: 1.1096324e+00, 1.56: 1.0783058e+00, 1.58: 1.0487589e+00,
    1.60: 1.0208390e+00, 1.62: 9.9441094e-01, 1.64: 9.6935444e-01,
    1.66: 9.4556217e-01, 1.68: 9.2293801e-01, 1.70: 9.0139555e-01,
    1.72: 8.8085691e-01, 1.74: 8.6125168e-01, 1.76: 8.4251603e-01,
    1.78: 8.2459196e-01, 1.80: 8.0742662e-01, 1.82: 7.9097174e-01,
    1.84: 7.7518314e-01, 1.86: 7.6002032e-01, 1.88: 7.4544605e-01,
    1.90: 7.3142605e-01, 1.92: 7.1792869e-01, 1.94: 7.0492476e-01,
    1.96: 6.9238723e-01, 1.98: 6.8029102e-01, 2.00: 6.6861287e-01,
}

# Benchmark 3.4, Table 3.4.4: critical scalar flux distributions phi(r)/phi(0)
# for selected values of c. Values are (r/R, flux).
_CRITICAL_FLUX_POSITIONS = (0.0, 0.25, 0.50, 0.75, 0.85, 0.91, 0.95, 0.98,
                            1.0)
CRITICAL_CYLINDER_FLUX = {
    c: list(zip(_CRITICAL_FLUX_POSITIONS, values)) for c, values in {
        1.05: (1.0, 9.29851e-01, 7.33990e-01, 4.52168e-01, 3.26662e-01,
               2.49166e-01, 1.95805e-01, 1.53086e-01, 1.17908e-01),
        1.1: (1.0, 9.36052e-01, 7.56084e-01, 4.92189e-01, 3.71791e-01,
              2.96040e-01, 2.43013e-01, 1.99922e-01, 1.64122e-01),
        1.6: (1.0, 9.54996e-01, 8.24845e-01, 6.21823e-01, 5.22344e-01,
              4.56627e-01, 4.08837e-01, 3.68806e-01, 3.35065e-01),
        2.0: (1.0, 9.59783e-01, 8.42634e-01, 6.56963e-01, 5.64397e-01,
              5.02561e-01, 4.57218e-01, 4.18990e-01, 3.86650e-01),
    }.items()
}

# Benchmark 3.4, Tables 3.4.1 and 3.4.2: surrogate flux F(r) = 1 - (1 - c)
# phi(r)/Q for a uniform isotropic source Q in a bare infinite cylinder of
# radius R. Keys are (R, c); values are (r, F).
_R1 = [0.1*i for i in range(11)]
_R10 = [float(i) for i in range(11)]
SURROGATE_FLUX = {
    (1.0, 0.1): list(zip(_R1, (
        2.99256e-01, 3.01163e-01, 3.06966e-01, 3.16920e-01, 3.31495e-01,
        3.51447e-01, 3.77974e-01, 4.13044e-01, 4.60197e-01, 5.27241e-01,
        6.52895e-01))),
    (1.0, 0.15): list(zip(_R1, (
        3.13617e-01, 3.15533e-01, 3.21362e-01, 3.31351e-01, 3.45954e-01,
        3.65905e-01, 3.92362e-01, 4.27224e-01, 4.73899e-01, 5.39894e-01,
        6.62411e-01))),
    (1.0, 0.2): list(zip(_R1, (
        3.29165e-01, 3.31087e-01, 3.36930e-01, 3.46935e-01, 3.61539e-01,
        3.81450e-01, 4.07785e-01, 4.42369e-01, 4.88471e-01, 5.53285e-01,
        6.72459e-01))),
    (1.0, 0.25): list(zip(_R1, (
        3.46040e-01, 3.47964e-01, 3.53809e-01, 3.63807e-01, 3.78379e-01,
        3.98205e-01, 4.24356e-01, 4.58580e-01, 5.04003e-01, 5.67491e-01,
        6.83091e-01))),
    (1.0, 0.3): list(zip(_R1, (
        3.64405e-01, 3.66325e-01, 3.72157e-01, 3.82121e-01, 3.96622e-01,
        4.16309e-01, 4.42206e-01, 4.75976e-01, 5.20597e-01, 5.82596e-01,
        6.94369e-01))),
    (1.0, 0.5): list(zip(_R1, (
        4.57065e-01, 4.58897e-01, 4.64449e-01, 4.73895e-01, 4.87552e-01,
        5.05927e-01, 5.29816e-01, 5.60510e-01, 6.00306e-01, 6.54245e-01,
        7.47538e-01))),
    (1.0, 0.7): list(zip(_R1, (
        5.96607e-01, 5.98138e-01, 6.02766e-01, 6.10603e-01, 6.21852e-01,
        6.36839e-01, 6.56077e-01, 6.80397e-01, 7.11290e-01, 7.52048e-01,
        8.19470e-01))),
    (1.0, 0.9): list(zip(_R1, (
        8.24677e-01, 8.25431e-01, 8.27704e-01, 8.31533e-01, 8.36986e-01,
        8.44173e-01, 8.53271e-01, 8.64571e-01, 8.78607e-01, 8.96589e-01,
        9.24929e-01))),
    (10.0, 0.1): list(zip(_R10, (
        2.18198e-05, 2.86832e-05, 5.60906e-05, 1.32442e-04, 3.42077e-04,
        9.33715e-04, 2.66667e-03, 7.98668e-03, 2.54732e-02, 9.07754e-02,
        5.26092e-01))),
    (10.0, 0.15): list(zip(_R10, (
        2.56864e-05, 3.36496e-05, 6.53170e-05, 1.52913e-04, 3.91281e-04,
        1.05700e-03, 2.98339e-03, 8.81266e-03, 2.76356e-02, 9.62760e-02,
        5.33447e-01))),
    (10.0, 0.2): list(zip(_R10, (
        3.06790e-05, 4.00294e-05, 7.70364e-05, 1.78578e-04, 4.52101e-04,
        1.20704e-03, 3.36212e-03, 9.78051e-03, 3.01085e-02, 1.02371e-01,
        5.41241e-01))),
    (10.0, 0.25): list(zip(_R10, (
        3.72505e-05, 4.83781e-05, 9.21787e-05, 2.11254e-04, 5.28318e-04,
        1.39184e-03, 3.81968e-03, 1.09243e-02, 3.29550e-02, 1.09156e-01,
        5.49525e-01))),
    (10.0, 0.3): list(zip(_R10, (
        4.60882e-05, 5.95324e-05, 1.12121e-04, 2.53579e-04, 6.25301e-04,
        1.62250e-03, 4.37880e-03, 1.22887e-02, 3.62559e-02, 1.16748e-01,
        5.58360e-01))),
    (10.0, 0.5): list(zip(_R10, (
        1.38859e-04, 1.73680e-04, 3.05085e-04, 6.37744e-04, 1.44881e-03,
        3.44884e-03, 8.48704e-03, 2.15284e-02, 5.66444e-02, 1.58763e-01,
        6.00996e-01))),
    (10.0, 0.7): list(zip(_R10, (
        8.01829e-04, 9.46293e-04, 1.45897e-03, 2.62685e-03, 5.12494e-03,
        1.04520e-02, 2.19360e-02, 4.71009e-02, 1.03616e-01, 2.37138e-01,
        6.63330e-01))),
    (10.0, 0.9): list(zip(_R10, (
        2.01898e-02, 2.16079e-02, 2.61624e-02, 3.48252e-02, 4.94721e-02,
        7.33358e-02, 1.11811e-01, 1.73872e-01, 2.74762e-01, 4.42664e-01,
        7.81243e-01))),
}

# Benchmark 5.1, Table 5.1.1: searchlight problem. Scalar flux at depth
# z = 2 in a half-space with c = 0.99 from a beam of unit strength normally
# incident at a point on the surface. Values are (rho, flux).
SEARCHLIGHT = {
    'c': 0.99,
    'z': 2.0,
    'flux': [
        (0.20, 3.0111e-01),
        (1.18, 8.0206e-02),
        (2.16, 3.7865e-02),
        (3.14, 1.9927e-02),
        (4.12, 1.1239e-02),
        (5.10, 6.6777e-03),
        (6.08, 4.1345e-03),
        (7.06, 2.6470e-03),
        (8.04, 1.7423e-03),
        (9.02, 1.1738e-03),
        (10.0, 8.0654e-04),
    ],
}


@dataclass
class Problem:
    """One-group analytic transport benchmark.

    Parameters
    ----------
    name : str
        Unique name of the problem
    kind : str
        One of 'critical-cylinder', 'surrogate-cylinder', 'infinite-beam',
        'half-space-beam', 'slab-beam', or 'searchlight'
    c : float
        Mean number of secondaries per collision
    size : float
        Characteristic dimension in mean free paths (cylinder radius or slab
        thickness), if applicable
    keff : float
        Reference eigenvalue for criticality problems
    flux : list of tuple
        Reference (position, value) pairs. The meaning of the values depends
        on the kind of problem.
    reflected : list of tuple
        For 'slab-beam', reference (|mu|, angular flux) exiting the left
        surface
    transmitted : list of tuple
        For 'slab-beam', reference (|mu|, angular flux) exiting the right
        surface

    """
    name: str
    kind: str
    c: float
    size: float = None
    keff: float = None
    flux: list = field(default_factory=list)
    reflected: list = field(default_factory=list)
    transmitted: list = field(default_factory=list)


def _build_problems():
    problems = [
        Problem('infinite-beam', 'infinite-beam', INFINITE_MEDIUM_BEAM['c'],
                flux=INFINITE_MEDIUM_BEAM['flux']),
        Problem('half-space-beam', 'half-space-beam', HALF_SPACE_BEAM['c'],
                flux=HALF_SPACE_BEAM['flux']),
        Problem('slab-beam', 'slab-beam', SLAB_BEAM['c'],
                SLAB_BEAM['thickness'], reflected=SLAB_BEAM['reflected'],
                transmitted=SLAB_BEAM['transmitted']),
        Problem('searchlight', 'searchlight', SEARCHLIGHT['c'],
                SEARCHLIGHT['z'], flux=SEARCHLIGHT['flux']),
    ]
    for c, radius in CRITICAL_CYLINDER_RADII.items():
        flux = CRITICAL_CYLINDER_FLUX.get(c, [])
        problems.append(Problem(f'critical-cylinder-c{c:.2f}',
                                'critical-cylinder', c, radius, keff=1.0,
                                flux=flux))
    for (radius, c), flux in SURROGATE_FLUX.items():
        problems.append(Problem(f'surrogate-cylinder-R{radius:g}-c{c:g}',
                                'surrogate-cylinder', c, radius, flux=flux))
    return problems


PROBLEMS = _build_problems()
PROBLEMS_BY_NAME = {p.name: p for p in PROBLEMS}
