"""Consistency checks of the transcribed analytic benchmark data.

These tests do not run OpenMC. They verify the cross sections, dimensions,
and reference results against the redundant information given in LA-13511 and
NEA/DB/DOC(2008)1 (balance of the cross sections, secondaries per collision,
k-infinity, flux ratios, and critical dimensions quoted in both centimeters
and mean free paths).
"""

import numpy as np
import pytest

from validation.analytic import ganapol, sood
from validation.analytic.sood import MATERIALS, PROBLEMS


# Number of secondaries per collision quoted for each one-group cross section
# set in Tables 5, 12, 16, 20, 24, 26, and 28 of LA-13511
SECONDARIES = {
    'PUa': 1.50, 'PUb': 1.40, 'Ua': 1.30, 'Ub': 1.3194202, 'Uc': 1.3014616,
    'Ud': 1.2958396, 'H2O': 0.90, 'UD2O': 1.02, 'H2O-UD2O': 0.90,
    'Ue': 1.230, 'Fe': 0.9980, 'Na': 1.00, 'PUa-P1': 1.40, 'PUa-P2': 1.40,
    'PUb-P1': 1.40, 'PUb-P2': 1.40, 'Ua-P1': 1.30, 'Ub-P1': 1.30,
    'UD2Oa-P1': 1.0308381, 'UD2Ob-P1': 1.0341086, 'UD2Oc-P1': 1.01964,
}


@pytest.mark.parametrize('material', MATERIALS.values(), ids=lambda m: m.name)
def test_cross_section_balance(material):
    """Total = capture + fission + outscatter in every group."""
    assert material.balance() == pytest.approx(0., abs=1e-6)


@pytest.mark.parametrize('material', MATERIALS.values(), ids=lambda m: m.name)
def test_fission_spectrum(material):
    if material.fissionable:
        assert material.chi.sum() == pytest.approx(1.0)
    else:
        assert not material.chi.any()


@pytest.mark.parametrize('name', SECONDARIES)
def test_secondaries_per_collision(name):
    assert MATERIALS[name].c()[0] == pytest.approx(SECONDARIES[name],
                                                   rel=1e-6)


def test_problem_numbers():
    assert [p.number for p in PROBLEMS] == list(range(1, 76))
    assert len(sood.PROBLEMS_BY_NAME) == 75


_INFINITE = [p for p in PROBLEMS if p.geometry == 'infinite']


@pytest.mark.parametrize('problem', _INFINITE, ids=lambda p: p.name)
def test_k_infinity(problem):
    """The cross sections reproduce the published k-infinity and flux
    ratios."""
    material = problem.materials[0]
    k, flux = material.k_inf()
    # The quoted nu values of U-235 (c) and (d) are rounded to seven digits
    assert k == pytest.approx(problem.keff, rel=5e-6)
    G = material.num_groups
    for (a, b), ratio in problem.flux_ratios.items():
        assert flux[G - a]/flux[G - b] == pytest.approx(ratio, rel=1e-6)


def _dimensions(problem):
    """Region dimensions [cm] and the total cross section used to convert
    them to mean free paths, in the order quoted by LA-13511."""
    regions = problem.regions
    if problem.geometry in ('cylinder', 'sphere'):
        inner = [0.] + [r for _, r in regions[:-1]]
        return [(r - r0, sood.MATERIALS[m].total[0])
                for (m, r), r0 in zip(regions, inner)]
    dims = [(t, sood.MATERIALS[m].total[0]) for m, t in regions]
    if problem.boundary == ('reflective', 'reflective'):
        # Infinite slab lattice cell: half of the fuel and moderator
        return dims
    if problem.number == 30:
        return dims
    # Critical half-thickness of the fuel followed by the reflector thickness
    fuel = [d for d, (m, _) in zip(dims, regions)
            if sood.MATERIALS[m].fissionable][0]
    reflector = [d for d, (m, _) in zip(dims, regions)
                 if not sood.MATERIALS[m].fissionable]
    return [(fuel[0]/2, fuel[1])] + reflector[:1]


_FINITE = [p for p in PROBLEMS if p.geometry != 'infinite']


@pytest.mark.parametrize('problem', _FINITE, ids=lambda p: p.name)
def test_dimensions(problem):
    """Dimensions in cm are consistent with those in mean free paths."""
    dims = _dimensions(problem)
    assert len(problem.mfp) >= 1
    for (cm, sigma_t), mfp in zip(dims, problem.mfp):
        assert cm*sigma_t == pytest.approx(mfp, rel=2e-4)


@pytest.mark.parametrize('problem', [p for p in _FINITE if p.flux_shape],
                         ids=lambda p: p.name)
def test_flux_positions_inside_problem(problem):
    if problem.geometry == 'slab':
        extent = sum(t for _, t in problem.regions)
    else:
        extent = problem.regions[-1][1]
    for position, _, value in problem.flux_shape:
        assert 0. <= position <= extent + 1e-9
        assert value > 0.


@pytest.mark.parametrize('number, c', [(7, 1.40), (13, 1.30), (23, 1.02)])
def test_cylinder_radii_agree_between_references(number, c):
    """Critical cylinder radii of LA-13511 agree with Ganapol Table 3.4.3."""
    problem = sood.PROBLEMS_BY_NUMBER[number]
    assert problem.mfp[0] == pytest.approx(
        ganapol.CRITICAL_CYLINDER_RADII[c], rel=2e-6)


def test_critical_radii_decrease_with_c():
    c = np.array(list(ganapol.CRITICAL_CYLINDER_RADII))
    radius = np.array(list(ganapol.CRITICAL_CYLINDER_RADII.values()))
    order = np.argsort(c)
    assert np.all(np.diff(radius[order]) < 0.)


@pytest.mark.parametrize('c', ganapol.CRITICAL_CYLINDER_FLUX)
def test_critical_flux_decreases_with_radius(c):
    values = [v for _, v in ganapol.CRITICAL_CYLINDER_FLUX[c]]
    assert values[0] == 1.0
    assert np.all(np.diff(values) < 0.)


@pytest.mark.parametrize('key', ganapol.SURROGATE_FLUX,
                         ids=lambda k: f'R{k[0]:g}-c{k[1]:g}')
def test_surrogate_flux_increases_with_radius(key):
    values = [F for _, F in ganapol.SURROGATE_FLUX[key]]
    assert len(values) == 11
    assert np.all(np.diff(values) > 0.)
    assert all(0. < F < 1. for F in values)


def test_ganapol_problem_names_unique():
    names = [p.name for p in ganapol.PROBLEMS]
    assert len(names) == len(set(names))
