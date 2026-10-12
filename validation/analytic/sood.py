"""Analytical benchmark test set for criticality code verification.

The 75 problems, cross sections, and reference results are taken from

    A. Sood, R. A. Forster, and D. K. Parsons, "Analytical Benchmark Test Set
    for Criticality Code Verification," Los Alamos National Laboratory report
    LA-13511 (1999). Summarized in LA-UR-99-2567.

All finite problems are exactly critical (keff = 1) at the quoted critical
dimensions. Infinite medium problems have the quoted k-infinity. The reference
documents number energy groups such that the highest group number is the
fastest group; the cross sections below are entered using that convention and
reordered fastest-first by :meth:`Material.from_reference`.

"""

from dataclasses import dataclass, field

from .xs import Material


# One-group isotropic cross sections (Tables 5, 12, 16, and 20)
def _one_group_fissile(name, nu, fission, capture, scatter, total):
    return Material.one_group(name, total, scatter, capture, fission, nu)


MATERIALS = {}


def _add(material):
    MATERIALS[material.name] = material


_add(_one_group_fissile('PUa', 3.24, 0.081600, 0.019584, 0.225216, 0.32640))
_add(_one_group_fissile('PUb', 2.84, 0.081600, 0.019584, 0.225216, 0.32640))
_add(_one_group_fissile('Ua', 2.70, 0.065280, 0.013056, 0.248064, 0.32640))
_add(_one_group_fissile('Ub', 2.797101, 0.065280, 0.013056, 0.248064, 0.32640))
_add(_one_group_fissile('Uc', 2.707308, 0.065280, 0.013056, 0.248064, 0.32640))
_add(_one_group_fissile('Ud', 2.679198, 0.065280, 0.013056, 0.248064, 0.32640))
_add(Material.one_group('H2O', 0.32640, 0.293760, 0.032640))
_add(_one_group_fissile('UD2O', 1.70, 0.054628, 0.027314, 0.464338, 0.54628))
_add(Material.one_group('H2O-UD2O', 0.54628, 0.491652, 0.054628))
_add(_one_group_fissile('Ue', 2.50, 0.06922744, 0.01013756, 0.328042,
                        0.407407))
_add(Material.one_group('Fe', 0.23256, 0.23209488, 0.00046512))
_add(Material.one_group('Na', 0.086368032, 0.086368032))

# One-group anisotropic cross sections (Tables 24, 26, and 28)
_add(_one_group_fissile('PUa-P1', 2.5, 0.266667, 0.0, [0.733333, 0.20], 1.0))
_add(_one_group_fissile('PUa-P2', 2.5, 0.266667, 0.0,
                        [0.733333, 0.20, 0.075], 1.0))
_add(_one_group_fissile('PUb-P1', 2.5, 0.266667, 0.0,
                        [0.733333, 0.333333], 1.0))
_add(_one_group_fissile('PUb-P2', 2.5, 0.266667, 0.0,
                        [0.733333, 0.333333, 0.125], 1.0))
_add(_one_group_fissile('Ua-P1', 2.70, 0.065280, 0.013056,
                        [0.248064, 0.042432], 0.32640))
_add(_one_group_fissile('Ub-P1', 2.70, 0.065280, 0.013056,
                        [0.248064, 0.212160], 0.32640))
_add(_one_group_fissile('UD2Oa-P1', 1.808381, 0.054628, 0.027314,
                        [0.464338, 0.056312624], 0.54628))
_add(_one_group_fissile('UD2Ob-P1', 1.841086, 0.054628, 0.027314,
                        [0.464338, 0.112982569], 0.54628))
_add(_one_group_fissile('UD2Oc-P1', 1.6964, 0.054628, 0.027314,
                        [0.464338, -0.27850447], 0.54628))


# Two-group cross sections (Tables 30-31, 33-34, 36-37, 39-40, 43-44, 46-47,
# 49-50, 52-53, and 56-57). Group 2 is the fast group.
def _two_group(name, fast, slow, s22, s21, s11, s12=0.):
    """Two-group set where s21 is scattering from group 2 to group 1."""
    return Material.from_reference(
        name, {2: fast, 1: slow},
        {(2, 2): s22, (2, 1): s21, (1, 1): s11, (1, 2): s12})


_add(_two_group(
    'PU-2G',
    dict(nu=3.10, fission=0.0936, capture=0.00480, total=0.2208, chi=0.575),
    dict(nu=2.93, fission=0.08544, capture=0.0144, total=0.3360, chi=0.425),
    s22=0.0792, s21=0.0432, s11=0.23616))
_add(_two_group(
    'U-2G',
    dict(nu=2.70, fission=0.06192, capture=0.00384, total=0.2160, chi=0.575),
    dict(nu=2.50, fission=0.06912, capture=0.01344, total=0.3456, chi=0.425),
    s22=0.078240, s21=0.0720, s11=0.26304))
_add(_two_group(
    'UAL-2G',
    dict(capture=0.000222, total=0.26817, chi=1.0),
    dict(nu=2.83, fission=0.06070636042, capture=0.00314363958,
         total=1.27698),
    s22=0.247516, s21=0.020432, s11=1.21313))
_add(_two_group(
    'URRa-2G',
    dict(nu=2.50, fission=0.0010484, capture=0.0010046, total=0.65696,
         chi=1.0),
    dict(nu=2.50, fission=0.050632, capture=0.025788, total=2.52025),
    s22=0.62568, s21=0.029227, s11=2.44383))
_add(_two_group(
    'URRa-2G-P1',
    dict(nu=2.50, fission=0.0010484, capture=0.0010046, total=0.65696,
         chi=1.0),
    dict(nu=2.50, fission=0.050632, capture=0.025788, total=2.52025),
    s22=[0.62568, 0.27459], s21=[0.029227, 0.0075737],
    s11=[2.44383, 0.83318]))
_add(_two_group(
    'URRb-2G',
    dict(nu=2.50, fission=0.000836, capture=0.001104, total=0.88721, chi=1.0),
    dict(nu=2.50, fission=0.029564, capture=0.024069, total=2.9727),
    s22=0.83892, s21=0.04635, s11=2.9183, s12=0.000767))
_add(_two_group(
    'URRc-2G',
    dict(nu=2.50, fission=0.001648, capture=0.001472, total=0.88655, chi=1.0),
    dict(nu=2.50, fission=0.057296, capture=0.029244, total=2.9628),
    s22=0.83807, s21=0.04536, s11=2.8751, s12=0.00116))
_add(_two_group(
    'H2Oa-2G',
    dict(capture=0.00074, total=0.88798),
    dict(capture=0.018564, total=2.9865),
    s22=0.83975, s21=0.04749, s11=2.9676, s12=0.000336))
_add(_two_group(
    'URRd-2G',
    dict(nu=1.004, fission=0.61475, capture=0.0019662, total=0.650917,
         chi=1.0),
    dict(nu=2.50, fission=0.045704, capture=0.023496, total=2.13800),
    s22=0.0, s21=0.0342008, s11=2.06880))
_add(_two_group(
    'H2Ob-2G',
    dict(capture=8.480293e-6, total=0.1106832906),
    dict(capture=0.00016, total=0.36355),
    s22=0.1096742149, s21=0.001000595707, s11=0.36339))
_add(_two_group(
    'H2Oc-2G',
    dict(capture=4.97229e-4, total=1.331518007),
    dict(capture=0.0188, total=4.37350),
    s22=1.226381244, s21=0.1046395340, s11=4.35470))
_add(_two_group(
    'UD2O-2G',
    dict(nu=2.50, fission=0.002817, capture=0.008708, total=0.33588, chi=1.0),
    dict(nu=2.50, fission=0.097, capture=0.02518, total=0.54628),
    s22=0.31980, s21=0.004555, s11=0.42410))
_add(_two_group(
    'UD2O-2G-P1',
    dict(nu=2.50, fission=0.002817, capture=0.008708, total=0.33588, chi=1.0),
    dict(nu=2.50, fission=0.097, capture=0.02518, total=0.54628),
    s22=[0.31980, 0.06694], s21=[0.004555, -0.0003972],
    s11=[0.42410, 0.05439]))

# Three-group cross sections (Tables 59-61). Group 3 is the fast group.
_add(Material.from_reference(
    'URR-3G',
    {3: dict(nu=3.0, fission=0.006, capture=0.006, total=0.240, chi=0.96),
     2: dict(nu=2.50, fission=0.060, capture=0.040, total=0.975, chi=0.04),
     1: dict(nu=2.0, fission=0.90, capture=0.20, total=3.10)},
    {(3, 3): 0.024, (3, 2): 0.171, (3, 1): 0.033,
     (2, 2): 0.60, (2, 1): 0.275, (1, 1): 2.0}))

# Six-group cross sections (Tables 62-67). Groups 6 and 1, 5 and 2, and 4 and
# 3 are equivalent; groups 1-3 scatter "up" into each other.
_add(Material.from_reference(
    'URR-6G',
    {6: dict(nu=3.0, fission=0.006, capture=0.006, total=0.240, chi=0.48),
     5: dict(nu=2.50, fission=0.060, capture=0.040, total=0.975, chi=0.02),
     4: dict(nu=2.0, fission=0.90, capture=0.20, total=3.10),
     3: dict(nu=2.0, fission=0.90, capture=0.20, total=3.10),
     2: dict(nu=2.50, fission=0.060, capture=0.040, total=0.975, chi=0.02),
     1: dict(nu=3.0, fission=0.006, capture=0.006, total=0.240, chi=0.48)},
    {(6, 6): 0.024, (6, 5): 0.171, (6, 4): 0.033,
     (5, 5): 0.60, (5, 4): 0.275, (4, 4): 2.0,
     (1, 1): 0.024, (1, 2): 0.171, (1, 3): 0.033,
     (2, 2): 0.60, (2, 3): 0.275, (3, 3): 2.0}))


@dataclass
class Problem:
    """Analytic criticality benchmark problem.

    Parameters
    ----------
    number : int
        Problem number in LA-13511
    name : str
        Problem identifier, e.g. 'PUa-1-0-SL'
    geometry : {'infinite', 'slab', 'cylinder', 'sphere'}
        Geometry type
    regions : list of tuple
        For a slab, (material, thickness [cm]) pairs from left to right. For a
        cylinder or sphere, (material, outer radius [cm]) pairs from the
        center outward. For an infinite medium, a single (material, None)
        pair.
    keff : float
        Reference eigenvalue (1.0 for finite critical problems, k-infinity for
        infinite media).
    boundary : tuple of str
        Boundary conditions of a slab on the left and right
    mfp : list of float
        Region dimensions in mean free paths as quoted in the reference, used
        to check the transcription of the dimensions in cm
    flux_ratios : dict
        For multigroup infinite medium problems, mapping of (group, group) in
        reference numbering to the ratio of their scalar fluxes
    flux_shape : list of tuple
        (position [cm], reference group, normalized scalar flux) giving
        spatial flux values. Positions are measured from the left boundary
        for slabs and from the center for cylinders and spheres.
    negative_pdf : bool
        Whether the Legendre scattering expansion becomes negative over part
        of the cosine interval. Such distributions cannot be sampled exactly
        by a Monte Carlo code.

    """
    number: int
    name: str
    geometry: str
    regions: list
    keff: float = 1.0
    boundary: tuple = ('vacuum', 'vacuum')
    mfp: list = field(default_factory=list)
    flux_ratios: dict = field(default_factory=dict)
    flux_shape: list = field(default_factory=list)
    negative_pdf: bool = False

    @property
    def materials(self):
        names = dict.fromkeys(m for m, _ in self.regions)
        return [MATERIALS[m] for m in names]

    @property
    def num_groups(self):
        return self.materials[0].num_groups


def _infinite(number, name, material, keff, flux_ratios=None,
              negative_pdf=False):
    return Problem(number, name, 'infinite', [(material, None)], keff,
                   flux_ratios=flux_ratios or {}, negative_pdf=negative_pdf)


def _bare_slab(number, name, material, rc, mfp, shape=(), **kwargs):
    """Bare slab of half-thickness rc. Flux values are at r/rc positions
    measured from the center of the slab."""
    flux = [(rc + f*rc, g, v) for f, g, v in shape]
    return Problem(number, name, 'slab', [(material, 2*rc)], mfp=[mfp],
                   flux_shape=flux, **kwargs)


def _reflected_slab(number, name, fuel, rc, refl, t, mfp):
    """Slab of half-thickness rc with a reflector of thickness t on both
    sides."""
    return Problem(number, name, 'slab', [(refl, t), (fuel, 2*rc), (refl, t)],
                   mfp=mfp)


def _round(number, name, geometry, regions, mfp, shape=(), **kwargs):
    """Cylinder or sphere. regions are (material, thickness) pairs."""
    radii = []
    r = 0.
    for material, thickness in regions:
        r += thickness
        radii.append((material, r))
    rc = regions[0][1]
    flux = [(f*rc, g, v) for f, g, v in shape]
    return Problem(number, name, geometry, radii, mfp=mfp, flux_shape=flux,
                   **kwargs)


def _profile(fractions, values, group=1):
    return [(f, group, v) for f, v in zip(fractions, values)]


_QUARTERS = (0.25, 0.5, 0.75, 1.0)

# Thickness of a reflector treated as infinite [cm]. This is 100 fast-group
# mean free paths of H2O (a), far more than needed for convergence of keff.
_INFINITE_H2O_A = 100/0.88798

PROBLEMS = [
    # One-group isotropic Pu-239 (Tables 5-11)
    _infinite(1, 'PUa-1-0-IN', 'PUa', 2.612903),
    _bare_slab(2, 'PUa-1-0-SL', 'PUa', 1.853722, 0.605055),
    Problem(3, 'PUa-H2O(1)-1-0-SL', 'slab',
            [('PUa', 2*1.478401), ('H2O', 3.063725)], mfp=[0.48255, 1]),
    _reflected_slab(4, 'PUa-H2O(0.5)-1-0-SL', 'PUa', 1.317831, 'H2O',
                    1.531863, [0.43014, 0.5]),
    _infinite(5, 'PUb-1-0-IN', 'PUb', 2.290323),
    _bare_slab(6, 'PUb-1-0-SL', 'PUb', 2.256751, 0.73660355, _profile(
        _QUARTERS, (0.9701734, 0.8810540, 0.7318131, 0.4902592))),
    _round(7, 'PUb-1-0-CY', 'cylinder', [('PUb', 4.279960)], [1.396979],
           _profile((0.5, 1.0), (0.8093, 0.2926))),
    _round(8, 'PUb-1-0-SP', 'sphere', [('PUb', 6.082547)], [1.9853434324],
           _profile(_QUARTERS,
                    (0.93538006, 0.75575352, 0.49884364, 0.19222603))),
    _round(9, 'PUb-H2O(1)-1-0-CY', 'cylinder',
           [('PUb', 3.397610), ('H2O', 3.063725)], [1.10898, 1]),
    _round(10, 'PUb-H2O(10)-1-0-CY', 'cylinder',
           [('PUb', 3.077574), ('H2O', 30.637255)], [1.00452, 10]),

    # One-group isotropic U-235 (Tables 12-15)
    _infinite(11, 'Ua-1-0-IN', 'Ua', 2.25),
    _bare_slab(12, 'Ua-1-0-SL', 'Ua', 2.872934, 0.93772556, _profile(
        _QUARTERS, (0.9669506, 0.8686259, 0.7055218, 0.4461912))),
    _round(13, 'Ua-1-0-CY', 'cylinder', [('Ua', 5.284935)], [1.72500292]),
    _round(14, 'Ua-1-0-SP', 'sphere', [('Ua', 7.428998)], [2.4248249802],
           _profile(_QUARTERS,
                    (0.93244907, 0.74553332, 0.48095413, 0.17177706))),
    _infinite(15, 'Ub-1-0-IN', 'Ub', 2.330917),
    _round(16, 'Ub-H2O(1)-1-0-SP', 'sphere',
           [('Ub', 6.12745), ('H2O', 3.063725)], [2, 1]),
    _infinite(17, 'Uc-1-0-IN', 'Uc', 2.256083),
    _round(18, 'Uc-H2O(2)-1-0-SP', 'sphere',
           [('Uc', 6.12745), ('H2O', 6.12745)], [2, 2]),
    _infinite(19, 'Ud-1-0-IN', 'Ud', 2.232667),
    _round(20, 'Ud-H2O(3)-1-0-SP', 'sphere',
           [('Ud', 6.12745), ('H2O', 9.191176)], [2, 3]),

    # One-group isotropic U-D2O reactor (Tables 16-19)
    _infinite(21, 'UD2O-1-0-IN', 'UD2O', 1.133333),
    _bare_slab(22, 'UD2O-1-0-SL', 'UD2O', 10.371065, 5.6655054562, _profile(
        _QUARTERS, (0.93945236, 0.76504084, 0.49690627, 0.13893858))),
    _round(23, 'UD2O-1-0-CY', 'cylinder', [('UD2O', 16.554249)], [9.043255]),
    _round(24, 'UD2O-1-0-SP', 'sphere', [('UD2O', 22.017156)],
           [12.0275320980], _profile(
               _QUARTERS, (0.91063756, 0.67099621, 0.35561622, 0.04678614))),
    _reflected_slab(25, 'UD2O-H2O(1)-1-0-SL', 'UD2O', 9.214139, 'H2O-UD2O',
                    1.830563, [5.0335, 1]),
    _reflected_slab(26, 'UD2O-H2O(10)-1-0-SL', 'UD2O', 8.428096, 'H2O-UD2O',
                    18.30563, [4.6041, 10]),
    _round(27, 'UD2O-H2O(1)-1-0-CY', 'cylinder',
           [('UD2O', 15.396916), ('H2O-UD2O', 1.830563)], [8.411027, 1]),
    _round(28, 'UD2O-H2O(10)-1-0-CY', 'cylinder',
           [('UD2O', 14.606658), ('H2O-UD2O', 18.30563)], [7.979325, 10]),

    # One-group isotropic U-235 reactor with Fe and Na (Tables 20-23). The
    # fluxes are at the material interfaces and the right boundary,
    # normalized to the flux at the left boundary.
    _infinite(29, 'Ue-1-0-IN', 'Ue', 2.1806667),
    Problem(30, 'Ue-Fe-Na-1-0-SL', 'slab',
            [('Fe', 0.317337461), ('Ue', 5.119720083), ('Fe', 0.317337461),
             ('Na', 2.002771002)],
            mfp=[0.0738, 2.0858098, 0.0738, 0.173],
            flux_shape=[(0.0, 1, 1.0), (0.317337461, 1, 1.229538),
                        (5.437057544, 1, 1.49712), (5.754395005, 1, 1.324899),
                        (7.757166007, 1, 0.912273)]),

    # One-group anisotropic Pu-239 (Tables 24-25); Sigma_t = 1 so cm = mfp
    _infinite(31, 'PU-1-1-IN', 'PUa-P1', 2.5),
    _bare_slab(32, 'PUa-1-1-SL', 'PUa-P1', 0.77032, 0.77032),
    _bare_slab(33, 'PUa-1-2-SL', 'PUa-P2', 0.76378, 0.76378),
    _bare_slab(34, 'PUb-1-1-SL', 'PUb-P1', 0.79606, 0.79606,
               negative_pdf=True),
    _bare_slab(35, 'PUb-1-2-SL', 'PUb-P2', 0.78396, 0.78396),

    # One-group anisotropic U-235 (Tables 26-27)
    _round(36, 'Ua-1-1-CY', 'cylinder', [('Ua-P1', 5.514296811)],
           [1.799866479]),
    _round(37, 'Ub-1-1-CY', 'cylinder', [('Ub-P1', 6.940205668)],
           [2.265283130], negative_pdf=True),

    # One-group anisotropic U-D2O (Tables 28-29)
    _infinite(38, 'UD2Oa-1-1-IN', 'UD2Oa-P1', 1.205587),
    _round(39, 'UD2Oa-1-1-SP', 'sphere', [('UD2Oa-P1', 18.30563081)], [10]),
    _infinite(40, 'UD2Ob-1-1-IN', 'UD2Ob-P1', 1.227391),
    _round(41, 'UD2Ob-1-1-SP', 'sphere', [('UD2Ob-P1', 18.30563081)], [10]),
    _infinite(42, 'UD2Oc-1-1-IN', 'UD2Oc-P1', 1.130933, negative_pdf=True),
    _round(43, 'UD2Oc-1-1-SP', 'sphere', [('UD2Oc-P1', 18.30563081)], [10],
           negative_pdf=True),

    # Two-group isotropic Pu-239 (Tables 30-32)
    _infinite(44, 'PU-2-0-IN', 'PU-2G', 2.683767, {(2, 1): 0.675229}),
    _bare_slab(45, 'PU-2-0-SL', 'PU-2G', 1.795602, 0.396469),
    _round(46, 'PU-2-0-SP', 'sphere', [('PU-2G', 5.231567)], [1.15513]),

    # Two-group isotropic U-235 (Tables 33-35)
    _infinite(47, 'U-2-0-IN', 'U-2G', 2.216349, {(2, 1): 0.474967}),
    _bare_slab(48, 'U-2-0-SL', 'U-2G', 3.006375, 0.649377),
    _round(49, 'U-2-0-SP', 'sphere', [('U-2G', 7.909444)], [1.70844]),

    # Two-group isotropic U-Al-water assembly (Tables 36-38)
    _infinite(50, 'UAL-2-0-IN', 'UAL-2G', 2.661745, {(2, 1): 3.1250}),
    _bare_slab(51, 'UAL-2-0-SL', 'UAL-2G', 7.830630, 2.09994),
    _round(52, 'UAL-2-0-SP', 'sphere', [('UAL-2G', 17.66738)], [4.73786]),

    # Two-group isotropic research reactor (Tables 39-48). The flux profile
    # is normalized to the fast flux at the center.
    _infinite(53, 'URRa-2-0-IN', 'URRa-2G', 1.631452, {(2, 1): 2.614706}),
    _bare_slab(54, 'URRa-2-0-SL', 'URRa-2G', 7.566853, 4.97112,
               [(0.0, 2, 1.0)] +
               _profile((0.241394, 0.502905, 0.744300, 1.0),
                        (0.943363, 0.761973, 0.504012, 0.147598), 2) +
               _profile((0.241394, 0.502905, 0.744300, 1.0),
                        (0.340124, 0.273056, 0.173845, 0.0212324), 1)),
    _round(55, 'URRa-2-0-SP', 'sphere', [('URRa-2G', 16.049836)], [10.5441]),
    _infinite(56, 'URRb-2-0-IN', 'URRb-2G', 1.365821, {(2, 1): 1.173679}),
    _infinite(57, 'URRc-2-0-IN', 'URRc-2G', 1.633380, {(2, 1): 1.933422}),
    _reflected_slab(58, 'URRb-H2Oa(1)-2-0-SL', 'URRb-2G', 6.696802,
                    'H2Oa-2G', 1.126152, [5.94147, 1]),
    _reflected_slab(59, 'URRb-H2Oa(5)-2-0-SL', 'URRb-2G', 4.863392,
                    'H2Oa-2G', 5.630757, [4.31485, 5]),
    _reflected_slab(60, 'URRb-H2Oa(IN)-2-0-SL', 'URRb-2G', 4.686230,
                    'H2Oa-2G', _INFINITE_H2O_A, [4.15767]),
    _reflected_slab(61, 'URRc-H2Oa(IN)-2-0-SL', 'URRc-2G', 2.461903,
                    'H2Oa-2G', _INFINITE_H2O_A, [2.1826]),
    _infinite(62, 'URRd-2-0-IN', 'URRd-2G', 1.034970, {(2, 1): 2.023344}),

    # Infinite slab lattice cells: half of the fuel and half of the moderator
    # between two reflective boundaries.
    Problem(63, 'URRd-H2Ob(1)-2-0-ISLC', 'slab',
            [('URRd-2G', 0.0329074), ('H2Ob-2G', 9.034787)],
            boundary=('reflective', 'reflective'), mfp=[0.02142, 1]),
    Problem(64, 'URRd-H2Ob(10)-2-0-ISLC', 'slab',
            [('URRd-2G', 0.460135), ('H2Ob-2G', 90.347875)],
            boundary=('reflective', 'reflective'), mfp=[0.29951, 10]),
    Problem(65, 'URRd-H2Oc(1)-2-0-ISLC', 'slab',
            [('URRd-2G', 0.341011), ('H2Oc-2G', 0.751023)],
            boundary=('reflective', 'reflective'), mfp=[0.22197, 1]),
    Problem(66, 'URRd-H2Oc(10)-2-0-ISLC', 'slab',
            [('URRd-2G', 2.719087), ('H2Oc-2G', 7.510225)],
            boundary=('reflective', 'reflective'), mfp=[1.7699, 10]),

    # Two-group isotropic U-D2O reactor (Tables 49-51)
    _infinite(67, 'UD2O-2-0-IN', 'UD2O-2G', 1.000196, {(2, 1): 26.823271}),
    _bare_slab(68, 'UD2O-2-0-SL', 'UD2O-2G', 846.632726, 284.367),
    _round(69, 'UD2O-2-0-SP', 'sphere', [('UD2O-2G', 1695.337621)],
           [569.430]),

    # Two-group linearly anisotropic research reactor (Tables 52-55)
    _infinite(70, 'URRa-2-1-IN', 'URRa-2G-P1', 1.631452, {(2, 1): 2.614706}),
    _bare_slab(71, 'URRa-2-1-SL', 'URRa-2G-P1', 9.491600, 6.2356,
               [(0.0, 2, 1.0)] +
               _profile((0.2, 0.5, 0.8, 1.0),
                        (0.963873, 0.781389, 0.472787, 0.189578), 2) +
               _profile((0.2, 0.5, 0.8, 1.0),
                        (0.349006, 0.280870, 0.157376, 0.0277639), 1),
               negative_pdf=True),

    # Two-group linearly anisotropic U-D2O reactor (Tables 56-58)
    _infinite(72, 'UD2O-2-1-IN', 'UD2O-2G-P1', 1.000196, {(2, 1): 26.823271}),
    _bare_slab(73, 'UD2O-2-1-SL', 'UD2O-2G-P1', 1000.506133, 336.05),

    # Three- and six-group infinite media (Tables 59-67)
    _infinite(74, 'URR-3-0-IN', 'URR-3G', 1.60,
              {(2, 3): 0.480, (1, 2): 0.3125, (1, 3): 0.150}),
    _infinite(75, 'URR-6-0-IN', 'URR-6G', 1.60,
              {(5, 6): 0.480, (2, 1): 0.480, (4, 5): 0.3125, (3, 2): 0.3125,
               (4, 6): 0.150, (3, 1): 0.150}),
]

PROBLEMS_BY_NAME = {p.name: p for p in PROBLEMS}
PROBLEMS_BY_NUMBER = {p.number: p for p in PROBLEMS}
