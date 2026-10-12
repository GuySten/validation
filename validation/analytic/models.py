"""Build and evaluate OpenMC multigroup models of the analytic benchmarks."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.special import j0

import openmc
import openmc.mgxs

from .xs import Material


# Linear extrapolation distance used only to shape the initial fission source
_EXTRAPOLATION = 0.7104

# Width of the y and z extent of slab problems, which have reflective
# boundaries in those directions [cm]
_TRANSVERSE = 1.0


@dataclass
class Comparison:
    """Comparison of a calculated quantity with its reference value.

    Parameters
    ----------
    quantity : str
        Description of the quantity
    reference : float
        Reference value
    mean : float
        Calculated value
    std_dev : float
        Standard deviation of the calculated value
    rtol : float
        Relative tolerance added to the statistical tolerance to account for
        discretization of the tallies
    atol : float
        Absolute tolerance added to the statistical tolerance to account for
        rounding of the reference value

    """
    quantity: str
    reference: float
    mean: float
    std_dev: float
    rtol: float = 0.0
    atol: float = 0.0

    def __post_init__(self):
        # Some estimates have no variance at all (e.g. k-infinity of a
        # one-group infinite medium), for which OpenMC may report NaN
        if np.isnan(self.std_dev):
            self.std_dev = 0.0

    @property
    def difference(self):
        return self.mean - self.reference

    @property
    def sigma(self):
        """Difference in units of the standard deviation."""
        if self.std_dev == 0.:
            return np.inf if self.difference != 0. else 0.
        return self.difference / self.std_dev

    def passed(self, n_sigma=4.0):
        tolerance = (n_sigma*self.std_dev + self.rtol*abs(self.reference) +
                     self.atol)
        return abs(self.difference) <= tolerance


@dataclass
class Results:
    """Results of running a benchmark problem."""
    name: str
    comparisons: list = field(default_factory=list)

    def passed(self, n_sigma=4.0):
        return all(c.passed(n_sigma) for c in self.comparisons)


def energy_groups(num_groups):
    """Arbitrary group structure with the requested number of groups."""
    return openmc.mgxs.EnergyGroups(np.logspace(-5, 7, num_groups + 1))


def export_library(materials, path):
    """Write the cross sections to an OpenMC MGXS library.

    Parameters
    ----------
    materials : iterable of Material
        Cross section sets with the same number of groups
    path : str or pathlib.Path
        HDF5 file to write

    """
    materials = list(materials)
    groups = energy_groups(materials[0].num_groups)
    library = openmc.MGXSLibrary(groups)
    for mat in materials:
        xsdata = openmc.XSdata(mat.name, groups)
        xsdata.order = mat.order
        xsdata.set_total(mat.total)
        xsdata.set_absorption(mat.absorption)
        # OpenMC expects [group in, group out, Legendre order]
        xsdata.set_scatter_matrix(np.transpose(mat.scatter, (1, 2, 0)))
        if mat.fissionable:
            xsdata.set_fission(mat.fission)
            xsdata.set_nu_fission(mat.nu_fission)
            xsdata.set_chi(mat.chi)
        library.add_xsdata(xsdata)
    library.export_to_hdf5(str(path))


def _materials(materials, directory):
    """Create OpenMC materials referencing a library written to directory."""
    path = Path(directory).resolve() / 'mgxs.h5'
    export_library(materials, path)
    omc_materials = {}
    for mat in materials:
        omc = openmc.Material(name=mat.name)
        omc.set_density('macro', 1.0)
        omc.add_macroscopic(mat.name)
        omc_materials[mat.name] = omc
    collection = openmc.Materials(omc_materials.values())
    collection.cross_sections = path
    return omc_materials, collection


def _settings(anisotropic):
    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.output = {'summary': False}
    if anisotropic:
        # Sample the Legendre expansion directly rather than converting it
        # to a tabular distribution
        settings.tabular_legendre = {'enable': False}
    return settings


def _source_energy(num_groups):
    """Discrete source energy in the fastest group."""
    edges = energy_groups(num_groups).group_edges
    return openmc.stats.Discrete([np.sqrt(edges[-1]*edges[-2])], [1.0])


def _histogram(edges, pdf):
    """Histogram distribution from bin edges and the pdf at bin centers."""
    p = np.append(pdf, 0.)
    return openmc.stats.Tabular(edges, p, 'histogram')


def _slab_planes(x_edges, boundary):
    planes = [openmc.XPlane(x) for x in x_edges]
    planes[0].boundary_type = boundary[0]
    planes[-1].boundary_type = boundary[1]
    h = _TRANSVERSE/2
    y0 = openmc.YPlane(-h, boundary_type='reflective')
    y1 = openmc.YPlane(h, boundary_type='reflective')
    z0 = openmc.ZPlane(-h, boundary_type='reflective')
    z1 = openmc.ZPlane(h, boundary_type='reflective')
    return planes, +y0 & -y1 & +z0 & -z1


def _z_slice(height=_TRANSVERSE):
    h = height/2
    z0 = openmc.ZPlane(-h, boundary_type='reflective')
    z1 = openmc.ZPlane(h, boundary_type='reflective')
    return +z0 & -z1


# -----------------------------------------------------------------------------
# LA-13511 criticality problems
# -----------------------------------------------------------------------------

def sood_model(problem, directory, particles=10000, batches=150, inactive=50):
    """Build an OpenMC model of a LA-13511 criticality problem.

    Parameters
    ----------
    problem : validation.analytic.sood.Problem
        Problem to model
    directory : str or pathlib.Path
        Directory where the cross section library is written
    particles : int
        Number of particles per batch
    batches : int
        Total number of batches
    inactive : int
        Number of inactive batches

    Returns
    -------
    openmc.Model
        Model of the problem

    """
    openmc.reset_auto_ids()
    mats, collection = _materials(problem.materials, directory)
    num_groups = problem.num_groups
    anisotropic = any(m.order > 0 for m in problem.materials)

    settings = _settings(anisotropic)
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    geometry, space, tallies = _SOOD_BUILDERS[problem.geometry](
        problem, mats, num_groups)

    settings.source = openmc.IndependentSource(
        space=space, energy=_source_energy(num_groups),
        constraints={'fissionable': True})

    return openmc.Model(geometry, collection, settings,
                        openmc.Tallies(tallies))


def _group_filter(num_groups):
    return openmc.EnergyFilter(energy_groups(num_groups).group_edges)


def _sood_infinite(problem, mats, num_groups):
    material = mats[problem.regions[0][0]]
    box = openmc.model.RectangularParallelepiped(
        -1., 1., -1., 1., -1., 1., boundary_type='reflective')
    cell = openmc.Cell(fill=material, region=-box)
    geometry = openmc.Geometry([cell])
    space = openmc.stats.Box((-1., -1., -1.), (1., 1., 1.))

    tallies = []
    if num_groups > 1:
        tally = openmc.Tally(name='group flux')
        tally.filters = [openmc.CellFilter([cell]),
                         _group_filter(num_groups)]
        tally.scores = ['flux']
        tallies.append(tally)
    return geometry, space, tallies


def _flux_bins(positions, lower, upper, width, lower_width=None):
    """Mesh grid with a bin of the given width at each position.

    A position on the lower or upper boundary gets a one-sided bin. The
    width of a one-sided bin at the lower boundary can be set separately,
    which is useful at the center of a cylinder or sphere where the flux is
    flat and a narrow bin has a tiny volume.

    """
    if lower_width is None:
        lower_width = width/2
    edges = set()
    for p in positions:
        if p - width/2 <= lower:
            edges.update((lower, lower + lower_width))
        elif p + width/2 >= upper:
            edges.update((upper - width/2, upper))
        else:
            edges.update((p - width/2, p + width/2))
    return np.array(sorted(edges))


def _bin_index(grid, position):
    """Index of the bin built around a position by _flux_bins."""
    i = np.searchsorted(grid, position, side='right') - 1
    return int(min(max(i, 0), grid.size - 2))


def _point_flux_tallies(positions, boundaries, mesh_factory, lower, upper,
                        width, num_groups=1, lower_width=None):
    """Tallies giving the scalar flux at a set of positions.

    Positions on a vacuum boundary use the surface-crossing flux estimator.
    All other positions use a track-length estimate in a narrow mesh bin
    built around the position.

    Parameters
    ----------
    positions : iterable of float
        Positions (x for slabs, r for cylinders and spheres)
    boundaries : dict
        Mapping of a name to a (position, surface, area) tuple for each
        vacuum boundary
    mesh_factory : callable
        Function that returns a mesh given the grid of bin edges
    lower, upper : float
        Extent of the problem along the mesh coordinate
    width : float
        Width of the mesh bins
    num_groups : int
        Number of energy groups
    lower_width : float, optional
        Width of a one-sided bin at the lower bound

    """
    tallies = []
    interior = sorted({p for p in positions
                       if _boundary_name(boundaries, p) is None})
    if interior:
        grid = _flux_bins(interior, lower, upper, width, lower_width)
        tally = openmc.Tally(name='flux')
        tally.filters = [openmc.MeshFilter(mesh_factory(grid))]
        tallies.append(tally)
    for name, (position, surface, _) in boundaries.items():
        if any(np.isclose(p, position) for p in positions):
            tally = openmc.Tally(name=f'surface flux {name}')
            tally.filters = [openmc.SurfaceFilter([surface])]
            tallies.append(tally)
    for tally in tallies:
        if num_groups > 1:
            tally.filters.append(_group_filter(num_groups))
        tally.scores = ['flux']
    return tallies


def _boundary_name(boundaries, position):
    for name, (p, _, _) in boundaries.items():
        if np.isclose(p, position):
            return name
    return None


def _point_fluxes(sp, positions, boundaries):
    """Read the scalar flux at each position from a statepoint.

    Returns a dictionary mapping each position to arrays of the mean and
    standard deviation of the group fluxes, ordered fastest group first.

    """
    fluxes = {}
    mesh_tally = None
    for p in positions:
        name = _boundary_name(boundaries, p)
        if name is not None:
            tally = sp.get_tally(name=f'surface flux {name}')
            area = boundaries[name][2]
            mean = tally.mean.ravel()/area
            std = tally.std_dev.ravel()/area
        else:
            if mesh_tally is None:
                mesh_tally = sp.get_tally(name='flux')
                mesh = mesh_tally.find_filter(openmc.MeshFilter).mesh
                volumes = _mesh_volumes(mesh)
                grid = _mesh_grid(mesh)
                n = volumes.size
                all_mean = mesh_tally.mean.reshape(n, -1)/volumes[:, None]
                all_std = mesh_tally.std_dev.reshape(n, -1)/volumes[:, None]
            i = _bin_index(grid, p)
            mean, std = all_mean[i], all_std[i]
        # Energy filter bins are in increasing energy
        fluxes[p] = (np.asarray(mean)[::-1], np.asarray(std)[::-1])
    return fluxes


def _slab_mesh(grid):
    h = _TRANSVERSE/2
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = grid
    mesh.y_grid = [-h, h]
    mesh.z_grid = [-h, h]
    return mesh


def _cylinder_mesh(grid):
    return openmc.CylindricalMesh(
        r_grid=grid, z_grid=[-_TRANSVERSE/2, _TRANSVERSE/2])


def _sphere_mesh(grid):
    return openmc.SphericalMesh(r_grid=grid)


def _sood_boundaries(problem, surfaces=None):
    """Vacuum boundaries of a finite LA-13511 problem as a mapping of name to
    (position, surface, area)."""
    if surfaces is None:
        surfaces = {}
    if problem.geometry == 'slab':
        width = sum(t for _, t in problem.regions)
        sides = {'left': 0., 'right': width}
        return {name: (x, surfaces.get(name), _TRANSVERSE**2)
                for (name, x), bc in zip(sides.items(), problem.boundary)
                if bc == 'vacuum'}
    radius = problem.regions[-1][1]
    if problem.geometry == 'cylinder':
        area = 2*np.pi*radius*_TRANSVERSE
    else:
        area = 4*np.pi*radius**2
    return {'outer': (radius, surfaces.get('outer'), area)}


def _sood_slab(problem, mats, num_groups):
    thicknesses = [t for _, t in problem.regions]
    x_edges = np.concatenate(([0.], np.cumsum(thicknesses)))
    planes, transverse = _slab_planes(x_edges, problem.boundary)
    cells = []
    for i, (name, _) in enumerate(problem.regions):
        region = +planes[i] & -planes[i + 1] & transverse
        cells.append(openmc.Cell(fill=mats[name], region=region))
    geometry = openmc.Geometry(cells)

    # Initial source shaped like the fundamental mode of a bare slab
    width = x_edges[-1]
    vacuum = 'vacuum' in problem.boundary
    sigma = min(m.total.min() for m in problem.materials)
    extrapolated = width + 2*_EXTRAPOLATION/sigma
    xs = {m.name: m for m in problem.materials}
    edges, pdf = [], []
    for i, (name, _) in enumerate(problem.regions):
        x = np.linspace(x_edges[i], x_edges[i + 1], 101)
        centers = 0.5*(x[1:] + x[:-1])
        if not xs[name].fissionable:
            values = np.zeros_like(centers)
        elif vacuum:
            values = np.cos(np.pi*(centers - width/2)/extrapolated)
        else:
            values = np.ones_like(centers)
        edges.append(x[:-1])
        pdf.append(values)
    edges = np.append(np.concatenate(edges), x_edges[-1])
    h = _TRANSVERSE/2
    space = openmc.stats.CartesianIndependent(
        _histogram(edges, np.concatenate(pdf)),
        openmc.stats.Uniform(-h, h), openmc.stats.Uniform(-h, h))

    tallies = []
    if problem.flux_shape:
        boundaries = _sood_boundaries(
            problem, {'left': planes[0], 'right': planes[-1]})
        positions = [p for p, _, _ in problem.flux_shape]
        tallies = _point_flux_tallies(positions, boundaries, _slab_mesh, 0.,
                                      width, 0.004*width, num_groups)
    return geometry, space, tallies


def _sood_round(problem, mats, num_groups):
    cylinder = problem.geometry == 'cylinder'
    radii = [r for _, r in problem.regions]
    if cylinder:
        surfaces = [openmc.ZCylinder(r=r) for r in radii]
        axial = _z_slice()
    else:
        surfaces = [openmc.Sphere(r=r) for r in radii]
    surfaces[-1].boundary_type = 'vacuum'

    cells = []
    inner = None
    for (name, _), surface in zip(problem.regions, surfaces):
        region = -surface if inner is None else +inner & -surface
        if cylinder:
            region &= axial
        cells.append(openmc.Cell(fill=mats[name], region=region))
        inner = surface
    geometry = openmc.Geometry(cells)

    # Initial source shaped like the fundamental mode of a bare cylinder or
    # sphere, restricted to the (central) fissile region
    sigma = min(m.total.min() for m in problem.materials)
    extrapolated = radii[-1] + _EXTRAPOLATION/sigma
    r = np.linspace(0., radii[0], 201)
    centers = 0.5*(r[1:] + r[:-1])
    if cylinder:
        pdf = centers*j0(2.404826*centers/extrapolated)
        space = openmc.stats.CylindricalIndependent(
            _histogram(r, pdf), openmc.stats.Uniform(0., 2*np.pi),
            openmc.stats.Uniform(-_TRANSVERSE/2, _TRANSVERSE/2))
    else:
        pdf = centers*np.sin(np.pi*centers/extrapolated)
        space = openmc.stats.SphericalIndependent(
            _histogram(r, pdf), openmc.stats.Uniform(-1., 1.),
            openmc.stats.Uniform(0., 2*np.pi))

    tallies = []
    if problem.flux_shape:
        boundaries = _sood_boundaries(problem, {'outer': surfaces[-1]})
        positions = [p for p, _, _ in problem.flux_shape]
        mesh_factory = _cylinder_mesh if cylinder else _sphere_mesh
        tallies = _point_flux_tallies(positions, boundaries, mesh_factory,
                                      0., radii[-1], 0.004*radii[0],
                                      num_groups)
    return geometry, space, tallies


_SOOD_BUILDERS = {
    'infinite': _sood_infinite,
    'slab': _sood_slab,
    'cylinder': _sood_round,
    'sphere': _sood_round,
}


def _mesh_volumes(mesh):
    return np.asarray(mesh.volumes).ravel(order='F')


def _ratio(a, sa, b, sb):
    r = a/b
    return r, abs(r)*np.sqrt((sa/a)**2 + (sb/b)**2)


def sood_results(problem, statepoint):
    """Compare the results of a LA-13511 problem with the reference.

    Parameters
    ----------
    problem : validation.analytic.sood.Problem
        Problem that was run
    statepoint : str or pathlib.Path
        Statepoint file from the run

    Returns
    -------
    Results
        Comparisons of keff and any flux ratios with the reference values

    """
    results = Results(problem.name)
    G = problem.num_groups
    with openmc.StatePoint(statepoint) as sp:
        # The reference eigenvalues are quoted to at least six decimal places
        keff = sp.keff
        results.comparisons.append(Comparison(
            'keff', problem.keff, keff.nominal_value, keff.std_dev,
            atol=1e-6))

        if problem.geometry == 'infinite' and problem.flux_ratios:
            tally = sp.get_tally(name='group flux')
            # Energy filter bins are in increasing energy
            mean = tally.mean.ravel()[::-1]
            std = tally.std_dev.ravel()[::-1]
            for (a, b), ref in problem.flux_ratios.items():
                ia, ib = G - a, G - b
                r, s = _ratio(mean[ia], std[ia], mean[ib], std[ib])
                results.comparisons.append(Comparison(
                    f'flux ratio phi{a}/phi{b}', ref, r, s, atol=1e-6))

        if problem.flux_shape:
            positions = [p for p, _, _ in problem.flux_shape]
            fluxes = _point_fluxes(sp, positions, _sood_boundaries(problem))

            # Normalize to the first reference point
            p0, g0, v0 = problem.flux_shape[0]
            mean0, std0 = fluxes[p0]
            for p, g, v in problem.flux_shape[1:]:
                mean, std = fluxes[p]
                r, s = _ratio(mean[G - g], std[G - g], mean0[G - g0],
                              std0[G - g0])
                results.comparisons.append(Comparison(
                    f'flux(x={p:.6g}, g={g})/flux(x={p0:.6g}, g={g0})',
                    v/v0, r, s, rtol=1e-3))
    return results


def _mesh_grid(mesh):
    if isinstance(mesh, openmc.RectilinearMesh):
        return np.asarray(mesh.x_grid)
    return np.asarray(mesh.r_grid)


# -----------------------------------------------------------------------------
# Ganapol one-group benchmarks
# -----------------------------------------------------------------------------

def _one_group_material(c, fissile=False):
    """One-group material with a total cross section of 1 cm^-1."""
    if fissile:
        # Half of the collisions are scattering and half are fission
        return Material.one_group('fuel', 1.0, 0.5, 0.0, 0.5, (c - 0.5)/0.5)
    return Material.one_group('medium', 1.0, c, 1.0 - c)


def ganapol_model(problem, directory, particles=100000, batches=100,
                  inactive=50):
    """Build an OpenMC model of one of the Ganapol benchmarks.

    Parameters
    ----------
    problem : validation.analytic.ganapol.Problem
        Problem to model
    directory : str or pathlib.Path
        Directory where the cross section library is written
    particles : int
        Number of particles per batch
    batches : int
        Number of batches
    inactive : int
        Number of inactive batches (criticality problems only)

    Returns
    -------
    openmc.Model
        Model of the problem

    """
    openmc.reset_auto_ids()
    critical = problem.kind == 'critical-cylinder'
    material = _one_group_material(problem.c, fissile=critical)
    mats, collection = _materials([material], directory)
    mat = mats[material.name]

    settings = _settings(False)
    settings.particles = particles
    settings.batches = batches
    if critical:
        settings.run_mode = 'eigenvalue'
        settings.inactive = inactive
    else:
        settings.run_mode = 'fixed source'

    builder = _GANAPOL_BUILDERS[problem.kind]
    geometry, source, tallies = builder(problem, mat)
    source.energy = _source_energy(1)
    settings.source = source
    return openmc.Model(geometry, collection, settings,
                        openmc.Tallies(tallies))


def _beam_slab(problem, mat, x_min, x_max, positions, width,
               surface_point=None):
    planes, transverse = _slab_planes([x_min, x_max], ('vacuum', 'vacuum'))
    cell = openmc.Cell(fill=mat, region=+planes[0] & -planes[1] & transverse)
    geometry = openmc.Geometry([cell])

    h = _TRANSVERSE/2
    x0 = max(x_min, 0.) + 1e-10
    source = openmc.IndependentSource(
        space=openmc.stats.CartesianIndependent(
            openmc.stats.Discrete([x0], [1.0]),
            openmc.stats.Uniform(-h, h), openmc.stats.Uniform(-h, h)),
        angle=openmc.stats.Monodirectional((1., 0., 0.)))

    tallies = []
    interior = [x for x in positions if x != surface_point]
    tally = openmc.Tally(name='flux')
    tally.filters = [openmc.MeshFilter(_slab_mesh(
        _point_bins(interior, width)))]
    tally.scores = ['flux']
    tallies.append(tally)
    if surface_point is not None:
        # Flux of the neutrons leaving through the free surface
        tally = openmc.Tally(name='surface flux left')
        tally.filters = [openmc.SurfaceFilter([planes[0]])]
        tally.scores = ['flux']
        tallies.append(tally)
    return geometry, source, tallies


def _point_bins(positions, width):
    """Grid with a bin of the given width at each position. A position
    within width/2 of zero gets a one-sided bin ending at zero so that bins
    do not straddle a source plane."""
    edges = set()
    for p in positions:
        if 0. <= p < width/2:
            edges.update((0., width))
        elif -width/2 < p < 0.:
            edges.update((-width, 0.))
        else:
            edges.update((p - width/2, p + width/2))
    return np.array(sorted(edges))


def _infinite_beam(problem, mat):
    # Neutrons that travel 60 mfp from the source plane have a negligible
    # chance of returning to the region of interest for c = 0.9
    positions = [x for x, _ in problem.flux]
    return _beam_slab(problem, mat, -60., 60., positions, 0.02)


def _half_space_beam(problem, mat):
    # The reference point closest to the free surface (x = 1e-5) is
    # evaluated on the surface itself with the surface flux estimator
    positions = [x for x, _ in problem.flux]
    return _beam_slab(problem, mat, 0., 300., positions, 0.02,
                      surface_point=min(positions))


def _slab_beam(problem, mat):
    planes, transverse = _slab_planes([0., problem.size],
                                      ('vacuum', 'vacuum'))
    cell = openmc.Cell(fill=mat, region=+planes[0] & -planes[1] & transverse)
    geometry = openmc.Geometry([cell])

    h = _TRANSVERSE/2
    source = openmc.IndependentSource(
        space=openmc.stats.CartesianIndependent(
            openmc.stats.Discrete([1e-10], [1.0]),
            openmc.stats.Uniform(-h, h), openmc.stats.Uniform(-h, h)),
        angle=openmc.stats.Monodirectional((1., 0., 0.)))

    # Cosine bins of width 0.02 centered on each reference |mu|. OpenMC
    # measures the cosine relative to the direction of crossing, so it is
    # positive on both surfaces. The transmitted bins stop just short of
    # mu = 1 to exclude the uncollided beam, which the reference omits.
    mus = sorted({mu for mu, _ in problem.reflected + problem.transmitted})
    edges = np.unique(np.clip(_point_bins(mus, 0.02), 0., 1.))
    tallies = []
    for name, surface in (('reflected', planes[0]),
                          ('transmitted', planes[1])):
        bins = edges.copy()
        if name == 'transmitted':
            bins[-1] = 1. - 1e-9
        tally = openmc.Tally(name=name)
        tally.filters = [openmc.SurfaceFilter([surface]),
                         openmc.MuSurfaceFilter(bins)]
        tally.scores = ['current']
        tallies.append(tally)
    return geometry, source, tallies


def _cylinder(radius, mat):
    cyl = openmc.ZCylinder(r=radius, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-cyl & _z_slice())
    return openmc.Geometry([cell]), cyl


def _cylinder_boundaries(radius, surface=None):
    return {'outer': (radius, surface, 2*np.pi*radius*_TRANSVERSE)}


def _critical_cylinder(problem, mat):
    R = problem.size
    geometry, cyl = _cylinder(R, mat)
    r = np.linspace(0., R, 201)
    centers = 0.5*(r[1:] + r[:-1])
    pdf = centers*j0(2.404826*centers/(R + _EXTRAPOLATION))
    source = openmc.IndependentSource(
        space=openmc.stats.CylindricalIndependent(
            _histogram(r, pdf), openmc.stats.Uniform(0., 2*np.pi),
            openmc.stats.Uniform(-_TRANSVERSE/2, _TRANSVERSE/2)))
    tallies = []
    if problem.flux:
        # The bin at the center has a tiny volume and is not compared
        positions = [f*R for f, _ in problem.flux if f > 0.]
        tallies = _point_flux_tallies(
            positions, _cylinder_boundaries(R, cyl), _cylinder_mesh, 0., R,
            0.004*R)
    return geometry, source, tallies


def _surrogate_cylinder(problem, mat):
    R = problem.size
    geometry, cyl = _cylinder(R, mat)
    source = openmc.IndependentSource(
        space=openmc.stats.CylindricalIndependent(
            openmc.stats.PowerLaw(0., R, 1), openmc.stats.Uniform(0., 2*np.pi),
            openmc.stats.Uniform(-_TRANSVERSE/2, _TRANSVERSE/2)),
        angle=openmc.stats.Isotropic())
    # The flux is flat near the axis, so the bin at the center is wider
    positions = [r for r, _ in problem.flux]
    tallies = _point_flux_tallies(
        positions, _cylinder_boundaries(R, cyl), _cylinder_mesh, 0., R,
        0.01*R, lower_width=0.05*R)
    return geometry, source, tallies


def _searchlight(problem, mat):
    # The half-space is truncated well beyond the diffusion length for c =
    # 0.99 (about 6 mfp)
    depth = 150.
    radius = 150.
    cyl = openmc.ZCylinder(r=radius, boundary_type='vacuum')
    bottom = openmc.ZPlane(0., boundary_type='vacuum')
    top = openmc.ZPlane(depth, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-cyl & +bottom & -top)
    geometry = openmc.Geometry([cell])

    source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 1e-10)),
        angle=openmc.stats.Monodirectional((0., 0., 1.)))

    z = problem.size
    positions = [rho for rho, _ in problem.flux]
    grid = _flux_bins(positions, 0., radius, 0.04)
    mesh = openmc.CylindricalMesh(r_grid=grid, z_grid=[z - 0.02, z + 0.02])
    tally = openmc.Tally(name='flux')
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ['flux']
    return geometry, source, [tally]


_GANAPOL_BUILDERS = {
    'infinite-beam': _infinite_beam,
    'half-space-beam': _half_space_beam,
    'slab-beam': _slab_beam,
    'critical-cylinder': _critical_cylinder,
    'surrogate-cylinder': _surrogate_cylinder,
    'searchlight': _searchlight,
}


def ganapol_results(problem, statepoint):
    """Compare the results of a Ganapol benchmark with the reference.

    Parameters
    ----------
    problem : validation.analytic.ganapol.Problem
        Problem that was run
    statepoint : str or pathlib.Path
        Statepoint file from the run

    Returns
    -------
    Results
        Comparisons with the reference values

    """
    results = Results(problem.name)
    with openmc.StatePoint(statepoint) as sp:
        if problem.kind == 'critical-cylinder':
            keff = sp.keff
            results.comparisons.append(Comparison(
                'keff', problem.keff, keff.nominal_value, keff.std_dev))

        if problem.kind == 'slab-beam':
            for name in ('reflected', 'transmitted'):
                tally = sp.get_tally(name=name)
                bins = tally.find_filter(openmc.MuSurfaceFilter).bins
                mean = np.abs(tally.mean.ravel())
                std = tally.std_dev.ravel()
                for mu, ref in getattr(problem, name):
                    i = _containing_bin(bins, min(mu, bins[-1][1]))
                    lo, hi = bins[i]
                    # Partial current in a cosine bin is the integral of
                    # |mu| psi(mu) over the bin
                    weight = abs(hi**2 - lo**2)/2
                    results.comparisons.append(Comparison(
                        f'{name} angular flux(|mu|={mu:g})', ref,
                        mean[i]/weight, std[i]/weight, rtol=5e-3))
            return results

        if not problem.flux:
            return results

        if problem.kind == 'critical-cylinder':
            # Normalize to the flux at r/R = 0.25 since the bin at the
            # center has a small volume
            R = problem.size
            points = [(f, v) for f, v in problem.flux if f > 0.]
            fluxes = _point_fluxes(sp, [f*R for f, _ in points],
                                   _cylinder_boundaries(R))
            f0, v0 = points[0]
            mean0, std0 = fluxes[f0*R]
            for f, v in points[1:]:
                mean, std = fluxes[f*R]
                r, s = _ratio(mean[0], std[0], mean0[0], std0[0])
                results.comparisons.append(Comparison(
                    f'flux(r/R={f:g})/flux(r/R={f0:g})', v/v0, r, s,
                    rtol=1e-3))
        elif problem.kind == 'surrogate-cylinder':
            # The source density is one particle per unit volume of the
            # cylinder, so (1 - c) phi / Q = 1 - F
            R = problem.size
            q = 1.0/(np.pi*R**2*_TRANSVERSE)
            fluxes = _point_fluxes(sp, [r for r, _ in problem.flux],
                                   _cylinder_boundaries(R))
            for r, F in problem.flux:
                mean, std = fluxes[r]
                results.comparisons.append(Comparison(
                    f'1 - F(r={r:g})', 1.0 - F, (1 - problem.c)*mean[0]/q,
                    (1 - problem.c)*std[0]/q, rtol=1e-3))
        else:
            # Beam problems: the source is one particle per unit area
            # (slabs) or a unit point beam (searchlight)
            boundaries = {}
            if problem.kind == 'half-space-beam':
                x_surface = min(x for x, _ in problem.flux)
                boundaries['left'] = (x_surface, None, _TRANSVERSE**2)
            fluxes = _point_fluxes(sp, [x for x, _ in problem.flux],
                                   boundaries)
            coordinate = 'rho' if problem.kind == 'searchlight' else 'x'
            for x, ref in problem.flux:
                mean, std = fluxes[x]
                if 'left' in boundaries and x == boundaries['left'][0]:
                    # Add the incident beam of unit strength (mu = 1)
                    mean = mean + 1.0
                results.comparisons.append(Comparison(
                    f'flux({coordinate}={x:g})', ref, mean[0], std[0],
                    rtol=2e-3))
    return results


def _containing_bin(bins, value):
    for i, (lo, hi) in enumerate(bins):
        if lo <= value <= hi:
            return i
    raise ValueError(f'No bin contains {value}')
