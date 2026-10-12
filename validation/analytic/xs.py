"""Macroscopic multigroup cross sections for analytic benchmark problems."""

from dataclasses import dataclass

import numpy as np


@dataclass
class Material:
    """Macroscopic multigroup cross section set.

    Groups are ordered from highest to lowest energy (index 0 is the fastest
    group), which is the ordering used by OpenMC. Note that many of the
    analytic benchmark references number the groups in the reverse order.

    Parameters
    ----------
    name : str
        Name of the cross section set
    total : numpy.ndarray
        Total cross section in each group [cm^-1]
    scatter : numpy.ndarray
        Legendre moments of the scattering matrix with shape
        (order + 1, groups, groups) where ``scatter[l, g, h]`` is the l-th
        moment for scattering from group g to group h [cm^-1]. The moments do
        not include the (2l + 1) factor.
    capture : numpy.ndarray
        Capture (absorption without neutron emission) cross section [cm^-1]
    fission : numpy.ndarray
        Fission cross section [cm^-1]
    nu : numpy.ndarray
        Average number of neutrons emitted per fission
    chi : numpy.ndarray
        Fission spectrum

    """
    name: str
    total: np.ndarray
    scatter: np.ndarray
    capture: np.ndarray
    fission: np.ndarray
    nu: np.ndarray
    chi: np.ndarray

    @classmethod
    def one_group(cls, name, total, scatter, capture=0., fission=0., nu=0.):
        """Create a one-group cross section set.

        Parameters
        ----------
        name : str
            Name of the cross section set
        total : float
            Total cross section [cm^-1]
        scatter : float or iterable of float
            Scattering cross section or its Legendre moments [cm^-1]
        capture : float
            Capture cross section [cm^-1]
        fission : float
            Fission cross section [cm^-1]
        nu : float
            Average number of neutrons emitted per fission

        """
        moments = np.atleast_1d(np.asarray(scatter, dtype=float))
        return cls(
            name=name,
            total=np.array([total], dtype=float),
            scatter=moments.reshape(-1, 1, 1),
            capture=np.array([capture], dtype=float),
            fission=np.array([fission], dtype=float),
            nu=np.array([nu], dtype=float),
            chi=np.array([1.0 if fission > 0. else 0.0]),
        )

    @classmethod
    def from_reference(cls, name, groups, scatter):
        """Create a cross section set from data that uses the reference
        group numbering, where the highest group number is the fastest group.

        Parameters
        ----------
        name : str
            Name of the cross section set
        groups : dict of int to dict
            Mapping of the group number to a dictionary with keys 'total',
            'capture', and optionally 'fission', 'nu', and 'chi'.
        scatter : dict of tuple to float or iterable of float
            Mapping of (from group, to group) to the scattering cross section
            or its Legendre moments [cm^-1]

        """
        numbers = sorted(groups, reverse=True)
        index = {g: i for i, g in enumerate(numbers)}
        n = len(numbers)
        order = max(len(np.atleast_1d(v)) for v in scatter.values()) - 1
        matrix = np.zeros((order + 1, n, n))
        for (g_from, g_to), value in scatter.items():
            moments = np.atleast_1d(np.asarray(value, dtype=float))
            matrix[:moments.size, index[g_from], index[g_to]] = moments

        def get(key):
            return np.array([groups[g].get(key, 0.) for g in numbers],
                            dtype=float)

        return cls(name=name, total=get('total'), scatter=matrix,
                   capture=get('capture'), fission=get('fission'),
                   nu=get('nu'), chi=get('chi'))

    @property
    def num_groups(self):
        return self.total.size

    @property
    def order(self):
        """Legendre order of the scattering matrix."""
        return self.scatter.shape[0] - 1

    @property
    def absorption(self):
        return self.capture + self.fission

    @property
    def nu_fission(self):
        return self.nu * self.fission

    @property
    def fissionable(self):
        return bool(np.any(self.nu_fission > 0.))

    def truncate(self, order):
        """Return a copy keeping only the first ``order`` Legendre moments
        of the scattering matrix."""
        return Material(self.name, self.total, self.scatter[:order + 1],
                        self.capture, self.fission, self.nu, self.chi)

    def balance(self):
        """Difference between the total cross section and the sum of the
        capture, fission, and outscatter cross sections in each group."""
        return self.total - (self.capture + self.fission +
                             self.scatter[0].sum(axis=1))

    def c(self):
        """Mean number of secondaries per collision in a one-group set."""
        return (self.scatter[0].sum(axis=1) + self.nu_fission) / self.total

    def k_inf(self):
        """Infinite medium multiplication factor and fundamental mode flux.

        Returns
        -------
        float
            Infinite medium multiplication factor
        numpy.ndarray
            Group fluxes of the fundamental mode normalized to sum to unity

        """
        removal = np.diag(self.total) - self.scatter[0].T
        production = np.outer(self.chi, self.nu_fission)
        eigvals, eigvecs = np.linalg.eig(np.linalg.solve(removal, production))
        i = np.argmax(eigvals.real)
        flux = np.abs(eigvecs[:, i].real)
        return eigvals[i].real, flux / flux.sum()
