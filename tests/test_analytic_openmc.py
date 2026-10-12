"""Run the analytic benchmark problems with OpenMC and compare with the
reference solutions.

These tests require the ``openmc`` executable and Python API, version 0.16.0
or later (for the surface flux estimator). Each test runs a full OpenMC
calculation, so the whole module takes a while; use ``-k`` to select
problems and the ``--particles``/``--batches`` options to change the
statistics.
"""

import shutil

import pytest

openmc = pytest.importorskip('openmc')

from validation.analytic import ganapol, sood
from validation.analytic.run import run_ganapol, run_sood


def _version(text):
    return tuple(int(x) for x in text.split('.')[:2])


pytestmark = [
    pytest.mark.openmc,
    pytest.mark.skipif(shutil.which('openmc') is None,
                       reason='openmc executable not found'),
    pytest.mark.skipif(_version(openmc.__version__) < (0, 16),
                       reason='OpenMC 0.16.0 or later is required'),
]


def _settings(request):
    names = ('particles', 'batches', 'inactive', 'threads')
    return {name: request.config.getoption(name) for name in names
            if request.config.getoption(name) is not None}


def _check(results, n_sigma):
    failed = [c for c in results.comparisons if not c.passed(n_sigma)]
    message = '\n'.join(
        f'{c.quantity}: reference {c.reference:.7g}, calculated '
        f'{c.mean:.7g} +/- {c.std_dev:.3g} ({c.sigma:+.2f} sigma)'
        for c in failed)
    assert not failed, message


@pytest.mark.parametrize('problem', sood.PROBLEMS,
                         ids=lambda p: f'{p.number:02d}-{p.name}')
def test_sood(problem, request, tmp_path):
    results = run_sood(problem, tmp_path, **_settings(request))
    _check(results, request.config.getoption('n_sigma'))


@pytest.mark.parametrize('problem', ganapol.PROBLEMS, ids=lambda p: p.name)
def test_ganapol(problem, request, tmp_path):
    results = run_ganapol(problem, tmp_path, **_settings(request))
    _check(results, request.config.getoption('n_sigma'))
