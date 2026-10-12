def pytest_addoption(parser):
    group = parser.getgroup('analytic benchmarks')
    group.addoption('--particles', type=int, default=None,
                    help='Particles per batch for the OpenMC runs')
    group.addoption('--batches', type=int, default=None,
                    help='Number of batches for the OpenMC runs')
    group.addoption('--inactive', type=int, default=None,
                    help='Number of inactive batches for eigenvalue runs')
    group.addoption('--threads', type=int, default=None,
                    help='Number of OpenMP threads for the OpenMC runs')
    group.addoption('--n-sigma', type=float, default=4.0,
                    help='Allowed difference from the reference values in '
                    'standard deviations')
