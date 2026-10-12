# Validation

This repository contains a collection of validation scripts, notebooks, results, etc. that have been used in the preparation of reports and articles. The following scripts are currently available:

- `openmc-run-benchmarks` -- Run a collection of ICSBEP benchmarks with OpenMC or MCNP6
- `openmc-plot-benchmarks` -- Plot results from `openmc-run-benchmarks`
- `openmc-report-benchmarks` -- Create .tex file comparing results from `openmc-run-benchmarks` to experimental values
- `openmc-validate-neutron-physics` -- Used for validating neutron transport in OpenMC. The script generates an infinite medium problem of a single nuclide with a point source at a single energy and compares the resulting neutron energy spectra from OpenMC and either MCNP6 or Serpent 2. This script requires that you have MCNP6 or Serpent 2 installed (`mcnp6` or `sss2` executable available).
- `openmc-validate-photon-physics` -- Used for validating photon transport in OpenMC. The script generates an infinite medium problem of a single element with a point source at a single energy and compares the resulting photon energy spectra from OpenMC and either MCNP6 or Serpent 2. This script requires that you have MCNP6 or Serpent 2 installed (`mcnp6` or `sss2` executable available).
- `openmc-validate-photon-production` -- Used for validating photon production in OpenMC. This script generates a broomstick problem where monoenergetic neutrons are shot down a thin, infinitely long cylinder. The energy spectra of photons produced from neutron reactions is tallied along the surface of the cylinder. Comparisons are made between OpenMC and either MCNP6 or Serpent 2. Again, this script requires that MCNP6 or Serpent 2 is properly installed.

- `openmc-run-analytic-benchmarks` -- Run multigroup problems with analytic solutions and compare OpenMC against them (see below).

## Analytic benchmarks

The `validation.analytic` package contains problems with exact (semi-)analytic solutions of the transport equation, which verify the OpenMC multigroup transport algorithms independently of any nuclear data:

- `sood` -- All 75 criticality problems from A. Sood, R. A. Forster, and D. K. Parsons, "Analytical Benchmark Test Set for Criticality Code Verification," LA-13511 (1999): infinite media, slabs, cylinders, spheres, reflected systems, and infinite slab lattice cells with one, two, three, and six groups and isotropic, P1, and P2 scattering. Each problem checks keff (or k-infinity), and, where the report gives them, group flux ratios and normalized spatial flux distributions.
- `ganapol` -- One-group problems from B. D. Ganapol, "Analytical Benchmarks for Nuclear Engineering Applications," NEA/DB/DOC(2008)1: the plane beam in an infinite medium (Benchmark 3.1), the half-space albedo problem (3.2), reflected and transmitted angular fluxes from a slab (3.3), critical radii, critical flux distributions, and fixed-source fluxes of infinite cylinders (3.4), and the searchlight problem (5.1).

Run the problems with, for example:

```bash
openmc-run-analytic-benchmarks sood                 # all LA-13511 problems
openmc-run-analytic-benchmarks sood 12 44 'URR*'    # by number or name
openmc-run-analytic-benchmarks ganapol searchlight -p 1000000
```

A problem passes when every calculated quantity is within `--n-sigma` (default 4) standard deviations of the reference value. Spatial fluxes add a small relative tolerance for the width of the tally bins. LA-13511 problems 34, 37, 43, and 71 use a Legendre scattering expansion that is negative over part of the cosine interval. Such distributions cannot be sampled exactly by a Monte Carlo code, so they are skipped unless `--include-negative-pdf` is given. Fluxes on vacuum boundaries use the surface flux tally added in OpenMC 0.16.0, so OpenMC 0.16.0 or later is required.

## Tests

The `tests` directory contains pytest tests. `tests/test_analytic_data.py` checks the transcribed benchmark data without running OpenMC: it confirms that the cross sections reproduce the published k-infinity and flux ratios and that the dimensions in centimeters match those in mean free paths. `tests/test_analytic_openmc.py` runs every analytic benchmark with OpenMC and is skipped if OpenMC is not available:

```bash
python -m pip install .[test]
pytest tests/test_analytic_data.py              # fast
pytest tests/test_analytic_openmc.py -k sood    # runs OpenMC
pytest tests/test_analytic_openmc.py --particles 100000 --batches 300
```

## Installation

The validation package can be installed via:

```bash
python -m pip install .
```

