#!/usr/bin/env python3
"""Run analytic benchmark problems with OpenMC and compare with the
reference solutions."""

import argparse
from pathlib import Path
import shutil
import sys

from . import ganapol, sood
from .models import (ganapol_model, ganapol_results, sood_exact, sood_model,
                     sood_results)


def run_sood(problem, directory, particles=10000, batches=150, inactive=50,
             threads=None, mpi_args=None):
    """Run a LA-13511 problem and compare with the reference.

    Parameters
    ----------
    problem : validation.analytic.sood.Problem
        Problem to run
    directory : str or pathlib.Path
        Working directory for the run
    particles : int
        Number of particles per batch
    batches : int
        Total number of batches
    inactive : int
        Number of inactive batches
    threads : int, optional
        Number of OpenMP threads
    mpi_args : list of str, optional
        MPI execute command and arguments

    Returns
    -------
    validation.analytic.models.Results
        Comparisons with the reference values

    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    model = sood_model(problem, directory, particles, batches, inactive)
    statepoint = model.run(cwd=directory, threads=threads, mpi_args=mpi_args,
                           output=False)
    return sood_results(problem, statepoint)


def run_ganapol(problem, directory, particles=None, batches=None,
                inactive=50, threads=None, mpi_args=None):
    """Run a Ganapol benchmark and compare with the reference.

    Parameters are the same as for :func:`run_sood`. By default, criticality
    problems use 10,000 particles and 150 batches and fixed-source problems
    use 100,000 particles and 100 batches.

    """
    critical = problem.kind == 'critical-cylinder'
    if particles is None:
        particles = 10000 if critical else 100000
    if batches is None:
        batches = 150 if critical else 100
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    model = ganapol_model(problem, directory, particles, batches, inactive)
    statepoint = model.run(cwd=directory, threads=threads, mpi_args=mpi_args,
                           output=False)
    return ganapol_results(problem, statepoint)


def _select(problems, names):
    if not names:
        return list(problems)
    selected = []
    for name in names:
        matches = [p for p in problems if name == p.name or
                   name == str(getattr(p, 'number', '')) or
                   (name.endswith('*') and p.name.startswith(name[:-1]))]
        if not matches:
            raise ValueError(f'Unknown problem {name}')
        selected.extend(matches)
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('suite', choices=['sood', 'ganapol'],
                        help='Benchmark suite: "sood" for the LA-13511 '
                        'criticality problems or "ganapol" for the problems '
                        'from NEA/DB/DOC(2008)1')
    parser.add_argument('problems', nargs='*',
                        help='Problem names or LA-13511 problem numbers to '
                        'run. A trailing "*" matches a prefix. Defaults to '
                        'all problems.')
    parser.add_argument('-p', '--particles', type=int,
                        help='Number of particles per batch')
    parser.add_argument('-b', '--batches', type=int,
                        help='Number of batches')
    parser.add_argument('-i', '--inactive', type=int,
                        help='Number of inactive batches')
    parser.add_argument('-n', '--n-sigma', type=float, default=4.0,
                        help='Number of standard deviations allowed between '
                        'the calculated and reference values')
    parser.add_argument('-o', '--output-dir', type=Path,
                        default=Path('analytic'),
                        help='Directory for the OpenMC runs')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of OpenMP threads')
    parser.add_argument('--mpi-args', default='',
                        help='MPI execute command and any additional MPI '
                        'arguments')
    parser.add_argument('--include-negative-pdf', action='store_true',
                        help='Include LA-13511 problems whose Legendre '
                        'scattering expansion is negative over part of the '
                        'cosine interval and cannot be made nonnegative')
    args = parser.parse_args()

    if not shutil.which('openmc'):
        raise IOError('Unable to locate executable openmc in path.')

    if args.suite == 'sood':
        problems = _select(sood.PROBLEMS, args.problems)
        if not args.include_negative_pdf and not args.problems:
            problems = [p for p in problems if sood_exact(p)]
        run = run_sood
    else:
        problems = _select(ganapol.PROBLEMS, args.problems)
        run = run_ganapol

    kwargs = {k: getattr(args, k) for k in ('particles', 'batches', 'inactive')
              if getattr(args, k) is not None}
    mpi_args = args.mpi_args.split() or None

    failures = []
    for problem in problems:
        directory = args.output_dir / args.suite / problem.name
        results = run(problem, directory, threads=args.threads,
                      mpi_args=mpi_args, **kwargs)
        status = 'PASS' if results.passed(args.n_sigma) else 'FAIL'
        print(f'{status} {problem.name}', flush=True)
        for c in results.comparisons:
            mark = ' ' if c.passed(args.n_sigma) else '*'
            print(f'  {mark} {c.quantity:45s} ref={c.reference:<12.7g} '
                  f'calc={c.mean:<12.7g} +/- {c.std_dev:<10.3g} '
                  f'({c.sigma:+.2f} sigma)', flush=True)
        if status == 'FAIL':
            failures.append(problem.name)

    print(f'{len(problems) - len(failures)} of {len(problems)} problems '
          'passed')
    if failures:
        print('Failed: ' + ', '.join(failures))
        sys.exit(1)


if __name__ == '__main__':
    main()
