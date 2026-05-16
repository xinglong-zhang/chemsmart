#!/usr/bin/env python
import os
os.environ['OMP_NUM_THREADS'] = '1'

from chemsmart.cli.run import run

def run_job():
    run(['--server', 'cu_dev', '--no-fake', '--no-delete-scratch', '--no-debug', '--stream', 'gaussian', '--filename', '2a.xyz', '--project', 'pka_thf', '--no-forces', '--no-remove-solvent', 'pka', '--reference', '1a.xyz', '--reference-proton-index', '10', '--reference-charge', '0', '--reference-multiplicity', '1', '--skip-completed', '--scheme', 'proton exchange', '--delta-g-proton', '-265.9', '--temperature', '298.15', '--concentration', '1.0', '--cutoff-entropy-grimme', '100.0', '--cutoff-enthalpy', '100.0', 'batch', '--skip-completed', '-pi', '9', '-c', '0', '-m', '1'])

if __name__ == '__main__':
    run_job()