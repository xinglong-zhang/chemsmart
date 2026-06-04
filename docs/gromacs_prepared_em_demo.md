# GROMACS prepared EM demo

This document describes the first real GROMACS validation case for the
ChemSmart GROMACS integration.

## Goal

The goal is to validate the minimal real execution chain:

```text
project.yaml -> GromacsProjectSettings -> GromacsEMJob -> GromacsJobRunner -> grompp -> mdrun