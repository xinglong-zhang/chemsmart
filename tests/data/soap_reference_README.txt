soap_reference.npz — golden SOAP vectors for regression tests.

Provenance: generated from DScribe 2.1.2 (GTO, average="off", compression="off")
using tests/data/generate_soap_reference.py in the chemsmart conda environment.

Package versions used for the current file:
  dscribe 2.1.2
  ase 3.24.0
  scipy 1.15.2
  numpy 1.26.4

Keys / geometries / hyperparameters:

  water_sigma_0_5
    water O/H/H positions from tests/test_soap.py water_molecule fixture
    n_max=4, l_max=2, sigma=0.5, r_cut=6.0, species=["H","O"]

  methanol_centers_2_1
    methanol C/O/H4 from methanol_molecule fixture
    centers=[2,1] (1-based) → rows for O then C
    n_max=4, l_max=2, sigma=1.0, r_cut=6.0, species=["C","H","O"]

  water_species_HCON
    water geometry; species=["H","C","O","N"]
    n_max=4, l_max=2, sigma=1.0, r_cut=6.0

  sulfur_species_HOS
    S/H/H at water coordinates; species=["H","O","S"]
    n_max=4, l_max=2, sigma=1.0, r_cut=6.0

  water_defaults
    water; public defaults n_max=8, l_max=6, sigma=1.0, r_cut=6.0
    species=["H","O"]

  water_lmax4_sigma1
    water; n_max=4, l_max=4, sigma=1.0, r_cut=6.0
    covers solid-harmonic / lpmv path beyond l=2

  methanol_sigma1
    methanol; n_max=4, l_max=2, sigma=1.0, r_cut=6.0
    species=["C","H","O"]

  padding_shell_He
    water + He at (8.0, 0.0, 0.11779) Å from O — inside DScribe's
    r_cut + sigma*sqrt(-2*ln(1e-3)) neighbor shell (~9.72 Å for sigma=1)
    n_max=4, l_max=2, sigma=1.0, r_cut=6.0, species=["H","He","O"]

Regenerate only when intentionally changing SOAP conventions:

  python tests/data/generate_soap_reference.py

Then re-verify parity against DScribe 2.1.2 or update tolerances accordingly.
