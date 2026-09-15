
---
# CHEMSMART - Chemistry Simulation and Modeling Automation Toolkit

[![codecov](https://codecov.io/gh/xinglong-zhang/chemsmart/graph/badge.svg?token=L01DRALL5E)](https://codecov.io/gh/xinglong-zhang/chemsmart)
[![CI](https://github.com/xinglong-zhang/chemsmart/actions/workflows/main.yml/badge.svg)](https://github.com/xinglong-zhang/chemsmart/actions/workflows/main.yml)

---
<p align="center">
  <img src="docs/source/_static/chemsmart_logo.png" alt="CHEMSMART Logo" width="600"/>
</p>

---
## Notice: 

**If you have cloned this package before and find something that did not work, updating this repo via `git pull` will likely fix it. If you need additional features, please do not hesitate to get in touch!**

GitHub does not natively support forcing links in a *README.md* file to open in a new tab or window. If needed, please use one of the following methods:

Mac: Hold CMD + click the link.

Windows / Linux: Hold CTRL + click the link.

Mouse: Click the link using the middle scroll wheel.


## Introduction

CHEMSMART is a Python-based toolkit for the automatic creation of input and submission script files, the submission and the analysis of quantum chemistry simulation jobs.

It uses the same submission command regardless of the queueing systems (SLURM, Torque or SLF) used by any High Performance Computing (HPC) cluster. 

Users can customize their own HPC server settings and project settings to run different jobs, without modifying the codes in this package.

For completed introduction videos, please refer to [YouTube](https://www.youtube.com/watch?v=hAQOqj9c_4w) and [Bilibili](https://www.bilibili.com/video/BV1d59LBiEb5/?spm_id_from=333.1387.homepage.video_card.click&vd_source=fd103f8baa97dc76be9c7402dceee380).


## Tutorial

Full tutorials are available on [Read the Docs](https://chemsmart.readthedocs.io/en/latest/).

**Getting Started**

First, users can select the appropriate installation method for their operating system.

- [Installation (Linux and macOS)](https://chemsmart.readthedocs.io/en/latest/installation-linux-macos.html)
- [Installation (Windows WSL/Ubuntu)](https://chemsmart.readthedocs.io/en/latest/installation-windows-wsl.html)
- [Installation (Windows Git Bash)](https://chemsmart.readthedocs.io/en/latest/installation-windows-gitbash.html)
- [Installation (Windows Anaconda PowerShell)](https://chemsmart.readthedocs.io/en/latest/installation-windows-powershell.html)
- [Installation (HPC Cluster)](https://chemsmart.readthedocs.io/en/latest/installation-hpc-cluster.html)

For Windows users, it is recommended to install using Windows WSL for full functionality support. 

The corresponding tutorial video for installation is available on [YouTube](https://www.youtube.com/watch?v=VtY7fhvfptQ) and [Bilibili](https://www.bilibili.com/video/BV1abQtBZEKf/?spm_id_from=333.1387.homepage.video_card.click&vd_source=fd103f8baa97dc76be9c7402dceee380).

**Configuration**

After installation, users can configure their user settings, server settings, and project settings according to their needs.

- [Configuration Overview](https://chemsmart.readthedocs.io/en/latest/configuration-overview.html)
- [User Settings](https://chemsmart.readthedocs.io/en/latest/configuration-user-settings.html)
- [Server Settings](https://chemsmart.readthedocs.io/en/latest/configuration-server-settings.html)
- [Project Settings](https://chemsmart.readthedocs.io/en/latest/configuration-project-settings.html)

The corresponding tutorial video for configuration is available on [YouTube](https://www.youtube.com/watch?v=orvZhi7FFBs) and [Bilibili](https://www.bilibili.com/video/BV1dDRMBQEQt/?spm_id_from=333.1387.homepage.video_card.click&vd_source=fd103f8baa97dc76be9c7402dceee380).

**CLI Reference**

Before starting computational workflows with CHEMSMART, we recommend taking a few minutes to familiarize yourself with the basic command-line interface, molecular input options, file conversion, and ChemDraw-based structure preparation.

- [Command Line Interface Overview](https://chemsmart.readthedocs.io/en/latest/cli-overview.html)
  Tutorial video: [YouTube](https://www.youtube.com/watch?v=ydZFesfYO3k) and [Bilibili](https://www.bilibili.com/video/BV1NjVw67Ejw/?spm_id_from=333.1387.homepage.video_card.click&vd_source=fd103f8baa97dc76be9c7402dceee380).

- [Molecule Input Formats](https://chemsmart.readthedocs.io/en/latest/molecule-input-formats.html)
- [Convert CLI Options](https://chemsmart.readthedocs.io/en/latest/convert-cli-options.html)
- [ChemDraw Files](https://chemsmart.readthedocs.io/en/latest/chemdraw-organometallic.html)

**Gaussian Jobs**

This section provides practical guides for setting up and running Gaussian calculations with CHEMSMART, covering geometry optimization, transition-state searches, conformational sampling, QRC calculations, electronic-structure analyses, QM/MM ONIOM calculations, and other common workflows.

- [Gaussian CLI Options](https://chemsmart.readthedocs.io/en/latest/gaussian-cli-options.html)
- [Structure Optimization (Gaussian)](https://chemsmart.readthedocs.io/en/latest/gaussian-structure-optimization.html)
- [Transition State Search (Gaussian)](https://chemsmart.readthedocs.io/en/latest/gaussian-transition-state.html)
- [Conformational Sampling & Dynamics (Gaussian)](https://chemsmart.readthedocs.io/en/latest/gaussian-conformational-sampling.html)
- [Quick Reaction Coordinate (QRC) Jobs](https://chemsmart.readthedocs.io/en/latest/gaussian-qrc.html)
- [Electronic Structure Properties & Analyses (Gaussian)](https://chemsmart.readthedocs.io/en/latest/gaussian-electronic-structure.html)
- [Gaussian QM/MM ONIOM Calculations Guide](https://chemsmart.readthedocs.io/en/latest/gaussian-qmmm-jobs.html)
- [Other Gaussian Jobs](https://chemsmart.readthedocs.io/en/latest/gaussian-other-jobs.html)

Two Gaussian tutorial videos are available:
Gaussian Jobs part I: [YouTube](https://www.youtube.com/watch?v=NIlXvSMUXYo) | [Bilibili](https://www.bilibili.com/video/BV1hSVB68EtA/?spm_id_from=333.1387.homepage.video_card.click)
Gaussian Jobs part II: [YouTube](https://www.youtube.com/watch?v=F8dKqplusgA) | [Bilibili](https://www.bilibili.com/video/BV1wGEb6YEqh/?spm_id_from=333.1387.homepage.video_card.click&vd_source=fd103f8baa97dc76be9c7402dceee380)

**ORCA Jobs**

This section provides practical guides for setting up and running ORCA calculations with CHEMSMART, including geometry optimization, single-point calculations, transition-state and reaction-path searches, direct ORCA input, and multiscale QM/MM calculations.

- [ORCA CLI Options](https://chemsmart.readthedocs.io/en/latest/orca-cli-options.html)
- [Structure Optimization (ORCA)](https://chemsmart.readthedocs.io/en/latest/orca-structure-optimization.html)
- [Transition State Search (ORCA)](https://chemsmart.readthedocs.io/en/latest/orca-transition-state.html)
- [Direct ORCA Input](https://chemsmart.readthedocs.io/en/latest/orca-direct-input.html)
- [ORCA QM/MM Multiscale Calculations Guide](https://chemsmart.readthedocs.io/en/latest/orca-multiscale-calculations.html)

**xTB Jobs**

This section provides guides for running xTB calculations with CHEMSMART, including geometry optimization, single-point calculations, Hessian and frequency calculations, and related CLI options.

- [xTB CLI Options](https://chemsmart.readthedocs.io/en/latest/xtb-cli-options.html)
- [Structure Optimization (xTB)](https://chemsmart.readthedocs.io/en/latest/xtb-structure-optimization.html)

**Crest Jobs**

This section provides guides for performing conformational sampling with CREST through CHEMSMART, including general CLI options, free conformational searches, and constrained searches for applications such as transition-state conformer sampling.

- [CREST CLI Options](https://chemsmart.readthedocs.io/en/latest/crest-cli-options.html)
- [Conformational Search (CREST)](https://chemsmart.readthedocs.io/en/latest/crest-conformational-search.html)

**pKa Calculations**

This section describes how to set up, run, and analyze pKa calculations using CHEMSMART, including single-system and batch workflows.

- [pKa Calculations](https://chemsmart.readthedocs.io/en/latest/pka-calculations.html)

**Thermochemistry**

This section introduces thermochemistry analysis for Gaussian, ORCA, and xTB calculations, including single-file and batch-processing workflows.

- [Thermochemistry Analysis](https://chemsmart.readthedocs.io/en/latest/thermochemistry-analysis.html)

**PyMOL Visualization**

This section provides guides for molecular visualization and analysis with PyMOL, including basic visualization, reaction analysis, electronic-structure analysis, and noncovalent interaction analysis.

- [PyMOL CLI Options](https://chemsmart.readthedocs.io/en/latest/pymol-cli-options.html)
- [Basic Visualization (PyMOL)](https://chemsmart.readthedocs.io/en/latest/pymol-visualization.html)
- [Reaction Analysis (PyMOL)](https://chemsmart.readthedocs.io/en/latest/pymol-reaction-analysis.html)
- [Electronic Structure Analysis (PyMOL)](https://chemsmart.readthedocs.io/en/latest/pymol-electronic-structure.html)
- [Interaction Analysis (PyMOL)](https://chemsmart.readthedocs.io/en/latest/pymol-interaction-analysis.html)

**Grouper Tool**

This section introduces the CHEMSMART Grouper tool for clustering and selecting molecular structures using different strategies, together with workflows for CREST conformer ensembles and molecular trajectories.

- [Grouper CLI Options](https://chemsmart.readthedocs.io/en/latest/grouper-cli-options.html)
- [Grouping Strategies](https://chemsmart.readthedocs.io/en/latest/grouper-strategies.html)
- [CREST and Trajectory Workflows](https://chemsmart.readthedocs.io/en/latest/grouper-crest-or-traj-workflow.html)

**NCIPLOT**

This section provides a tutorial for performing and visualizing noncovalent interaction analysis with NCIPLOT.

- [NCIPLOT Tutorial](https://chemsmart.readthedocs.io/en/latest/nciplot-tutorial.html)

**Auxiliary Scripts**

This section introduces auxiliary scripts for common data-processing tasks, including file management and electronic-structure analysis.

- [Auxiliary Scripts Overview](https://chemsmart.readthedocs.io/en/latest/scripts-overview.html)
- [File Management](https://chemsmart.readthedocs.io/en/latest/scripts-data-management.html)
- [Electronic Structure Analysis](https://chemsmart.readthedocs.io/en/latest/scripts-electronic-analysis.html)

**API Reference**

This section provides the API reference for CHEMSMART modules, classes, and functions.

- [CHEMSMART Modules](https://chemsmart.readthedocs.io/en/latest/modules.html)

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.

## 📖 Citations

If you use **CHEMSMART** in your work, please follow good scholarly practice and kindly cite our work: [https://arxiv.org/abs/2508.20042](https://arxiv.org/abs/2508.20042). 

### Plain Text (ACS Style)

Zhang, X.; Tan, H.; Liu, J.; Li, Z.; Wang, L.; Chen, B. W. J. CHEMSMART: Chemistry Simulation and Modeling Automation Toolkit for High-Efficiency Computational Chemistry Workflows. arXiv 2025, arXiv:2508.20042. https://doi.org/10.48550/arXiv.2508.20042.


### BibTeX

```bibtex
@misc{zhang2025chemsmartchemistrysimulationmodeling,
  title        = {CHEMSMART: Chemistry Simulation and Modeling Automation Toolkit for High-Efficiency Computational Chemistry Workflows},
  author       = {Xinglong Zhang and Huiwen Tan and Jingyi Liu and Zihan Li and Lewen Wang and Benjamin W. J. Chen},
  year         = {2025},
  eprint       = {2508.20042},
  archivePrefix= {arXiv},
  primaryClass = {physics.chem-ph},
  url          = {https://arxiv.org/abs/2508.20042}
}
```


---
In addition, if you use **ASE** Atoms object in **CHEMSMART**, please cite:
### Plain Text (ACS Style)

Ask Hjorth Larsen et al The atomic simulation environment—a Python library for working with atoms. J. Phys.: Condens. Matter, 2017, 29, 273002.

### BibTeX
```bibtex
@article{Hjorth Larsen_2017,
doi = {10.1088/1361-648X/aa680e},
url = {https://dx.doi.org/10.1088/1361-648X/aa680e},
year = {2017},
month = {jun},
publisher = {IOP Publishing},
volume = {29},
number = {27},
pages = {273002},
author = {Hjorth Larsen, Ask and Jørgen Mortensen, Jens and Blomqvist, Jakob and Castelli, Ivano E and Christensen, Rune and Dułak, Marcin and Friis, Jesper and Groves, Michael N and Hammer, Bjørk and Hargus, Cory and Hermes, Eric D and Jennings, Paul C and Bjerre Jensen, Peter and Kermode, James and Kitchin, John R and Leonhard Kolsbjerg, Esben and Kubal, Joseph and Kaasbjerg, Kristen and Lysgaard, Steen and Bergmann Maronsson, Jón and Maxson, Tristan and Olsen, Thomas and Pastewka, Lars and Peterson, Andrew and Rostgaard, Carsten and Schiøtz, Jakob and Schütt, Ole and Strange, Mikkel and Thygesen, Kristian S and Vegge, Tejs and Vilhelmsen, Lasse and Walter, Michael and Zeng, Zhenhua and Jacobsen, Karsten W},
title = {The atomic simulation environment—a Python library for working with atoms},
journal = {Journal of Physics: Condensed Matter},
abstract = {The atomic simulation environment (ASE) is a software package written in the Python programming language with the aim of setting up, steering, and analyzing atomistic simulations. In ASE, tasks are fully scripted in Python. The powerful syntax of Python combined with the NumPy array library make it possible to perform very complex simulation tasks. For example, a sequence of calculations may be performed with the use of a simple ‘for-loop’ construction. Calculations of energy, forces, stresses and other quantities are performed through interfaces to many external electronic structure codes or force fields using a uniform interface. On top of this calculator interface, ASE provides modules for performing many standard simulation tasks such as structure optimization, molecular dynamics, handling of constraints and performing nudged elastic band calculations.}
}
```
---
If you use RDKit functionalities in **CHEMSMART**, please cite:

### Plain Text (ACS Style)

ARDKit: Open-source cheminformatics. https://www.rdkit.org

### BibTeX
```bibtex
@article{Landrum2016RDKit2016_09_4,
  added-at = {2017-04-11T06:11:47.000+0200},
  author = {Landrum, Greg},
  biburl = {https://www.bibsonomy.org/bibtex/28d01fceeccd6bf2486e47d7c4207b108/salotz},
  description = {Release 2016_09_4 (Q3 2016) Release · rdkit/rdkit},
  interhash = {ee9a4ddeff3121aa622cf35709fa6e21},
  intrahash = {8d01fceeccd6bf2486e47d7c4207b108},
  keywords = {chemoinformatics drug-design pharmacophores software},
  timestamp = {2017-04-11T06:11:47.000+0200},
  title = {RDKit: Open-Source Cheminformatics Software},
  url = {https://github.com/rdkit/rdkit/releases/tag/Release_2016_09_4},
  year = 2016
}
```
---
Our package has minimal dependencies on **pymatgen**, but if you convert **CHEMSMART** molecule into pymatgen **AseAtomsAdaptor**, please cite:

### Plain Text (ACS Style)
A. Jain, S.P. Ong, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson
The Materials Project: A materials genome approach to accelerating materials innovation.
*APL Materials*, 2013, 1(1), 011002.

### BibTeX
```bibtex
@article{Jain2013,
author = {Jain, Anubhav and Ong, Shyue Ping and Hautier, Geoffroy and Chen, Wei and Richards, William Davidson and Dacek, Stephen and Cholia, Shreyas and Gunter, Dan and Skinner, David and Ceder, Gerbrand and Persson, Kristin a.},
doi = {10.1063/1.4812323},
issn = {2166532X},
journal = {APL Materials},
number = {1},
pages = {011002},
title = {{The Materials Project: A materials genome approach to accelerating materials innovation}},
url = {http://link.aip.org/link/AMPADS/v1/i1/p011002/s1\&Agg=doi},
volume = {1},
year = {2013}
}
```

---
If you use **scikit-learn**, please cite

### Plain Text (ACS Style)

Pedregosa et al., Scikit-learn: Machine Learning in Python, *J. Mach. Learn. Res* 2011, 12, 2825-2830.

### BibTeX
```bibtex
@article{scikit-learn,
  title={Scikit-learn: Machine Learning in {P}ython},
  author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
          and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
          and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
          Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
  journal={Journal of Machine Learning Research},
  volume={12},
  pages={2825--2830},
  year={2011}
}
```

---
**Please also cite other relavant software (e.g., Gaussian, ORCA, NCIPLOT, PyMOL) and DFT functionals and basis sets you use in your research accordingly.**
