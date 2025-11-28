
---
# chemsmart - Chemistry Simulation and Modeling Automation Toolkit

[![codecov](https://codecov.io/gh/xinglong-zhang/chemsmart/branch/main/graph/badge.svg?token=chemsmart_token_here)](https://codecov.io/gh/xinglong-zhang/chemsmart)
[![CI](https://github.com/xinglong-zhang/chemsmart/actions/workflows/main.yml/badge.svg)](https://github.com/xinglong-zhang/chemsmart/actions/workflows/main.yml)

---
<p align="center">
  <img src="docs/source/_static/chemsmart_logo.png" alt="CHEMSMART Logo" width="600"/>
</p>

---
Notice: If you have cloned this package before and find something that did not work, updating this repo via `git pull` will likely fix it. If you need additional features, please do not hesitate to get in touch!

Chemsmart is a python-based toolkit for the automatic creation of input and submission script files, the submission and the analysis of quantum chemistry simulation jobs.

It uses the same submission command regardless of the queueing systems (SLURM, Torque or SLF) used by any High Performance Computing (HPC) cluster. 

Users can customize their own HPC server settings and project settings to run different jobs, without modifying the codes in this package.

## Installation

Users are recommended to install conda environments to mange the packages required by this software toolkit. Either Anaconda3 or Miniconda3 may be installed; see more information on conda installation process at [Conda Installation Page](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html).

Once conda has been installed successfully, one can clone this package to a local directory via

```bash
git clone https://github.com/xinglong-zhang/chemsmart.git
```
First, cd into chemsmart folder. For linux and MacOS systems which support `make`, users can run 

```bash
make env
```
to create a running environment.

By default, this will create a conda environment named `chemsmart`, which installs all the required python packages for this toolkit.

If conda is not installed, one can run

```bash
make env USE_CONDA=false
```

or 

```bash
make virtualenv
```
to install using virtualenv. It is however recommanded that `conda` be used.

Help options are available by typing `make help`.

After the virtual conda environment is created and activated via `conda activate chemsmart`, one can run

```bash
make install
```
which installs the packages and dependencies required for `chemsmart` package.

For developers, one may run

```bash
make install-dev
```
which installs additoinal packages and dependencies (dev, test, docs dependencies in pyproject.toml) required for developing `chemsmart` package.

Next, run
```bash
make configure
```
to sets up the user-specific directory `~/.chemsmart` automatically. You will be prompt to enter the paths to g16 and ORCA software, which will then be added automatically. The correct `conda` path for the user will also be updated.

The configuration also adds the environment variables for chemsmart to the user `~/.bashrc` file.

-----
The `~/.chemsmart/usersettings.yaml` file contains informations such as project number or account number that are required in a typical submission script that specifies the account for use at some HPC servers. It can also contain options specifying user's email to inform user of the job start and job end once a job is submitted. If more features are needed, please submit a request via `Issues`. A typical `~/.chemsmart/usersettings.yaml` file looks like this:
```text
PROJECT: 1234567  # alias ACCOUNT FOR SLURM
EMAIL: abc@gmail.com
```
-----
The `~/.chemsmart/server/` directory contains files related to server setup for a particular HPC cluster that the user is using. For example, we can specify a SLURM based server setting as `~/.chemsmart/server/shared.yaml` with the following information:

```text
SERVER:
    SCHEDULER: SLURM
    QUEUE_NAME: RM-shared
    NUM_HOURS: 48
    MEM_GB: 100
    NUM_CORES: 64
    NUM_GPUS: Null
    NUM_THREADS: 64
    SUBMIT_COMMAND: sbatch
    ##PROJECT: 13003611
    ##PROJECT: 13002374
    SCRATCH_DIR: null
    USE_HOSTS: true
    EXTRA_COMMANDS: |
        export PATH=$HOME/bin/chemsmart:$PATH
        export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
        export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
        export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
GAUSSIAN:
    EXEFOLDER: ~/bin/g16 
    LOCAL_RUN: True 
    SCRATCH: True  # set scratch to True to run in scratch folder
    CONDA_ENV: |   # program-specific conda env
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate chemsmart
    MODULES: |
        module purge
        # module load craype-x86-rome
        # module load libfabric/1.11.0.4.125
    SCRIPTS: |
        tcsh -c "source ~/programs/g16/bsd/g16.login"
    ENVARS: |
        export SCRATCH=/tmp # required if scratch is true
        export GAUSS_EXEDIR=~/bin/g16
        export g16root=~/bin/g16
           
ORCA:
    EXEFOLDER: ~/bin/orca_6_0_1
    LOCAL_RUN: False 
    ENVARS: |
        export PATH=$HOME/bin/openmpi-4.1.6/build/bin:$PATH
        export LD_LIBRARY_PATH=$HOME/bin/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH
```
This file can be customized by user for different submission systems. This file contains the server configuration information that is needed for chemsmart to automatically write the submission script for each job.

-----
The `~/.chemsmart/gaussian/` directory contains files related to gaussian project settings, which contain DFT functional and basis set etc, that is required to write the input file for running a gaussian job. For example, we can specify a test project settings in `~/.chemsmart/gaussian/test.yaml` with the following information:

```test
gas:
  functional: m062x  # quotes required for string with spaces
  basis: def2svp
  solvent_model: smd
  solvent_id: dichloroethane
solv: 
  functional: m062x
  basis: def2tzvp
  freq: False
  solvent_model: smd
  solvent_id: dichloroethane
td:
  functional: cam-b3lyp
  basis: genecp
  heavy_elements: ['I']
  heavy_elements_basis: def2-SVPD
  light_elements_basis: def2SVP
  freq: False
  ##solvent_model: smd
  ##solvent_id: DiethylEther
```
By default, the `gas` phase settings are used for all jobs such as geometry optimization, transition state search etc, and the `solv` settings are used for single point calculations; the `td` settings are used to run TD-DFT calculations. One can specify additional project settings in `~/.chemsmart/gaussian/` in a similar way to adapt to each project that one wishes to run. If setting
```text
gas: Null
```
Then all jobs will use settings specified in `solv`, i.e., all calculations will be run in implicit solvation model.

---
The `~/.chemsmart/orca/` directory contains files related to ORCA project settings, which contain DFT functional and basis set etc, that is required to write the input file for running an ORCA job. For example, we can specify a test project settings in `~/.chemsmart/orca/test.yaml` with the following information:

```text
gas:
  functional: M062X
  basis: def2-SVP
solv:
  ab_initio: DLPNO-CCSD(T)
  functional: Null
  basis: Extrapolate(2/3,cc)
  aux_basis: AutoAux
  defgrid: DEFGRID3
  freq: False
  scf_tol: TightSCF
  scf_algorithm: KDIIS
  scf_maxiter: 500
  mdci_cutoff: Normal
  mdci_density: None    
  dipole: False
  solvent_model: SMD
  solvent_id: "toluene"
```
This will run jobs in the gas phase (geometry and TS opt etc) using M062X/def2-SVP method and run single point with solvent correction using DLPNO-CCSD(T)/CBS with cc-pVDZ/cc-pVTZ extrapolation in SMD(toluene), for example. Again, users can customize different settings in different `~/.chemsmart/orca/*project_settings*.yaml` files to adapt to different project requirements.

---
Although `make configure` would set up `~/.chemsmart` mostly correctly, a user should check the contents in `~/.chemsmart` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g., modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one may copy e.g., `~/.chemsmart/server/SLURM.yaml` to your own customised server `~/.chemsmart/server/custom.yaml` and modify it accordingly, such that the submission becomes `chemsmart sub -s custom <other commands>`.

One also need to set up scratch directories where scratch jobs may be run (for Gaussian and ORCA jobs, by default, these are run in scratch folder), one may do `ls -s /path/to/scratch/ ~/scratch`.

Note also that a user can modify the contents in `~/.chemsmart` files freely without affecting or needing to know the `chemsmart` source code.

The `make configure` will also add the required paths to the user `~/.bashrc` file. User may need to do 

```bash
source ~/.bashrc
```

to effect the changes.

<!-- ---
Once `make configure` is done, one can optionally run 
```bash
make fmt
```
and

```bash
make lint
```
to format and lint the codes (this should have been handled by the developers). Also optionally, one can run 

```bash
make test
```

to make sure that all tests in chemsmart pass.

---
Finally one can clean up by running

```bash
make clean
``` -->

## Testing Installations

Installations is deemed successfully if the commands `make install` and `make configure` do not return any errors. Installation will also create a `~/.chemsmart` containing the required files. In addition, the paths for chemsmart packages should be correctly added to the user `~/.bashrc` file. Finally, one should be able to run 

```bash
chemsmart --help
```
to get the options for running chemsmart package.

## Usage

With setup completed, one is able to run different Gaussian jobs via command-line interface (CLI).

To submit (and run) a geometry optimization job, do:

```py
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt
```
where `<server_name>` is the one of the servers specified in `~/.chemsmart/server/*.yaml` files, without `.yaml` extension; `<project>` is one of the project settings specified in `~/.chemsmart/gaussian/*.yaml` files, without `.yaml` extension; and `<input_file>` is an input file the user wishes to run job on. Note that this input file can be any format, such as `.xyz`, Gaussian `.com`, `.gjf` or `.log` file or ORCA `.inp` or `.out` file.

---

If one wants to submit geometry optimization with frozen atoms in the molecule (such as https://www.researchgate.net/post/Freezing-atoms-in-gaussian-how-to-do-it), one can do:

```bash 
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt -f <indices_of_atoms_to_freeze>
```
For example, to submit the geometry optimization job with atoms numbered 1 to 10 frozen, one can do

```bash 
chemsmart sub -s shared gaussian -p test -f input.com opt -f 1-10
```
Note that 1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for moleculare are 1-indexed.

---
To submit transition state modredundant job (frozen coordinates optimization), do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> modred -c <list_of_coords_to_constraint>
```
For example, to submit a modredundant job with constraints on bond between atom 4 and atom 17 and on bond between atom 9 and atom 10, do:
```bash
chemsmart sub -s shared gaussian -p test -f input.com modred -c [[4,17],[9,10]]
```

---
To submit transition state search job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> ts
```

---
To submit intrinsic reaction coordinate (IRC) job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> irc
```

---
To submit relaxed potential energy surface (PES) scan, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> scan -c <list_of_coords_to_constraint> -s <scan_step_size> -n <num_scan_steps>
```

For example, to submit the PES scan job with along bond between atom 4 and atom 17 for 10 steps with 0.1Å increment per step:

```bash
chemsmart sub -s shared gaussian -p test -f input.com scan -c [[4,17]] -s 0.1 -n 10
```

---
To submit single point job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp
```

For single-point job that user wants to test which uses different solvent model and id from that specified in `<project>`, one can do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp -sm <user_solvent_model> -si <user_solvent_id>
```
to specify a different solvent model `<user_solvent_model>` and solvent `<user_solvent_id>`.

---
To submit non-covalent interaction job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> nci
```

---
To submit RESP charges fitting job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> resp
```

Note that this creates an input file with fix route for RESP job:

```text
HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)
```

To run opt or modred or ts conformers from crest run output, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt
```
or

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j modred -c [1,4]
```

or 

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j ts
```
respectively

This optimizes all the conformers available in the `<input_file>`. Typically, the `<input_file>` is a list of all conformers obtained by CREST program and named `crest_conformers.xyz`.

To optimize a fixed number of lowest energy conformers, `n_conformers_to_opt`, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt -n <n_conformers_to_opt>
```

If the job terminates before `<n_conformers_to_opt>` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all `<n_conformers_to_opt>`are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .xyz file. In fact, whenever .xyz file is used as input, the charge and multiplicity should be specified via `-c <charge> -m <multiplicity` via CLI.

---

To optimize unique structure from an md trajectory file, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> traj 
```
This optimizes all the unique structures available in the md trajectory `<input_file>`. Typically, the `<input_file>` is a list of all structures on an md trajectory obtained by ASE md run and named `md.traj`. (TODO: this method is not properly implemented in chemsmart yet.)

To optimize a fixed number of lowest energy structures, `<num_structures_to_opt>`, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> traj -n <n_conformers_to_opt>
```
If the job terminates before `<n_conformers_to_opt>` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all `<n_conformers_to_opt>`are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .traj file. 

Two grouper types for determining/clustering unique structures are available from CLI option `-g`:

1) Sequential grouper (default), selected by option value of seq, which sequentially checks for unique structures in a given list of md structures, and 

2) Self-consistent grouper, selected by option value of sc, which self-consistently checks for unique structures in a given list of md structures using the reverse Cuthill–McKee algorithm for structure clustering. By default, only the last 0.1 proportion of the structures of the md.traj file is considered. This can be changed via cli option `-x <proportion_structures_to_use>`.

For example, to consider the last 20% of the structures in md.traj trajectory file, then uses Sequential grouper to group those structures into unique structures and run the 10 lowest energy structures from the list of unique structures found by the grouper:

```bash
chemsmart sub -s shared gaussian -p test -f imd.traj traj -x 0.2 -n 10 -g seq
```

---
To run distortion-interaction/activation-strain (DI-AS) job, do:

```bash
chemsmart sub -s <server_name> gaussian -p <project> -f <irc_output_file_for_dias> dias -i <indices_of_any_one_fragment> -n <number_of_every_n_step_along_irc_to_run>
```

For example to run DI-AS job for fragment 1 with atoms numbered from 5-17 at every 10 steps along the irc.log file:

```bash
chemsmart sub -s shared gaussian -p test -f irc.log dias -i 5-17 -n 10
```

---
If a user wants to run a job with pre-prepared Gaussian input file directly, one can run the job directly using:

```bash
chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> com
```

---
Generally, if a user wants to run job that is currently not present in our package, one can run custom job using:

```bash
chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> userjob -r <user_defined_gaussian_route> -a <appending_information_as_string_at_the_end_of_input_file_after_coordinates_specification>
```

For example, to create an input file named `user_defined_job.com` with user-specified route `mnr functional/basis solvent` etc and `B 1 2 F\nA 1 2 3 F` at the end of the input file after the specification of coordinates, run

```bash
chemsmart sub -s shared gaussian -p test -f test.com -l user_defined_job userjob -r 'mnr functional/basis solvent etc' -a 'B 1 2 F\nA 1 2 3 F'
```

---

### General options available to all jobs:

Users can specify the name of the file to be created for the job, without file extension, they want to run by using the option `-l`, e.g.:

```bash
chemsmart sub -s shared gaussian -p test -f test.com -l custom_job_name opt
```

will create input file named `custom_job_name.com` instead of the default `test_opt.com`.

Users can also simply append a string to the base name of the filename supplied, e.g.:


```bash
chemsmart sub -s shared gaussian -p test -f test.com -a append_string ts
```

will create input file named `test_append_string.com` instead of the default `test_ts.com`.

---
Users can also modify the charge and multiplicity from the CLI, e.g.:

Modify the charge in ``test.com`` to charge of +1 in the newly created input file ``test_charge.com`` via:

```bash
chemsmart sub -s shared gaussian -p test -f test.com -c 1 -a charge opt
```

Modify the multiplicity in ``test.com`` to multiplicity of 3 in the newly created input file ``test_multiplicity.com`` via

```bash
chemsmart sub -s shared gaussian -p test -f test.com -m 3 -a multiplicity opt
```

Modify the charge to +1 and multiplicity to 2 in the newly created input file ``test_charge_multiplicity.com`` via:
```bash
chemsmart sub -s shared gaussian -p test -f test.com -c 1 -m 2 -l test_charge_multiplicity opt
```

This can be useful when, e.g., using optimized structure of a neutral closed-shell (charge 0, multiplicity 1) system to run a charged radical ion (e.g., charge +1 and multiplicity 2 in radical cation).

---
Users can also modify the functional and basis from the CLI to differ from those in project settings, e.g.:

Modify the functional to `b3lyp` in the newly created input file ``test_functional.com`` via

```bash
chemsmart sub -s shared gaussian -p test -f test.com -x b3lyp -a functional opt
```

Modify the basis to `6-31G*` in the newly created input file ``test_basis.com`` via


```bash
chemsmart sub -s shared gaussian -p test -f test.com -b "6-31G*" -a basis opt
```

Users can also specify additional optimization options for `opt=()` in the route, for example,

```bash
chemsmart sub -s shared gaussian -p test -f test.com -o maxstep=8,maxsize=12 -a opt_options opt
```

will create `opt=(maxstep=8,maxsize=12)` as part of the route in the newly created input file `test_opt_options.com`.

Users can also add in additional parameters used in the route, e.g.,

```bash
chemsmart sub -s shared gaussian -p test -f test.com --r nosymm -a route_params opt
```

will add in `nosymm` as part of the route in the newly created input file `test_route_params.com`.

If one has more than one structure in the supplied file for input preparation, one can select the particular structure to perform job on by using the `-i/--index` option, e.g.:

```bash
chemsmart sub -s shared gaussian -p test -f small.db -i 5 -c 0 -m 1 opt
```
will take the 5th structure (1-indexed, as in chemsmart) from ase database file, `small.db`, to create the input file for geometry optimization.

---
Similar commands exists for ORCA job submssions. One can run 

```bash
chemsmart sub orca --help
```
to find out more.

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
