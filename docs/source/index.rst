#########################
 CHEMSMART Documentation
#########################

.. image:: _static/chemsmart_logo.png
   :width: 400
   :align: center

CHEMSMART is the canonical command-line and project-YAML hub for automating computational-chemistry workflows across
Gaussian, ORCA, PySCF, xTB, and related tools. It centralizes scientific configuration, backend input generation, local
or scheduled execution, and result inspection so users and agents can operate supported programs through one transparent
interface.

The provider-neutral ``chemsmart agent`` pipeline -- ``plan``, ``review``, ``run``, and the ``tui`` terminal -- uses the
same project loaders and command compiler. Models propose typed scientific intent; CHEMSMART owns executable commands,
approvals, execution state, and deterministic validation.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   introduction
   installation-ubuntu-cpu-server
   installation-linux-macos
   installation-windows-wsl
   installation-windows-gitbash
   installation-windows-powershell
   installation-hpc-cluster

.. toctree::
   :maxdepth: 2
   :caption: Configuration

   configuration-overview
   configuration-test
   configuration-user-settings
   configuration-server-settings
   configuration-project-settings

.. toctree::
   :maxdepth: 2
   :caption: CLI Reference

   cli-overview
   agent-workflows
   agent-tui
   molecule-input-formats
   chemdraw-organometallic

.. toctree::
   :maxdepth: 2
   :caption: Gaussian Jobs

   gaussian-cli-options
   gaussian-structure-optimization
   gaussian-transition-state
   gaussian-conformational-sampling
   gaussian-qrc
   gaussian-electronic-structure
   gaussian-qmmm-jobs
   gaussian-other-jobs

.. toctree::
   :maxdepth: 2
   :caption: ORCA Jobs

   orca-cli-options
   orca-structure-optimization
   orca-transition-state
   orca-direct-input
   orca-multiscale-calculations

.. toctree::
   :maxdepth: 2
   :caption: PySCF and xTB Jobs

   pyscf-cli-options
   xtb-cli-options

.. toctree::
   :maxdepth: 2
   :caption: pKa Calculations

   pka-calculations

.. toctree::
   :maxdepth: 2
   :caption: Thermochemistry

   thermochemistry-analysis

.. toctree::
   :maxdepth: 2
   :caption: Molecular Database

   database-overview
   database-workflow
   database-assemble
   database-inspect
   database-query
   database-export

.. toctree::
   :maxdepth: 2
   :caption: ITERATE Workflows

   iterate-cli-options
   iterate-structure-generation

.. toctree::
   :maxdepth: 2
   :caption: PyMOL Visualization

   pymol-cli-options
   pymol-visualization
   pymol-reaction-analysis
   pymol-electronic-structure
   pymol-interaction-analysis

.. toctree::
   :maxdepth: 2
   :caption: Grouper Tool

   grouper-cli-options
   grouper-strategies
   grouper-crest-or-traj-workflow

.. toctree::
   :maxdepth: 2
   :caption: NCIPLOT

   nciplot-tutorial

.. toctree::
   :maxdepth: 2
   :caption: Auxiliary Scripts

   scripts-overview
   scripts-data-management
   scripts-electronic-analysis

********************
 Indices and Tables
********************

-  :ref:`search`
