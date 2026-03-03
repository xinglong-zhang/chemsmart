#################################################
 Installation for Windows Using Conda PowerShell
#################################################

This guide covers installing Chemsmart on Windows using Conda PowerShell.

********************
 Create Environment
********************

#. Install the required software:

   -  **Conda**: Install either Anaconda3 or Miniconda3 for Windows from the `conda installation guide
      <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_.

   -  **Git**: Install Git from https://git-scm.com/downloads to enable *git* commands.

   -  **Visual Studio Build Tools**: Install from https://visualstudio.microsoft.com/visual-cpp-build-tools/ . During
      installation, select the **Desktop development with C++** workload to ensure the necessary compilers and build
      tools are installed.

.. note::

   Ensure the environment variables for Conda are added correctly.

2. Open Anaconda PowerShell Prompt (Administrator) and install **Chocolatey** following the instructions on the
   Chocolatey website https://chocolatey.org/install .

#. Install ``make``:

   .. code:: bash

      choco install make

#. Clone the repository:

   .. code:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

#. Initialize conda for PowerShell:

   .. code:: bash

      conda init bash

#. Create the conda environment:

   .. code:: bash

      cd chemsmart
      conda env create -f environment.yml

#. Activate the environment:

   .. code:: bash

      conda activate chemsmart

.. note::

   You may need to run ``conda init`` first on Windows.

*******************
 Make Installation
*******************

#. Install the package and dependencies:

   .. code:: bash

      conda run -n chemsmart python -m pip install -e .

#. For developers, need to install additional packages:

   .. code:: bash

      make install-dev
