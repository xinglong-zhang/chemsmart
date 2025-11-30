#########################################
 Installation for Windows Using Git Bash
#########################################

This guide covers installing Chemsmart on Windows using Git Bash.

********************
 Create Environment
********************

#. Install the required software:

   -  **Conda**: Install either Anaconda3 or Miniconda3 for Windows from
      the `conda installation guide
      <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_.

   -  **Git**: Install Git from https://git-scm.com/downloads (Git Bash
      is included).

   -  **Chocolatey**: Install from https://chocolatey.org/install to
      enable ``make`` support.

.. note::

   Ensure the environment variables for Git and Conda are added
   correctly.

2. Open Git Bash and install ``make``:

   .. code:: bash

      choco install make

#. Clone the repository:

   .. code:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

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

      make install

#. For developers, install additional packages:

   .. code:: bash

      make install-dev
