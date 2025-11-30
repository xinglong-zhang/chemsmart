#########################################
 Installation for Windows Using Git Bash
#########################################

This guide covers installing Chemsmart on Windows using Git Bash.

********************
 Create Environment
********************

1. Install the required software:

   - **Conda**: Install either Anaconda3 or Miniconda3 for Windows from the
     `conda installation guide <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_.

   - **Git**: Install Git from https://git-scm.com/downloads (Git Bash is
     included).

   - **Chocolatey**: Install from https://chocolatey.org/install to enable
     ``make`` support.

.. note::

   Ensure the environment variables for Git and Conda are added correctly.

2. Open Git Bash and install ``make``:

   .. code-block:: bash

      choco install make

3. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

4. Create the conda environment:

   .. code-block:: bash

      cd chemsmart
      conda env create -f environment.yml

5. Activate the environment:

   .. code-block:: bash

      conda activate chemsmart

.. note::

   You may need to run ``conda init`` first on Windows.

*******************
 Make Installation
*******************

1. Install the package and dependencies:

   .. code-block:: bash

      make install

2. For developers, install additional packages:

   .. code-block:: bash

      make install-dev
