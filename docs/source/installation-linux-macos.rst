##################################
 Installation for Linux and macOS
##################################

This guide covers installing Chemsmart on Linux and macOS systems.

********************
 Create Environment
********************

We recommend using conda to manage the packages required by Chemsmart.
Either Anaconda3 or Miniconda3 may be installed. See the
`conda installation guide <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_
for more information.

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

2. Change to the chemsmart directory and create the environment:

   .. code-block:: bash

      cd chemsmart
      make env

   This creates a conda environment named ``chemsmart`` with all required
   Python packages.

   If conda is not installed, you can use virtualenv instead:

   .. code-block:: bash

      make env USE_CONDA=false

   or:

   .. code-block:: bash

      make virtualenv

   However, using conda is recommended.

.. tip::

   Run ``make help`` to see all available make targets.

3. Activate the conda environment:

   .. code-block:: bash

      conda activate chemsmart

*******************
 Make Installation
*******************

1. Install the package and dependencies:

   .. code-block:: bash

      make install

2. For developers, install additional packages (dev, test, docs):

   .. code-block:: bash

      make install-dev
