##################################
 Installation for Linux and macOS
##################################

This guide covers installing Chemsmart on Linux and macOS systems.

********************
 Create Environment
********************

We recommend using conda to manage the packages required by Chemsmart. Either Anaconda3 or Miniconda3 may be installed.
See the `conda installation guide <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_ for
more information.

#. Clone the repository:

   .. code:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

#. Change to the chemsmart directory and create the environment:

   .. code:: bash

      cd chemsmart
      make env

   This creates a conda environment named ``chemsmart`` with all required Python packages.

   If conda is not installed, you can use virtualenv instead:

   .. code:: bash

      make env USE_CONDA=false

   or:

   .. code:: bash

      make virtualenv

   However, using conda is recommended.

.. tip::

   Run ``make help`` to see all available make targets.

3. Activate the conda environment:

   .. code:: bash

      conda activate chemsmart

*******************
 Make Installation
*******************

#. Install the package and dependencies:

   .. code:: bash

      make install

#. For developers, install additional packages (dev, test, docs):

   .. code:: bash

      make install-dev
