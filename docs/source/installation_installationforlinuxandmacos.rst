Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the
corresponding installation steps.

If you encounter any problems, please feel free to contact us.

##################################
 Installation for Linux and MacOS
##################################

********************
 Create Environment
********************

Users are recommended to install conda environments to mange the packages required by this software toolkit. Either
Anaconda3 or Miniconda3 may be installed. See more information on conda installation process at
https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

-  Once conda has been installed successfully, one can clone this package to a local directory via：

   .. code:: console

      git clone https://github.com/xinglong-zhang/chemsmart.git

-  Cd into chemsmart folder. For linux and MacOS systems which support ``make``, users can run：

   .. code:: console

      make env

   to create a running environment.

   By default, this will create a conda environment named ``chemsmart``, which installs all the required python packages
   for this toolkit.

   If conda is not installed, one can run：

   .. code:: console

      make env USE_CONDA=false

   or

   .. code:: console

      make virtualenv

   to install using virtualenv. It is however recommanded that ``conda`` be used.

.. tip::

   Help options are available by typing ``make help``.

-  Virtual conda environment can be activated via：

   .. code:: console

      conda activate chemsmart

*******************
 Make Installation
*******************

-  One can install the packages and dependencies required for ``chemsmart`` package via：

   .. code:: console

      make install

-  For developers, one needs to install additional packages and dependencies (dev, test, docs targets in pyproject.toml)
   required for ``chemsmart`` package via：

   .. code:: console

      make install-dev

********************
 Make Configuration
********************

-  Next, one can run:

   .. code:: console

      make configure

   to sets up the user-specific directory ``~/.chemsmart`` automatically. You will be prompt to enter the paths to g16
   and ORCA software, which will then be added automatically. The correct ``conda`` path for the user will also be
   updated.

   **The configuration also adds the environment variables for chemsmart to the user ``~/.bashrc`` file.**

.. warning::

   ``make configure`` would set up ``~/.chemsmart`` mostly correctly, a user should check the contents in
   ``~/.chemsmart`` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g.,
   modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one
   may copy e.g., ``~/.chemsmart/server/SLURM.yaml`` to your own customised server ``~/.chemsmart/server/custom.yaml``
   and modify it accordingly, such that the submission becomes ``chemsmart sub -s custom <other commands>``.

-  The ``make configure`` will also add the required paths to the user ``~/.bashrc`` file. User may need to do

   .. code:: console

      source ~/.bashrc

   to effect the changes.

..
   -  Once ``make configure`` is done, one can optionally run

..
   .. code:: console

..
   make fmt

..
   and

..
   .. code:: console

..
   make lint

..
   to format and lint the codes (this should have been handled by the developers). Also optionally, one can run

..
   .. code:: console

..
   make test

..
   to make sure that all tests in chemsmart pass.

..
   -  Finally one can clean up by running

..
   .. code:: console

..
   make clean
