Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the
corresponding installation steps.

If you encounter any problems, please feel free to contact us.

#########################################
 Installation for Windows Using Git Bash
#########################################

It is recommended that Windows users set up chemsmart using a Bash-based terminal application, such as Git Bash.

********************
 Create Environment
********************

-  Conda is also recommended on Windows to manage the packages required by this software toolkit. Either Anaconda3 or
   Miniconda3 *for Windows* may be installed. See more information on conda installation process at
   https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

-  Git is essential for installation. Users can visit https://git-scm.com/downloads to install it (Git Bash will be
   installed at the same time).

-  To support ``make`` in windows, chocolatey is required to be installed from https://chocolatey.org/install#generic in
   advance.

.. note::

   Please make sure the environmental variables for Git and Conda are added correctly.

-  Once three apps are installed successfully, one can open Git Bash to install ``make`` via:

   .. code:: console

      choco install make

-  Next, one can clone chemsmart package to a local directory via:

   .. code:: console

      git clone https://github.com/xinglong-zhang/chemsmart.git

-  Cd into chemsmart folder. Users can run:

   .. code:: console

      conda env create -f environment.yml

   to create a running environment from environment.yml file.

   By default, this will create a conda environment named chemsmart, which installs all the required python packages for
   this toolkit.

-  Virtual conda environment can be activated via：

   .. code:: console

      conda activate chemsmart

.. note::

   For windows system users, the ``conda init`` command may need to be run first.

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
