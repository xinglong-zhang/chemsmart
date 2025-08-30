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

********************
 Make Configuration
********************

-  Since the Windows system does not come with built-in .zshrc files, users need to run:

   .. code:: console

      touch ~/.zshrc

   to create the ``~/.zshrc`` file first.

-  Next, one can run:

   .. code:: console

      make configure

   to sets up the user-specific directory ``~/.chemsmart`` automatically. You will be prompt to enter the paths to g16
   and ORCA software, which will then be added automatically. The correct ``conda`` path for the user will also be
   updated.

   **The configuration also adds the environment variables for chemsmart to the user ``~/.zshrc`` file.**

.. warning::

   ``make configure`` would set up ``~/.chemsmart`` mostly correctly, a user should check the contents in
   ``~/.chemsmart`` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g.,
   modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one
   may copy e.g., ``~/.chemsmart/server/SLURM.yaml`` to your own customised server ``~/.chemsmart/server/custom.yaml``
   and modify it accordingly, such that the submission becomes ``chemsmart sub -s custom <other commands>``.

-  The ``make configure`` will also add the required paths to the user ``~/.zshrc`` file. User may need to do

   .. code:: console

      source ~/.zshrc

   to effect the changes.
