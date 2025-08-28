Installation
============

Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the corresponding installation steps.

If you encounter any problems, please feel free to contact us.

Installation for Linux and MacOS
-------------------------------

Create Environment
^^^^^^^^^^^^^^^^^^^^^^^^

Users are recommended to install conda environments to mange the packages required by this software toolkit. Either Anaconda3 or Miniconda3 may be installed. See more information on conda installation process at https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

*   Once conda has been installed successfully, one can clone this package to a local directory via：

    .. code-block:: console

        git clone https://github.com/xinglong-zhang/chemsmart.git

*   Cd into chemsmart folder. For linux and MacOS systems which support ``make``, users can run：

    .. code-block:: console

        make env

    to create a running environment.

    By default, this will create a conda environment named ``chemsmart``, which installs all the required python packages for this toolkit.

    If conda is not installed, one can run：

    .. code-block:: console

        make env USE_CONDA=false

    or

    .. code-block:: console

        make virtualenv

    to install using virtualenv. It is however recommanded that ``conda`` be used.

.. tip::

    Help options are available by typing ``make help``.

*   Virtual conda environment can be activated via：

    .. code-block:: console

        conda activate chemsmart

Make Installation
^^^^^^^^^^^^^^^^^^^^^^^^

*   One can install the packages and dependencies required for ``chemsmart`` package via：

    .. code-block:: console

        make install

Make Configuration
^^^^^^^^^^^^^^^^^^^^^^^^

*   Next, one can run:

    .. code-block:: console

        make configure

    to sets up the user-specific directory ``~/.chemsmart`` automatically. You will be prompt to enter the paths to g16 and ORCA software, which will then be added automatically. The correct ``conda`` path for the user will also be updated.

    **The configuration also adds the environment variables for chemsmart to the user ``~/.bashrc`` file.**

.. warning::

    ``make configure`` would set up ``~/.chemsmart`` mostly correctly, a user should check the contents in ``~/.chemsmart`` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g., modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one may copy e.g., ``~/.chemsmart/server/SLURM.yaml`` to your own customised server ``~/.chemsmart/server/custom.yaml`` and modify it accordingly, such that the submission becomes ``chemsmart sub -s custom <other commands>``.

*   The ``make configure`` will also add the required paths to the user ``~/.bashrc`` file. User may need to do

    .. code-block:: console

        source ~/.bashrc

    to effect the changes.


*   Once ``make configure`` is done, one can optionally run

    .. code-block:: console

        make fmt

    and

    .. code-block:: console

        make lint

    to format and lint the codes (this should have been handled by the developers). Also optionally, one can run

    .. code-block:: console

        make test

    to make sure that all tests in chemsmart pass.


*   Finally one can clean up by running

    .. code-block:: console

        make clean


Installation for Windows Using Git Bash
-------------------------------

It is recommended that Windows users set up chemsmart using a Bash-based terminal application, such as Git Bash.

Create Environment
^^^^^^^^^^^^^^^^^^^^^^^^

*   Conda is also recommended on Windows to manage the packages required by this software toolkit. Either Anaconda3 or Miniconda3 *for Windows* may be installed. See more information on conda installation process at https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

*   Git is essential for installation. Users can visit https://git-scm.com/downloads to install it (Git Bash will be installed at the same time).

*   To support ``make`` in windows, chocolatey is required to be installed from https://chocolatey.org/install#generic in advance.

.. note::

    Please make sure the environmental variables for Git and Conda are added correctly.

*   Once three apps are installed successfully, one can open Git Bash to install ``make`` via:

    .. code-block:: console

        choco install make

*   Next, one can clone chemsmart package to a local directory via:

    .. code-block:: console

        git clone https://github.com/xinglong-zhang/chemsmart.git

*   Cd into chemsmart folder. Users can run:

    .. code-block:: console

        conda env create -f environment.yml

    to create a running environment from environment.yml file.

    By default, this will create a conda environment named chemsmart, which installs all the required python packages for this toolkit.

*   Virtual conda environment can be activated via：

    .. code-block:: console

        conda activate chemsmart

.. note::

    For windows system users, the ``conda init`` command may need to be run first.


Make Installation
^^^^^^^^^^^^^^^^^^^^^^^^
*   One can install the packages and dependencies required for ``chemsmart`` package via：

    .. code-block:: console

        make install

Make Configuration
^^^^^^^^^^^^^^^^^^^^^^^^

*   Since the Windows system does not come with built-in .zshrc files, users need to run:

    .. code-block:: console

        touch ~/.zshrc

    to create the ``~/.zshrc`` file first.

*   Next, one can run:

    .. code-block:: console

        make configure

    to sets up the user-specific directory ``~/.chemsmart`` automatically. You will be prompt to enter the paths to g16 and ORCA software, which will then be added automatically. The correct ``conda`` path for the user will also be updated.

    **The configuration also adds the environment variables for chemsmart to the user ``~/.zshrc`` file.**

.. warning::

    ``make configure`` would set up ``~/.chemsmart`` mostly correctly, a user should check the contents in ``~/.chemsmart`` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g., modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one may copy e.g., ``~/.chemsmart/server/SLURM.yaml`` to your own customised server ``~/.chemsmart/server/custom.yaml`` and modify it accordingly, such that the submission becomes ``chemsmart sub -s custom <other commands>``.

*   The ``make configure`` will also add the required paths to the user ``~/.zshrc`` file. User may need to do

    .. code-block:: console

        source ~/.zshrc

    to effect the changes.


*   Once ``make configure`` is done, one can optionally run

    .. code-block:: console

        make fmt

    and

    .. code-block:: console

        make lint

    to format and lint the codes (this should have been handled by the developers). Also optionally, one can run

    .. code-block:: console

        make test

    to make sure that all tests in chemsmart pass.


*   Finally one can clean up by running

    .. code-block:: console

        make clean

Installation for Windows Using Ubuntu
-------------------------------

Alternatively, Windows users may install chemsmart via Ubuntu (Windows Subsystem for Linux).

*   Ubuntu can be accessed and downloaded from https://ubuntu.com, and is also available directly on the Microsoft Store.

*   Once Ubuntu is installed, one can proceed to install either Anaconda3 or Miniconda3 *for Linux* via https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

*   Next, the ``git`` and ``make`` can be installed in Ubuntu via:

    .. code-block:: console

        sudo apt install git

    and

    .. code-block:: console

        sudo apt install make

As Ubuntu is fundamentally based on the Linux kernel, one can follow the instructions in **Installation - Installation for Linux and macOS** to complete the following setup of chemsmart.

Installation for HPC Cluster
-------------------------------
As a powerful toolkit, chemsmart can work on any High Performance Computing (HPC) cluster, using same commands as on the local machine to accomplish the same tasks. If you only need to use chemsmart on a HPC cluster, you can install chemsmart solely on the cluster.

*   Before starting the installation, please consult the HPC cluster administrator to confirm the pre-installed software on cluster, the support for required software, and your installation permissions in the specified directory.

*   Since most clusters are based on the Linux system, one can refer to **Installation - Installation for Linux and MacOS** to complete the installation and configuration of chemsmart on your cluster.



Test Installation
-------------------------------

For users of any operating system, installations is deemed successfully if the commands ``make install`` and ``make configure`` do not return any errors. Installation will also create a ``~/.chemsmart`` containing the required files. In addition, the paths for chemsmart packages should be correctly added to the user ``~/.bashrc`` or ``~/.zshrc`` file. Finally, one should be able to run

.. code-block:: console

    chemsmart --help

to get the options for running chemsmart package.