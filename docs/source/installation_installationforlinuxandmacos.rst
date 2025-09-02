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
