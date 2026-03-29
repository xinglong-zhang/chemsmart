#########################################
 Installation for Windows Using Git Bash
#########################################

This guide covers installing Chemsmart on Windows using Git Bash.

********************
 Create Environment
********************

#. Install the required software:

   -  **Conda**: Install either Anaconda3 or Miniconda3 for Windows from the `conda installation guide
      <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_.
   -  **Git**: Install Git from https://git-scm.com/downloads (Git Bash is included).
   -  **Chocolatey**: Install from https://chocolatey.org/install to enable ``make`` support.

.. note::

   Ensure the environment variables for Git and Conda are added correctly.

2. Open Git Bash and install ``make``:

   .. code:: bash

      choco install make

#. Clone the repository:

   .. code:: bash

      git clone https://github.com/xinglong-zhang/chemsmart.git

#. Initialize conda for Git Bash in the Anaconda prompt:

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

      make install

#. For developers, install additional packages:

   .. code:: bash

      make install-dev

********************
 Configure Chemsmart
********************

Run the ``make configure`` command to set up the ``~/.chemsmart`` templates and register the
``chemsmart`` command in your Git Bash environment:

.. code:: bash

   make configure

What ``make configure`` does on Git Bash:

1. **Copies templates** — copies the bundled ``.chemsmart`` configuration templates to
   ``~/.chemsmart``.
2. **Updates** ``~/.bashrc`` — appends ``export PATH=...`` and ``export PYTHONPATH=...`` lines so
   that the ``chemsmart`` command is available in new Git Bash sessions.
3. **Configures the conda path** — auto-detects your conda installation via ``which conda`` and
   updates the ``~/.chemsmart/server/*.yaml`` files with the correct conda path for your remote HPC
   cluster.  You can override the detected path with the ``--conda-path`` flag:

   .. code:: bash

      chemsmart config server --conda-path ~/miniconda3

To apply the updated ``PATH`` in your **current** Git Bash session without restarting, run:

.. code:: bash

   source ~/.bashrc

After this, you can verify the installation:

.. code:: bash

   chemsmart --version

.. note::

   If ``~/.bashrc`` already contains a chemsmart section (i.e. ``make configure`` has been run
   before), it will *not* be modified again to avoid duplicate entries.

.. tip::

   If you prefer to use the Anaconda or Miniconda PowerShell Prompt instead of Git Bash, see
   :doc:`installation-windows-powershell`.
