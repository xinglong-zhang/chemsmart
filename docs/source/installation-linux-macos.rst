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

*********************
 Configure Chemsmart
*********************

Run the ``make configure`` command to set up the ``~/.chemsmart`` templates and register the ``chemsmart`` command in
your shell environment:

.. code:: bash

   make configure

What ``make configure`` does on Linux and macOS:

#. **Copies templates** — copies the bundled ``.chemsmart`` configuration templates to ``~/.chemsmart``.

#. **Updates your shell rc file** — appends ``export PATH=...`` and ``export PYTHONPATH=...`` lines to ``~/.bashrc``
   (bash) or ``~/.zshrc`` (zsh) so that the ``chemsmart`` command is available in new terminal sessions.

#. **Prompts for optional software paths** — ``chemsmart config`` asks for the installation folders of Gaussian g16,
   ORCA, and NCIPLOT and writes them into ``~/.chemsmart/server/*.yaml``. Each prompt can be skipped by pressing Enter.

   .. note::

      When ``make configure`` runs ``chemsmart config`` via ``conda run``, stdin is not connected to your terminal and
      the path prompts are **automatically skipped** with an informational message. To configure these paths
      interactively, activate the ``chemsmart`` conda environment and run ``chemsmart config`` directly in your terminal
      after ``make configure`` completes:

      .. code:: bash

         conda activate chemsmart
         chemsmart config

#. **Configures the conda path** — auto-detects your conda installation via ``which conda`` and updates the
   ``~/.chemsmart/server/*.yaml`` files with the correct conda path for your remote HPC cluster. If conda is not found
   in PATH, a message is logged — add conda to your PATH and re-run ``chemsmart config server``.

#. **Automatically sources your shell config** — after writing the ``export`` lines, ``make configure`` sources every
   shell rc file that exists (``~/.bashrc``, ``~/.zshrc``, ``~/.profile``) so that ``chemsmart`` is active for the rest
   of the current make session, regardless of which file chemsmart wrote to.

.. note::

   ``make configure`` sources the rc files inside its own sub-shell process. Because a child process cannot export
   environment changes back to its parent, your **current** terminal session is *not* automatically updated. You must
   run the ``source`` command below in your own terminal, or open a new terminal window.

After ``make configure`` completes, the ``chemsmart`` command is available immediately in any **new** terminal. To
activate it in your **current** terminal without opening a new one, run the command for your shell:

.. code:: bash

   source ~/.bashrc   # bash
   # or
   source ~/.zshrc    # zsh
   # or
   source ~/.profile  # sh / other

Then verify the installation:

.. code:: bash

   chemsmart --version

.. note::

   If ``~/.bashrc`` (or ``~/.zshrc``) already contains a chemsmart section (i.e. ``make configure`` has been run
   before), it will *not* be modified again to avoid duplicate entries.
