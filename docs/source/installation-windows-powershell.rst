###################################################
 Installation for Windows Using Anaconda PowerShell
###################################################

This guide covers installing Chemsmart on Windows using the Anaconda or Miniconda PowerShell prompt.

********************
 Create Environment
********************

#. Install the required software:

   -  **Conda**: Install either `Anaconda3 <https://www.anaconda.com/products/distribution>`_ or
      `Miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_ for Windows.  The
      Anaconda / Miniconda installer registers the *Anaconda PowerShell Prompt* shortcut in the
      Start Menu.
   -  **Chocolatey**: Install from https://chocolatey.org/install to enable ``make`` support.

.. note::

   Ensure the Anaconda / Miniconda environment variables are added correctly during installation
   (select *"Add Anaconda to my PATH environment variable"* or run the installer with default
   settings and use the dedicated Anaconda Prompt).

2. Open the **Anaconda PowerShell Prompt** (or **Miniconda PowerShell Prompt**) from the Start Menu
   and install ``make``:

   .. code:: powershell

      choco install make

#. Clone the repository:

   .. code:: powershell

      git clone https://github.com/xinglong-zhang/chemsmart.git

#. Create the conda environment:

   .. code:: powershell

      cd chemsmart
      conda env create -f environment.yml

#. Activate the environment:

   .. code:: powershell

      conda activate chemsmart

.. note::

   You may need to run ``conda init powershell`` once to allow conda environments to be activated
   inside PowerShell.  If prompted, close and reopen the Anaconda PowerShell Prompt afterwards.

*******************
 Make Installation
*******************

#. Install the package and dependencies:

   .. code:: powershell

      make install

#. For developers, install additional packages:

   .. code:: powershell

      make install-dev

********************
 Configure Chemsmart
********************

Run the ``make configure`` command to set up the ``~/.chemsmart`` templates and register the
``chemsmart`` command in your PowerShell environment:

.. code:: powershell

   make configure

What ``make configure`` does on Anaconda / Miniconda PowerShell:

1. **Copies templates** — copies the bundled ``.chemsmart`` configuration templates to
   ``~\.chemsmart`` (``%USERPROFILE%\.chemsmart``).
2. **Updates PowerShell profiles** — appends ``$env:PATH`` and ``$env:PYTHONPATH`` entries to:

   - ``~/Documents/WindowsPowerShell/Microsoft.PowerShell_profile.ps1`` (Windows PowerShell 5.x,
     used by Anaconda / Miniconda prompts)
   - ``~/Documents/PowerShell/Microsoft.PowerShell_profile.ps1`` (PowerShell 7+, if installed)

3. **Configures the conda path** — auto-detects your conda installation via ``which conda`` and
   updates the ``~/.chemsmart/server/*.yaml`` files with the correct conda path for your remote HPC
   cluster.  You can override the detected path with the ``--conda-path`` flag:

   .. code:: powershell

      python -m chemsmart config server --conda-path ~/miniconda3

To apply the updated PATH in your **current** PowerShell session without restarting, run:

.. code:: powershell

   . $PROFILE

After this, you can verify the installation:

.. code:: powershell

   chemsmart --version

.. note::

   If the profile scripts already contain a chemsmart section (i.e. ``make configure`` has been run
   before), the profiles will *not* be modified again to avoid duplicate entries.

*****************************
 Troubleshooting
*****************************

**PowerShell execution policy error**

If you see an error like ``cannot be loaded because running scripts is disabled on this system``,
you need to allow local scripts to run.  Open PowerShell **as Administrator** and run:

.. code:: powershell

   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

**conda not found after make configure**

If ``conda`` is not in your PATH when ``make configure`` runs, the server YAML files will not be
updated automatically.  You can configure the conda path manually afterwards:

.. code:: powershell

   chemsmart config server --conda-path ~/miniconda3

Replace ``~/miniconda3`` with the actual path to your conda installation on the **remote HPC
cluster** (e.g. ``/home/username/miniconda3``).

**chemsmart command not found after make configure**

Make sure you have reloaded the PowerShell profile:

.. code:: powershell

   . $PROFILE

Or open a new Anaconda PowerShell Prompt — the profile is sourced automatically on startup.
