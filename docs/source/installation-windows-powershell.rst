####################################################
 Installation for Windows Using Anaconda PowerShell
####################################################

This guide covers installing Chemsmart on Windows using the Anaconda or Miniconda PowerShell prompt.

********************
 Create Environment
********************

#. Install the required software:

   -  **Conda**: Install either `Anaconda3 <https://www.anaconda.com/products/distribution>`_ or `Miniconda3
      <https://docs.conda.io/en/latest/miniconda.html>`_ for Windows. The Anaconda / Miniconda installer registers the
      *Anaconda PowerShell Prompt* shortcut in the Start Menu.

   -  **Chocolatey**: Install from https://chocolatey.org/install to enable ``make`` support.

.. note::

   Ensure the Anaconda / Miniconda environment variables are added correctly during installation (select *"Add Anaconda
   to my PATH environment variable"* or run the installer with default settings and use the dedicated Anaconda Prompt).

2. Open the **Anaconda PowerShell Prompt** (or **Miniconda PowerShell Prompt**) from the Start Menu and install
   ``make``:

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

   You may need to run ``conda init powershell`` once to allow conda environments to be activated inside PowerShell. If
   prompted, close and reopen the Anaconda PowerShell Prompt afterwards.

*******************
 Make Installation
*******************

#. Install the package and dependencies:

   .. code:: powershell

      make install

#. For developers, install additional packages:

   .. code:: powershell

      make install-dev

*********************
 Configure Chemsmart
*********************

Run the ``make configure`` command to set up the ``~/.chemsmart`` templates and register the ``chemsmart`` command in
your PowerShell environment:

.. code:: powershell

   make configure

What ``make configure`` does on Anaconda / Miniconda PowerShell:

#. **Copies templates** — copies the bundled ``.chemsmart`` configuration templates to ``~\.chemsmart``
   (``%USERPROFILE%\.chemsmart``).

#. **Updates PowerShell profiles** — writes a ``$env:PYTHONPATH`` entry and a ``Set-Alias -Name chemsmart -Value
   chemsmart.exe`` declaration to:

   -  ``~/Documents/WindowsPowerShell/Microsoft.PowerShell_profile.ps1`` (Windows PowerShell 5.x, used by Anaconda /
      Miniconda prompts)
   -  ``~/Documents/PowerShell/Microsoft.PowerShell_profile.ps1`` (PowerShell 7+, if installed)

   .. note::

      A **``Set-Alias``** declaration (rather than a bare ``$env:PATH`` addition) is used because PowerShell resolves
      aliases **before** external commands in PATH. This means that even if a bare POSIX-style ``chemsmart`` script (no
      ``.exe`` extension) is present somewhere on PATH — for example in ``chemsmart/cli/`` in a development install —
      the alias wins and ``chemsmart.exe`` (which pip places in the conda ``Scripts/`` directory) is called directly.
      This avoids the Windows *"Open with"* popup.

      If you previously ran ``make configure`` and your profile contained old ``$env:PATH`` or ``function chemsmart``
      entries, re-running ``make configure`` will automatically remove the old block and write the new ``Set-Alias``
      declaration in its place.

#. **Prompts for optional software paths** — ``chemsmart config`` interactively asks for the installation folders of
   Gaussian g16, ORCA, and NCIPLOT and writes them into ``~/.chemsmart/server/*.yaml``. Each prompt can be skipped by
   pressing Enter. These prompts are handled by Python (``click.prompt``) and work reliably on all platforms —
   including Conda PowerShell — without any shell-specific syntax.

#. **Configures the conda path** — auto-detects your conda installation via ``conda`` in PATH and updates the
   ``~/.chemsmart/server/*.yaml`` files with the correct conda path for your remote HPC cluster. If conda is not found
   in PATH, a message is logged — add conda to your PATH and re-run ``chemsmart config server``.

#. **Outputs activation instructions** — at the end of ``make configure``, the Makefile prints a reminder to run ``.
   $PROFILE`` so that ``chemsmart`` is immediately available in the current session.

.. note::

   On Windows, ``make configure`` does **not** execute ``. $PROFILE`` automatically — it only prints the instruction.
   You must run ``. $PROFILE`` manually in your Anaconda/Miniconda PowerShell Prompt to load the ``chemsmart`` alias
   into your current session. Alternatively, simply open a new Anaconda PowerShell Prompt.

After ``make configure`` completes, the ``chemsmart`` command is available immediately in any **new** Anaconda
PowerShell Prompt. To activate it in your **current** session without opening a new one, run:

.. code:: powershell

   . $PROFILE

After this, you can verify the installation:

.. code:: powershell

   chemsmart --version

.. note::

   Each time ``make configure`` is run, the installer replaces the previous chemsmart block in the PS profile with the
   latest settings. Re-running configure is therefore safe and idempotent.

*****************
 Troubleshooting
*****************

**PowerShell execution policy error**

If you see an error like ``cannot be loaded because running scripts is disabled on this system``, you need to allow
local scripts to run. Open PowerShell **as Administrator** and run:

.. code:: powershell

   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

**conda not found after make configure**

If ``conda`` is not in your PATH when ``make configure`` runs, the server YAML files will not be updated automatically.
Add conda to your PATH (or activate the conda base environment) and then re-run:

.. code:: powershell

   chemsmart config server

**chemsmart command not found after make configure**

Make sure you have reloaded the PowerShell profile:

.. code:: powershell

   . $PROFILE

Or open a new Anaconda PowerShell Prompt — the profile is sourced automatically on startup.

**chemsmart opens an "Open with" popup instead of running**

This means Windows found a bare POSIX script named ``chemsmart`` (no ``.exe`` extension) in ``$env:PATH`` before
resolving the alias. Re-run ``make configure`` to update the profile with the ``Set-Alias`` declaration that takes
precedence over external files, then reload the profile:

.. code:: powershell

   make configure
   . $PROFILE
