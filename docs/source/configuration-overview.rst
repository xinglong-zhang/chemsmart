########################
 Configuration Overview
########################

After installation, run simple commands to configure Chemsmart on your local machine or cluster.

****************************************
 For Linux, macOS, Ubuntu, and Clusters
****************************************

Run the configuration command:

.. code:: bash

   make configure

This command:

-  Sets up the user-specific directory ``~/.chemsmart`` automatically.
-  Prompts for paths to Gaussian (g16) and ORCA software.
-  Updates the correct conda path for your user.
-  Adds environment variables to your ``~/.bashrc`` file.

After configuration, reload your shell:

.. code:: bash

   source ~/.bashrc

.. warning::

   ``make configure`` sets up ``~/.chemsmart`` with reasonable defaults. You should review the contents to ensure they
   match your server configuration (modules, scratch directories, etc.).

   Depending on your queue system (SLURM, Torque, etc.), copy and customize the server configuration:

   .. code:: bash

      cp ~/.chemsmart/server/SLURM.yaml ~/.chemsmart/server/myserver.yaml

   Then use: ``chemsmart sub -s myserver <other commands>``

*******************
 For Windows Users
*******************

Since Windows does not have a ``.zshrc`` file by default, create it first:

.. code:: bash

   touch ~/.zshrc

Then run:

.. code:: bash

   make configure

This sets up ``~/.chemsmart`` and adds environment variables to ``~/.zshrc``.

After configuration, reload your shell:

.. code:: bash

   source ~/.zshrc

.. warning::

   Review ``~/.chemsmart`` to ensure the settings match your server configuration. Customize server YAML files as
   needed.
