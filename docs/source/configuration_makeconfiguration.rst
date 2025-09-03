Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the
corresponding configuration steps.

If you encounter any problems, please feel free to contact us.

####################
 Make Configuration
####################

After installation, one can run simple commands to configure chemsmart in local machine or cluster.

********************************************
 For Linux, MacOS, Ubuntu and Cluster users
********************************************

-  one can run:

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

*******************
 For Windows users
*******************

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
