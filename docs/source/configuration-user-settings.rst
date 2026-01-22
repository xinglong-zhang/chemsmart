###############
 User Settings
###############

Configure user-specific settings in the ``~/.chemsmart`` directory.

.. note::

   You can modify files in ``~/.chemsmart`` freely without affecting or needing to understand the Chemsmart source code.

********************
 User Settings File
********************

The ``~/.chemsmart/usersettings.yaml`` file contains:

-  Project or account numbers required for HPC job submissions.
-  Email address for job notifications.

Example configuration:

.. code:: yaml

   PROJECT: 1234567  # Account number for SLURM
   EMAIL: user@example.com

To request additional features, submit an issue on the `GitHub repository
<https://github.com/xinglong-zhang/chemsmart/issues>`_.
