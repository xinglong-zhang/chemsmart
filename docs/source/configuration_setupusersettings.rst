Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the
corresponding configuration steps.

If you encounter any problems, please feel free to contact us.

######################
 Set up User Settings
######################

Once the test is passed, users can modify the ``~/.chemsmart`` file to set up custom settings.

.. note::

   A user can modify the contents in ``~/.chemsmart`` files freely without affecting or needing to know the
   ``chemsmart`` source code.

-  The ``~/.chemsmart/usersettings.yaml`` file contains informations such as project number or account number that are
   required in a typical submission script that specifies the account for use at some HPC servers. It can also contain
   options specifying user's email to inform user of the job start and job end once a job is submitted. If more features
   are needed, please submit a request via `Issues`. A typical ``~/.chemsmart/usersettings.yaml`` file looks like this:

      .. code:: console

         PROJECT: 1234567  # alias ACCOUNT FOR SLURM
         EMAIL: abc@gmail.com
