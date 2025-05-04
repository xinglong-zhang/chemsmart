
Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.

Submit a Energy and Thermodynamic Analysis Job Using Gaussian
-------------------------------------------------------------

Single Point Job
^^^^^^^^^^^^^^^^

*   To submit single point job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp

    For single-point job that user wants to test which uses different solvent model and id from that specified in ``<project>``, one can do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp -sm <user_solvent_model> -si <user_solvent_id>

    to specify a different solvent model ``<user_solvent_model>`` and solvent ``<user_solvent_id>``.

Distortion-Interaction/Activation-Strain (DI-AS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To run distortion-interaction/activation-strain (DI-AS) job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <irc_output_file_for_dias> dias -i <indices_of_any_one_fragment> -n <number_of_every_n_step_along_irc_to_run>

    For example to run DI-AS job for fragment 1 with atoms numbered from 5-17 at every 10 steps along the irc.log file:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f irc.log dias -i 5-17 -n 10
