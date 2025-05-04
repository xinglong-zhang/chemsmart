
Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.


Submit a Transition State and Reaction Pathway Analysis Job Using Gaussian
--------------------------------------------------------------------------

Transition State Modredundant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit transition state modredundant job (frozen coordinates optimization), do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> modred -c <list_of_coords_to_constraint>

    For example, to submit a modredundant job with constraints on bond between atom 4 and atom 17 and on bond between atom 9 and atom 10, do:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com modred -c [[4,17],[9,10]]

Transition State Search
^^^^^^^^^^^^^^^^^^^^^^^

*   To submit transition state search job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> ts

Intrinsic Reaction Coordinate (IRC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit intrinsic reaction coordinate (IRC) job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> irc

Relaxed Potential Energy Surface (PES) Scan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit relaxed potential energy surface (PES) scan, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> scan -c <list_of_coords_to_constraint> -s <scan_step_size> -n <num_scan_steps>

*   For example, to submit the PES scan job with along bond between atom 4 and atom 17 for 10 steps with 0.1Å increment per step:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com scan -c [[4,17]] -s 0.1 -n 10
