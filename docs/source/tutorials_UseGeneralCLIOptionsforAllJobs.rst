
Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.

Use General CLI Options for All Jobs
------------------------------------

Specify Name of the File
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can specify the name of the file to be created for the job, without file extension, they want to run by using the option ``-l``, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -l custom_job_name opt

    will create input file named ``custom_job_name.com`` instead of the default ``test_opt.com``.

Append String to Input File Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also simply append a string to the base name of the filename supplied, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -a append_string ts

    will create input file named ``test_append_string.com`` instead of the default ``test_ts.com``.

Modify the Charge and Multiplicity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also modify the charge and multiplicity from the CLI, e.g.:

    Modify the charge in ``test.com`` to charge of +1 in the newly created input file ``test_charge.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -c 1 -a charge opt

    Modify the multiplicity in ``test.com`` to multiplicity of 3 in the newly created input file ``test_multiplicity.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -m 3 -a multiplicity opt

    Modify the charge to +1 and multiplicity to 2 in the newly created input ``file test_charge_multiplicity.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -c 1 -m 2 -l test_charge_multiplicity opt

.. tip::

    This can be useful when, e.g., using optimized structure of a neutral closed-shell (charge 0, multiplicity 1) system to run a charged radical ion (e.g., charge +1 and multiplicity 2 in radical cation).

Modify the Functional and Basis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also modify the functional and basis from the CLI to differ from those in project settings, e.g.:

    Modify the functional to ``b3lyp`` in the newly created input file ``test_functional.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -x b3lyp -a functional opt

    Modify the basis to ``6-31G*`` in the newly created input file ``test_basis.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -b "6-31G*" -a basis opt

Specify Additional Optimization Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also specify additional optimization options for opt=() in the route, for example,

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -o maxstep=8,maxsize=12 -a opt_options opt

    will create ``opt=(maxstep=8,maxsize=12)`` as part of the route in the newly created input file ``test_opt_options.com``.

Add in Additional Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also add in additional parameters used in the route, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com --r nosymm -a route_params opt

    will add in ``nosymm`` as part of the route in the newly created input file ``test_route_params.com``.

Select the Particular Structure in file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*   If one has more than one structure in the supplied file for input preparation, one can select the particular structure to perform job on by using the ``-i/--index`` option, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f small.db -i 5 -c 0 -m 1 opt

    will take the 5th structure (1-indexed, as in chemsmart) from ase database file, ``small.db``, to create the input file for geometry optimization.
