
Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.

Submit a Geometry Optimization Job Using Gaussian
-------------------------------------------------
With setup completed, one is able to run different Gaussian jobs via command-line interface (CLI).

Basic Geometry Optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit (and run) a geometry optimization job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt

    where ``<server_name>`` is the one of the servers specified in ``~/.chemsmart/server/*.yaml`` files, without ``.yaml`` extension; ``<project>`` is one of the project settings specified in ``~/.chemsmart/gaussian/*.yaml`` files, without ``.yaml`` extension; and ``<input_file>`` is an input file the user wishes to run job on. Note that this input file can be any format, such as ``.xyz``, Gaussian ``.com``, ``.gjf`` or ``.log`` file or ORCA ``.inp`` or ``.out`` file.

Geometry Optimization Job with Frozen Atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   If one wants to submit geometry optimization with frozen atoms in the molecule (such as https://www.researchgate.net/post/Freezing-atoms-in-gaussian-how-to-do-it), one can do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt -f <indices_of_atoms_to_freeze>

    or example, to submit the geometry optimization job with atoms numbered 1 to 10 frozen, one can do

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com opt -f 1-10

.. Note::

    1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for moleculare are 1-indexed.

Optimize a Fixed Number of Lowest Energy Conformers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize a fixed number of lowest energy conformers, ``n_conformers_to_opt``, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt -n <n_conformers_to_opt>

.. note::

    If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .xyz file. In fact, whenever .xyz file is used as input, the charge and multiplicity should be specified via ``-c <charge> -m <multiplicity`` via CLI.

Optimize Unique Structure from an MD Trajectory File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize unique structure from an MD trajectory file, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> saopt

.. note::

    This optimizes all the unique structures available in the MD trajectory ``<input_file>``. Typically, the ``<input_file>`` is a list of all structures on an MD trajectory obtained by ASE MD run and named ``md.traj``. **(TODO: this method is not properly implemented in chemsmart yet.)**

Optimize a Fixed Number of Lowest Energy Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize a fixed number of lowest energy structures, ``<num_structures_to_opt>``, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> saopt -n <n_conformers_to_opt>


    If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .traj file.

    Two grouper types for determining/clustering unique structures are available from CLI option ``-g``:

    1. Sequential grouper (default), selected by option value of seq, which sequentially checks for unique structures in a given list of md structures, and

    2. Self-consistent grouper, selected by option value of sc, which self-consistently checks for unique structures in a given list of md structures using the reverse Cuthill–McKee algorithm for structure clustering. By default, only the last 0.1 proportion of the structures of the md.traj file is considered. This can be changed via cli option ``-x <proportion_structures_to_use>``.

    For example, to consider the last 20% of the structures in md.traj trajectory file, then uses Sequential grouper to group those structures into unique structures and run the 10 lowest energy structures from the list of unique structures found by the grouper:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f imd.traj saopt -x 0.2 -n 10 -g seq

run opt or modred or ts conformers from CREST run output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To run opt or modred or ts conformers from CREST run output, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt

    or

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j modred -c [1,4]

    or

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j ts

    respectively

.. note::

    This optimizes all the conformers available in the ``<input_file>``. Typically, the ``<input_file>`` is a list of all conformers obtained by CREST program and named ``crest_conformers.xyz``.

