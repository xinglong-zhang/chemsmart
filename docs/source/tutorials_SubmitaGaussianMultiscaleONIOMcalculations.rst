Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.


Submit a Multiscale (ONIOM) Job Using Gaussian
--------------------------------------------------------------------------

*   Chemsmart helps users easily submit multiscale ONIOM calculations with Gaussian via：

    .. code-block:: console

        chemsmart run gaussian -p <project> -f <input_file> qmmm <option> <value>

    The available options and corresponding meaning are listed in the table below.

            .. list-table:: Available Options
                :header-rows: 1
                :widths: 15 55 30

                *   - options
                    - meaning
                    - example of value
                *   - -hx
                    - high-level functional
                    - mn15
                *   - -hb
                    - high-level basis
                    - def2svp
                *   - -hf
                    - high-level forcefield **(do not recommended)**
                    - --
                *   - -mx
                    - medium-level functional
                    - pbe
                *   - -mb
                    - medium-level basis
                    - 6-31g(d,p)/auto
                *   - -mf
                    - medium-level forcefield
                    - amber
                *   - -lx
                    - low-level functional
                    - --
                *   - -lb
                    - low-level basis
                    - --
                *   - -lf
                    - low-level forcefield
                    - amber=hardfirst
                *   - -cr
                    - charge of real system
                    - 0
                *   - -cm
                    - multiplicity of real system
                    - 1
                *   - -ci
                    - charge of intermediate layer
                    - 0
                *   - -mi
                    - multiplicity of intermediate layer
                    - 1
                *   - -cm
                    - charge of model system
                    - 0
                *   - -mm
                    - multiplicity of model system
                    - 1
                *   - -ha
                    - indices of high-level atoms
                    - '[5,6,7]'
                *   - -ma
                    - indices of medium-level atoms
                    - '[1,2,3]'
                *   - -la
                    - indices of low-level atoms
                    - '[4,8,10,12,11]'
                *   - -ba
                    - indices of bonded atoms
                    - '[(4,5), (1,4)]'
                *   - -s
                    - value of scale factors
                    - 0.738 (default 1.0)

    This is an example of calling command:

    .. code-block:: console

        chemsmart run gaussian -p test_qmmm -f methanol_in_water.com qmmm -hx mn15 -hb def2svp -mx hf -mb def2svp -lf uff -cr 0 -mr 1 -ci 0 -mi 1 -cm 0 -mm 1 -ha '[5,6,7]' -ma 4 -la '[1,2,3]' -ba '[(4,5), (1,4)]'

.. note::

    1.	For charge and multiplicity, ``-cr`` and ``-cm`` options are mandatory, while ``-ci``, ``-mi``, ``-cm``, and ``-mm`` are optional. If the charge or multiplicity of a layer is omitted by user, the program will be assigned them in the default way. Please be careful when assigning charge and multiplicity, as the model system is embedded in the intermediate layer.

    2.  For ``-ha``, ``-ma``, and ``-la``, please specify them in string format instead of lists, i.e., using '[1,2,3]' instead of [1,2,3]. Here ``-ha`` is mandatory, while ``-ma`` and ``-la`` are optional. If only ``-ha`` is given, the program will assign the rest atoms as low-level atoms automatically.

    3.  ``-ba`` (bonded atoms) is mandatory, which indicates the covalent bond to be cut for dividing different partition layers. Please specify your bonded atoms as the string format of a list of tuples, where the number in each tuple specifying the indices of atoms forming the bond.
