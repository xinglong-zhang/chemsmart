
Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.

Submit a Electronic Structure Analysis Job Using Gaussian
---------------------------------------------------------

RESP Charges Fitting
^^^^^^^^^^^^^^^^^^^^

*   To submit RESP charges fitting job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> resp

.. note::

    This creates an input file with fix route for RESP job:
    ``HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)``

Non-covalent Interaction Job
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit non-covalent interaction job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> nci




