

Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.


Submit a Custom Job Using Gaussian
-------------------------------------------

Run a Job with Pre-prepared Gaussian Input File Directly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*   If a user wants to run a job with pre-prepared Gaussian input file directly, one can run the job directly using:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> com

Run a Custom Job
^^^^^^^^^^^^^^^^

*   Generally, if a user wants to run job that is currently not present in our package, one can run custom job using:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> userjob -r <user_defined_gaussian_route> -a <appending_information_as_string_at_the_end_of_input_file_after_coordinates_specification>

    For example, to create an input file named ``user_defined_job.com`` with user-specified route ``mnr functional/basis solvent`` etc and ``B 1 2 F\nA 1 2 3 F`` at the end of the input file after the specification of coordinates, run

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -l user_defined_job userjob -r 'mnr functional/basis solvent etc' -a 'B 1 2 F\nA 1 2 3 F'
