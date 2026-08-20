############################
 Auxiliary Scripts Overview
############################

CHEMSMART includes auxiliary Python scripts for data management and electronic structure analysis. These scripts are
standalone utilities that complement the main CLI. They ship as package
modules without console entry points, so run them as
``python -m chemsmart.scripts.<name>`` (or directly from a source
checkout); a plain ``pip install`` does not place them on ``PATH``.

Available script categories:

-  **File Management**: Organize and convert computational chemistry files
-  **Electronic Structure Analysis**: Extract and analyze properties from output files

See the following pages for detailed documentation:

-  :doc:`scripts-data-management`
-  :doc:`scripts-electronic-analysis`
