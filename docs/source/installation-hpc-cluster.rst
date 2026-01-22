##############################
 Installation for HPC Cluster
##############################

Chemsmart can run on any High Performance Computing (HPC) cluster, using the same commands as on a local machine. If you
only need to use Chemsmart on an HPC cluster, you can install it directly there.

***************
 Prerequisites
***************

Before starting:

-  Consult your HPC cluster administrator about pre-installed software and your installation permissions.
-  Confirm the supported queue system (SLURM, Torque, PBS, etc.).
-  Ensure that a suitable C/C++ compiler is available (e.g. via environment modules such as ``gcc``), as some
   dependencies require compilation.

**************
 Installation
**************

Since most HPC clusters run Linux, follow the instructions in :doc:`installation-linux-macos` to install Chemsmart on
your cluster.

Before running ``make env``, you may need to load an appropriate compiler module (e.g. ``gcc``) to ensure that C++
extensions such as ``pyvoro`` can be built successfully.
