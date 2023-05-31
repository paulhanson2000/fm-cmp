#!/bin/bash

# Must be run using source: "source this_script.sh"
# Otherwise, the script will execute in a sub-shell, which will load modules in its own ~sub-environment~ and then just exit.
# We want the modules to be loaded to ~our~ environment.

module purge

module load r

module load gsl # For MsCAVIAR
module try-load flexiblas # For MsCAVIAR on Compute Canada
#module try-load liblas # For MsCaviar on BU SCC # TODO, automatically loaded by LAPACK maybe?
module try-load lapack # For MsCaviar on BU SCC
