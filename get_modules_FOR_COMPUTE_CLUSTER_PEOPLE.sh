#!/bin/bash

# Must be run using source: "source this_script.sh"
# Otherwise, the script will execute in a sub-shell, which will load modules in its own ~sub-environment~ and then just exit.
# We want the modules to be loaded to ~our~ environment.

module purge

module load r

module load scipy-stack            # For SuSiEx TODO eventually calc ld myself and remove this
module load plink/1.9b_6.21-x86_64 # For SuSiEx. module load plink/2.00a3.6` (plink 2.0) doesn't seem to work. TODO eventually calc ld myself and remove this

module load gsl           # For MsCAVIAR
module try-load flexiblas # For MsCAVIAR on Compute Canada
#module try-load liblas   # For MsCaviar on BU SCC # TODO, automatically loaded by LAPACK maybe so no need?
module try-load lapack    # For MsCaviar on BU SCC
