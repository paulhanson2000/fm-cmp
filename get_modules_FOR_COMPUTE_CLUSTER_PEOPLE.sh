#!/bin/bash

# Must be run using source: "source this_script.sh"
# Otherwise, the script will execute in a sub-shell, which will load modules in its own ~sub-environment~ and then just exit.
# We want the modules to be loaded to ~our~ environment.

module purge

module load r/4.2 # Compute Canada
module try-load python/3.11.2  # Compute Canada version
module try-load python3/3.10.5 # BU SCC version

module load plink/1.9b_6.21-x86_64 # For SuSiEx. module load plink/2.00a3.6` (plink 2.0) doesn't seem to work. TODO eventually calc ld myself and remove this
module load htslib

module load boost         # For fGWAS
module load gsl           # For fGWAS and MsCAVIAR
module try-load flexiblas # For MsCAVIAR on Compute Canada
#module try-load liblas   # For MsCaviar on BU SCC # TODO, automatically loaded by LAPACK maybe so no need?
module try-load lapack    # For MsCaviar on BU SCC
#module load cuda          # For PulyFun dependency arrow -> pyarrow # TODO: "module spider arrow" shows that this seems optional. Is having this module activated important in compiling important GPU-accel stuff for pyarrow?
module load gcc           # For PolyFun dependency arrow -> pyarrow
module try-load arrow     # For PolyFun dependency pyarrow on Compute Canada
module try-load thrift    # For PolyFun dependency pyarrow on Compute Canada
module try-load arrow_cpp # For PolyFun dependency pyarrow on BU SCC # TODO: future self trying to install pyarrow on BU SCC: only pyarrow version 11 worked seemingly b/c I was using ComputeCanada's arrow/11.0. BU SCC arrow_cpp versions only go up to 8, so may only work w/ a different pyarrow version.
#module try-load rpy2      # For Polyfun on BU SCC # TODO module present on BU might prevent installation, test, see if just using this module works, or install again w/ pip --ignore-installedd
#module try-load java/1.8 # For Hail TODO: r/4.2.2 requires higher Java version, so this is a conflict. Oh dear. 
