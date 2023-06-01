if(!require("BiocManager")) install.packages("BiocManager")
if(!require("remotes"))     install.packages("remotes")
if(!require("quarto"))      install.packages("quarto")

if(!require("SeqArray"))  BiocManager::install("SeqArray")
if(!require("SNPRelate")) BiocManager::install("SNPRelate")

if(!require("hJAM"))      remotes::install_github("USCbiostats/hJAM")
if(!require("caviarbf"))  remotes::install_bitbucket("Wenan/caviarbf", subdir="caviarbf-r-package/caviarbf/")
