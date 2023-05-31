if(!require("BiocManager")) install.packages("BiocManager")
if(!require("remotes")) install.packages("remotes")
if(!require("quarto")) install.packages("quarto")

if(!require("BiocManager")) install.packages("SeqArray")
if(!require("BiocManager")) install.packages("SNPRelate")
if(!require("hJAM")) remotes::install_github("USCbiostats/hJAM")
if(!require("caviarbf")) remotes::install_bitbucket("Wenan/caviarbf", subdir="caviarbf-r-package/caviarbf/")
