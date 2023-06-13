install.packages(c("BiocManager", "remotes", "quarto", "data.table", "R.utils"), repos="https://cloud.r-project.org")

BiocManager::install(c("SeqArray", "SNPRelate"),
                     update=F, ask=F)

remotes::install_github("USCbiostats/hJAM")
remotes::install_bitbucket("Wenan/caviarbf", subdir="caviarbf-r-package/caviarbf/")
