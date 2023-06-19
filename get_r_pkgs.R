install.packages(c("BiocManager", "remotes", "quarto", "data.table", "R.utils", "susieR", "Ckmeans.1d.dp"), repos="https://cloud.r-project.org") # TODO: is Ckmeans.1d.dp (recommended but optional by PolyFun) actually needed?


BiocManager::install(c("SeqArray", "SNPRelate"),
                     update=F, ask=F)

remotes::install_github("USCbiostats/hJAM")
remotes::install_bitbucket("Wenan/caviarbf", subdir="caviarbf-r-package/caviarbf/")
remotes::install_version("susieR", version = "0.11.92", repos = "http://cran.us.r-project.org") # susieR version >0.12 breaking change removal of susie_suff_stat(), breaks PolyFun. So use older version.
