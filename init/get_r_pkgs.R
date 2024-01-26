install.packages(c("BiocManager", "remotes", "quarto", "kableExtra", "yaml", "jsonlite", "R.utils", "Ckmeans.1d.dp"), repos="https://cloud.r-project.org") # TODO: is Ckmeans.1d.dp (recommended but optional by PolyFun) actually needed?

BiocManager::install(c("SeqArray", "SNPRelate"), update=F, ask=F, force=T)

remotes::install_github("paulhanson2000/refpanelutil")

remotes::install_github("mailund/tailr")
remotes::install_github("USCbiostats/hJAM")
remotes::install_bitbucket("Wenan/caviarbf", subdir="caviarbf-r-package/caviarbf/")
remotes::install_version("susieR", version = "0.11.92", repos = "http://cran.us.r-project.org") # susieR version >0.12 breaking change removal of susie_suff_stat(), breaks PolyFun. So use older version.

# Multithread support for data.table on Mac requires fiddling
is_mac <- "Darwin" == Sys.info()["sysname"]
is_arm <- "arm64"  == Sys.info()["machine"]
if(is_mac) {
  if(!dir.exists("/opt/homebrew/opt/gcc") || !dir.exists("/opt/homebrew/opt/libomp")) {
    warning(
      "Your computer is a Mac. The data.table R package requries extra steps to use multithreading on Macs.\n",
      "Please install homebrew by entering the following command in a shell:\n",
      "/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"\n",
      "Then, also enter \"brew install gcc\" and \"brew install libomp\"."
    )
    response <- readline("\nAlternatively, you could install the slower, single-threaded version. Would you prefer that? [y/N]: ")
    if(tolower(response) %in% c("y","yes")) {
      install.packages("data.table", repos="https://cloud.r-project.org")
    } else cat("Please run this script again once you have installed the dependencies mentioned above.\n")

  } else { # User has dependencies

    if(is_arm) {
      writeLines(c("LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp", "CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp"), "Makevars")
    } else { # x86_64
      writeLines(c("CPPFLAGS += -Xclang -fopenmp", "LDFLAGS += -lomp"), "Makevars")
    }
    Sys.setenv(R_MAKEVARS_USER=paste0(getwd(), "/Makevars"))

    if(system.file(package="data.table")!="") {remove.packages("data.table")} # Start fresh
    install.packages("data.table", type="source", repos="https://cloud.r-project.org") # Install from source. Uses the Makevars file we just made. 
  }
}
