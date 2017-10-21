# Install R pacakges
#
# This is to be replaced by conda later.
# a few of the packages need to be build for bioconda

# install core bioconductor packages
# Bioc 3.5
# this can take a while...
source("https://bioconductor.org/biocLite.R")
biocLite()

# install core dependencies
biocLite(c("Biobase", "edgeR", "DESeq2", "limma", "SCnorm",
           "scran", "BPSC", "monocle", "ROTS", "scde", "MAST"))

# scde:
#  There were issues installing scde current release
#  need to install scde1.99.2 along with an earlier version of flexmix (2.13-13)
#  link address to scde1.99.2:; https://github.com/hms-dbmi/scde/archive/1.99.2.tar.gz
#  require(devtools)
#  install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

# scran:
# this version of scran is dependent on singleCellExperiment objects
# introduced by scater 1.5
# since Bioc 3.5 has scater 1.4, we hold off on running scran

# BPSC and SCnorm not available for R3.4
# so install from github
devtools::install_github("BPSC","nghiavtr")
devtools::install_github("rhondabacher/SCnorm")

# ashbun package
# install from github repo
# developing locally
devtools::install_github("jhsiao999/ashbun", ref = "benchmarkMethods")

# install data packages
devtools::install_github("jhsiao999/singleCellRNASeqMouseZeiselBrain")
devtools::install_github("jhsiao999/singleCellRNASeqHumanTungiPSC")
devtools::install_github("jhsiao999/singleCellRNASeqMouseKleinESC")

# simulated data and output are stored in /scratch/midway2/joycehsiao
