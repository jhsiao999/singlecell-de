# Download lists of mitochrondrial genes
#
# Download lists of mitochrondrial genes from here: https://www.broadinstitute.org/scientific-community/science/programs/metabolic-disease-program/publications/mitocarta/mitocarta-in-0

library(gdata)

mouse <- read.xls("data/Mouse.MitoCarta2.0.xls",
                  sheet = 2, header = TRUE, stringsAsFactors=FALSE)
saveRDS(with(mouse, data.frame(Symbol, EnsemblGeneID)), file = "data/mitogenes.mouse.rds")

human <- read.xls("data/Human.MitoCarta2.0.xls",
                   sheet = 2, header = TRUE, stringsAsFactors=FALSE)
saveRDS(with(human, data.frame(Symbol, EnsemblGeneID)), file = "data/mitogenes.human.rds")
