# This script parses data for plots to put them into a more manageable form

library(glow)
library(dplyr)
library(data.table)
library(arrow)
library(qs)
library(stringr)

# COVID example
library(sf)
library(USAboundaries) # install.packages("USAboundaries", repos = "https://ropensci.r-universe.dev", type = "source")
library(USAboundariesData) # install.packages("USAboundariesData", repos = "https://ropensci.r-universe.dev", type = "source")

# Volcano
library(minfi) # BiocManager::install(c("minfi", "methylationArrayAnalysis"))
library(methylationArrayAnalysis)

tempdir <- tempdir()
outdir <- "plot_data"
dir.create(outdir, showWarnings=FALSE)

nt <- 8

# Galaxy #######################################################################

# Some URIs take longer than 60 seconds 
options(timeout=600)
files <- c("http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_1584380076484244352_2200921635402776448.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_2200921875920933120_3650804325670415744.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_2851858288640_1584379458008952960.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_3650805523966057472_4475721411269270528.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_4475722064104327936_5502601461277677696.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_5502601873595430784_5933051501826387072.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_5933051914143228928_6714230117939284352.csv.gz",
           "http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_6714230465835878784_6917528443525529728.csv.gz")

stars <- lapply(files, function(f) {
  data.table::fread(f, select = c("parallax", "dec", "ra", "astrometric_pseudo_colour", "lum_val"))
}) %>% rbindlist %>% as.data.frame
stars2 <- stars[complete.cases(stars),]
stars2 <- stars2 %>% arrange_all

qsave(stars2, "plot_data/gaia_stars.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads = 8)


# COVID #######################################################################

cov_cases <- fread("https://usafactsstatic.blob.core.windows.net/public/data/covid-19/covid_confirmed_usafacts.csv", data.table=F)

cov_cases <- cov_cases %>% 
  dplyr::transmute(state_abbr = State, name = stringr::str_trim(gsub(" County", "" ,`County Name`)), total = `2021-08-16`)

plot_states <- setdiff(USAboundaries::us_states()$state_name, c("Hawaii", "Northern Mariana Islands", "Puerto Rico", "Alaska", "American Samoa", "Guam", "Virgin Islands"))
county <- us_counties(map_date = NULL, resolution = "high", states = plot_states)
state <- us_states(map_date = NULL, resolution = "high", states = plot_states)

centroids <- st_centroid(county)
centroids <- cbind(centroids, st_coordinates(centroids$geometry))
centroids <- left_join(centroids, cov_cases, by = c("name", "state_abbr"))
centroids$total[is.na(centroids$total)] <- 0
centroids <- filter(centroids, total > 0)

qsave(list(centroids=centroids, county=county, state=state), "plot_data/covid_confirmed_usafacts.qs", preset = "custom", algorithm = "zstd", compress_level = 22)


# Volcano ##########################################

# https://www.bioconductor.org/packages/release/workflows/tests/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#differential-methylation-analysis

dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
rgSet <- read.metharray.exp(targets=targets)
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
detP <- detectionP(rgSet)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]
mSetSq <- preprocessQuantile(rgSet, quantileNormalize=FALSE)
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mVals <- getM(mSetSqFlt)
cellType <- factor(targets$Sample_Group)
individual <- factor(targets$Sample_Source) 
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(naive-rTreg, levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)


DMPs <- DMPs %>% dplyr::select(logFC, P.Value, adj.P.Val)
rownames(DMPs) <- NULL

qsave(DMPs, file = "plot_data/methylation_data.qs", preset = "custom", algorithm = "zstd", compress_level = 22)

# Airline ######################################################################
# https://www.r-bloggers.com/visualize-large-data-sets-with-the-bigvis-package/

# revolutionanalytics SSL cert expired
# Some versions of zip don't work on < 4GB zip files
system(sprintf("wget -v --no-check-certificate https://packages.revolutionanalytics.com/datasets/AirOnTime87to12/AirOnTimeCSV.zip -P %s", tempdir))
system(sprintf("7za x %s/AirOnTimeCSV.zip -o%s/AirOnTimeCSV", tempdir, tempdir))

files <- list.files(paste0(tempdir, "/AirOnTimeCSV/AirOnTimeCSV"), full.names=T)
air <- lapply(1:length(files), function(i) {
  print(i)
  z <- fread(files[i], select = c("DEP_DELAY", "ARR_DELAY")) %>% 
    filter(complete.cases(.))
  return(z)
})
qsave(air, file="plot_data/AirOnTime.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads=nt)

# GPS traces ######################################################################
# https://planet.osm.org/gps/

system(sprintf("wget -v https://planet.osm.org/gps/simple-gps-points-120312.txt.xz -P %s", tempdir))
system(sprintf("unxz %s/simple-gps-points-120312.txt.xz", tempdir))

gps <- data.table::fread(sprintf("%s/simple-gps-points-120312.txt", tempdir), header=FALSE, data.table=FALSE)
gps <- list(latitude = gps[[1]], longitude = gps[[2]]) # Long dataframes not supported in R as of 4.2.3

qsave(gps, file="plot_data/simple-gps-points.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads=nt)


