# Glow plot examples
# 
# Milky way galaxy
# U.S. COVID cases
# Volcano plot
# Iris dataset
# Fractal (TO DO)
# Spiral

library(glow)
library(dplyr)
library(data.table)
library(arrow)
library(duckdb)
library(qs)
library(stringr)

# COVID example
library(sf)
library(USAboundaries) # install.packages("USAboundariesData", repos = "https://ropensci.r-universe.dev", type = "source")
library(USAboundariesData)

# Volcano
library(minfi)
library(methylationArrayAnalysis)

tempdir <- tempdir()
outdir <- "plot_data"
dir.create(outdir, showWarnings=FALSE)

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
  arrow::read_csv_arrow(f, col_select = c("parallax", "dec", "ra", "astrometric_pseudo_colour", "lum_val"))
}) %>% rbindlist %>% as.data.frame
stars2 <- stars[complete.cases(stars),]
stars2 <- stars2 %>% arrange_all

qsave(stars2, "plot_data/gaia_stars.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads = 8)


# COVID #######################################################################

cov_cases <- fread("https://usafactsstatic.blob.core.windows.net/public/data/covid-19/covid_confirmed_usafacts.csv", data.table=F)

cov_cases <- cov_cases %>% 
  dplyr::transmute(state_abbr = State, name = stringr::str_trim(gsub(" County", "" ,`County Name`)), total = `2020-08-21`)

plot_states <- setdiff(USAboundaries::us_states()$state_name, c("Hawaii", "Northern Mariana Islands", "Puerto Rico", "Alaska", "American Samoa", "Guam", "Virgin Islands"))
county <- us_counties(map_date = NULL, resolution = "high", states = plot_states)
state <- us_states(map_date = NULL, resolution = "high", states = plot_states)

centroids <- st_centroid(county)
centroids <- cbind(centroids, st_coordinates(centroids$geometry))
centroids <- left_join(centroids, cov_cases, by = "name")
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
mSetSq <- preprocessQuantile(rgSet)
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
colnames(DMPs) <- NULL

qsave(DMPs, file = "plot_data/methylation_data.qs", preset = "custom", algorithm = "zstd", compress_level = 22)




# Airline ######################################################################
# https://www.r-bloggers.com/visualize-large-data-sets-with-the-bigvis-package/

# revolutionanalytics SSL cert expired
system(sprintf("wget -v --no-check-certificate https://packages.revolutionanalytics.com/datasets/AirOnTime87to12/AirOnTimeCSV.zip %s", tempdir))


files <- list.files(paste0(tempdir, "/AirOnTimeCSV/"), full.names=T)
air <- lapply(1:length(files), function(i) {
  print(i)
  z <- fread(files[i], select = c("DEP_DELAY", "ARR_DELAY")) %>% 
    filter(complete.cases(.))
  file.remove(files[i])
  return(z)
})
qsave(air, file="~/N/R_stuff/AirOnTime.qs", preset = "custom", algorithm = "zstd", compress_level = 12, nthreads=nt)

air <- qread("~/N/R_stuff/AirOnTime.qs", nthreads=nt)

temp <- rbindlist(air)
qlo1 <- temp$ARR_DELAY %>% quantile(0.0025)
qhi1 <- temp$ARR_DELAY %>% quantile(0.9975)
qlo2 <- temp$DEP_DELAY %>% quantile(0.0025)
qhi2 <- temp$DEP_DELAY %>% quantile(0.9975)
rm(temp)

gm <- GlowMapper$new(xdim = 2000, ydim = 2000, blend_mode = "screen", nthreads=nt)


# GPS traces ######################################################################

# https://planet.osm.org/gps/

# Attribution must be to OpenStreetMap.
# 
# Attribution must also make it clear that the data is available under the Open 
# Database License. This may be done by making the text OpenStreetMap a link to 
# openstreetmap.org/copyright, which has information about OpenStreetMaps data 
# sources (which OpenStreetMap needs to credit) as well as the ODbL.
