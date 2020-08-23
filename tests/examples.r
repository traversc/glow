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
library(qs)
library(ggplot2)
library(patchwork)
library(magick)
library(viridisLite)

# COVID example
library(sf)
library(USAboundaries)
library(USAboundariesData)

# Volcano
library(minfi)
library(methylationArrayAnalysis)

nt <- 8

# result images are compressed with https://tinypng.com/ for vignettes

# Galaxy #######################################################

# run once
if(F) {
cols <- c("parallax", "dec", "ra", "astrometric_pseudo_colour", "lum_val")
files <- c("http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_1584380076484244352_2200921635402776448.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_2200921875920933120_3650804325670415744.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_2851858288640_1584379458008952960.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_3650805523966057472_4475721411269270528.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_4475722064104327936_5502601461277677696.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_5502601873595430784_5933051501826387072.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_5933051914143228928_6714230117939284352.csv.gz",
"http://cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source_with_rv/csv/GaiaSource_6714230465835878784_6917528443525529728.csv.gz")

stars <- lapply(files, function(f) {
  fread(f, select = cols, colClasses = "numeric")
}) %>% rbindlist %>% as.data.frame
stars2 <- stars[complete.cases(stars),]
stars2 <- stars2 %>% arrange_all
qsave(stars2, file="~/N/R_stuff/GAIA_plot_data.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads=8)
}

stars <- qread("~/N/R_stuff/GAIA_plot_data.qs")

# Transform to galactic coordinates
# https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html
Ag <- cbind(c(-0.0548755604162154,+0.4941094278755837,-0.8676661490190047),
            c(-0.8734370902348850,-0.4448296299600112,-0.1980763734312015),
            c(-0.4838350155487132,+0.7469822444972189,+0.4559837761750669))

cos_alpha <- cos(pi/180 * stars$ra)
cos_delta <- cos(pi/180 * stars$dec)
sin_alpha <- sin(pi/180 * stars$ra)
sin_delta <- sin(pi/180 * stars$dec)

X <- cos_alpha * cos_delta
Y <- sin_alpha * cos_delta
Z <- sin_delta

r_icrs <- rbind(X, Y, Z)
r_gal <- t(Ag %*% r_icrs)

l <- atan2(r_gal[,2], r_gal[,1])
b <- atan2(r_gal[,3], sqrt(r_gal[,1]^2 + r_gal[,2]^2))

# Transform to mollweide projection
proj <- mollweide_projection(latitude = b, longitude = l, meridian = 0)

# old monochrome plot
if(F) {
gm <- GlowMapper$new(xdim=3000, ydim = 1500, blend_mode = "additive", nthreads=nt)
gm$map(x=proj$x, y=proj$y, intensity = stars$lum_val, radius = .05, distance_exponent = 2)
pd <- gm$output_dataframe(saturation = quantile(gm$output, 0.95))
g <- ggplot(pd, aes(x = x, y = y, fill = value)) + geom_raster(show.legend = F) +
  scale_fill_gradientn(colors=additive_alpha(c("black", "blue", "white"))) +
  coord_fixed(ratio = gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = T) + 
  theme_night(bgcolor = "black") + labs(x = "Right ascension", y = "Declination")

outfile <- "~/GoogleDrive/glow/tests/GAIA_galaxy4.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)
}

ramp <- colorRamp(c("red", "#555555", "blue"))
color <- 1 / stars$astrometric_pseudo_colour
q5 <- quantile(color, 0.05)
q95 <- quantile(color, 0.95)
color <- pmin(pmax(color, q5), q95)
color <- color - min(color)
color <- ramp(color / max(color))/255

gm <- GlowMapper4$new(xdim=3000, ydim = 1500, blend_mode = "additive", nthreads=nt)
gm$map(x=proj$x, y=proj$y, r = stars$lum_val*color[,1]/1e3,
       g = stars$lum_val*color[,2]/1e3, b = stars$lum_val*color[,3]/1e3, radius = .05, distance_exponent = 2)
pd <- gm$output_dataframe(saturation = 1)
g <- ggplot(pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) + geom_raster(show.legend = F) +
  scale_fill_identity() +
  coord_fixed(ratio = gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = T) + 
  theme_night(bgcolor = "black") + labs(x = "Right ascension", y = "Declination")

outfile <- "~/GoogleDrive/glow/tests/GAIA_galaxy_pseudocolor.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300, quality = 95)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)


# COVID #######################################################


cov_cases <- fread("https://usafactsstatic.blob.core.windows.net/public/data/covid-19/covid_confirmed_usafacts.csv")
cov_cases <- cov_cases %>% 
  dplyr::transmute(name = gsub(" County", "" ,`County Name`), total = `8/21/20`)

county <- us_boundaries(map_date = NULL, type = "county", resolution = "high", states = NULL)
county <- county %>% filter(!state_name %in% c("Hawaii", "Northern Mariana Islands", "Puerto Rico", "Alaska", "American Samoa", "Guam", "Virgin Islands"))

state <- us_boundaries(map_date = NULL, type = "state", resolution = "high", states = NULL) %>% filter(!state_name %in% c("Hawaii", "Northern Mariana Islands", "Puerto Rico", "Alaska", "American Samoa", "Guam", "Virgin Islands"))

centroids <- st_centroid(county)
centroids <- cbind(centroids, st_coordinates(centroids$geometry))
centroids <- left_join(centroids, cov_cases, by = "name")
centroids$total[is.na(centroids$total)] <- 0
centroids <- filter(centroids, total > 0)


aspect <- 994/533
gm <- GlowMapper$new(xdim=round(2000*aspect), ydim = 2000, blend_mode = "additive", nthreads=nt)
gm$map(x=centroids$X, y=centroids$Y, intensity=centroids$total, radius = .5, #log10(centroids$total)/3
       distance_exponent = 2)
pd <- gm$output_dataframe(saturation = quantile(gm$output, 0.995))

g <- ggplot() + 
  geom_sf(data = county, fill = "black") + coord_sf(expand = F) + 
  geom_sf(data = state, fill = "#00000000", color = "#005555FF") + 
  coord_sf(expand = F) + 
  # geom_point(data = centroids, aes(x = X, y = Y, color = log10(total))) + 
  # scale_color_viridis_c() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(magma(12))) +
  theme_night() + 
  labs(x = "Longitude", y = "latitude")

outfile <- "~/GoogleDrive/glow/tests/US_coronavirus_8_19_2020.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)

# Spiral ######################################

N <- 5000
seq <- seq(0,sqrt(100), length.out=N)^2
x <- numeric(length=length(seq))
y <- numeric(length=length(seq))
radius <- numeric(length=length(seq))
for(i in 1:length(seq)) {
  t <- seq[i]
  xy <- exp(1i * t/2 - t/20)
  # xy <- exp( 1i * t/2 - t/30)
  x[i] <- Re(xy)
  y[i] <- Im(xy)
  radius[i] <- (x[i]^2 + y[i]^2)
}
df <- data.frame(x, y, radius)
ramp <- colorRamp(additive_alpha(viridis(12)))
color <- ramp(seq(0,1, length.out=N))/255/5
df$r <- color[,1]
df$g <- color[,2]
df$b <- color[,3]

gm <- GlowMapper4$new(xdim=2000, ydim=2000, background_color = "#00000000", blend_mode = "screen", nthreads=nt)
gm$map(x=df$x, y=df$y, radius=df$radius/8+0.1, 
       r=df$r, g=df$g, b=df$b, 
       xlimits=c(-0.9,1.2), ylimits = c(-0.9,1.2))

pd <- gm$output_dataframe(saturation = 1, saturation_mode = "overflow")
g <- ggplot(pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) + 
  geom_raster() +
  scale_fill_identity() +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = F) + 
  theme_night()

outfile <- "~/GoogleDrive/glow/tests/glow_spiral2.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)

# IRIS ##########################################

df <- iris
df <- df %>% mutate(
  color = case_when(iris$Species == "setosa" ~ "green",
                    iris$Species == "versicolor" ~ "red",
                    iris$Species == "virginica" ~ "blue"),
  exponent = case_when(iris$Species == "setosa" ~ 1.2,
                    iris$Species == "versicolor" ~ 2,
                    iris$Species == "virginica" ~ 2.8),
  radius = case_when(iris$Species == "setosa" ~ 0.11,
                    iris$Species == "versicolor" ~ 0.13,
                    iris$Species == "virginica" ~ 0.13)
)


lg <- df %>% dplyr::select(Species, color, exponent, radius) %>%
  distinct %>%
  arrange(Species) %>%
  mutate(x = 2, y = c(4.4, 4.3, 4.2), 
         radius = radius * 1.2,
         text = c("Setosa", "Versicolor", "Virginica"))

gm <- GlowMapper4$new(xdim=2000, ydim=1600, background_color = c(0,0,0,0), blend_mode = "screen", nthreads=nt)
gm$map(x=df$Petal.Width, y=df$Sepal.Width, 
       color = df$color, radius = df$radius, distance_exponent = df$exponent)
pd <- gm$output_dataframe(saturation = 1)

l <- GlowMapper4$new(xdim=2000, ydim=1600, background_color = c(0,0,0,0), blend_mode = "screen", nthreads=nt)
l$map(x=lg$x, y=lg$y, 
       color = lg$color, radius = lg$radius, distance_exponent = lg$exponent, xlimits = gm$xlim(), ylimits = gm$ylim())
ld <- l$output_dataframe(saturation = 1)

g <- ggplot() + 
  geom_raster(data=pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) +
  geom_raster(data=ld, aes(x = x, y = y, fill = rgb(r,g,b,a))) +
  geom_text(data=lg, aes(x = x+0.05, y = y, label = text), color = "gray80", hjust = 0, size = 3) +
  scale_fill_identity() +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = F) + 
  labs(x = "Petal Width", y = "Sepal Width") +
  theme_night()

# doesn't look good enough currently
if(F) {
outfile <- "~/GoogleDrive/glow/tests/iris_example.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)
}

# Diamonds ##########################################
gm <- GlowMapper$new(xdim=2000, ydim = 2000, blend_mode = "screen", nthreads=nt)
gm$map(x=diamonds$carat, y=diamonds$price, intensity=0.5, radius = .1)
pd <- gm$output_dataframe(saturation = 1)
gglow <- ggplot() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(magma(12))) +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "carat", y = "price") + 
  theme_night(bgcolor = magma(1))

# ghex <-  ggplot(diamonds, aes(carat, price)) + 
#   geom_hex(bins=50, show.legend=F) + 
#   scale_fill_gradientn(colors = magma(12)[3:12]) +
#   coord_fixed(ratio=gm$aspect()) + 
#   theme_night(bgcolor = magma(1))

# g <- ghex + gglow + plot_layout(ncol = 2)

outfile <- "~/GoogleDrive/glow/tests/diamonds.png"
ggsave(gglow, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)

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

adj_pval_threshold <- DMPs %>% filter(adj.P.Val < 0.05) %>%
  pull(P.Value) %>% max

# ggplot(DMPs, aes(x=logFC, y = -log10(P.Value))) +
#   geom_point(shape = ".", color = "white") +
#   geom_hline(yintercept = -log10(adj_pval_threshold), lty=2, color = "red") +
#   labs(y = "-log10 P-value") +
#   theme_night()

gm <- GlowMapper$new(xdim = 2000, ydim = 2000, blend_mode = "screen", nthreads=nt)
gm$map(x = DMPs$logFC, y = -log10(DMPs$P.Value), radius = 0.1, intensity = .5)
pd <- gm$output_dataframe()

g <- ggplot(pd, aes(x=x, y=y, fill = value)) + 
  geom_raster(show.legend=F) + 
  geom_hline(yintercept = -log10(adj_pval_threshold), lty=2, color = "red") + 
  labs(x = "log Fold Change", y = "-log10 P-value") +
  scale_fill_gradientn(colors = additive_alpha(magma(12))) + 
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  theme_night(bgcolor = magma(12)[1])

outfile <- "~/GoogleDrive/glow/tests/volcano.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)


# Airline ###################################
# https://www.r-bloggers.com/visualize-large-data-sets-with-the-bigvis-package/
## wget https://packages.revolutionanalytics.com/datasets/AirOnTime87to12/AirOnTimeCSV.zip .

# run once
if(F) {
files <- list.files("~/N/R_stuff/AirOnTimeCSV/", full.names=T)
air <- lapply(1:length(files), function(i) {
  print(i)
  z <- fread(files[i], select = c("DEP_DELAY", "ARR_DELAY")) %>% 
    filter(complete.cases(.))
  file.remove(files[i])
  return(z)
})
qsave(air, file="~/N/R_stuff/AirOnTime.qs", preset = "custom", algorithm = "zstd", compress_level = 22, nthreads=8)
}

air <- qread("~/N/R_stuff/AirOnTime.qs", nthreads=8)

temp <- rbindlist(air)
qlo1 <- temp$ARR_DELAY %>% quantile(0.0025)
qhi1 <- temp$ARR_DELAY %>% quantile(0.9975)
qlo2 <- temp$DEP_DELAY %>% quantile(0.0025)
qhi2 <- temp$DEP_DELAY %>% quantile(0.9975)
rm(temp)

gm <- GlowMapper$new(xdim = 2000, ydim = 2000, blend_mode = "screen", nthreads=nt)

# We can load each CSV and append plotting data to the existing raster
# We've loaded the entire dataset as a list here since we have a lot of memory :)
for(i in 1:length(air)) {
  print(i)
  a <- air[[i]]
  a <- filter(a, ARR_DELAY < qhi1, ARR_DELAY > qlo1, DEP_DELAY < qhi1, DEP_DELAY > qlo1)
  gm$map(x = a$DEP_DELAY, y = a$ARR_DELAY, radius = .5, intensity = .05, append=T)
}
pd <- gm$output_dataframe()

g <- ggplot(pd, aes(x=x, y=y, fill = value)) + 
  geom_raster(show.legend=F) + 
  scale_fill_gradientn(colors = additive_alpha(magma(12))) + 
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "Departure Delay (minutes)", y = "Arrival Delay (minutes)") +
  theme_night(bgcolor = magma(12)[1])

outfile <- "~/GoogleDrive/glow/tests/airline_mt.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
img <- magick::image_read(outfile)
img <- magick::image_trim(img, fuzz = 50)
magick::image_write(img, outfile)

# Other interesting datasets ###################################

# https://hackernoon.com/drawing-2-7-billion-points-in-10s-ecc8c85ca8fa
# https://stackoverflow.com/questions/51122970/efficiently-plotting-hundreds-of-millions-of-points-in-r
