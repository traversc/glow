# Glow plot examples
# 
# Milky way galaxy
# U.S. COVID cases
# Volcano plot
# Iris dataset
# Spiral

# All examples in this script can be easily run on a desktop computer
# Just downscale the output resolution

# See examples_parse_data.R for data pre-processing

library(dplyr)
library(data.table)
library(qs2)

library(glow)
library(ggplot2)
library(sf)
library(EBImage)
library(viridisLite)

nt <- 8 # number of threads to use

outdir <- "plots"
dir.create(outdir, showWarnings=FALSE)

create_EBImage <- function(raw, bgcolor = "black", colormode = "Color") {
  bgcolor_rgb <- as.vector(col2rgb(bgcolor))
  if(colormode == "Color") {
    image_array <- array(0, dim=c(dim(raw[[1]]), 3))
    print(dim(image_array))
    image_array[,,1] <- pd[[1]]*pd[[4]] + bgcolor_rgb[1] * (1-pd[[4]])
    image_array[,,2] <- pd[[2]]*pd[[4]] + bgcolor_rgb[2] * (1-pd[[4]])
    image_array[,,3] <- pd[[3]]*pd[[4]] + bgcolor_rgb[3] * (1-pd[[4]])
  } else if(colormode == "Grayscale") {
    image_array <- array(0, dim=c(dim(raw), 1))
    image_array <- raw
  } else {
    stop("colormode must be Color or Grayscale")
  }
  EBImage::Image(image_array, colormode=colormode)
}

# ugly hack to trim plots with coord_cartesian
trim_image <- function(file, bgcolor) {
  # normalize input color
  col <- col2rgb(bgcolor)
  col <- rgb(red=col[1]/255,green=col[2]/255,blue=col[3]/255)
  
  cmd <- sprintf('convert %s -trim -bordercolor "%s" -border 20 %s', file, col, file)
  system(cmd)
}

# Gaia Stars Galaxy #######################################################

output_width = 1920*4
output_height = 1080*4
outfile <- "plots/GAIA_galaxy_pseudocolor.png"

stars <- qs_read("plot_data/gaia_stars.qs2")

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

ramp <- colorRamp(c("red", "#555555", "blue"))
color <- 1 / stars$astrometric_pseudo_colour
q5 <- quantile(color, 0.05)
q95 <- quantile(color, 0.95)
color <- pmin(pmax(color, q5), q95)
color <- color - min(color)
color <- ramp(color / max(color))/255

gm <- GlowMapper4$new(xdim=output_width, ydim=output_height, blend_mode = "additive", nthreads=nt)
gm$map(x=proj$x, y=proj$y, r = stars$lum_val*color[,1]/1e3,
       g = stars$lum_val*color[,2]/1e3, b = stars$lum_val*color[,3]/1e3, radius = .05, distance_exponent = 2)

pd <- gm$output_raw(saturation = 1)
img <- create_EBImage(pd, "black")
writeImage(img, "plots/GAIA_galaxy_pseudocolor.png")

# pd <- gm$output_dataframe(saturation = 1)
# g <- ggplot(pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) + geom_raster(show.legend = F) +
#   scale_fill_identity() +
#   coord_fixed(ratio = gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = T) + 
#   theme_night(bgcolor = "black") + labs(x = "Right ascension", y = "Declination")
# ggsave(g, file=outfile, width=10, dpi=300)
# img <- magick::image_read(outfile)
# img <- magick::image_trim(img, fuzz = 50)
# magick::image_write(img, outfile)


# COVID ########################################################################

outfile <- "plots/US_coronavirus_2021.png"

cov_cases <- qs_read("plot_data/covid_confirmed_usafacts.qs2")
centroids <- cov_cases$centroids
state <- cov_cases$state
county <- cov_cases$county

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
  scale_fill_gradientn(colors = additive_alpha(cividis(12))) +
  theme_night() + 
  labs(x = "Longitude", y = "latitude")

ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile, "black")

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

outfile <- "plots/glow_spiral2.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile, "black")

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

gm <- GlowMapper4$new(xdim=1920, ydim=1440, background_color = c(0,0,0,0), blend_mode = "screen", nthreads=nt)
gm$map(x=df$Petal.Width, y=df$Sepal.Width, 
       color = df$color, radius = df$radius, distance_exponent = df$exponent)
pd <- gm$output_dataframe(saturation = 1)

l <- GlowMapper4$new(xdim=1920, ydim=1440, background_color = c(0,0,0,0), blend_mode = "screen", nthreads=nt)
l$map(x=lg$x, y=lg$y, 
       color = lg$color, radius = lg$radius, distance_exponent = lg$exponent, xlimits = gm$xlim(), ylimits = gm$ylim())
ld <- l$output_dataframe(saturation = 1)

g <- ggplot() + 
  geom_raster(data=pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) +
  geom_rect(aes(xmin=1.9+0.025, xmax = 2.35-0.025, ymin = 4.1+0.025, ymax = 4.5-0.025), fill = "white") + 
  geom_raster(data=ld, aes(x = x, y = y, fill = rgb(r,g,b,a))) +
  geom_text(data=lg, aes(x = x+0.05, y = y, label = text), color = "black", hjust = 0, size = 3) +
  scale_fill_identity() +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim(), expand = F) + 
  labs(x = "Petal Width", y = "Sepal Width") +
  theme_night(bgcolor = "whitesmoke") + 
  theme(panel.grid.major=element_line(colour="darkgray"),
        panel.grid.minor=element_line(colour="darkgray"), axis.title = element_text(color = "darkgray"))

outfile <- "plots/iris_example.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile, "whitesmoke")

# Diamonds ##########################################
gm <- GlowMapper$new(xdim=2000, ydim = 2000, blend_mode = "screen", nthreads=nt)
gm$map(x=diamonds$carat, y=diamonds$price, intensity=0.5, radius = .1)
pd <- gm$output_dataframe(saturation = 1)
g <- ggplot() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(cividis(12))) +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "carat", y = "price") + 
  theme_night(bgcolor = cividis(12)[1])


outfile <- "plots/diamonds.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile, cividis(12)[1])

# Diamonds vignette examples ##########################################
library(glow)
library(ggplot2)
library(viridisLite) # Magma color scale

data(diamonds)
gm <- GlowMapper$new(xdim=800, ydim = 640, blend_mode = "screen", nthreads=nt)
gm$map(x=diamonds$carat, y=diamonds$price, intensity=1, radius = rely(0.002))
pd <- gm$output_dataframe(saturation = 1)

# Dark color theme
g <- ggplot() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(magma(12))) +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "carat", y = "price") + 
  theme_night(bgcolor =  magma(12)[1])

outfile <- "plots/diamonds_vignette_dark.png"
ggsave(g, file=outfile, width=10, height=4, dpi=96)
trim_image(outfile, magma(12)[1])


# light color theme
light_colors <- glow::light_heat_colors(144)
g <- ggplot() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(light_colors)) +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "carat", y = "price") + 
  theme_bw(base_size = 14)

outfile <- "plots/diamonds_vignette_light.png"
ggsave(g, file=outfile, width=10, height=4, dpi=96)
trim_image(outfile, "white")

# light color theme with cool colors
light_colors <- light_cool_colors(144)
g <- ggplot() + 
  geom_raster(data = pd, aes(x = pd$x, y = pd$y, fill = pd$value), show.legend = F) +
  scale_fill_gradientn(colors = additive_alpha(light_colors)) +
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "carat", y = "price") + 
  theme_bw(base_size = 14)

outfile <- "plots/diamonds_vignette_cool.png"
ggsave(g, file=outfile, width=10, height=4, dpi=96)
trim_image(outfile, "white")


# Volcano ##########################################

DMPs <- qs_read("plot_data/methylation_data.qs2")

adj_pval_threshold <- DMPs %>% filter(adj.P.Val < 0.05) %>%
  pull(P.Value) %>% max

gm <- GlowMapper$new(xdim = 2000, ydim = 2000, blend_mode = "screen", nthreads=nt)
gm$map(x = DMPs$logFC, y = -log10(DMPs$P.Value), radius = 0.1, intensity = .5)
pd <- gm$output_dataframe()

g <- ggplot(pd, aes(x=x, y=y, fill = value)) + 
  geom_raster(show.legend=F) + 
  geom_hline(yintercept = -log10(adj_pval_threshold), lty=2, color = "red") + 
  labs(x = "log Fold Change", y = "-log10 P-value") +
  scale_fill_gradientn(colors = additive_alpha(light_heat_colors(12))) + 
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  theme_bw(base_size=14)

outfile <- "plots/volcano_white_bg.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile,  "white")


# Airline ###################################
# https://www.r-bloggers.com/visualize-large-data-sets-with-the-bigvis-package/
## wget https://packages.revolutionanalytics.com/datasets/AirOnTime87to12/AirOnTimeCSV.zip .

air <- qs_read("plot_data/AirOnTime.qs2", nthreads=nt)

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
  scale_fill_gradientn(colors = additive_alpha(mako(12))) + 
  coord_fixed(gm$aspect(), xlim = gm$xlim(), ylim = gm$ylim()) + 
  labs(x = "Departure Delay (minutes)", y = "Arrival Delay (minutes)") +
  theme_night(bgcolor = mako(12)[1])

outfile <- "plots/airline_mt.png"
ggsave(g, file=outfile, width=10, height=4, dpi=300)
trim_image(outfile,  mako(12)[1])

# Other interesting datasets ###################################

# ??