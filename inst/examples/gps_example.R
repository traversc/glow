library(qs)
library(glow)
library(EBImage)

# Note: This example takes a decent amount of resources to run, especially memory
# This could be re-worked for low memory if the data were chunked up front

output_width = 1920*4
output_height = 1080*4
radius <- 0.2 # Increase this value if plotting fewer points
intensity <- 0.6 # Increase this value if plotting fewer points

gps <- qread("plot_data/simple-gps-points.qs", nthreads=32)

length(gps$latitude) # 2770233904

gps$latitude <- gps$latitude/1e7
gps$longitude <- gps$longitude/1e7
gps$latitude <- -gps$latitude # need to flip Y-axis for some reason...?

xlimits <- range(gps$longitude)
ylimits <- range(gps$latitude)

# we need to chunk the data due to data.frame nrow limits used within the R6 class...
N <- length(gps$latitude)
starts <- seq(1,N,1e7)
stops <- c(starts[2:length(starts)]-1, N)

time <- Sys.time()

gm <- GlowMapper$new(xdim=output_width, ydim = output_height, blend_mode = "additive", nthreads=32)
for(i in 1:length(starts)) {
  print(i)
  gm$map(x=gps$longitude[starts[i]:stops[i]], y=gps$latitude[starts[i]:stops[i]], radius=radius, intensity=intensity, xlimits=xlimits, ylimits=ylimits, append=TRUE)
  if(i %% 10 == 0) print(range(gm$output_raw()))
}

pd <- gm$output_raw(saturation=1)
image_array <- array(1, dim=c(output_width, output_height, 1))
image_array[,,1] <- pd
img <- Image(image_array, colormode='Grayscale')
writeImage(img, "plots/gps_trace_8K.png")

image_array[,,1] <- 1-pd
img <- Image(image_array, colormode='Grayscale')
writeImage(img, "plots/gps_trace_invert_8K.png")


print(Sys.time() - time)



