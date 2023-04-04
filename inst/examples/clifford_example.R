library(glow)
library(viridisLite)
library(EBImage)

# Note: This example takes a LOT of resources to run
# Decrease the parameter scaling below to make it run in a reasonable amount of time
# parameters
output_width = 1920*4
output_height = 1080*4
N <- 20e6 # plot N * M points
M <- 50
radius <- 0.005 # Increase this value if plotting fewer points

polar_rainbow <- function(x, min_limit=NULL, max_limit=NULL, pal_function = rainbow, invert=FALSE) {
  if(is.null(min_limit)) min_limit <- min(x)
  if(is.null(max_limit)) max_limit <- max(x)
  pal <- pal_function(1024, alpha = 0.5)
  if(invert) pal <- rev(pal)
  pal <- c(pal, rev(pal))
  pal[findInterval(x, seq(min_limit, max_limit, length.out = length(pal) + 1), all.inside = TRUE)]
}

time <- Sys.time()
gm <- GlowMapper4$new(xdim=output_width, ydim = output_height, blend_mode = "additive", nthreads=32)
x0 <- 0.1
y0 <- 0
for(i in 1:M) {
  print(i)
  cliff_points <- clifford_attractor(N, 1.886,-2.357,-0.328, 0.918, x0, y0)
  # color <- polar_rainbow(cliff_points$angle, -pi, pi, viridis)
  color <- polar_rainbow(cliff_points$distance, pal_function=viridis, invert=TRUE)
  gm$map(x=cliff_points$x[-1], y=cliff_points$y[-1], radius = radius, color = color[-1], append = TRUE) # radius for 1 billion
  # gm$map(x=cliff_points$x[-1], y=cliff_points$y[-1], radius = radius, color = color[-1], append = TRUE) # radius for 10 million
  x0 <- tail(cliff_points$x, 1)
  y0 <- tail(cliff_points$y, 1)
  
  pd <- gm$output_raw(saturation = 1)
  image_array <- array(1, dim=c(output_width, output_height, 3))
  image_array[,,1] <- pd[[1]]*pd[[4]]
  image_array[,,2] <- pd[[2]]*pd[[4]]
  image_array[,,3] <- pd[[3]]*pd[[4]]
  img = Image(image_array, colormode='Color')
  writeImage(img, sprintf("plots/clifford_distance_inverse_viridis_%s.png", i))
}
print(Sys.time() - time)

################################################################################
# example for vignette

cliff_points <- clifford_attractor(1e6, 1.886,-2.357,-0.328, 0.918, 0.1, 0)
color_pal <- circular_palette(n=144, pal_function=rainbow, alpha=0.5)
cliff_points$color <- map_colors(color_pal, cliff_points$angle, min_limit=-pi, max_limit=pi)

gm <- GlowMapper4$new(xdim=960, ydim = 540, blend_mode = "additive", nthreads=4)
gm$map(x=cliff_points$x, y=cliff_points$y, radius=0.02, color=cliff_points$color)
pd <- gm$output_raw(saturation = 1)

image_array <- array(1, dim=c(960, 540, 3))
image_array[,,1] <- pd[[1]]*pd[[4]]
image_array[,,2] <- pd[[2]]*pd[[4]]
image_array[,,3] <- pd[[3]]*pd[[4]]
img <- EBImage::Image(image_array, colormode='Color')
writeImage(img, "plots/clifford_vignette.png")



