library(glow)
library(viridisLite)
library(Rcpp)
library(EBImage)

# output parameters
output_width = 1920/2
output_height = 1080/2
N <- 1e5
M <- 10

# C++ function to rapidly generate points
sourceCpp(code=
  ' #include <Rcpp.h>
    using namespace Rcpp;
    // [[Rcpp::export(rng = false)]]
    DataFrame clifford_attractor(size_t n_iter, double A, double B, double C, double D, double x0, double y0) {
    NumericVector x(n_iter);
    NumericVector y(n_iter);
    NumericVector angle(n_iter);
    double * xp = REAL(x);
    double * yp = REAL(y);
    double * anglep = REAL(angle);
    xp[0] = x0;
    yp[0] = y0;
    anglep[0] = 0;
    for (size_t i=1; i < n_iter; ++i) {
      double diffx = (sin(A*yp[i-1]) + C*cos(A*xp[i-1]));
      double diffy = (sin(B*xp[i-1]) + D*cos(B*yp[i-1]));
      anglep[i] = atan2(diffy, diffx);
      xp[i] += diffx;
      yp[i] += diffy;
    }
    return DataFrame::create(_["x"]=x, _["y"]=y, _["angle"]=angle);
    }'
)

polar_viridis <- function(angle) {
  pal <- viridis(1024, alpha = 0.5)
  pal <- c(pal, rev(pal))
  pal[findInterval(angle, seq(-pi, pi, length.out = length(pal) + 1), all.inside = TRUE)]
}

time <- Sys.time()
gm <- GlowMapper4$new(xdim=output_width, ydim = output_height, blend_mode = "additive", nthreads=16)

x0 <- 0.1
y0 <- 0
for(i in 1:M) {
  print(i)
  cliff_points <- clifford_attractor(N, 1.886,-2.357,-0.328, 0.918, x0, y0)
  cliff_points$angle <- polar_viridis(cliff_points$angle)
  gm$map(x=cliff_points$x[-1], y=cliff_points$y[-1], radius = 0.005, color = cliff_points$angle[-1], append = TRUE)
  x0 <- tail(cliff_points$x, 1)
  y0 <- tail(cliff_points$y, 1)
}

pd <- gm$output_raw(saturation = 1)
image_array <- array(1, dim=c(output_width, output_height, 3))
image_array[,,1] <- pd[[1]]*pd[[4]]
image_array[,,2] <- pd[[2]]*pd[[4]]
image_array[,,3] <- pd[[3]]*pd[[4]]
img = Image(image_array, colormode='Color')
writeImage(img, "plots/clifford.png")
print(Sys.time() - time)



