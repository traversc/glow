# theme_night ###################################################

#' theme_night
#' 
#' A dark ggplot2 theme with a default black background intended to be used with the glow package.
#' @usage theme_night(bgcolor = "black", base_size = 14, base_family = "")
#' @param bgcolor Background color, default black. Generally you want to match the background with the lowest color value on a color scale. 
#' @param base_size Base default font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme.
#' @details The theme is heavily modified from the minimal ggplot theme. It is intended to be use with dark background colors and should not be used with white or light backgrounds. 
#' @examples
#' ggplot(mtcars, aes(x = mpg, y = wt)) + 
#'   geom_point(color = "white") +
#'   theme_night(bgcolor = "black")
#' @name theme_night
NULL

# additive_alpha ###################################################

#' additive_alpha
#' 
#' Simulates additive blending on a dark to light color scale.
#' @usage additive_alpha(colors)
#' @param colors colors
#' @return A simulated additive scale of the input colors.
#' @details In R plotting (both ggplot and base R) blending is performed by alpha blending, which is an averaging effect. When combining light and glow effects, additive blending is more appropriate. 
#' 
#' This function simulates additive blending by increasing color on a color scale to compensate for the averaging effect of alpha blending. 
#' 
#' Note: this function is only appropriate for dark to light color scales.  
#' @examples
#' m_solid <- viridisLite::magma(12)
#' m_additive <- additive_alpha(m_solid)
#' @name additive_alpha
NULL

# mollweide_projection ###################################################

#' mollweide_projection
#' 
#' Performs a cartographic mollweide projection from polar coordinates (latitude/longitude) to X-Y map coordinates
#' @usage mollweide_projection(latitude, longitude, meridian)
#' @param latitude Latitude (aka declination) of points
#' @param longitude Longitude (aka right ascension) of points
#' @param meridian The x=0 center of the plot
#' @return X/Y coordinates
#' @details 
#' This function uses the "Newton-Raphson with fast convergence everywhere" algorithm. 
#' 
#' Latitude and longitude should be in units of radians not degrees. Latitude ranges from +/- pi/2 and longitude ranges from +/- pi. 
#' @seealso https://en.wikipedia.org/wiki/Talk:Mollweide_projection
#' @examples
#' longitude <- pi / 4
#' latitude <- pi / 4
#' mollweide_projection(longitude, latitude, meridian = 0)
#' @name mollweide_projection
NULL

# GlowMapper ###################################################

#' GlowMapper
#'
#' This class provides a framework for creating scatter plots based on a glow simulation. Points are mapped with a gaussian gradient to a raster with specified dimensions and properties. 
#'
#' @section Usage:
#' \preformatted{m <- GlowMapper$new(xdim=1000, ydim=800, blend_mode = "screen", contrast_limit = 1e5, nthreads = 1)
#' 
#' m$map(x, y, radius, intensity = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE)
#' 
#' m$output_raw(saturation = NA_real_)
#' 
#' m$output_dataframe(saturation = NA_real_)
#' 
#' m$aspect()
#' 
#' m$xlim()
#' 
#' m$ylim()
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{xdim}{- The first dimension of the output matrix raster.}
#'   \item{ydim}{- The second dimension of the output matrix raster.}
#'   \item{blend_mode}{- Either screen or additive blending mode. See details.}
#'   \item{contrast_limit}{- Determines the distance to search from a point. You shouldn't need to change this unless you have a lot of points in the plot stacked on top of each other.}
#'   \item{nthreads}{- Number of threads to use.}
#'   \item{x}{- X coordinate of points.}
#'   \item{y}{- Y coordinate of points.}
#'   \item{radius}{- Relative spread of glow intensity. The radius should be proportional to the x and y-ranges of the plot. Values between 1/10 to 1/100 of the range of the plot generally produce good results.}
#'   \item{intensity}{- Maximum intensity at the center of a point.}
#'   \item{distance_exponent}{- Exponent of the distance calculation when calculating intensities. A value of 2 corresponds to euclidean distance; a value of 1 corresponds to manhattan distance.}
#'   \item{xlimits}{- The x-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{ylimits}{- The y-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{append}{ - Whether to add to the existing output or overwrite.}
#'   \item{saturation}{- When retrieving the output with \code{$output_raw} or \code{$output_dataframe}, maximum intensity values are capped at the given value. This is often useful when using additive blend mode to increase contrast.}
#' }
#'
#' @section Details:
#' \code{$new()} creates a new GlowMapper object, which holds parameters, plotting data, and the output (a matrix of glow intensities). Creates a canvas to plot point data. With additive blending, the intensities of each point are added arithmetically, which is how light intensities are added in the physical world. This is equivalent to an fast/approximate un-normalized 2D kernel density estimate. 
#' 
#' With "screen" blending,  two intensities are added according to the formula: \code{I_out = 1 - (1-I_a)*(1-I_b)}. Both additive blending and screen blending are commutative operations, meaning the order of points in a plot does not affect the output.
#' 
#' Screen blending can often improve contrast in a plot and is the default. 
#'
#' \code{$map()} maps points to the canvas. 
#' 
#' \code{$output_raw()} output raw matrix rasters. Useful for plotting in base R.
#' 
#' \code{$output_dataframe()} output the raster as a dataframe with XY coordinates. This is meant to pipe directly into ggplot.
#' 
#' \code{$aspect()}, \code{$xlim()}, \code{$ylim()} return the aspect ratio, x-limits and y-limits of the raster. These functions are intended to be used with plotting functions (e.g. \code{ggplot2::coord_fixed()}) so that the output raster is not distorted. See example below. 
#'
#' @examples
#' # Plot Data: x,y,r
#' x <- numeric(length=50)
#' y <- numeric(length=50)
#' r <- numeric(length=50)
#' for(t in 1:50) {
#'   xy <- exp(1i * t/2 - t/12)
#'   x[t] <- Re(xy)
#'   y[t] <- Im(xy)
#'   r[t] <- sqrt(x[t]^2 + y[t]^2)
#' }
#' 
#' # New class object
#' m <- GlowMapper$new(xdim=500, ydim = 400, blend_mode = "screen")
#' 
#' # Map data on to raster
#' m$map(x=x, y=y, intensity = 1, radius = r/4 + 0.1, distance_exponent = 2)
#' 
#' # Output raster data as a dataframe
#' pd <- m$output_dataframe(saturation = 1)
#' 
#' # Plot with ggplot
#' ggplot(pd, aes(x = x, y = y, fill = value)) + 
#'   geom_raster(show.legend = FALSE) +
#'   scale_fill_gradientn(colors=additive_alpha(c("black", "purple", "white"))) +
#'   coord_fixed(ratio = m$aspect(), xlim = m$xlim(), ylim = m$ylim(), expand = FALSE) + 
#'   theme_night(bgcolor = "black")
#' @name GlowMapper
NULL


# GlowMapper4 ###################################################

#' GlowMapper4
#'
#'This class provides a framework for creating scatter plots based on a glow simulation with explicit color intensities. Points are mapped with a gaussian gradient to a raster with specified dimensions and properties. 
#'
#' @section Usage:
#' \preformatted{m <- GlowMapper$new(xdim=1000, ydim=800, blend_mode = "additive", background_color = "#00000000", contrast_limit = 1e5, nthreads = 1)
#' 
#' m$map(x, y, radius, color = NULL, r=NULL, g=NULL, b=NULL, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE)
#' 
#' m$output_raw(saturation = 1, saturation_mode = "overflow")
#' 
#' m$output_dataframe(saturation = 1, saturation_mode = "overflow")
#' 
#' m$aspect()
#' 
#' m$xlim()
#' 
#' m$ylim()
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{xdim}{- The first dimension of the output matrix raster.}
#'   \item{ydim}{- The second dimension of the output matrix raster.}
#'   \item{blend_mode}{- Either screen or additive blending mode. See details.}
#'   \item{background_color}{- A color that can be coerced to RGBA with `col2rgb`, or a vector of four values between 0 and 1.}
#'   \item{contrast_limit}{- Determines the distance to search from a point. You shouldn't need to change this unless you have a lot of points in the plot stacked on top of each other.}
#'   \item{nthreads}{- Number of threads to use.}
#'   \item{x}{- X coordinate of points.}
#'   \item{y}{- Y coordinate of points.}
#'   \item{radius}{- Relative spread of glow intensity. The radius should be proportional to the x and y-ranges of the plot. Values between 1/10 to 1/100 of the range of the plot generally produce good results.}
#'   \item{color}{- Color of points. If NULL, r, g, and b parameters must be defined (and vice versa).}
#'   \item{r}{- Red intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{g}{- Green intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{b}{- Blue intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{distance_exponent}{- Exponent of the distance calculation when calculating intensities. A value of 2 corresponds to euclidean distance; a value of 1 corresponds to manhattan distance.}
#'   \item{xlimits}{- The x-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{ylimits}{- The y-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{append}{ - Whether to add to the existing output or overwrite.}
#'   \item{saturation}{- When retrieving the output with \code{$output_raw} or \code{$output_dataframe}, maximum intensity values are capped at the given value. This is often useful when using additive blend mode to increase contrast.}
#'   \item{saturation_mode}{- When intensity values are above the saturation threshold, values can be overflowed into other color channels ("overflow") or simply clipped at the threshold ("clip"). "Overflow" always produces a gradient to white for intensities above the threshold, which may produce artistically better results.}
#' }
#'
#' @section Details:
#' 
#' This `GlowMapper4` class is similar to the `GlowMapper` class, but instead of a single intensity matrix output, color is specified explicitly.
#' 
#' \code{$new()} creates a new GlowMapper object, which holds parameters, plotting data, and the output (a matrix of glow intensities). Creates a canvas to plot point data. With additive blending, the intensities of each point are added arithmetically, which is how light intensities are added in the physical world. This is equivalent to an fast/approximate un-normalized 2D kernel density estimate. 
#' 
#' With "screen" blending,  two intensities are added according to the formula: \code{I_out = 1 - (1-I_a)*(1-I_b)}. Both additive blending and screen blending are commutative operations, meaning the order of points in a plot does not affect the output.
#' 
#' Screen blending can often improve contrast in a plot and is the default. 
#'
#' \code{$map()} maps points to the canvas. 
#' 
#' \code{$output_raw()} output raw matrix rasters (a list of four matrices, one for each RGBA channel). Useful for plotting in base R.
#' 
#' \code{$output_dataframe()} output the raster as a dataframe with XY coordinates. This is meant to pipe directly into ggplot.
#' 
#' \code{$aspect()}, \code{$xlim()}, \code{$ylim()} return the aspect ratio, x-limits and y-limits of the raster. These functions are intended to be used with plotting functions (e.g. \code{ggplot2::coord_fixed()}) so that the output raster is not distorted. See example below. 
#'
#' @examples
#' # Plot Data: x,y,r
#' x <- numeric(length=50)
#' y <- numeric(length=50)
#' r <- numeric(length=50)
#' color <- character(length=50)
#' for(t in 1:50) {
#'   xy <- exp(1i * t/2 - t/12)
#'   x[t] <- Re(xy)
#'   y[t] <- Im(xy)
#'   r[t] <- sqrt(x[t]^2 + y[t]^2)
#'   color[t] <- rgb(t/50,0,1-t/50)
#' }
#' 
#' # New class object
#' m <- GlowMapper4$new(xdim=500, ydim = 400, blend_mode = "additive")
#' 
#' # Map data on to raster
#' m$map(x=x, y=y, color = color, radius = r/4 + 0.1, distance_exponent = 2)
#' 
#' # Output raster data as a dataframe
#' pd <- m$output_dataframe(saturation = 1, saturation_mode = "overflow")
#' 
#' # Plot with ggplot
#' ggplot(pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) + 
#'   geom_raster(show.legend = FALSE) +
#'   scale_fill_identity() +
#'   coord_fixed(ratio = m$aspect(), xlim = m$xlim(), ylim = m$ylim(), expand = FALSE) + 
#'   theme_night(bgcolor = "black")
#' @name GlowMapper4
NULL

# LightMapper ###################################################

#' LightMapper
#'
#' This class provides a framework for creating scatter plots based on a point light simulation. Points are mapped with a inverse power gradient to a raster with specified dimensions and properties. 
#'
#' @section Usage:
#' \preformatted{m <- LightMapper$new(xdim=1000, ydim=800, blend_mode = "screen", contrast_limit = 1e5, nthreads = 1)
#' 
#' m$map(x, y, radius, intensity = 1, radius, falloff_exponent = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE)
#' 
#' m$output_raw(saturation = NA_real_)
#' 
#' m$output_dataframe(saturation = NA_real_)
#' 
#' m$aspect()
#' 
#' m$xlim()
#' 
#' m$ylim()
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{xdim}{- The first dimension of the output matrix raster.}
#'   \item{ydim}{- The second dimension of the output matrix raster.}
#'   \item{blend_mode}{- Either screen or additive blending mode. See details.}
#'   \item{nthreads}{- Number of threads to use.}
#'   \item{x}{- X coordinate of points.}
#'   \item{y}{- Y coordinate of points.}
#'   \item{radius}{- Relative spread of glow intensity. The radius should be proportional to the x and y-ranges of the plot. Values between 1/10 to 1/100 of the range of the plot generally produce good results.}
#'   \item{intensity}{- Maximum intensity at the center of a point.}
#'   \item{falloff_exponent}{- Exponent to determine how fast light intensity decreases from the point origin. A value of 0.5 corresponds to a linear falloff; a value of 2 corresponds to an inverse square. Generally you want this value to be high, otherwise you'll flood your plot with light.}
#'   \item{distance_exponent}{- Exponent of the distance calculation when calculating intensities. A value of 2 corresponds to euclidean distance; a value of 1 corresponds to manhattan distance.}
#'   \item{xlimits}{- The x-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{ylimits}{- The y-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{append}{ - Whether to add to the existing output or overwrite.}
#'   \item{saturation}{- When retrieving the output with \code{$output_raw} or \code{$output_dataframe}, maximum intensity values are capped at the given value. This is often useful when using additive blend mode to increase contrast.}
#' }
#'
#' @section Details:
#' \code{$new()} creates a new LightMapper object, which holds parameters, plotting data, and the output (a matrix of glow intensities). Creates a canvas to plot point data. With additive blending, the intensities of each point are added arithmetically, which is how light intensities are added in the physical world. This is equivalent to an fast/approximate un-normalized 2D kernel density estimate. 
#' 
#' #' With "screen" blending, two intensities are added according to the formula: \code{I_out = 1 - (1-I_a)*(1-I_b)}. Both additive blending and screen blending are commutative operations, meaning the order of points in a plot does not affect the output.
#' 
#' Note: Mapping "lights" (inverse power intensity) is much slower than "glow" effects (gaussian intensities) for various reasons. Plotting more than a few hundred points with \code{LightMapper} or \code{LightMapper4} may be computationally prohibitive. 
#' 
#' Screen blending can often improve contrast in a plot and is the default. 
#'
#' \code{$map()} maps points to the canvas. 
#' 
#' \code{$output_raw()} output raw matrix rasters. Useful for plotting in base R.
#' 
#' \code{$output_dataframe()} output the raster as a dataframe with XY coordinates. This is meant to pipe directly into ggplot.
#' 
#' \code{$aspect()}, \code{$xlim()}, \code{$ylim()} return the aspect ratio, x-limits and y-limits of the raster. These functions are intended to be used with plotting functions (e.g. \code{ggplot2::coord_fixed()}) so that the output raster is not distorted. See example below. 
#'
#' @examples
#' # Plot Data: x,y,r
#' x <- numeric(length=50)
#' y <- numeric(length=50)
#' r <- numeric(length=50)
#' for(t in 1:50) {
#'   xy <- exp(1i * t/2 - t/12)
#'   x[t] <- Re(xy)
#'   y[t] <- Im(xy)
#'   r[t] <- sqrt(x[t]^2 + y[t]^2)
#' }
#' 
#' # New class object
#' m <- LightMapper$new(xdim=500, ydim = 400, blend_mode = "screen")
#' 
#' # Map data on to raster
#' m$map(x=x, y=y, intensity = 1, radius = r/100, falloff_exponent = 0.5, distance_exponent = 2)
#' 
#' # Output raster data as a dataframe
#' pd <- m$output_dataframe(saturation = 1)
#' 
#' # Plot with ggplot
#' ggplot(pd, aes(x = x, y = y, fill = value)) + 
#'   geom_raster(show.legend = FALSE) +
#'   scale_fill_gradientn(colors=additive_alpha(c("black", "purple", "white"))) +
#'   coord_fixed(ratio = m$aspect(), xlim = m$xlim(), ylim = m$ylim(), expand = FALSE) + 
#'   theme_night(bgcolor = "black")
#' @name LightMapper
NULL

# LightMapper4 ###################################################

#' LightMapper4
#'
#'This class provides a framework for creating scatter plots based on a point light simulation  with explicit color intensities. Points are mapped with a inverse power gradient to a raster with specified dimensions and properties. 
#'
#' @section Usage:
#' \preformatted{m <- GlowMapper$new(xdim=1000, ydim=800, blend_mode = "additive", background_color = "#00000000", nthreads = 1)
#' 
#' m$map(x, y, radius, color = NULL, r=NULL, g=NULL, b=NULL, falloff_exponent = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE)
#' 
#' m$output_raw(saturation = 1, saturation_mode = "overflow")
#' 
#' m$output_dataframe(saturation = 1, saturation_mode = "overflow")
#' 
#' m$aspect()
#' 
#' m$xlim()
#' 
#' m$ylim()
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{xdim}{- The first dimension of the output matrix raster.}
#'   \item{ydim}{- The second dimension of the output matrix raster.}
#'   \item{blend_mode}{- Either screen or additive blending mode. See details.}
#'   \item{background_color}{- A color that can be coerced to RGBA with `col2rgb`, or a vector of four values between 0 and 1.}
#'   \item{nthreads}{- Number of threads to use.}
#'   \item{x}{- X coordinate of points.}
#'   \item{y}{- Y coordinate of points.}
#'   \item{radius}{- Relative spread of glow intensity. The radius should be proportional to the x and y-ranges of the plot. Values between 1/10 to 1/100 of the range of the plot generally produce good results.}
#'   \item{color}{- Color of points. If NULL, r, g, and b parameters must be defined (and vice versa).}
#'   \item{r}{- Red intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{g}{- Green intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{b}{- Blue intensity of points. Must be between 0 and 1 if using screen blending.}
#'   \item{falloff_exponent}{- Exponent to determine how fast light intensity decreases from the point origin. A value of 0.5 corresponds to a linear falloff; a value of 2 corresponds to an inverse square. Generally you want this value to be high, otherwise you'll flood your plot with light.}
#'   \item{distance_exponent}{- Exponent of the distance calculation when calculating intensities. A value of 2 corresponds to euclidean distance; a value of 1 corresponds to manhattan distance.}
#'   \item{xlimits}{- The x-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{ylimits}{- The y-limits of the output plot. If NA, the limits are +/- 5\% of the maximum and minimum points.}
#'   \item{append}{ - Whether to add to the existing output or overwrite.}
#'   \item{saturation}{- When retrieving the output with \code{$output_raw} or \code{$output_dataframe}, maximum intensity values are capped at the given value. This is often useful when using additive blend mode to increase contrast.}
#'   \item{saturation_mode}{- When intensity values are above the saturation threshold, values can be overflowed into other color channels ("overflow") or simply clipped at the threshold ("clip"). "Overflow" always produces a gradient to white for intensities above the threshold, which may produce artistically better results.}
#' }
#'
#' @section Details:
#' 
#' This `LightMapper4` class is similar to the `LightMapper` class, but instead of a single intensity matrix output, color is specified explicitly.
#' 
#' \code{$new()} creates a new LightMapper4 object, which holds parameters, plotting data, and the output (a matrix of glow intensities). Creates a canvas to plot point data. With additive blending, the intensities of each point are added arithmetically, which is how light intensities are added in the physical world. This is equivalent to an fast/approximate un-normalized 2D kernel density estimate. 
#' 
#' With "screen" blending,  two intensities are added according to the formula: \code{I_out = 1 - (1-I_a)*(1-I_b)}. Both additive blending and screen blending are commutative operations, meaning the order of points in a plot does not affect the output.
#' 
#' Screen blending can often improve contrast in a plot and is the default. 
#'
#' \code{$map()} maps points to the canvas. 
#' 
#' \code{$output_raw()} output raw matrix rasters (a list of four matrices, one for each RGBA channel). Useful for plotting in base R.
#' 
#' \code{$output_dataframe()} output the raster as a dataframe with XY coordinates. This is meant to pipe directly into ggplot.
#' 
#' \code{$aspect()}, \code{$xlim()}, \code{$ylim()} return the aspect ratio, x-limits and y-limits of the raster. These functions are intended to be used with plotting functions (e.g. \code{ggplot2::coord_fixed()}) so that the output raster is not distorted. See example below. 
#'
#' @examples
#' # Plot Data: x,y,r
#' x <- numeric(length=50)
#' y <- numeric(length=50)
#' r <- numeric(length=50)
#' color <- character(length=50)
#' for(t in 1:50) {
#'   xy <- exp(1i * t/2 - t/12)
#'   x[t] <- Re(xy)
#'   y[t] <- Im(xy)
#'   r[t] <- sqrt(x[t]^2 + y[t]^2)
#'   color[t] <- rgb(t/50,0,1-t/50)
#' }
#' 
#' # New class object
#' m <- LightMapper4$new(xdim=500, ydim = 400, blend_mode = "additive")
#' 
#' # Map data on to raster
#' m$map(x=x, y=y, color = color, radius = r/30+0.01, falloff_exponent = 1, distance_exponent = 2)
#' 
#' # Output raster data as a dataframe
#' pd <- m$output_dataframe(saturation = 1)
#' 
#' # Plot with ggplot
#' ggplot(pd, aes(x = x, y = y, fill = rgb(r,g,b,a))) + 
#'   geom_raster(show.legend = FALSE) +
#'   scale_fill_identity() +
#'   coord_fixed(ratio = m$aspect(), xlim = m$xlim(), ylim = m$ylim(), expand = FALSE) + 
#'   theme_night(bgcolor = "black")
#' @name LightMapper4
NULL


