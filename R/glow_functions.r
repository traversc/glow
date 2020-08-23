theme_night <- function (bgcolor = "black", base_size = 14, base_family = ""){
  theme_minimal() + 
    theme(
      line = element_line(color = "white", size = 0.5, linetype = 1, lineend = "butt"), 
      rect = element_rect(fill = "white", color = "white", size = 0.5, linetype = 1), 
      text = element_text(family = base_family, face = "plain", color = "white", size = base_size),
      plot.background = element_rect(color = bgcolor, fill = bgcolor),
      panel.background = element_rect(color = bgcolor, fill = bgcolor),
      panel.border = element_rect(fill = NA, color = "white"), 
      axis.text = element_text(color = "grey50"),
      panel.grid.major = element_line(color = "grey10", size = 0.3), 
      panel.grid.minor = element_line(color = "grey5", size = 0.3),
      strip.background = element_rect(fill = "grey30", color = "grey30")
    )
}

additive_alpha <- function(colors) {
  s <- seq(2, length(colors))
  x <- t(col2rgb(colors, alpha=F))/255
  cmax <- apply(x, 1, max)
  for(i in s) {
    x[i,] <- x[i,] - x[1,] * (1-cmax[i])
  }
  x[s,] <- x[s,] / cmax[s]
  x <- pmax(x, 0)
  x <- cbind(x, alpha = c(0,cmax[-1]))
  rgb(x[,1], x[,2], x[,3], x[,4])
}


# GlowMapper ###################################
GlowMapper <- R6Class("GlowMapper", list(
  # scene variables
  xdim = NULL,
  ydim = NULL,
  blend_mode = NULL,
  
  # plot data etc
  plot_data = NULL,
  xmin = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  x_aspect_ratio = NULL,
  y_aspect_ratio = NULL,
  
  # options
  contrast_limit = NULL,
  nthreads = NULL,
  
  output = NULL,
  
  initialize = function(xdim=1000, ydim=800, blend_mode = "screen", contrast_limit = 1e5, nthreads = 1) {
    stopifnot(is.numeric(xdim), length(xdim) == 1)
    stopifnot(is.numeric(ydim), length(ydim) == 1)
    stopifnot(is.character(blend_mode), length(blend_mode) == 1)
    stopifnot(blend_mode %in% c("screen", "additive"))
    stopifnot(is.numeric(contrast_limit), length(contrast_limit) == 1)
    stopifnot(is.numeric(nthreads), length(nthreads) == 1)
    self$xdim <- xdim
    self$ydim <- ydim
    self$blend_mode <- blend_mode
    self$contrast_limit <- contrast_limit
    self$nthreads <- nthreads
  },
  map = function(x, y, radius, intensity = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE) {
    stopifnot(is.numeric(x), length(x) >= 1)
    stopifnot(is.numeric(y), length(y) >= 1)
    stopifnot(is.numeric(intensity), length(intensity) >= 1)
    if(self$blend_mode == "screen") {
      if(any(intensity > 1 | intensity < 0)) stop("intensity should be between 0 and 1 if using screen blending")
    }
    stopifnot(is.numeric(radius), length(radius) >= 1)
    stopifnot(is.numeric(distance_exponent), length(distance_exponent) >= 1)
    
    stopifnot(is.numeric(xlimits), length(xlimits) == 2)
    stopifnot(is.numeric(ylimits), length(ylimits) == 2)
    
    self$plot_data <- data.frame(x, y, intensity, radius, distance_exponent)
    
    xdiff <- diff(range(x))
    ydiff <- diff(range(y))
    default_margin <- 0.05
    self$xmin <- ifelse(is.na(xlimits[1]), min(x) - xdiff * default_margin, xlimits[1])
    self$xmax <- ifelse(is.na(xlimits[2]), max(x) + xdiff * default_margin, xlimits[2])
    self$ymin <- ifelse(is.na(ylimits[1]), min(y) - ydiff * default_margin, ylimits[1])
    self$ymax <- ifelse(is.na(ylimits[2]), max(y) + ydiff * default_margin, ylimits[2])
    
    xincrement <- (self$xmax - self$xmin) / self$xdim
    yincrement <- (self$ymax - self$ymin) / self$ydim
    
    self$x_aspect_ratio <- max(xincrement / yincrement,1)
    self$y_aspect_ratio <- max(yincrement / xincrement,1)
    
    out <- c_map_glow(self$plot_data$x/self$x_aspect_ratio, 
                      self$plot_data$y/self$y_aspect_ratio, 
                      self$plot_data$intensity, self$plot_data$radius, self$plot_data$distance_exponent,
                      self$xdim, self$ydim, 
                      self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                      self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                      0, self$blend_mode, self$contrast_limit, self$nthreads)
    
    if(append == F || is.null(self$output)) {
      self$output <- out
    } else {
      if(self$blend_mode == "screen") {
        self$output <- 1 - (1-self$output)*(1-out)
      } else {
        self$output <- self$output + out
      }
    }
  },
  output_raw = function(saturation = NA_real_) {
    if(!is.na(saturation)) {
      out <- self$output
      out[] <- ifelse(out > saturation, saturation, out)
    } else {
      out <- self$output
    }
  },
  output_dataframe = function(saturation = NA_real_) {
    if(!is.na(saturation)) {
      out <- self$output
      out[] <- ifelse(out > saturation, saturation, out)
    } else {
      out <- self$output
    }
    
    xincrement <- (self$xmax - self$xmin)/self$xdim
    yincrement <- (self$ymax - self$ymin)/self$ydim
    
    # matrices are stored by column
    x <- xincrement * (seq(0,self$xdim-1)) + self$xmin
    y <- yincrement * (seq(0,self$ydim-1)) + self$ymin
    df <- expand.grid(x=x, y=y)
    df$value <- as.vector(out)
    df
  },
  aspect = function() {
    abs(self$xmax - self$xmin) / abs(self$ymax - self$ymin) * self$ydim / self$xdim
  },
  xlim = function() {
    c(self$xmin, self$xmax)
  },
  ylim = function() {
    c(self$ymin, self$ymax)
  }
))

# GlowMapper4 ###################################
GlowMapper4 <- R6Class("GlowMapper4", list(
  # scene variables
  xdim = NULL,
  ydim = NULL,
  blend_mode = NULL,
  
  # plot data etc
  background_color = NULL,
  plot_data = NULL,
  xmin = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  x_aspect_ratio = NULL,
  y_aspect_ratio = NULL,
  
  # options
  contrast_limit = NULL,
  nthreads = NULL,
  
  #output -- list of matrices per channel
  outputR = NULL,
  outputG = NULL,
  outputB = NULL,
  
  initialize = function(xdim=1000, ydim=800, blend_mode = "additive", background_color = "#00000000", contrast_limit = 1e5, nthreads = 1) {
    stopifnot(is.numeric(xdim), length(xdim) == 1)
    stopifnot(is.numeric(ydim), length(ydim) == 1)
    stopifnot(is.character(blend_mode), length(blend_mode) == 1)
    stopifnot(blend_mode %in% c("screen", "additive"))
    stopifnot(is.numeric(contrast_limit), length(contrast_limit) == 1)
    stopifnot(is.numeric(nthreads), length(nthreads) == 1)
    self$xdim <- xdim
    self$ydim <- ydim
    self$blend_mode <- blend_mode
    self$contrast_limit <- contrast_limit
    self$nthreads <- nthreads
    
    if(is.numeric(background_color)) {
      stopifnot(length(background_color) == 4)
      if(any(background_color > 1 | background_color < 0)) stop("background_color should be values between 0 and 1")
      self$background_color <- background_color
    } else {
      stopifnot(is.character(background_color))
      rgba <- grDevices::col2rgb(background_color, alpha=T)
      rgba[] <- rgba/255
      self$background_color <- rgba[,1]
    }
  },
  map = function(x, y, radius, color=NULL, r=NULL, g=NULL, b=NULL, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE) {
    stopifnot(is.numeric(x), length(x) >= 1)
    stopifnot(is.numeric(y), length(y) >= 1)
    stopifnot(is.numeric(radius), length(radius) >= 1)
    stopifnot(is.numeric(distance_exponent), length(distance_exponent) >= 1)
    
    if(is.null(color)) {
      if(is.null(r) || !is.numeric(r)) stop("must define color or rgb values")
      if(is.null(g) || !is.numeric(g)) stop("must define color or rgb values")
      if(is.null(b) || !is.numeric(b)) stop("must define color or rgb values")
      if(self$blend_mode == "screen") {
        if(any(r > 1 | r < 0)) stop("color should be between 0 and 1 if using screen blending")
        if(any(g > 1 | g < 0)) stop("color should be between 0 and 1 if using screen blending")
        if(any(b > 1 | b < 0)) stop("color should be between 0 and 1 if using screen blending")
      }
    } else {
      stopifnot(is.character(color))
      rgb <- grDevices::col2rgb(color, alpha=F)
      rgb[] <- rgb/255
      r <- rgb[1,]
      g <- rgb[2,]
      b <- rgb[3,]
    }
    
    stopifnot(is.numeric(xlimits), length(xlimits) == 2)
    stopifnot(is.numeric(ylimits), length(ylimits) == 2)
    
    self$plot_data <- data.frame(x, y, r,g,b, radius, distance_exponent)
    
    default_margin <- 0.05
    xdiff <- diff(range(x))
    ydiff <- diff(range(y))
    # if(xdiff == 0 || ydiff == 0) stop("set plot xlimits and ylimits")
    self$xmin <- ifelse(is.na(xlimits[1]), min(x) - xdiff * default_margin, xlimits[1])
    self$xmax <- ifelse(is.na(xlimits[2]), max(x) + xdiff * default_margin, xlimits[2])
    self$ymin <- ifelse(is.na(ylimits[1]), min(y) - ydiff * default_margin, ylimits[1])
    self$ymax <- ifelse(is.na(ylimits[2]), max(y) + ydiff * default_margin, ylimits[2])
    
    xincrement <- (self$xmax - self$xmin) / self$xdim
    yincrement <- (self$ymax - self$ymin) / self$ydim
    
    self$x_aspect_ratio <- max(xincrement / yincrement,1)
    self$y_aspect_ratio <- max(yincrement / xincrement,1)
    
    outR <- c_map_glow(self$plot_data$x/self$x_aspect_ratio, 
                       self$plot_data$y/self$y_aspect_ratio, 
                       self$plot_data$r, self$plot_data$radius, self$plot_data$distance_exponent,
                       self$xdim, self$ydim, 
                       self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                       self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                       self$background_color[1], self$blend_mode, self$contrast_limit, self$nthreads)
    outG <- c_map_glow(self$plot_data$x/self$x_aspect_ratio, 
                       self$plot_data$y/self$y_aspect_ratio, 
                       self$plot_data$g, self$plot_data$radius, self$plot_data$distance_exponent,
                       self$xdim, self$ydim, 
                       self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                       self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                       self$background_color[2], self$blend_mode, self$contrast_limit, self$nthreads)
    outB <- c_map_glow(self$plot_data$x/self$x_aspect_ratio, 
                       self$plot_data$y/self$y_aspect_ratio, 
                       self$plot_data$b, self$plot_data$radius, self$plot_data$distance_exponent,
                       self$xdim, self$ydim, 
                       self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                       self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                       self$background_color[3], self$blend_mode, self$contrast_limit, self$nthreads)
    
    if(append == F || is.null(self$outputR)) {
      self$outputR <- outR
      self$outputG <- outG
      self$outputB <- outB
    } else {
      if(self$blend_mode == "screen") {
        self$outputR <- 1 - (1-self$outputR)*(1-outR)
        self$outputG <- 1 - (1-self$outputG)*(1-outG)
        self$outputB <- 1 - (1-self$outputB)*(1-outB)
      } else {
        self$outputR <- self$outputR + outR
        self$outputG <- self$outputG + outG
        self$outputB <- self$outputB + outB
      }
    }
  },
  output_raw = function(saturation = 1, saturation_mode = "overflow") {
    stopifnot(saturation_mode %in% c(c("overflow", "clip", "none")))
    outR <- self$outputR
    outG <- self$outputG
    outB <- self$outputB
    if(saturation_mode == "overflow") {
      outR[] <- outR + ifelse(self$outputG > saturation, (self$outputG-saturation)/2, 0) + 
        ifelse(self$outputB > saturation, (self$outputB-saturation)/2, 0)
      outR[] <- ifelse(outR > saturation, saturation, outR)
      outG[] <- outG + ifelse(self$outputR > saturation, (self$outputR-saturation)/2, 0) + 
        ifelse(self$outputB > saturation, (self$outputB-saturation)/2, 0)
      outG[] <- ifelse(outG > saturation, saturation, outG)
      outB[] <- outB + ifelse(self$outputR > saturation, (self$outputR-saturation)/2, 0) + 
        ifelse(self$outputG > saturation, (self$outputG-saturation)/2, 0)
      outB[] <- ifelse(outB > saturation, saturation, outB)
    } else if(saturation_mode == "clip") {
      outR <- self$outputR
      outR[] <- ifelse(outR > saturation, saturation, outR)
      outG <- self$outputG
      outG[] <- ifelse(outG > saturation, saturation, outG)
      outB <- self$outputB
      outB[] <- ifelse(outB > saturation, saturation, outB)
    } else{
      saturation <- max(outR, outG, outB)
    }
    
    # if background is completely opaque, we don't need to handle alpha
    if(self$background_color[4] == 1) {
      outA <- matrix(1.0, nrow=nrow(outR), ncol=ncol(outR))
      return(list(r=outR,g=outG,b=outB, a=outA))
    }
    cmax <- pmax(outR, outG, outB) # opacity ~ max channel
    satmat <- matrix(saturation, nrow=nrow(outR), ncol=ncol(outR))
    alpha <- cmax/satmat # scale to [0,1]
    
    # scale to minimum alpha: if bg[4] = 0 ~ alpha; if bg[4] = 1 ~ 1 (opaque)
    #                         if bg[4] = 0.5 ~ transparency is cut in half
    alpha <- 1 - (1-alpha) * (1-self$background_color[4])
    
    outR[] <- ifelse(alpha == 0, 0, (outR - self$background_color[1] * (1-alpha))/alpha)
    outR <- pmin(pmax(outR, 0), saturation)
    outG[] <- ifelse(alpha == 0, 0, (outG - self$background_color[2] * (1-alpha))/alpha)
    outG <- pmin(pmax(outG, 0), saturation)
    outB[] <- ifelse(alpha == 0, 0, (outB - self$background_color[3] * (1-alpha))/alpha)
    outB[] <- pmin(pmax(outB, 0), saturation)
    return(list(r=outR, g=outG, b=outB, a=alpha))
  },
  output_dataframe = function(saturation = 1, saturation_mode = "overflow") {
    out <- self$output_raw(saturation, saturation_mode)
    
    xincrement <- (self$xmax - self$xmin)/self$xdim
    yincrement <- (self$ymax - self$ymin)/self$ydim
    
    # matrices are stored by column
    x <- xincrement * (seq(0,self$xdim-1)) + self$xmin
    y <- yincrement * (seq(0,self$ydim-1)) + self$ymin
    df <- expand.grid(x=x, y=y)
    df$r <- as.vector(out$r)
    df$g <- as.vector(out$g)
    df$b <- as.vector(out$b)
    df$a <- as.vector(out$a)
    df
  },
  aspect = function() {
    abs(self$xmax - self$xmin) / abs(self$ymax - self$ymin) * self$ydim / self$xdim
  },
  xlim = function() {
    c(self$xmin, self$xmax)
  },
  ylim = function() {
    c(self$ymin, self$ymax)
  }
))


# Lightmapper #########################################

LightMapper <- R6Class("LightMapper", list(
  # scene variables
  xdim = NULL,
  ydim = NULL,
  blend_mode = NULL,
  
  # plot data etc
  plot_data = NULL,
  xmin = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  x_aspect_ratio = NULL,
  y_aspect_ratio = NULL,
  
  # options
  nthreads = NULL,
  
  output = NULL,
  
  initialize = function(xdim=1000, ydim=800, blend_mode = "screen", nthreads = 1) {
    stopifnot(is.numeric(xdim), length(xdim) == 1)
    stopifnot(is.numeric(ydim), length(ydim) == 1)
    stopifnot(is.character(blend_mode), length(blend_mode) == 1)
    stopifnot(blend_mode %in% c("screen", "additive"))
    stopifnot(is.numeric(nthreads), length(nthreads) == 1)
    self$xdim <- xdim
    self$ydim <- ydim
    self$blend_mode <- blend_mode
    self$nthreads <- nthreads
  },
  map = function(x, y, intensity, radius, falloff_exponent = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE) {
    stopifnot(is.numeric(x), length(x) >= 1)
    stopifnot(is.numeric(y), length(y) >= 1)
    stopifnot(is.numeric(intensity), length(intensity) >= 1)
    stopifnot(is.numeric(radius), length(radius) >= 1)
    stopifnot(is.numeric(falloff_exponent), length(falloff_exponent) >= 1)
    stopifnot(is.numeric(distance_exponent), length(distance_exponent) >= 1)
    
    stopifnot(is.numeric(xlimits), length(xlimits) == 2)
    stopifnot(is.numeric(ylimits), length(ylimits) == 2)
    
    self$plot_data <- data.frame(x, y, intensity, radius, falloff_exponent, distance_exponent)
    
    xdiff <- diff(range(x))
    ydiff <- diff(range(y))
    default_margin <- 0.05
    self$xmin <- ifelse(is.na(xlimits[1]), min(x) - xdiff * default_margin, xlimits[1])
    self$xmax <- ifelse(is.na(xlimits[2]), max(x) + xdiff * default_margin, xlimits[2])
    self$ymin <- ifelse(is.na(ylimits[1]), min(y) - ydiff * default_margin, ylimits[1])
    self$ymax <- ifelse(is.na(ylimits[2]), max(y) + ydiff * default_margin, ylimits[2])
    
    xincrement <- (self$xmax - self$xmin) / self$xdim
    yincrement <- (self$ymax - self$ymin) / self$ydim
    
    self$x_aspect_ratio <- max(xincrement / yincrement,1)
    self$y_aspect_ratio <- max(yincrement / xincrement,1)
    
    out <- c_map_light(self$plot_data$x/self$x_aspect_ratio, 
                       self$plot_data$y/self$y_aspect_ratio, 
                       self$plot_data$intensity, self$plot_data$radius, 
                       self$plot_data$falloff_exponent, self$plot_data$distance_exponent,
                       self$xdim, self$ydim, 
                       self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                       self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                       0, self$blend_mode, self$nthreads)
    
    if(append == F || is.null(self$output)) {
      self$output <- out
    } else {
      if(self$blend_mode == "screen") {
        self$output <- 1 - (1-self$output)*(1-out)
      } else {
        self$output <- self$output + out
      }
    }
  },
  output_raw = function(saturation = NA_real_) {
    if(!is.na(saturation)) {
      out <- self$output
      out[] <- ifelse(out > saturation, saturation, out)
    } else {
      out <- self$output
    }
  },
  output_dataframe = function(saturation = NA_real_) {
    if(!is.na(saturation)) {
      out <- self$output
      out[] <- ifelse(out > saturation, saturation, out)
    } else {
      out <- self$output
    }
    
    xincrement <- (self$xmax - self$xmin)/self$xdim
    yincrement <- (self$ymax - self$ymin)/self$ydim
    
    # matrices are stored by column
    x <- xincrement * (seq(0,self$xdim-1)) + self$xmin
    y <- yincrement * (seq(0,self$ydim-1)) + self$ymin
    df <- expand.grid(x=x, y=y)
    df$value <- as.vector(out)
    df
  },
  aspect = function() {
    abs(self$xmax - self$xmin) / abs(self$ymax - self$ymin) * self$ydim / self$xdim
  },
  xlim = function() {
    c(self$xmin, self$xmax)
  },
  ylim = function() {
    c(self$ymin, self$ymax)
  }
))

# LightMapper4 ##########################################

LightMapper4 <- R6Class("GlowMapper4", list(
  # scene variables
  xdim = NULL,
  ydim = NULL,
  blend_mode = NULL,
  
  # plot data etc
  background_color = NULL,
  plot_data = NULL,
  xmin = NULL,
  xmax = NULL,
  ymin = NULL,
  ymax = NULL,
  x_aspect_ratio = NULL,
  y_aspect_ratio = NULL,
  
  # options
  contrast_limit = NULL,
  nthreads = NULL,
  
  #output -- list of matrices per channel
  outputR = NULL,
  outputG = NULL,
  outputB = NULL,
  
  initialize = function(xdim=1000, ydim=800, blend_mode = "additive", background_color = "#00000000", nthreads = 1) {
    stopifnot(is.numeric(xdim), length(xdim) == 1)
    stopifnot(is.numeric(ydim), length(ydim) == 1)
    stopifnot(is.character(blend_mode), length(blend_mode) == 1)
    stopifnot(blend_mode %in% c("screen", "additive"))
    stopifnot(is.numeric(nthreads), length(nthreads) == 1)
    self$xdim <- xdim
    self$ydim <- ydim
    self$blend_mode <- blend_mode
    self$nthreads <- nthreads
    
    if(is.numeric(background_color)) {
      stopifnot(length(background_color) == 4)
      self$background_color <- background_color
    } else {
      stopifnot(is.character(background_color))
      rgba <- grDevices::col2rgb(background_color, alpha=T)
      rgba[] <- rgba/255
      self$background_color <- rgba[,1]
    }
  },
  map = function(x, y, radius, color=NULL, r=NULL, g=NULL, b=NULL, falloff_exponent = 1, distance_exponent = 2, xlimits = c(NA_real_, NA_real_), ylimits = c(NA_real_, NA_real_), append = FALSE) {
    stopifnot(is.numeric(x), length(x) >= 1)
    stopifnot(is.numeric(y), length(y) >= 1)
    stopifnot(is.numeric(radius), length(radius) >= 1)
    stopifnot(is.numeric(distance_exponent), length(distance_exponent) >= 1)
    stopifnot(is.numeric(falloff_exponent), length(falloff_exponent) >= 1)
    
    if(is.null(color)) {
      if(is.null(r) || !is.numeric(r)) stop("must define color or rgb values")
      if(is.null(g) || !is.numeric(g)) stop("must define color or rgb values")
      if(is.null(b) || !is.numeric(b)) stop("must define color or rgb values")
    } else {
      stopifnot(is.character(color))
      rgb <- grDevices::col2rgb(color, alpha=F)
      rgb[] <- rgb/255
      r <- rgb[1,]
      g <- rgb[2,]
      b <- rgb[3,]
    }
    
    stopifnot(is.numeric(xlimits), length(xlimits) == 2)
    stopifnot(is.numeric(ylimits), length(ylimits) == 2)
    
    self$plot_data <- data.frame(x, y, r,g,b, radius, distance_exponent, falloff_exponent)
    
    default_margin <- 0.05
    xdiff <- diff(range(x))
    ydiff <- diff(range(y))
    # if(xdiff == 0 || ydiff == 0) stop("set plot xlimits and ylimits")
    self$xmin <- ifelse(is.na(xlimits[1]), min(x) - xdiff * default_margin, xlimits[1])
    self$xmax <- ifelse(is.na(xlimits[2]), max(x) + xdiff * default_margin, xlimits[2])
    self$ymin <- ifelse(is.na(ylimits[1]), min(y) - ydiff * default_margin, ylimits[1])
    self$ymax <- ifelse(is.na(ylimits[2]), max(y) + ydiff * default_margin, ylimits[2])
    
    xincrement <- (self$xmax - self$xmin) / self$xdim
    yincrement <- (self$ymax - self$ymin) / self$ydim
    
    self$x_aspect_ratio <- max(xincrement / yincrement,1)
    self$y_aspect_ratio <- max(yincrement / xincrement,1)
    
    outR <- c_map_light(self$plot_data$x/self$x_aspect_ratio, 
                        self$plot_data$y/self$y_aspect_ratio, 
                        self$plot_data$r, self$plot_data$radius, self$plot_data$falloff_exponent, self$plot_data$distance_exponent,
                        self$xdim, self$ydim, 
                        self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                        self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                        self$background_color[1], self$blend_mode, self$nthreads)
    outG <- c_map_light(self$plot_data$x/self$x_aspect_ratio, 
                        self$plot_data$y/self$y_aspect_ratio, 
                        self$plot_data$g, self$plot_data$radius, self$plot_data$falloff_exponent, self$plot_data$distance_exponent,
                        self$xdim, self$ydim, 
                        self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                        self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                        self$background_color[2], self$blend_mode, self$nthreads)
    outB <- c_map_light(self$plot_data$x/self$x_aspect_ratio, 
                        self$plot_data$y/self$y_aspect_ratio, 
                        self$plot_data$b, self$plot_data$radius, self$plot_data$falloff_exponent, self$plot_data$distance_exponent,
                        self$xdim, self$ydim, 
                        self$xmin/self$x_aspect_ratio, self$xmax/self$x_aspect_ratio, 
                        self$ymin/self$y_aspect_ratio, self$ymax/self$y_aspect_ratio, 
                        self$background_color[3], self$blend_mode, self$nthreads)
    
    if(append == F || is.null(self$outputR)) {
      self$outputR <- outR
      self$outputG <- outG
      self$outputB <- outB
    } else {
      if(self$blend_mode == "screen") {
        self$outputR <- 1 - (1-self$outputR)*(1-outR)
        self$outputG <- 1 - (1-self$outputG)*(1-outG)
        self$outputB <- 1 - (1-self$outputB)*(1-outB)
      } else {
        self$outputR <- self$outputR + outR
        self$outputG <- self$outputG + outG
        self$outputB <- self$outputB + outB
      }
    }
  },
  output_raw = function(saturation = 1, saturation_mode = "overflow") {
    stopifnot(saturation_mode %in% c(c("overflow", "clip", "none")))
    outR <- self$outputR
    outG <- self$outputG
    outB <- self$outputB
    if(saturation_mode == "overflow") {
      outR[] <- outR + ifelse(self$outputG > saturation, (self$outputG-saturation)/2, 0) + 
        ifelse(self$outputB > saturation, (self$outputB-saturation)/2, 0)
      outR[] <- ifelse(outR > saturation, saturation, outR)
      outG[] <- outG + ifelse(self$outputR > saturation, (self$outputR-saturation)/2, 0) + 
        ifelse(self$outputB > saturation, (self$outputB-saturation)/2, 0)
      outG[] <- ifelse(outG > saturation, saturation, outG)
      outB[] <- outB + ifelse(self$outputR > saturation, (self$outputR-saturation)/2, 0) + 
        ifelse(self$outputG > saturation, (self$outputG-saturation)/2, 0)
      outB[] <- ifelse(outB > saturation, saturation, outB)
    } else if(saturation_mode == "clip") {
      outR <- self$outputR
      outR[] <- ifelse(outR > saturation, saturation, outR)
      outG <- self$outputG
      outG[] <- ifelse(outG > saturation, saturation, outG)
      outB <- self$outputB
      outB[] <- ifelse(outB > saturation, saturation, outB)
    } else{
      saturation <- max(outR, outG, outB)
    }
    
    # if background is completely opaque, we don't need to handle alpha
    if(self$background_color[4] == 1) {
      outA <- matrix(1.0, nrow=nrow(outR), ncol=ncol(outR))
      return(list(r=outR,g=outG,b=outB, a=outA))
    }
    cmax <- pmax(outR, outG, outB) # opacity ~ max channel
    satmat <- matrix(saturation, nrow=nrow(outR), ncol=ncol(outR))
    alpha <- cmax/satmat # scale to [0,1]
    
    # scale to minimum alpha: if bg[4] = 0 ~ alpha; if bg[4] = 1 ~ 1 (opaque)
    #                         if bg[4] = 0.5 ~ transparency is cut in half
    alpha <- 1 - (1-alpha) * (1-self$background_color[4])
    
    outR[] <- ifelse(alpha == 0, 0, (outR - self$background_color[1] * (1-alpha))/alpha)
    outR <- pmin(pmax(outR, 0), saturation)
    outG[] <- ifelse(alpha == 0, 0, (outG - self$background_color[2] * (1-alpha))/alpha)
    outG <- pmin(pmax(outG, 0), saturation)
    outB[] <- ifelse(alpha == 0, 0, (outB - self$background_color[3] * (1-alpha))/alpha)
    outB[] <- pmin(pmax(outB, 0), saturation)
    return(list(r=outR, g=outG, b=outB, a=alpha))
  },
  output_dataframe = function(saturation = 1, saturation_mode = "overflow") {
    out <- self$output_raw(saturation, saturation_mode)
    
    xincrement <- (self$xmax - self$xmin)/self$xdim
    yincrement <- (self$ymax - self$ymin)/self$ydim
    
    # matrices are stored by column
    x <- xincrement * (seq(0,self$xdim-1)) + self$xmin
    y <- yincrement * (seq(0,self$ydim-1)) + self$ymin
    df <- expand.grid(x=x, y=y)
    df$r <- as.vector(out$r)
    df$g <- as.vector(out$g)
    df$b <- as.vector(out$b)
    df$a <- as.vector(out$a)
    df
  },
  aspect = function() {
    abs(self$xmax - self$xmin) / abs(self$ymax - self$ymin) * self$ydim / self$xdim
  },
  xlim = function() {
    c(self$xmin, self$xmax)
  },
  ylim = function() {
    c(self$ymin, self$ymax)
  }
))
