#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#if RCPP_PARALLEL_USE_TBB
#include <tbb/task_arena.h>
#endif

using namespace Rcpp;
using namespace RcppParallel;

// https://stackoverflow.com/a/4609795/2723734
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export(rng = false)]]
bool is_tbb() {
#if RCPP_PARALLEL_USE_TBB
  return true;
#endif
  return false;
}

/////////////////////////////////////////////////////////

struct GlowMapper {
  const double xmin;
  const double ymin;
  const size_t xdim;
  const size_t ydim;
  const double xincrement;
  const double yincrement;
  const double contrast_limit;
  const Eigen::VectorXd xvec;
  const Eigen::RowVectorXd yvec;
  GlowMapper(const size_t xdim, const size_t ydim,
             const double xmin, const double xmax,
             const double ymin, const double ymax, const double contrast_limit) : 
    xmin(xmin), ymin(ymin), xdim(xdim), ydim(ydim),
    xincrement((xmax-xmin) / static_cast<double>(xdim-1)),
    yincrement((ymax-ymin) / static_cast<double>(ydim-1)),
    contrast_limit(log(contrast_limit)),
    xvec(transform_x(xdim, xmin, xmax)), yvec(transform_y(ydim, ymin, ymax)) {}
  
  inline static Eigen::VectorXd transform_x(const size_t xdim, const double xmin, const double xmax) {
    double xincrement = (xmax-xmin) / static_cast<double>(xdim-1);
    Eigen::VectorXd x(xdim);
    for(size_t i=0; i<xdim; ++i) {
      x(i) = xmin + xincrement * i;
    }
    return x;
  }
  inline static Eigen::RowVectorXd transform_y(const size_t ydim, const double ymin, const double ymax) {
    double yincrement = (ymax-ymin) / static_cast<double>(ydim-1);
    Eigen::RowVectorXd y(ydim);
    for(size_t i=0; i<ydim; ++i) {
      // y(i) = ymin + yincrement * (ydim - i - 1);
      y(i) = ymin + yincrement * i;
    }
    return y;
  }
  inline double nearest_x_idx(const double x) const {
    return  (x - xmin)/xincrement ; // this may be less than zero, since a point may be outside plotting range
  }
  inline double nearest_y_idx(const double y) const {
    return (y - ymin)/yincrement ;
  }
  inline size_t bound_x(long x) const {
    return std::min(std::max(x,0L),static_cast<long>(xdim)-1L);
  }
  inline size_t bound_y(long y) const {
    return std::min(std::max(y,0L),static_cast<long>(ydim)-1L);
  }
  void screen_update(Eigen::MatrixXd & mat, const double x, const double y, const double intensity, const double radius, const double exponent) const {
    double search_distance;
    if(exponent >= 1) {
      search_distance = pow( pow(radius, exponent) * (contrast_limit + log(intensity)), 1.0/exponent );
    } else {
      search_distance = 2*pow( pow(radius, exponent) * (contrast_limit + log(intensity)) / 2.0, 1.0/exponent );
    }
    double search_distance_x = search_distance / xincrement;
    double search_distance_y = search_distance / yincrement;
    size_t xblock_min = bound_x(lrint(nearest_x_idx(x) - search_distance_x));
    size_t xblock_size = bound_x(lrint(nearest_x_idx(x) + search_distance_x)) - xblock_min;
    size_t yblock_min = bound_y(lrint(nearest_y_idx(y) - search_distance_y));
    size_t yblock_size = bound_y(lrint(nearest_y_idx(y) + search_distance_y)) - yblock_min;
    
    // std::cout << search_distance << " xblock_min: " << xblock_min << " yblock_min: " <<
    //   yblock_min << " xblock_size: " << xblock_size << " yblock_size: " << yblock_size << std::endl;
    
    mat.block(xblock_min, yblock_min, xblock_size, yblock_size).array() *= 
      1.0-(
          (((xvec.segment(xblock_min, xblock_size).array() - x))/pow(radius,2.0)).abs().pow(exponent).exp().inverse().matrix() *
            (((yvec.segment(yblock_min, yblock_size).array() - y))/pow(radius,2.0)).abs().pow(exponent).exp().inverse().matrix() * intensity
      ).array();
  }
  void additive_update(Eigen::MatrixXd & mat, const double x, const double y, const double intensity, const double radius, const double exponent) const {
    double search_distance;
    if(exponent >= 1) {
      search_distance = pow( pow(radius, exponent) * (contrast_limit + log(intensity)), 1.0/exponent );
    } else {
      search_distance = 2*pow( pow(radius, exponent) * (contrast_limit + log(intensity)) / 2.0, 1.0/exponent );
    }
    double search_distance_x = search_distance / xincrement;
    double search_distance_y = search_distance / yincrement;
    size_t xblock_min = bound_x(lrint(nearest_x_idx(x) - search_distance_x));
    size_t xblock_size = bound_x(lrint(nearest_x_idx(x) + search_distance_x)) - xblock_min;
    size_t yblock_min = bound_y(lrint(nearest_y_idx(y) - search_distance_y));
    size_t yblock_size = bound_y(lrint(nearest_y_idx(y) + search_distance_y)) - yblock_min;
    
    // std::cout << search_distance << " xblock_min: " << xblock_min << " yblock_min: " <<
    //   yblock_min << " xblock_size: " << xblock_size << " yblock_size: " << yblock_size << std::endl;
    
    mat.block(xblock_min, yblock_min, xblock_size, yblock_size) += 
      (((xvec.segment(xblock_min, xblock_size).array() - x))/pow(radius,2.0)).abs().pow(exponent).exp().inverse().matrix() *
      (((yvec.segment(yblock_min, yblock_size).array() - y))/pow(radius,2.0)).abs().pow(exponent).exp().inverse().matrix() * intensity;
  }
};

#if RCPP_PARALLEL_USE_TBB
struct GlowWorker : public Worker {
  GlowMapper & gm;
  const std::string blend_mode; // we should probably not use a string here. Do blend mode with polymorphism?
  const double * const x;
  const double * const y;
  const double * const intensity;
  const double * const radius;
  const double * const exponent;
  Eigen::MatrixXd output;
  GlowWorker(GlowMapper & gm,
             const std::string blend_mode,
             const double * const x,
             const double * const y,
             const double * const intensity,
             const double * const radius,
             const double * const exponent) : 
    gm(gm), blend_mode(blend_mode), x(x), y(y), intensity(intensity), radius(radius), exponent(exponent), 
    output(Eigen::MatrixXd::Constant(gm.xdim, gm.ydim, blend_mode == "screen" ? 1.0 : 0.0)) {}
  GlowWorker(const GlowWorker & w, Split) : 
    gm(w.gm), blend_mode(w.blend_mode), x(w.x), y(w.y), intensity(w.intensity), radius(w.radius), exponent(w.exponent), 
    output(Eigen::MatrixXd::Constant(gm.xdim, gm.ydim, blend_mode == "screen" ? 1.0 : 0.0)) {}
  
  void join(const GlowWorker & rhs) {
    if(blend_mode == "screen") {
      output.array() *= rhs.output.array();
    } else {
      output.array() += rhs.output.array();
    }
  }
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      if(blend_mode == "screen") {
        gm.screen_update(output, x[i], y[i], intensity[i], radius[i], exponent[i]);
      } else {
        gm.additive_update(output, x[i], y[i], intensity[i], radius[i], exponent[i]);
      }
    }
  }
};
#endif

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd c_map_glow(NumericVector x, NumericVector y, NumericVector intensity, NumericVector radius, NumericVector exponent,
                    const size_t xdim, const size_t ydim,
                    const double xmin, const double xmax,
                    const double ymin, const double ymax,
                    const double background, const std::string blend_mode, const double contrast_limit, const int nthreads) {
  
  GlowMapper gm(xdim, ydim, xmin, xmax, ymin, ymax, contrast_limit);
  size_t len = Rf_xlength(x);
  
#if RCPP_PARALLEL_USE_TBB
  if(nthreads > 1) {
    GlowWorker w(gm, blend_mode, REAL(x), REAL(y), REAL(intensity), REAL(radius), REAL(exponent));
    parallelReduce(0, len, w, 100, nthreads);
    if(blend_mode == "screen") {
      w.output.array() = 1.0 - w.output.array() * (1.0-background);
    } else {
      w.output.array() += background;
    }
    return w.output;
  }
#endif
  
  Eigen::MatrixXd output;
  if(blend_mode == "screen") {
    output = Eigen::MatrixXd::Constant(xdim, ydim, 1-background);
    for(size_t i=0; i<len; ++i) {
      gm.screen_update(output, x[i], y[i], intensity[i], radius[i], exponent[i]);
    }
    output.array() = 1.0 - output.array();
  } else { // additive
    output = Eigen::MatrixXd::Constant(xdim, ydim, background);
    for(size_t i=0; i<len; ++i) {
      gm.additive_update(output, x[i], y[i], intensity[i], radius[i], exponent[i]);
    }
  }
  return output;
}

struct LightMapper {
  const double xmin;
  const double ymin;
  const size_t xdim;
  const size_t ydim;
  const double xincrement;
  const double yincrement;
  const Eigen::ArrayXXd xarr;
  const Eigen::ArrayXXd yarr;
  LightMapper(const size_t xdim, const size_t ydim,
              const double xmin, const double xmax,
              const double ymin, const double ymax) : 
    xmin(xmin), ymin(ymin), xdim(xdim), ydim(ydim),
    xincrement((xmax-xmin) / static_cast<double>(xdim-1)),
    yincrement((ymax-ymin) / static_cast<double>(ydim-1)),
    xarr(transform_x(xdim, ydim, xmin, xmax)), yarr(transform_y(xdim, ydim, ymin, ymax)) {}
  
  inline static Eigen::ArrayXXd transform_x(const size_t xdim, const size_t ydim, const double xmin, const double xmax) {
    double xincrement = (xmax-xmin) / static_cast<double>(xdim-1);
    Eigen::VectorXd x(xdim);
    for(size_t i=0; i<xdim; ++i) {
      x(i) = xmin + xincrement * i;
    }
    Eigen::ArrayXXd xa(xdim, ydim);
    xa.matrix().colwise() = x;
    return xa;
  }
  inline static Eigen::ArrayXXd transform_y(const size_t xdim, const size_t ydim, const double ymin, const double ymax) {
    double yincrement = (ymax-ymin) / static_cast<double>(ydim-1);
    Eigen::RowVectorXd y(ydim);
    for(size_t i=0; i<ydim; ++i) {
      y(i) = ymin + yincrement * i;
    }
    Eigen::ArrayXXd ya(xdim, ydim);
    ya.matrix().rowwise() = y;
    return ya;
  }
  void additive_update(Eigen::ArrayXXd & arr, const double x, const double y, const double intensity, const double radius, 
                       const double falloff_exponent, const double distance_exponent) const {
    arr += (((xarr-x)).abs().pow(distance_exponent)/(pow(radius,2.0)) + ((yarr-y)).abs().pow(distance_exponent)/(pow(radius,2.0)) + 1).pow(falloff_exponent).inverse() * intensity;
  }
  void screen_update(Eigen::ArrayXXd & arr, const double x, const double y, const double intensity, const double radius, 
                     const double falloff_exponent, const double distance_exponent) const {
    arr *= 1.0 - (
      (((xarr-x)).abs().pow(distance_exponent)/(pow(radius,2.0)) + ((yarr-y)).abs().pow(distance_exponent)/(pow(radius,2.0)) + 1).pow(falloff_exponent).inverse() * intensity
    );
  }
};

#if RCPP_PARALLEL_USE_TBB
struct LightWorker : public Worker {
  LightMapper & lm;
  const std::string blend_mode;
  const double * const x;
  const double * const y;
  const double * const intensity;
  const double * const radius;
  const double * const falloff_exponent;
  const double * const distance_exponent;
  Eigen::ArrayXXd output;
  LightWorker(LightMapper & lm,
             const std::string blend_mode,
             const double * const x,
             const double * const y,
             const double * const intensity,
             const double * const radius,
             const double * const falloff_exponent,
             const double * const distance_exponent) : 
    lm(lm), blend_mode(blend_mode), x(x), y(y), intensity(intensity), radius(radius), 
    falloff_exponent(falloff_exponent), distance_exponent(distance_exponent), 
    output(Eigen::ArrayXXd::Constant(lm.xdim, lm.ydim, blend_mode == "screen" ? 1.0 : 0.0)) {}
  LightWorker(const LightWorker & w, Split) : 
    lm(w.lm), blend_mode(w.blend_mode), x(w.x), y(w.y), intensity(w.intensity), radius(w.radius), 
    falloff_exponent(w.falloff_exponent), distance_exponent(w.distance_exponent), 
    output(Eigen::ArrayXXd::Constant(lm.xdim, lm.ydim, blend_mode == "screen" ? 1.0 : 0.0)) {}
  
  void join(const LightWorker & rhs) {
    if(blend_mode == "screen") {
      output *= rhs.output;
    } else {
      output += rhs.output;
    }
  }
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      if(blend_mode == "screen") {
        lm.screen_update(output, x[i], y[i], intensity[i], radius[i], falloff_exponent[i], distance_exponent[i]);
      } else {
        lm.additive_update(output, x[i], y[i], intensity[i], radius[i], falloff_exponent[i], distance_exponent[i]);
      }
    }
  }
};
#endif

// [[Rcpp::export(rng = false)]]
Eigen::ArrayXXd c_map_light(NumericVector x, NumericVector y, NumericVector intensity, NumericVector radius, 
                     NumericVector falloff_exponent, NumericVector distance_exponent,
                     const size_t xdim, const size_t ydim,
                     const double xmin, const double xmax,
                     const double ymin, const double ymax,
                     const double background, const std::string blend_mode, const int nthreads) {
  LightMapper lm(xdim, ydim, xmin, xmax, ymin, ymax);
  size_t len = Rf_xlength(x);
  
#if RCPP_PARALLEL_USE_TBB
  if(nthreads > 1) {
    LightWorker w(lm, blend_mode, REAL(x), REAL(y), REAL(intensity), REAL(radius), REAL(falloff_exponent), REAL(distance_exponent));
    parallelReduce(0, len, w, 100, nthreads);
    if(blend_mode == "screen") {
       w.output = 1.0 - w.output * (1.0-background);
    } else {
       w.output += background;
    }
    return w.output;
  }
#endif
  
  Eigen::ArrayXXd output;
  if(blend_mode == "screen") {
    output = Eigen::ArrayXXd::Constant(xdim, ydim, 1-background);
    for(size_t i=0; i<len; ++i) {
      lm.screen_update(output, x[i], y[i], intensity[i], radius[i], falloff_exponent[i], distance_exponent[i]);
    }
    output = 1.0 - output;
  } else {
    output = Eigen::ArrayXXd::Constant(xdim, ydim, background);
    for(size_t i=0; i<len; ++i) {
      lm.additive_update(output, x[i], y[i], intensity[i], radius[i], falloff_exponent[i], distance_exponent[i]);
    }
  }
  return output;
}

// https://en.wikipedia.org/wiki/Talk:Mollweide_projection
// This formula has no citation but it seems to work
inline double mollweide_newton_rhapson(const double lambda, const double phi, const size_t n_iter) {
  double theta = ( M_PI / 2.0 - pow(3.0 * M_PI / 8.0 * pow((M_PI/2.0 - abs(phi)),2.0), 1.0/3.0 ) ) * sgn(phi);
  if(phi > 1.570762) return theta;
  for(size_t i=0; i<n_iter; ++i) {
    theta -= (2.0 * theta + sin(2.0*theta) - M_PI * sin(phi)) / (2.0 + 2.0 * cos(2.0 * theta));
  }
  return theta;
}
// [[Rcpp::export(rng = false)]]
DataFrame mollweide_projection(NumericVector latitude, NumericVector longitude, const double meridian) {
  size_t N = Rf_xlength(latitude);
  NumericVector x(N);
  NumericVector y(N);
  for(size_t i=0; i<N; ++i) {
    double phi = latitude[i];
    double lambda = longitude[i];
    double theta = mollweide_newton_rhapson(lambda, phi, 3);
    x[i] = 2.0 * sqrt(2.0) / M_PI * (lambda-meridian) * cos(theta);
    y[i] = sqrt(2.0) * sin(theta);
  }
  return DataFrame::create(_["x"] = x, _["y"] = y);
}

// 1.886,-2.357,-0.328, 0.918
// [[Rcpp::export(rng = false)]]
DataFrame clifford_attractor(size_t n_iter, double A=1.886, double B=-2.357, double C=-0.328, double D=0.918, double x0=0.1, double y0=0) {
  NumericVector x(n_iter);
  NumericVector y(n_iter);
  NumericVector angle(n_iter);
  NumericVector distance(n_iter);
  double * xp = REAL(x);
  double * yp = REAL(y);
  double * anglep = REAL(angle);
  double * distancep = REAL(distance);
  xp[0] = x0;
  yp[0] = y0;
  anglep[0] = 0;
  distancep[0] = 0;
  for (size_t i=1; i < n_iter; ++i) {
    xp[i] = (sin(A*yp[i-1]) + C*cos(A*xp[i-1]));
    yp[i] = (sin(B*xp[i-1]) + D*cos(B*yp[i-1]));
    anglep[i] = atan2(yp[i], xp[i]);
    double ydiff = yp[i] - yp[i-1];
    double xdiff = xp[i] - xp[i-1];
    distancep[i] = sqrt(ydiff*ydiff + xdiff*xdiff);
  }
  return DataFrame::create(_["x"]=x, _["y"]=y, _["angle"]=angle, _["distance"]=distance);
}


