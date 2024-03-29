// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// is_tbb
bool is_tbb();
RcppExport SEXP _glow_is_tbb() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(is_tbb());
    return rcpp_result_gen;
END_RCPP
}
// c_map_glow
Eigen::MatrixXd c_map_glow(NumericVector x, NumericVector y, NumericVector intensity, NumericVector radius, NumericVector exponent, const size_t xdim, const size_t ydim, const double xmin, const double xmax, const double ymin, const double ymax, const double background, const std::string blend_mode, const double contrast_limit, const int nthreads);
RcppExport SEXP _glow_c_map_glow(SEXP xSEXP, SEXP ySEXP, SEXP intensitySEXP, SEXP radiusSEXP, SEXP exponentSEXP, SEXP xdimSEXP, SEXP ydimSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP yminSEXP, SEXP ymaxSEXP, SEXP backgroundSEXP, SEXP blend_modeSEXP, SEXP contrast_limitSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exponent(exponentSEXP);
    Rcpp::traits::input_parameter< const size_t >::type xdim(xdimSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ydim(ydimSEXP);
    Rcpp::traits::input_parameter< const double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const double >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const double >::type ymin(yminSEXP);
    Rcpp::traits::input_parameter< const double >::type ymax(ymaxSEXP);
    Rcpp::traits::input_parameter< const double >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< const std::string >::type blend_mode(blend_modeSEXP);
    Rcpp::traits::input_parameter< const double >::type contrast_limit(contrast_limitSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_map_glow(x, y, intensity, radius, exponent, xdim, ydim, xmin, xmax, ymin, ymax, background, blend_mode, contrast_limit, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// c_map_light
Eigen::ArrayXXd c_map_light(NumericVector x, NumericVector y, NumericVector intensity, NumericVector radius, NumericVector falloff_exponent, NumericVector distance_exponent, const size_t xdim, const size_t ydim, const double xmin, const double xmax, const double ymin, const double ymax, const double background, const std::string blend_mode, const int nthreads);
RcppExport SEXP _glow_c_map_light(SEXP xSEXP, SEXP ySEXP, SEXP intensitySEXP, SEXP radiusSEXP, SEXP falloff_exponentSEXP, SEXP distance_exponentSEXP, SEXP xdimSEXP, SEXP ydimSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP yminSEXP, SEXP ymaxSEXP, SEXP backgroundSEXP, SEXP blend_modeSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type falloff_exponent(falloff_exponentSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distance_exponent(distance_exponentSEXP);
    Rcpp::traits::input_parameter< const size_t >::type xdim(xdimSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ydim(ydimSEXP);
    Rcpp::traits::input_parameter< const double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const double >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const double >::type ymin(yminSEXP);
    Rcpp::traits::input_parameter< const double >::type ymax(ymaxSEXP);
    Rcpp::traits::input_parameter< const double >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< const std::string >::type blend_mode(blend_modeSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_map_light(x, y, intensity, radius, falloff_exponent, distance_exponent, xdim, ydim, xmin, xmax, ymin, ymax, background, blend_mode, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// mollweide_projection
DataFrame mollweide_projection(NumericVector latitude, NumericVector longitude, const double meridian);
RcppExport SEXP _glow_mollweide_projection(SEXP latitudeSEXP, SEXP longitudeSEXP, SEXP meridianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type longitude(longitudeSEXP);
    Rcpp::traits::input_parameter< const double >::type meridian(meridianSEXP);
    rcpp_result_gen = Rcpp::wrap(mollweide_projection(latitude, longitude, meridian));
    return rcpp_result_gen;
END_RCPP
}
// clifford_attractor
DataFrame clifford_attractor(size_t n_iter, double A, double B, double C, double D, double x0, double y0);
RcppExport SEXP _glow_clifford_attractor(SEXP n_iterSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP x0SEXP, SEXP y0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< size_t >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type y0(y0SEXP);
    rcpp_result_gen = Rcpp::wrap(clifford_attractor(n_iter, A, B, C, D, x0, y0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_glow_is_tbb", (DL_FUNC) &_glow_is_tbb, 0},
    {"_glow_c_map_glow", (DL_FUNC) &_glow_c_map_glow, 15},
    {"_glow_c_map_light", (DL_FUNC) &_glow_c_map_light, 15},
    {"_glow_mollweide_projection", (DL_FUNC) &_glow_mollweide_projection, 3},
    {"_glow_clifford_attractor", (DL_FUNC) &_glow_clifford_attractor, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_glow(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
