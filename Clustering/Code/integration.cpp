// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
#include <Rcpp.h>
#include <cmath>

double m_hat_cpp(double x, double b, Rcpp::NumericVector data_p, 
             Rcpp::NumericVector grid_p, double bw){
  double result = 0.0;
  double norm   = 0.0;
  int len = data_p.size();
  for (int i = 0; i < len; i++) {
    if (abs((grid_p[i] - x * b) / bw) <= 1.0){
      result += data_p[i];
      norm += 1.0;
    }
  }
  return (result / norm);
}
  
class Integrand: public Func
{
private:
  const double b;
  Rcpp::NumericVector data_p;
  Rcpp::NumericVector grid_p;
  const double bw;
public:
  Integrand(const double& b_, Rcpp::NumericVector& data_p_,
            Rcpp::NumericVector& grid_p_,
            const double& bw_) : b(b_), data_p(data_p_), grid_p(grid_p_), bw(bw_) {}
  
  double operator()(const double& x) const
  {
    return m_hat_cpp(x, b, data_p, grid_p, bw);
  }
};

class Integrand2: public Func
{
private:
  const double b;
  Rcpp::NumericVector data_p_1;
  Rcpp::NumericVector data_p_2;
  const double norm_1;
  const double norm_2;
  Rcpp::NumericVector grid_p;
  const double bw;
public:
  Integrand2(const double& b_, Rcpp::NumericVector& data_p_1_,
             Rcpp::NumericVector& data_p_2_, const double& norm_1_,
             const double& norm_2_, Rcpp::NumericVector& grid_p_, 
             const double& bw_) : b(b_), data_p_1(data_p_1_),
             data_p_2(data_p_2_), norm_1(norm_1_), norm_2(norm_2_),
             grid_p(grid_p_), bw(bw_) {}
  
  double operator()(const double& x) const
  {
    double tmp = m_hat_cpp(x, b, data_p_1, grid_p, bw) / (norm_1 * 1 / b) -
      m_hat_cpp(x, 1.0, data_p_2, grid_p, bw) / (norm_2 * 1 / b);
    return tmp * tmp;
  }
};

// [[Rcpp::export]]
Rcpp::List integrate1_cpp(double b, Rcpp::NumericVector data_points,
                          Rcpp::NumericVector grid_points, double bw, int subdiv)
{
  const double lower = 0.0, upper = 1 / b;

  Integrand f(b, data_points, grid_points, bw);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code, subdiv);
  return Rcpp::List::create(
    Rcpp::Named("res") = res,
    Rcpp::Named("err_est") = err_est,
    Rcpp::Named("err_code") = err_code
  );
}

// [[Rcpp::export]]
Rcpp::List integrate2_cpp(double b, Rcpp::NumericVector data_points_1,
                         Rcpp::NumericVector data_points_2,
                         double norm_1, double norm_2,
                         Rcpp::NumericVector grid_points, double bw, int subdiv)
{
  const double lower = 0.0, upper = 1 / b;
  
  Integrand2 f(b, data_points_1, data_points_2, norm_1, norm_2, grid_points, bw);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code, subdiv);
  return Rcpp::List::create(
    Rcpp::Named("res") = res,
    Rcpp::Named("err_est") = err_est,
    Rcpp::Named("err_code") = err_code
  );
}